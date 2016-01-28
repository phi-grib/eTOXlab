# -*- coding: utf-8 -*-

##    Description    eTOXlab model template
##                   
##    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
##
##    Copyright 2013 Manuel Pastor
##
##    This file is part of eTOXlab.
##
##    eTOXlab is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation version 3.
##
##    eTOXlab is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with eTOXlab.  If not, see <http://www.gnu.org/licenses/>.

from model import model

import os
import shutil
import subprocess
import cPickle as pickle
import time

from model import model
from utils import removefile
from utils import updateProgress
from utils import writeError

import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy import cluster
from scipy import spatial


class imodel(model):
    def __init__ (self, vpath):
        model.__init__(self, vpath)

        ##
        ## General settings
        ##
        self.buildable = True
        self.quantitative = True
        self.confidential = False
        self.identity = False
        self.experimental = False
        self.SDFileName = 'name'
        self.SDFileActivity = 'activity'
        self.SDFileExperimental = ''
        
        ##
        ## Normalization settings
        ##
        self.norm = True
        self.normStand = True
        self.normNeutr = True
        self.normNeutrMethod = 'moka'
        self.normNeutr_pH = 7.4
        self.norm3D = True

        ##
        ## Molecular descriptor settings
        ##

        #DEBUG SAVE-TIME TRICK
        self.MD = 'adriana'                       # 'padel'|'pentacle'|'adriana'
        self.padelMD = ['-2d']                    # '-2d'|'-3d'
        self.padelMaxRuntime = None
        self.padelDescriptor = None
        self.pentacleProbes = ['DRY','O','N1','TIP']
        self.pentacleOthers = []

        ##
        ## Modeling settings
        ##
        self.model = 'pls'
        self.modelLV = 2
        self.modelAutoscaling = True
        self.modelCutoff = 'auto'
        self.selVar = False
        #self.selVarMethod = GOLPE
        self.selVarLV = 2
        #self.selVarCV = 'LOO'
        self.selVarRun = 2
        self.selVarMask = None

        ##
        ## View settings
        ##
        self.viewType = 'property'    # 'pca' | 'property' | 'project' 
        self.viewBackground = False
        self.viewReferenceEndpoint = None
        self.viewReferenceVersion = 0

        self.plotPCAColor = 'red'
        self.plotPCAMarkerShape = 'D'
        self.plotPCAMarkerSize = 40
        self.plotPCAMarkerLine = 0
        
        self.plotPRPColor = 'red'
        self.plotPRPMarkerShape = 'D'
        self.plotPRPMarkerSize = 40
        self.plotPRPMarkerLine = 0
        
        self.plotPRJColor = 'DModX'    # DModX | [color] (e.g. red)
        self.plotPRJMarkerShape = 'o'
        self.plotPRJMarkerSize = 50
        self.plotPRJMarkerLine = 1
        
        self.plotBGColor = '#aaaaaa'
        self.plotBGMarkerShape = 'o'
        self.plotBGMarkerSize = 20
        self.plotBGMarkerLine = 0

        ##
        ## Path to external programs
        ##
        self.mokaPath = '/opt/blabber/blabber4eTOX/'
        self.padelPath = '/opt/padel/padel218ws/'
        self.padelURL = 'http://localhost:9000/computedescriptors?params=' 
        self.pentaclePath = '/opt/pentacle/pentacle106eTOX/'
        self.adrianaPath = '/opt/AdrianaCode/AdrianaCode226/'
        self.corinaPath = '/opt/corina/corina3494/'
        self.javaPath = '/usr/java/jdk1.7.0_51/'
        self.RPath = '/opt/R/R-3.0.2/'
        self.standardiserPath = '/opt/standardiser/standardise20140206/'

        ## FLEXY
        self.flexyData = []
        self.flexyMorganFingerprintRatius = 2
        self.flexyClusteringCutoff = 0.6
        self.flexyChemClassSDCutoff = 4
        self.flexyChemClassNCutoff = 5

    def buildFlexy (self):
 
        X,Y = self.getMatrices ()

        nobj, nvarx = np.shape(X)
        if (nobj==0) or (nvarx==0) : return (False, 'failed to extract activity or to generate MD')

        nobj = np.shape(Y)
        if (nobj==0) : return (False, 'no activity found')

        #FLEXY: process Muckro scaffolds removing entries with "", "c1ccccc1" and those withs more too wide activity range

        stdcutoff = np.std(Y) / self.flexyChemClassSDCutoff
        blacklist = ['','c1ccccc1']
        scaffolds = {}
        i=0
        for vi in self.flexyData:
            si = Chem.MolToSmiles(vi[0])  ## scaffold
            if si in blacklist:
                i = i+1
                continue
            
            if not si in scaffolds:
                scaffolds[si]=([],[])
                
            scaffolds[si][0].append(i)
            i = i+1

        for i in scaffolds:
            v = scaffolds[i][0]
            w = scaffolds[i][1]
            w.append(np.median(Y[v]))
            w.append(np.std(Y[v]))
            w.append(len(v) <= self.flexyChemClassNCutoff and w[1] < stdcutoff)

        print scaffolds
        
        #FLEXY: run clustering and obtain centroids
        # obtain a distance matrix calling

        nmols = len(self.flexyData)
        moldistance = np.zeros ((nmols,nmols),dtype=np.float64)

        # compute distance using Morgan fingerprints
        i=0
        for vi in self.flexyData:
            j=0
            for vj in self.flexyData:
                # only upper triangle
                if i<j :
                    # use the default distance metric: Tanimoto
                    moldistance [i,j] = 1.00 - DataStructs.FingerprintSimilarity(vi[1],vj[1])
                j=j+1
            i=i+1

        # transform to reduced distance matrix format needed by 
        cdistance = spatial.distance.squareform (moldistance, checks=False)
        
        Z = cluster.hierarchy.complete(cdistance)

        # extract clusters of compounds showing internal distances no higher than 10%
        # this is a very strict criteria, guaranteeing that clusters contains very similar compounds
        cindex = cluster.hierarchy.fcluster(Z, t=self.flexyClusteringCutoff, criterion='distance')

        print cindex

        # analyze custer index
        chemclusters = {}
        i=0
        for c in cindex:
            if not c in chemclusters:
                # include two lists:
                # [0] for the compound indexes
                # [1] for cluster description (median, sd, local-like)
                chemclusters[c]=([],[])
        
            chemclusters[c][0].append(i) # insert the index in the first list
            i=i+1

        print "numclusters is:", len(chemclusters), "stdcutoff is:", stdcutoff

        BASEDIR = '/home/modeler/soft/eTOXlab/src/'
        itag = self.vpath.split('/')[-2]

        nLocalModels = 0
        for i in chemclusters:
            v = chemclusters[i][0]
            w = chemclusters[i][1]
            w.append(np.median(Y[v]))
            w.append(np.std(Y[v]))
            w.append(len(v) > self.flexyChemClassNCutoff and w[1] > stdcutoff)

            #print chemclusters[i]

            if w[2] :
                # build a sdfile with compounds in cluster
                print 'local model for ', i, w[1]

                finp = open (self.vpath+'/training.sdf', 'r')
                fout = open (self.vpath+'/piece%0.4d.sdf' % i,'w') 
                mindex = 0
                for line in finp:
                    if mindex in v:               
                        fout.write(line)
        
                    if "$$$$" in line:
                        mindex=mindex+1

                finp.close()
                fout.close()
                        
                # build a localmodel (local%d,i)
                ndir = self.vpath+'/local%0.4d' % i
        
                if os.path.isdir (ndir):
                    shutil.rmtree (ndir, ignore_errors=True)
            
                os.mkdir (ndir)
            
                call =['/usr/bin/python', BASEDIR+'build.py','-e',itag,
                       '-f', self.vpath+'/piece%0.4d.sdf' % i, '-m', self.vpath+'/imodel.py',
                       '-s', str(i)]

                retcode = subprocess.call(call)

                if retcode !=0 : return (False, 'error in local model computation')

                nLocalModels = nLocalModels + 1
                print 'local model complete!'

        ## save everything needed for prediction
        try:
            f = open (self.vpath+'/flexydata.pkl','wb')
        except:
            return (False, 'error saving flexy data')

        morgans = []
        for vi in self.flexyData:
            morgans.append(vi[1])
                
        pickle.dump(morgans, f)            # Morgan fingerprints
        pickle.dump(cindex, f)             # cluster index
        pickle.dump(scaffolds,f)           # Muckro scaffolds (dictionary)
        pickle.dump(chemclusters,f)        # flat clusters (dictionary)
        f.close()

        self.infoNotes = []
        self.infoNotes.append ( ('flexy','%d scaffolds' %len(scaffolds)) )
        self.infoNotes.append ( ('flexy','%d clusters' %len(chemclusters)) )
        self.infoNotes.append ( ('flexy','%d local models' %nLocalModels) )
           
        return (True, 'success')


    def computeMurckoScaffold(self, molFile):
        suppl = Chem.SDMolSupplier(molFile)
        mol = suppl.next()
        try:
            MS = MurckoScaffold.GetScaffoldForMol(mol)
        except:
            MS = None

        #debug
        #print (Chem.MolToSmiles(MS))

        return ( (True, MS) )
                    
    def computeMorganFingerprint(self, molFile):
        suppl = Chem.SDMolSupplier(molFile)
        mol = suppl.next()
        try:
            # note: we are using a radius of two 2, which produces
            # results similar to those obtained with PP ECFP4
            MF = AllChem.GetMorganFingerprintAsBitVect(mol,2)
        except:
            MF = None
        
        return ( (True, MF) )

    def buildWorkflow(self, molecules):

        # if this is a submodel
        if self.vpath [-9:-4] == 'local':
            return (model.buildWorkflow (self, molecules))
        
        success, result = self.licenseTesting ()
        if not success: return (False, result)
        
        if not self.buildable:
            success, result = self.log ()
            if not success:
                return (False, result)
            return (False, 'Non-buildable model')

        # load data, if stored, or compute it from the provided SDFile

        dataReady = False

        if not molecules:

            # FLEXY: ***TODO*** add new loadData for extra arrays
            dataReady = self.loadData ()

            if not self.loadSeriesInfo ():
                self.setSeries ('training.sdf', len(self.tdata))

        if not dataReady: # datList was not completed because load failed or new series was set

            # estimate number of molecules inside the SDFile
            nmol = 0
            try:
                f = open (self.vpath+'/training.sdf','r')
            except:
                return (False,"Unable to open file %s" % molecules)
            for line in f:
                if '$$$$' in line: nmol+=1
            f.close()

            if not nmol:
                return (False,"No molecule found in %s:  SDFile format not recognized" % molecules)

            self.setSeries (molecules, nmol)

            print (molecules, nmol)  # DEBUG ONLY. REMOVE!!!!

            i = 0
            fout = None
            mol = ''

            # open SDFfile and iterate for every molecule
            f = open (self.vpath+'/training.sdf','r')

            # clean normalized structures
            removefile (self.vpath+'/tstruct.sdf')
            
            updateProgress (0.0)

            for line in f:
                if not fout or fout.closed:
                    i += 1
                    mol = 'm%0.10d.sdf' % i
                    fout = open(mol, 'w')

                fout.write(line)

                if '$$$$' in line:
                    fout.close()

                    ## workflow for molecule i (mol) ############
                    success, result = self.normalize (mol)
                    if not success:
                       writeError('error in normalize: '+result)
                       continue

                    molFile   = result[0]
                    molName   = result[1]
                    molCharge = result[2]
                    molPos    = self.saveNormalizedMol(molFile)

                    ## FLEXY: compute Muckro scaffolds
                    success, result = self.computeMurckoScaffold(molFile)
                    if not success:
                       writeError('error in Murcko scaffold computation')
                    molMS = result
                    
                    ## FLEXY: compute compute Morgan
                    success, result = self.computeMorganFingerprint(molFile)
                    if not success:
                       writeError('error in Morgan fingerprint')
                    molMF = result

                    self.flexyData.append ( (molMS, molMF) )

                    success, infN = self.extract (molFile,molName,molCharge,molPos)
                    if not success:
                       writeError('error in extract: '+ str(infN))
                       continue

                    updateProgress (float(i)/float(nmol))
                    ##############################################

                    removefile (mol)

            f.close()
            if fout :
                fout.close()

            self.saveData ()

        # build the model with the datList stored data

        success, result = self.buildFlexy ()  ## <- FLEXY
        if not success:
            return (False, result)

        print 'building regular model'
        success, result = self.build ()
        if not success:
            return (False, result)

        print 'log'
        success, result = self.log ()
        if not success:
            return (False, result)

        return (result)
    

    def predictFlexy (self, molFile, molName, molCharge, detail, morgans, cindex, scaffolds, clusters):

        #print molName
        
        predScaffold = ( (False,0.0), (False,0), (False,0.0) )
        predMorgan   = ( (False,0.0), (False,0), (False,0.0) )
        predLocal    = ( (False,0.0), (False,0), (False,0.0) )
        predGlobal   = ( (False,0.0), (False,0), (False,0.0) )

        # compute scaffold and Morgan fingerprints for query compound
        suppl = Chem.SDMolSupplier(molFile)
        mol = suppl.next()
        try:
            iscaffold = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
        except:
            return (False, 'error computing scaffolds')
        try:
            imorgan = AllChem.GetMorganFingerprintAsBitVect(mol,2)
        except:
            return (False, 'error computing fingerprint')

        # compare Muckro scaffolds
        if iscaffold in scaffolds:
            v = scaffolds[iscaffold][1]
            if v[2]:
                predScaffold = ( (True,v[0]), (True,0), (True,1.96*v[1]) )
        

        # compare Morgan distance
        dmin = 1.00
        imin = 1
        i = 0
        for m in morgans: 
            d = 1.0 - DataStructs.FingerprintSimilarity(imorgan,m)
            if d < dmin:
                dmin = d
                imin = i
            i=i+1

        icluster = cindex[imin]
        
        # find closest cluster
        
        v = clusters[icluster][1]

        #print v

        if dmin < self.flexyClusteringCutoff:  # if the point if farther than 0.6 it don't belongs to the class
            if v[2]:  # call local model

                # guess endpoint tag and version from vpath
                itag = self.vpath.split('/')[-2]
                ver = int (self.vpath.split('/')[-1][-4:])
                BASEDIR = '/home/modeler/soft/eTOXlab/src/'
                
                call =['/usr/bin/python', BASEDIR+'predict.py','-e',itag, '-v', str(ver),
                       '-f', molFile, '-s', str(icluster)]

                #print call

                retcode = subprocess.call(call)

                if retcode != 0 : return (False, 'prediction computation failed')
                
                f = open ('results.pkl','rb')
                iresult = pickle.load(f)
                f.close()

                if not iresult[0]:
                    return (False, iresult[1])
                else:
                    predLocal = iresult[1][0][1]

            else: # estimate using the median and sd      
                predMorgan = ( (True,v[0]), (True,0), (True,1.96*v[1]) )
            

        # global prediction

        predGlobal = self.predict (molFile, molName, molCharge, detail)

        # analyze results
        
        #print 'scaffold: ', predScaffold
        #print 'morgan:   ', predMorgan
        #print 'local:    ', predLocal
        #print 'global:   ', predGlobal

        table = []
        if predScaffold[2][0]: table.append (predScaffold[2][1])
        else: table.append (1.e10)
        if predMorgan[2][0]: table.append (predMorgan[2][1])
        else: table.append (1.e10)
        if predLocal[2][0]: table.append (predLocal[2][1])
        else: table.append (1.e10)
        if predGlobal[2][0]: table.append (predGlobal[2][1])
        else: table.append (1.e10)

        #print table

        tableMin = min(table)
        
        if table[0] == tableMin : return (predScaffold)
        elif table[1] == tableMin : return (predMorgan)
        elif table[2] == tableMin : return (predLocal)
        elif table[3] == tableMin : return (predGlobal)
      
        return (False, 'something went terribly wrong!')
    

    def predictWorkflow(self, molecules, detail, progress):

        # if this is a submodel
        if self.vpath [-9:-4] == 'local':
            return (model.predictWorkflow (self, molecules, detail, progress))

        #print 'prediction flexy', self.vpath

        success, result = self.licenseTesting ()
        if not success: return (False, result)
        
        datList = []
        datList = self.loadData ()
        
        i=0
        pred = []
        mol=''
        fout = None

        # open SDFfile and iterate for every molecule
        try:
            f = open (molecules,'r')
        except:
            return (False,"No molecule found in %s; SDFile format not recognized" % molecules)

        try:
            fflexy = open (self.vpath+'/flexydata.pkl','rb')
        except:
            return (False, 'error loading flexy data')
                
        morgans   = pickle.load(fflexy)   # Morgan fingerprints
        cindex    = pickle.load(fflexy)   # cluster index
        scaffolds = pickle.load(fflexy)   # Muckro scaffolds (dictionary)
        clusters  = pickle.load(fflexy)   # flat clusters (dictionary)
        fflexy.close()
        
        for line in f:
            if not fout or fout.closed:
                i += 1
                mol = 'm%0.10d.sdf' % i
                fout = open(mol, 'w')

            fout.write(line)

            if '$$$$' in line:
                fout.close()

                ## workflow for molecule i (mol) ###########
                success, result  = self.normalize (mol)
                if not success:
                    pred.append((False, result))
                    continue

                molFile   = result[0]
                molName   = result[1]
                molCharge = result[2]

                predN = self.predictFlexy (molFile, molName, molCharge, detail, morgans, cindex, scaffolds, clusters)

                pred.append((True, predN))
                ############################################

                if progress:
                    sys.stdout.write('completed: %d\n'%i)
                    sys.stdout.flush()

                removefile(mol)

        return (True, pred)
