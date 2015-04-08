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

from utils import wkd
from utils import updateProgress

import sys
import os
import shutil
import numpy as np
import subprocess
import cPickle as pickle
import time

from rdkit import Chem
from rdkit.Chem import Descriptors
from utils import removefile

class imodel(model):
    def __init__ (self, vpath):
        model.__init__(self, vpath)

        ##
        ## General settings
        ##
        self.buildable = True
        self.quantitative = True
        self.confidential = False
        self.identity = True
        self.experimental = True
        self.SDFileName = 'name'
        self.SDFileActivity = 'activity'
        self.SDFileExperimental = 'IKr'
        
        ##
        ## Normalization settings
        ##
        self.norm = True
        self.normStand = True
        self.normNeutr = True
        self.normNeutrMethod = 'moka'
        self.normNeutr_pH = 6.3
        self.norm3D = True

        ##
        ## Molecular descriptor settings
        ##
        self.MD = 'pentacle'                         # 'padel'|'pentacle'|'adriana'
        self.padelMD = ['-2d']                       # '-2d'|'-3d'
        self.padelMaxRuntime = None
        self.padelDescriptor = None
        self.pentacleProbes = ['DRY','O','N1','TIP']
        self.pentacleOthers = []

        ##
        ## Modeling settings
        ##
        self.model = 'pls'
        self.modelLV = 3
        self.modelAutoscaling = False
        self.modelCutoff = 'auto'
        self.selVar = True
        #self.selVarMethod = GOLPE
        self.selVarLV = 2
        #self.selVarCV = 'LOO'
        self.selVarRun = 2
        self.selVarMask = None

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
        self.standardiserPath = '/opt/standardise/standardise20140206/'
        
    def splitSet (self, molecules, nchunks):

        # open molecules and extract MW for every mol
        mw = []
        suppl = Chem.SDMolSupplier (molecules)
        for moli in suppl:
            mw.append(Descriptors.MolWt(moli))
            
        # compute cutoffs
        pstep = 100/nchunks
        vlow = []
        vhig = []
        fset = []
        for i in range(nchunks):
            vlow.append( np.percentile(mw, pstep*i) )
            vhig.append( np.percentile(mw, pstep*(i+1)) )
            fset.append( open (self.vpath+'/piece%0.4d.sdf' % i,'w') )

        fpickle = open (self.vpath+'/cutoffs.pkl','wb')
        pickle.dump (nchunks, fpickle)
        pickle.dump (vlow, fpickle)
        pickle.dump (vhig, fpickle)
        fpickle.close()
        
        # open molecules and, according to mw, dumpt to appropriate piece
        finp = open (molecules, 'r')

        i = 0
        assigned = False
        
        for line in finp:
            if not assigned:
                for j in range(nchunks):
                    # the top compound will not be selected, since mw == vhigh,
                    # but the it will be correctly assigned to the last piece nonetheless
                    if (   (mw[i] >= vlow[j]) & (mw[i] < vhig[j])   ):
                        break
                assigned=True
                
            fset[j].write(line)
            
            if "$$$$" in line:
                i=i+1
                assigned=False

        finp.close()
        for fi in fset: fi.close()
                         
        return (True,len(mw))


    def buildWorkflow (self, molecules):
        
        if not self.buildable:
            success, result = self.log ()
            if not success:
                return (False, result)
            return (result)

        if not molecules:
            molecules=self.vpath+'/training.sdf'

        # if this is a submodel
        if self.vpath [-9:-4] == 'local':
            return (model.buildWorkflow (self, molecules))

        # if this is the top model
        BASEDIR = '/home/modeler/soft/eTOXlab/src/'
        result = ''

        nchunks = 4
        
        itag = self.vpath.split('/')[-2]
                
        # split original SDFile and create 'nchunks' files called 'piece0000.sdf'...
        success, results = self.splitSet(molecules, nchunks)
        if not success:
            return (False, result)

        # make 'nchunk' calls to 'build' command using the respective pieces

        exper  = np.zeros(results, dtype=np.float64)
        recalc = np.zeros(results, dtype=np.float64)
        predic = np.zeros(results, dtype=np.float64)
        rcount = 0
        pcount = 0
        
        for ichunk in range (nchunks):

            print 'LOCAL MODEL %d' % ichunk

            ndir = self.vpath+'/local%0.4d' % ichunk
            
            if (molecules) :    
                if os.path.isdir (ndir):
                    shutil.rmtree (ndir, ignore_errors=True)
                os.mkdir (ndir)
                    
                call =['/usr/bin/python', BASEDIR+'build.py','-e',itag,
                       '-f', self.vpath+'/piece%0.4d.sdf' % ichunk, '-m', self.vpath+'/imodel.py',
                       '-s', str(ichunk)]
            else:
                if not os.path.isdir (ndir):
                    return (False, 'local models are empty, please provide training series again')
                
                ver = int (self.vpath.split('/')[-1][-4:])
                
                call =['/usr/bin/python', BASEDIR+'build.py','-e',itag,  
                       '-v', str(ver), '-m', self.vpath+'/imodel.py',
                       '-s', str(ichunk)]

            retcode = subprocess.call(call)

            if retcode !=0 : return (False, 'error in computation')

            # collect and accumulate

            f = open (ndir+'/pls-recalculated.txt','r')
            header = True
            for line in f:
                if header :
                    header = False
                    continue

                exper [rcount] = float(line.split()[1])
                recalc[rcount] = float(line.split()[-1])
                rcount+=1
                
            f.close()

            f = open (ndir+'/pls-predicted.txt','r')
            header = True
            for line in f:
                if header :
                    header = False
                    continue

                predic[pcount] = float(line.split()[-1])
                pcount+=1
                
            f.close()

        SSYp=np.sum(np.square(predic-exper))
        SSYr=np.sum(np.square(recalc-exper))

        emean = 0.00
        for i in range(pcount):
            emean+=exper[i]
        emean/=pcount

        SSY0 = 0.00
        for i in range(pcount):
            SSY0 += np.square(exper[i]-emean)
        
        R2 = 1.00-(SSYr/SSY0)
        Q2 = 1.00-(SSYp/SSY0)
        SDEP = np.sqrt(SSYp/pcount)

        print R2, Q2, SDEP

        self.infoResult = []    
        self.infoResult.append( ('nobj',pcount) )
        self.infoResult.append( ('R2','%5.3f' % R2) )
        self.infoResult.append( ('Q2','%5.3f' % Q2) )
        self.infoResult.append( ('SDEP','%5.3f' % SDEP) )

        
        success, result = self.log ()
        if not success:
            return (False, result)
        
        return (result)
    

    def splitQuery (self,molecules):

        fpickle = open (self.vpath+'/cutoffs.pkl','rb')       
        nchunks = pickle.load(fpickle)
        vlow    = pickle.load(fpickle)
        vhig    = pickle.load(fpickle)
        fpickle.close()

        fset  = []
        order = []
        for i in range(nchunks):
##            fset.append( open (self.vpath+'/query%0.4d.sdf' % i,'w') )
            fset.append( open ('/var/tmp/query%0.4d.sdf' % i,'w') )
            lista = []
            order.append([])

        i = 0
        suppl = Chem.SDMolSupplier (molecules)
        for moli in suppl:
            mwi = Descriptors.MolWt(moli)
            assigned=False
            for j in range(nchunks):
                if (   (mwi >= vlow[j]) & (mwi < vhig[j])   ):
                    assigned=True
                    break
            if not assigned:
                if mwi < vlow[0] :
                    j = 0
                else:
                    j = nchunks-1

            fset[j].write(Chem.MolToMolBlock(moli))

            if self.SDFileExperimental:
                if moli.HasProp(self.SDFileExperimental):
                    exp = moli.GetProp(self.SDFileExperimental)
                    fset[j].write('>  <'+self.SDFileExperimental+'>\n'+exp+'\n\n$$$$')
            
            fset[j].write('$$$$\n')
            order[j].append(i)
            i = i+1
            
        for fi in fset: fi.close()
        
        return (True, 'Split OK', order)
    

    def predictWorkflow (self, molecules, detail, progress):
        
        # if this is a submodel
        if self.vpath [-9:-4] == 'local':
            return (model.predictWorkflow (self, molecules, detail, progress))

        # if this is the top model
        BASEDIR = '/home/modeler/soft/eTOXlab/src/'
        result = ''

        fpickle = open (self.vpath+'/cutoffs.pkl','rb')       
        nchunks = pickle.load(fpickle)
        fpickle.close()

        # guess endpoint tag and version from vpath
        itag = self.vpath.split('/')[-2]
        ver = int (self.vpath.split('/')[-1][-4:])

        # split original SDFile and create 'nchunks' files called 'piece0000.sdf'...        
        success, results, resultsOrder = self.splitQuery(molecules)
        if not success:
            return (False, result)
        
        aggregatedResults = []
        sys.path.append ('/opt/RDKit/')
        sys.path.append ('/opt/standardiser/standardise20140206/')
        
        # make 'nchunk' calls to 'build' command using the respective pieces  
        for ichunk in range (nchunks):
            
            #if not os.path.isfile(self.vpath+'/query%0.4d.sdf' % ichunk):
            if not os.path.isfile('/var/tmp/query%0.4d.sdf' % ichunk):
                continue

            #print 'PREDICTION IN MODEL %d' % ichunk
                
##            call =['/usr/bin/python', BASEDIR+'predict.py','-e',itag, '-v', str(ver),
##                   '-f', self.vpath+'/query%0.4d.sdf' % ichunk, '-s', str(ichunk)]
            call =['/usr/bin/python', BASEDIR+'predict.py','-e',itag, '-v', str(ver),
                   '-f', '/var/tmp/query%0.4d.sdf' % ichunk, '-s', str(ichunk)]

            retcode = subprocess.call(call)

            removefile('/var/tmp/query%0.4d.sdf' % ichunk)
            
            if retcode != 0 : return (False, 'prediction computation failed')

            f = open ('results.pkl','rb')
            iresult = pickle.load(f)
            f.close()

            if not iresult[0]:
                return (False, iresult[1]) 
            
            for oi,iri in zip(resultsOrder[ichunk],iresult[1]):
                aggregatedResults.append( (oi,iri) )
                
        #process output and reorder the results
        
        aggregatedResults.sort()
        
        results = []
        for ri in aggregatedResults:
            results.append(ri[1])
        
        return (True, results)
    
    def log (self):
        """Stores information about the Series, the MD, the Model and the Results to the info.pkl file

           This information us used only for inform the user and not to recover the status of the model 
        """
        self.infoID = []
        self.infoID.append (('version', '*'))
        self.infoID.append (('date', time.asctime(time.localtime(time.time()))))

        self.infoID.append (('buildable', str(self.buildable) ))
        if self.quantitative:
            self.infoID.append (('dependent', 'quantitative'))
        else:
            self.infoID.append (('dependent', 'qualitative'))

        self.infoID.append (('confident', str(self.confidential)))
            
        self.infoMD = []

        if 'pentacle' in self.MD:
            self.infoMD.append( ('MD','Pentacle') )
            for probe in self.pentacleProbes:
                self.infoMD.append ( ('probe',probe) )
            for key in self.pentacleOthers:
                self.infoMD.append ( ('key',key) )
        elif 'padel' in self.MD:
            self.infoMD.append( ('MD','PaDEL') )
            if self.padelDescriptor:
                self.infoMD.append( ('descriptors', self.padelDescriptor) )
            if self.padelMaxRuntime:
                self.infoMD.append( ('max runtime', str(self.padelMaxRuntime)) )
        elif 'adriana' in self.MD:
            self.infoMD.append( ('MD','Adriana') )
            
        try:
            modelInfo = open (self.vpath+'/info.pkl','wb')
        except:
            return (False, 'Failed to write model log')
                
        pickle.dump(self.infoID, modelInfo)
        pickle.dump(self.infoSeries, modelInfo)
        pickle.dump(self.infoMD, modelInfo)
        pickle.dump(self.infoModel, modelInfo)
        pickle.dump(self.infoResult, modelInfo)

        modelInfo.close()
        
        return (True, "Model OK")
    
    def loadSeriesInfo (self):
        """Gets information about the series used to build the model (model.infoSeries) stored in file info.pkl
           This information is used only for informative purposes and does not change the model properties
        """
        try:
            modelInfo = open (self.vpath+'/info.pkl','rb')
        except:
            return (False)
        infoID          = pickle.load(modelInfo)
        self.infoSeries = pickle.load(modelInfo)
        modelInfo.close()
        return (True)
    
    def loadData (self):
        """Gets all model settings stored for the last model and compares with current settings.
           In case of any disagreement, a False is returned.

           Please note that this method does not intend to set up the model settings (__init__)
        """

        if self.confidential:
            return False
        
        if not os.path.isfile (self.vpath+'/tdata.pkl'):
            return False

        try:
            f = open (self.vpath+'/tdata.pkl','rb')
        except:
            return False
        
        norm = pickle.load(f)
        if norm != self.norm:
            return False

        if norm:
            normStand = pickle.load(f)
            if normStand != self.normStand:
                return False
            
            normNeutr = pickle.load(f)
            if normNeutr != self.normNeutr:
                return False

            if normNeutr:
                normNeutrMethod = pickle.load(f)
                if normNeutrMethod != self.normNeutrMethod:
                    return False

                normNeutr_pH = pickle.load(f)
                if normNeutr_pH != self.normNeutr_pH:
                    return False
                
            norm3D = pickle.load(f)
            if norm3D != self.norm3D:
                return False
        
        MD = pickle.load(f)
        if MD != self.MD:
            return False
        
        if 'pentacle' in MD:
            pentacleProbes = pickle.load(f)
            if pentacleProbes != self.pentacleProbes:
                return False

            pentacleOthers = pickle.load(f)
            if pentacleOthers != self.pentacleOthers:
                return False

        elif 'padel' in MD:
            padelMD = pickle.load(f)
            if padelMD != self.padelMD:
                return False
        
        self.tdata = pickle.load(f)
        f.close()
