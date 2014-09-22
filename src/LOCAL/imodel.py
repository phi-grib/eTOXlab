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
        self.identity = False
        self.SDFileName = 'name'
        self.SDFileActivity = 'activity'

        
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
        self.modelAutoscaling = True
        self.modelCutoff = 'auto'
        self.selVar = False
        #self.selVarMethod = GOLPE
        self.selVarLV = 2
        #self.selVarCV = 'LOO'
        self.selVarRun = 1
        self.selVarMask = None

        ##
        ## Path to external programs
        ##
        self.mokaPath = '/opt/blabber/blabber110/'
        self.padelPath = '/opt/padel/padel218ws/'
        self.padelURL = 'http://localhost:9000/computedescriptors?params=' 
        self.pentaclePath = '/opt/pentacle/pentacle106/'
        self.adrianaPath = '/opt/AdrianaCode/AdrianaCode226/'
        self.corinaPath = '/opt/corina/corina24/'
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
                         
        return (True,'splitting OK')


    def buildWorkflow (self, molecules):
        
        if not self.buildable:
            success, result = self.log ()
            if not success:
                return (False, result)
            return (result)

        # if this is a submodel
        if self.vpath [-9:-4] == 'local':
            return (model.buildWorkflow (self, molecules))

        # if this is the top model
        BASEDIR = '/home/modeler/soft/eTOXlab/src/'
        result = ''

        nchunks = 2
        
        itag = self.vpath.split('/')[-2]

        # split original SDFile and create 'nchunks' files called 'piece0000.sdf'...
        if (molecules) :
            success, results = self.splitSet(molecules, nchunks)
            if not success:
                return (False, result)

        # make 'nchunk' calls to 'build' command using the respective pieces  
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
            fset.append( open (self.vpath+'/query%0.4d.sdf' % i,'w') )
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
        
        # make 'nchunk' calls to 'build' command using the respective pieces  
        for ichunk in range (nchunks):

            if not os.path.isfile(self.vpath+'/query%0.4d.sdf' % ichunk):
                continue

            #print 'PREDICTION IN MODEL %d' % ichunk
                
            call =['/usr/bin/python', BASEDIR+'predict.py','-e',itag, '-v', str(ver),
                   '-f', self.vpath+'/query%0.4d.sdf' % ichunk, '-s', str(ichunk)]

            retcode = subprocess.call(call)
            
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
        
        return (results)
