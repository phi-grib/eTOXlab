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

from utils import removefile
##from rdkit.Chem import Descriptors
##from rdkit import Chem

from logp import computeLogP

class imodel(model):
    def __init__ (self, vpath):
        model.__init__(self, vpath)
        
        ##
        ## General settings
        ##
        self.buildable = False
        self.quantitative = False
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
        self.normNeutr_pH = 4.0
        self.norm3D = True

        ##
        ## Molecular descriptor settings
        ##
        self.MD = 'pentacle'                         # 'padel'|'pentacle'|'adriana'
        self.padelMD = ['-3d']                       # '-2d'|'-3d'
        self.padelMaxRuntime = None
        self.padelDescriptor = None
        self.pentacleProbes = ['DRY','O','N1','TIP']
        self.pentacleOthers = []

        ##
        ## Modeling settings
        ##
        self.model = 'pls'
        self.modelLV = 2
        self.modelAutoscaling = False
        self.modelCutoff = 'auto'
        self.selVar = False
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
        self.pentaclePath = '/opt/pentacle/pentacle106/'
        self.adrianaPath = '/opt/AdrianaCode/AdrianaCode226/'
        self.corinaPath = '/opt/corina/corina24/'
        self.javaPath = '/usr/java/jdk1.7.0_51/'
        self.RPath = '/opt/R/R-3.0.2/'
        self.standardiserPath = '/opt/standardise/standardise20140206/'


    def setSeries (self, molecules, numMol):

        self.infoSeries = []
        
    def log (self):
        
        self.infoModel = []
        self.infoModel.append( ('model','Decision tree') )
        
        result = model.log(self)
        return (result)

##  Example of how code can be moved to a new Python module (logp.py)
    
##    def computeLogP (self, mol):
##
##        lp = []
##        try:
##            suppl = Chem.SDMolSupplier(mol)
##            mi = suppl.next()
##
##            if mi is None:
##                return (False, 'wrong input format')
##
##            lp.append(Descriptors.MolLogP(mi))
##            
##        except:
##            return (False, 'wrong input format')
##
##        if len(lp) == 0:
##            return (False,'error in logP computation')
##        else:
##            return (True, lp)  


    def computePrediction (self, logP, charge):

        result = 'negative'

        if charge == 1:
            if logP[0] >= 1.61 :
                result = 'positive'
        elif charge == 2 :
            result = 'positive'
        elif charge > 2 :
            return (False, 'charge out of range')

        return (True, result)


    def predict (self, molFile, molName, molCharge, detail, clean=True):

        # default return values
        molPR=molCI=molAD=(False,0.0)

##        success, molMD = self.computeLogP (molFile)

        success, molMD = computeLogP (molFile)
        
        if not success: return (molPR,molAD,molCI)

        success, pr  = self.computePrediction (molMD,molCharge)
        molPR = (success, pr)
        if not success: return (molPR,molAD,molCI)

        if clean: removefile (molFile)
            
        return (molPR,molAD,molCI)
