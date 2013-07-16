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
from rdkit.Chem import Descriptors
from rdkit import Chem

class imodel(model):
    def __init__ (self, vpath):
        model.__init__(self, vpath)
        
        ##
        ## General settings
        ##
        self.buildable = False
        self.quantitative = False
        self.confidential = False
        
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
        self.MD = 'pentacle'                         # 'padel'|'pentacle'
        self.padelMD = ['-3d']                       # '-2d'|'-3d'
        self.pentacleProbes = ['DRY','O','N1','TIP']
        self.pentacleOthers = []

        ##
        ## Modeling settings
        ##
        self.model = 'pls'
        self.modelLV = 2
        self.modelAutoscaling = False


    def setSeries (self, molecules, numMol):

        self.infoSeries = []
        
    def log (self):
        
        self.infoModel = []
        self.infoModel.append( ('model','Decision tree') )
        
        result = model.log(self)
        return (result)

    def computeLogP (self, mol):

        lp = []
        try:
            suppl = Chem.SDMolSupplier(mol)
            mi = suppl.next()

            if mi is None:
                return (False, 'wrong input format')

            lp.append(Descriptors.MolLogP(mi))
            
        except:
            return (False, 'wrong input format')

        if len(lp) == 0:
            return (False,'error in logP computation')
        else:
            return (True, lp)  


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

                   
    def predict (self, molN, detail, clean=True):
        
        # alias
        mol, charge = molN[0], molN[1]

        # default return values
        pr=ri=ad=(False,0.0)

        md = self.computeLogP (mol)
        if not md[0]: return (pr,ad,ri)

        pr = self.computePrediction (md[1],charge)
        if not pr[0]: return (pr,ad,ri)

        if clean:
            removefile (mol)
            
        return (pr,ad,ri)
