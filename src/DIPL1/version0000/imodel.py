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
import numpy as np
from pls import pls
from qualit import *
from utils import removefile

class imodel(model):
    def __init__ (self, vpath):
        model.__init__(self, vpath)
        
        ##
        ## General settings
        ##
        self.buildable = True
        self.quantitative = False
        self.confidential = False
        self.identity = True
        
        ##
        ## Normalization settings
        ##
        self.norm = True
        self.normStand = False
        self.normNeutr = True           # in this model must be True
        self.normNeutrMethod = 'moka'
        self.normNeutr_pH = 4.8         
        self.norm3D = True

        ##
        ## Molecular descriptor settings
        ##
        self.MD = 'pentacle'                         # 'padel'|'pentacle'
        self.padelMD = ['-3d']                       # '-2d'|'-3d'
        self.padelMaxRuntime = None
        self.padelDescriptor = None
        self.pentacleProbes = ['DRY','O','N1']       # 'DRY','O','N1','TIP'
        self.pentacleOthers = ['macc2_window 1.6','step 1.3']

        ##
        ## Modeling settings
        ##
        self.model = 'pls'
        self.modelLV = 4
        self.modelAutoscaling = False
        self.modelCutoff = 'auto'


    def getMatrices (self):  
        ncol = 0
        xx = []
        yy = []
        
        # obtain X and Y
        for i in self.tdata:
            if i[3]<1 :  # only for neutral or positive compounds
                continue
            if len(i[2])>ncol: ncol = len(i[2])
            xx.append(i[2])
            yy.append(i[4])  # notice there is one more column!!!!

        nrow = len (xx)
        
        Y = np.array (yy)
        X = np.empty ((nrow,ncol),dtype=np.float64)
      
        i=0
        for row in xx:
            if 'pentacle' in self.MD:  
                X[i,:]=self.adjustPentacle(row,len(self.pentacleProbes),ncol)
            else:
                X[i,:]=np.array(row)
            i+=1

        return X, Y

    def diagnosePLS_DA (self, model):

        if 'auto' == self.modelCutoff:
            model.calcOptCutoff ()
        else:
            model.calcConfussion(self.modelCutoff)

        for a in range (self.modelLV):
            TP = model.TP[a]
            TN = model.TN[a]
            FP = model.FP[a]
            FN = model.FN[a]
            
            # correct TP and FN using data not included in the PLS model
            for i in self.tdata:
                if i[3]<1 :# negative
                    if i[4] < 0.5 :
                        TP += 1
                    else:
                        FN += 1

            sens = sensitivity(TP,FN)
            spec = specificity(TN,FP)
            mcc  = MCC(TP,TN,FP,FN)

            print "LV:%-2d cutoff:%4.2f TP:%3d TN:%3d FP:%3d FN:%3d spec:%5.3f sens:%5.3f MCC:%5.3f" % (a+1,
                    model.cutoff[a], TP, TN, FP, FN, spec, sens, mcc)

        self.infoResult = []    
        self.infoResult.append( ('nobj',model.nobj) )
        self.infoResult.append( ('cutoff',str(self.modelCutoff) ) )
        self.infoResult.append( ('sens','%5.3f' % sens ) )
        self.infoResult.append( ('spec','%5.3f' % spec ) )
        self.infoResult.append( ('MCC','%5.3f' % mcc ) )
        
       
    def predict (self, molFile, molName, molCharge, detail, clean=True):
        
        if molCharge < 1:
            if clean: removefile(molFile)
            return ((True,'negative'), (True, 0), (True, 0.0))

        return model.predict(self, molFile, molName, molCharge, detail, clean)

