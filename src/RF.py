# -*- coding: utf-8 -*-

##    Description    RF model classifier and regressor
##                   
##    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
##
##    Copyright 2017 Manuel Pastor
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


import numpy as np
import sys
import os

from scipy import stats
from scipy.stats import t
import matplotlib.pyplot as plt
from collections import OrderedDict

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.externals import joblib

from scale import center, scale
from qualit import *
from utils import updateProgress

class RF:

    def __init__ (self):

        self.X = None
        self.Y = None
        
        self.nobj = 0
        self.nvarx = 0

        self.quantitative = False
        self.autoscale = False
        self.estimators = 0
        self.features = ''
        self.random = False
        
        self.mux = None
        self.muy = None
        self.wgx = None
        
        self.TP = 0
        self.TN = 0
        self.FP = 0
        self.FN = 0
        
        self.SDEP = 0.00    # SD error of the predictions
        self.Q2   = 0.00    # cross-validated R2

        self.clf = None

    def saveModel(self,filename):
        """Saves the model to a binary file in pkl format

        """

        f = file(filename,'wb')

        np.save(f,self.nobj)
        np.save(f,self.nvarx)

        np.save(f,self.quantitative)
        np.save(f,self.autoscale)
        np.save(f,self.estimators)
        np.save(f,self.features)
        np.save(f,self.random)
        
        np.save(f,self.mux)
        np.save(f,self.muy)
        np.save(f,self.wgx)
        
        np.save(f,self.TP)
        np.save(f,self.TN)
        np.save(f,self.FP)
        np.save(f,self.FN)
        
        np.save(f,self.SDEP)
        np.save(f,self.Q2)
        
        f.close()

        # the classifier is not saved correctly using cpickl
        joblib.dump(self.clf, os.path.dirname(filename)+'/clasifier.pkl')
        

            
    def loadModel(self,filename):
        """Loads the model from two files in pkl format
        """

        f = file(filename,'rb')
        
        self.nobj = np.load(f)
        self.nvarx = np.load(f)

        self.quantitative = np.load(f)
        self.autoscale = np.load(f)
        self.estimators = np.load(f)
        self.features = np.load(f)
        self.random = np.load(f)
        
        self.mux = np.load(f)
        self.muy = np.load(f)
        self.wgx = np.load(f)
        
        self.TP = np.load(f)
        self.TN = np.load(f)
        self.FP = np.load(f)
        self.FN = np.load(f)

        self.SDEP = np.load(f)
        self.Q2   = np.load(f)
        
        f.close()

        # the classifier is not saved correctly using cpickl
        self.clf = joblib.load(os.path.dirname(filename)+'/clasifier.pkl')

        
    def build (self, X, Y, quantitative=False, autoscale=False, nestimators=0, features='', random=False, tune=False):
        """Build a new RF model with the X and Y numpy matrices

        """
        
        nobj, nvarx= np.shape(X)

        self.nobj  = nobj
        self.nvarx = nvarx

        self.quantitative = quantitative
        self.autoscale = autoscale
        self.estimators = nestimators
        self.features = features
        self.random = random
        
        self.X = X.copy()
        self.Y = Y.copy()

        if autoscale:
            self.X, self.mux = center(self.X)
            self.Y, self.muy = center(self.Y)
            self.X, self.wgx = scale(self.X, autoscale)

        if random :
            RANDOM_STATE = 1226 # no reason to pick this number...
        else:
            RANDOM_STATE = None

        if tune :
            self.estimators, self.features = self.optimize ()
            
        if self.quantitative:
            print "Building Quantitative RF model"
            self.clf = RandomForestRegressor(n_estimators = int(nestimators),
                                        warm_start=False,
                                        max_features=features,
                                        oob_score=True,
                                        random_state=RANDOM_STATE)
        else:
            print "Building Qualitative RF_model"
            self.clf = RandomForestClassifier(n_estimators = int(nestimators),
                                         warm_start=False,
                                         max_features=features,
                                         oob_score=True,
                                         random_state=RANDOM_STATE)
        self.clf.fit(self.X, self.Y)
            
        # Regenerate the X and Y, since they might have been centered/scaled
        self.X = X.copy()
        self.Y = Y.copy()


    def validate (self):
        """ Validates 

        """

        if self.X == None or self.clf == None:
            return 
        
        X = self.X.copy()
        Y = self.Y.copy()

        if self.autoscale:
            X = X-self.mux
            Y = Y-self.muy
            X = X*self.wgx

        Yp = self.clf.predict(X)
           
        if self.quantitative:
            Ym   = np.mean(Y)
            SSY0 = np.sum (np.square(Ym-Y))
            SSY  = np.sum (np.square(Yp-Y))
            
            self.SDEP = np.sqrt(SSY/j)
            self.Q2   = 1.00 - (SSY/SSY0)
            self.OOBe = 1.00 - self.clf.oob_score_

            print "SDEP:", self.SDEP, "Q2:", self.Q2, "OOB error:", self.OOBe

        else:

            # I think this is not needed.... by the characteristics of RF it allways shows perfect performance
            if len(Yp) != len(Y):
                return

            for i in range(len(Y)):
                print Y[i], Yp[i]
            
            TP=TN=FP=FN=0
            
            for i in range(len(Y)):
                    
                    if Y[i] == 1.0:
                        if Yp[i] == 1.0:
                            TP+=1
                        else:
                            FN+=1
                    else:
                        if Yp[i] == 1.0:
                            FP+=1
                        else:
                            TN+=1
                            
            if TP+TN+FP+FN == 0:
                #print 'no objects'
                return          

            self.TP = TP
            self.TN = TN
            self.FP = FP
            self.FN = FN
            
            sens = sensitivity (TP, FN)
            spec = specificity (TN, FP)
            mcc  = MCC (TP, TN, FP, FN)
            
            self.OOBe = 1.00 - self.clf.oob_score_

            print "sens:", sens, "spec:", spec, "MCC:", mcc, "OOB error:", self.OOBe

        return (Yp)

    
    def project (self, Xb):
        """ Validates 

        """

        if self.clf == None:
            print 'failed to load clasifier'
            return
        
        if self.autoscale:
            Xb = Xb-self.mux
            Xb = Xb*self.wgx
        
            
        Yp = self.clf.predict(Xb)
        
        return (Yp)
    
                 
    def optimize (self):

        optEstimators = 100
        optFeatures = 'sqrt'

        return (optEstimators, optFeatures)

  
