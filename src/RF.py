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
        self.wgx = None
        
        self.TP = 0
        self.TN = 0
        self.FP = 0
        self.FN = 0
        
        self.SDEC = 0.00    # SD error of the calculations
        self.R2   = 0.00    # determination coefficient
        
        self.OOBe = 0.00

        self.clf = None

    def saveModel(self,filename):
        """Saves the model to a binary file in numpy file and another in pkl format

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
        np.save(f,self.wgx)
        
        np.save(f,self.TP)
        np.save(f,self.TN)
        np.save(f,self.FP)
        np.save(f,self.FN)
        
        np.save(f,self.SDEC)
        np.save(f,self.R2)

        np.save(f,self.OOBe)
        
        f.close()

        # the classifier cannot be saved with numpy
        joblib.dump(self.clf, os.path.dirname(filename)+'/clasifier.pkl')
        
            
    def loadModel(self,filename):
        """Loads the model from two files, one in numpy and another in pkl format
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
        self.wgx = np.load(f)
        
        self.TP = np.load(f)
        self.TN = np.load(f)
        self.FP = np.load(f)
        self.FN = np.load(f)

        self.SDEC = np.load(f)
        self.R2   = np.load(f)

        self.OOBe = np.load(f)
        
        f.close()

        # the classifier cannot be loaded with numpy
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
            self.X, self.wgx = scale(self.X, autoscale)

        if random :
            RANDOM_STATE = None
        else:
            RANDOM_STATE = 1226 # no reason to pick this number

        if tune :
            self.estimators, self.features = self.optimize (X,Y)

            if self.features=='none':
                self.features = None
            
        if self.quantitative:
            print "Building Quantitative RF model"
            self.clf = RandomForestRegressor(n_estimators = int(self.estimators),
                                        warm_start=False,
                                        max_features=self.features,
                                        oob_score=True,
                                        random_state=RANDOM_STATE)
        else:
            print "Building Qualitative RF_model"
            self.clf = RandomForestClassifier(n_estimators = int(self.estimators),
                                         warm_start=False,
                                         max_features=self.features,
                                         oob_score=True,
                                         random_state=RANDOM_STATE)
        self.clf.fit(self.X, self.Y)
            
        # Regenerate the X and Y, since they might have been centered/scaled
        self.X = X.copy()
        self.Y = Y.copy()


    def validate (self):
        """ Validates the models and completes suitable scoring values

        """

        if self.X == None or self.clf == None:
            return 
        
        X = self.X.copy()
        Y = self.Y.copy()

        if self.autoscale:
            X = X-self.mux
            X = X*self.wgx
        
        Yp = self.clf.predict(X)
           
        if self.quantitative:
            Ym   = np.mean(Y)
            SSY0 = np.sum (np.square(Ym-Y))
            SSY  = np.sum (np.square(Yp-Y))
            
            self.SDEC = np.sqrt(SSY/self.nobj)
            self.R2   = 1.00 - (SSY/SSY0)
            self.OOBe = 1.00 - self.clf.oob_score_

            print "SDEC:", self.SDEC, "R2:", self.R2, "OOB error:", self.OOBe

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
        """ Uses the X matrix provided as argument to predict Y

        """

        if self.clf == None:
            print 'failed to load clasifier'
            return
        
        if self.autoscale:
            Xb = Xb-self.mux
            Xb = Xb*self.wgx
        
            
        Yp = self.clf.predict(Xb)
        
        return (Yp)
    
                 
    def optimize (self, X, Y ):
        """ Optimizes the number of trees (estimators) and max features used (features)
            and returns the best values, acording to the OOB criteria

            The results are shown in a diagnostic plot

            To avoid including many trees to produce tiny improvements, increments of OOB error
            below 0.01 are considered irrelevant
        """
                    
        RANDOM_STATE = 1226
        errors = {}
        features = ['sqrt','log2','none']

        if self.quantitative:
            tclf = {'sqrt': RandomForestRegressor(warm_start=False, oob_score=True, max_features="sqrt",random_state=RANDOM_STATE),
                    'log2': RandomForestRegressor(warm_start=False, oob_score=True, max_features="log2",random_state=RANDOM_STATE),
                    'none': RandomForestRegressor(warm_start=False, oob_score=True, max_features=None  ,random_state=RANDOM_STATE) }
        else:
            tclf = {'sqrt': RandomForestClassifier(warm_start=False, oob_score=True, max_features="sqrt",random_state=RANDOM_STATE),
                    'log2': RandomForestClassifier(warm_start=False, oob_score=True, max_features="log2",random_state=RANDOM_STATE),
                    'none': RandomForestClassifier(warm_start=False, oob_score=True, max_features=None  ,random_state=RANDOM_STATE) }

        # Range of `n_estimators` values to explore.
        min_estimators = 15
        max_estimators = 700
        stp_estimators = 100

        num_steps = int((max_estimators-min_estimators)/stp_estimators)

        print 'optimizing RF....'
        updateProgress (0.0)
        
        optValue = 1.0e10
        j = 0
        for fi in features:
            errors[fi] = []
            count = 0
            for i in range(min_estimators, max_estimators + 1,stp_estimators):
                clf = tclf[fi]
                clf.set_params(n_estimators=i)
                clf.fit(X,Y)
                oob_error = 1 - clf.oob_score_
                errors[fi].append((i,oob_error))
                if oob_error < optValue:
                    if np.abs(oob_error - optValue) > 0.01:
                        optValue = oob_error
                        optEstimators = i
                        optFeatures = fi

                updateProgress (float(count+(j*num_steps))/float(len(features)*num_steps))
                count = count+1
            j=j+1
            
        for ie in errors:
            xs, ys = zip (*errors[ie])
            plt.plot(xs, ys, label=ie)     

        plt.xlim(min_estimators, max_estimators)
        plt.xlabel("n_estimators (Trees)")
        plt.ylabel("OOB error rate")
        plt.legend(loc="upper right")
        plt.show()

        #plt.savefig(self.vpath+"/rf-OOB-parameter-tuning.png")
        plt.savefig("./rf-OOB-parameter-tuning.png")

        print 'optimum features:', optFeatures, 'optimum estimators:', optEstimators, 'best OOB:', optValue
        
        return (optEstimators, optFeatures)

  
