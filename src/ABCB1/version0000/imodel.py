# -*- coding: utf-8 -*-

##    Description    eTAM model template
##                   
##    Authors:       Manuel Pastor (manuel.pastor@upf.edu)
##                   Pau Carrio (pcarrio@imim.es)
##
##    Copyright 2013 Manuel Pastor
##
##    This file is part of eTAM.
##
##    eTAM is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation version 3.
##
##    eTAM is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with eTAM.  If not, see <http://www.gnu.org/licenses/>.

from model import model
import subprocess 
import numpy as np
from math import sqrt
from utils import randomName
import os
from qualit import *

class imodel(model):
    def __init__ (self, vpath):
        model.__init__(self, vpath)

        ##
        ## General settings
        ##
        self.buildable = True
        self.quantitative = False

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
        self.MD = 'padel'                         # 'padel'|'pentacle'
        self.padelMD = ['-3d']                       # '-2d'|'-3d'
        self.padelMaxRuntime = 12000
        self.padelDescriptor = '/opt/padel/descriptors_etam.xml'        
        self.pentacleProbes = ['DRY','O','N1']       # 'DRY','O','N1','TIP'
        self.pentacleOthers = ['macc2_window 1.6','step 1.3']

        ##
        ## Modeling settings
        ##
        self.model = 'myRF'
        self.modelLV = 3
        self.modelAutoscaling = True

            
    def build (self, data):
        print 'build RF'

        self.saveTraining (data)
 
        X,Y = self.getMatrices (data)

        ## new code 
        nrow, ncol = np.shape (X)

        np.savetxt( self.vpath+"/nvars.csv" , np.array([ncol]) , fmt='%5.5f' , delimiter=',')
        np.savetxt( self.vpath+"/tmpX.csv" , X , fmt='%5.5f' , delimiter=',')
        np.savetxt( self.vpath+"/tmpY.csv" , Y , fmt='%5.5f' , delimiter=',')
        # call R to build the model
        try:
            call = ['/usr/local/bin/Rscript',self.vpath+'/buildMyModel.R',self.vpath]
            retcode = subprocess.call(call)
        except:
            raise Exception('R model building produced errors')
        #os.remove(self.vpath+"/X.csv")
        #os.remove(self.vpath+"/Y.csv")
        
        # retrieve information from the model ( produced by buildMyModel.R )
        infoModelR = np.loadtxt( self.vpath+'/infoModelR.csv',delimiter=",",skiprows=1)
        TP = float(infoModelR[0])
        TN = float(infoModelR[1])
        FP = float(infoModelR[2])
        FN = float(infoModelR[3])

        sens = sensitivity (TP, FN)
        spec = specificity (TN, FP)
        mcc  = MCC (TP, TN, FP, FN)

        print  TP, TN, FP, FN, spec, sens, mcc

        self.infoModel = []
        self.infoModel.append( ('model','myRmodel') )

        self.infoResult = []    
        self.infoResult.append( ('nobj',nrow) )
        self.infoResult.append( ('sens','%6.3f' % sens ) )
        self.infoResult.append( ('spec','%6.3f' % spec ) )
        self.infoResult.append( ('MCC','%6.3f' % mcc ) )

        # so far we are not computing ADRI
        #success, result = self.ADRI (X,Y)

        return (True, 'Model OK')


    def computePR (self, md, charge):
        # adjust md
        nvars = np.loadtxt( self.vpath+"/nvars.csv" ,   delimiter=',' )
        mda = self.adjustPentacle(md,len(self.pentacleProbes),nvars)

        tmpfile = "./RF-"+randomName(20)+".csv"
        np.savetxt( tmpfile , mda , fmt='%5.5f' , delimiter=',')
        
        # call R to predict
        try:
            call = ['/usr/local/bin/Rscript',self.vpath+'/predictMyModel.R', self.vpath, tmpfile]
            retcode = subprocess.call(call)
        except:
            return (False,'R model predict produce errors')

        # prediction was stored in tmpfile 
        npr = np.loadtxt( tmpfile )

        os.remove(tmpfile)

        if npr < 0.5:
            return (True, 'negative')
        else:
            return (True, 'positive')
        

    def computeAD (self, md, pr, detail):

        return (True, 'not implemented for this model')

    def computeRI (self, ad):

        return (True, 'not implemented for this model')    

