# -*- coding: utf-8 -*-
#
#    Description    eTAM model class
#                
#
#    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
#
#    (c) PhI 2013

import os
import sys
import shutil
import subprocess

import openbabel as ob
##import rpy2.robjects as ro
import pybel
import numpy as np

from pls import pls
from StringIO import StringIO
from utils import removefile
from utils import opt

class model:
            
    def __init__ (self, vpath):

        self.vpath = vpath
        self.numLV = 5
        self.pH = 7.4

##        ro.r("""
##        computePR <- function( test ){
##        suppressPackageStartupMessages(library(pls))
##        
##        #load model
##        load(\""""+self.vpath+"""/model.Rdata")
##        
##        #adjust model and test dimensions
##        testMD <- mdl$pentacleAdjustOutputs( mdl$ncolTrain, test,
##        nrProbes = mdl$nrProbes)
##        
##        #perform prediction and return value for best LV
##        predict( mdl , newdata = as.data.frame(matrix(testMD,nrow=1)) , type="response")[,,mdl$bestncomp]
##        }""")
##        self.Rmodel = ro.globalenv['computePR']
##
##        ro.r("""
##        buildModel <- function( xx, yy, Rdatafile ,ncomp , nrProbes = 4 ){
##        suppressPackageStartupMessages(library(pls))
##        Y <- unlist(yy)
##        X <- as.data.frame(matrix(unlist(xx),nrow=nrow(xx),byrow=TRUE))
##        colnames(X) <- paste("V",1:(ncol(X)),sep="")
##        # build pls model
##        mdl <- plsr( Y ~ . , data=X , ncomp= ncomp , method="oscores")
##        mdl$ncolTrain = ncol(X)
##        mdl$bestncomp = ncomp
##        mdl$nrProbes = nrProbes
##        pentacleAdjustOutputs <- function(
##        ncolTrain ,i2 , nrProbes = mdl$nrProbes
##         ){
##            nrProbeBlocs = switch( nrProbes, 1 , 2 , 6 , 10 )
##            i2    <- unlist(i2)
##            delta <- ( ncolTrain - length(i2) )
##            if(delta == 0){
##            # nothing to do
##            }
##            if(delta < 0){
##                # i2 is bigger than already exisiting data
##                # remove columns
##                delta2 <- abs( delta / nrProbeBlocs )
##                icol <- ncolTrain / nrProbeBlocs
##                indx <- rep(c(rep(TRUE,icol),rep(FALSE, delta2)),nrProbeBlocs)
##                i2 <- matrix( i2[indx] ,nrow=1)
##            }
##            if(delta > 0 ){
##                # i2 is smaller than already exisiting data
##                # add columns  
##                delta2 <- delta / nrProbeBlocs
##                icol <- length(i2)/nrProbeBlocs
##                ni2 <- c(  i2[1:icol] , rep(0, delta2 ))
##                tmp <- i2[-c(1:icol)] 
##                while(length(tmp)>0){
##                    ni2 <- c( ni2 ,c(  tmp[1:icol],rep(0,delta2)))
##                    tmp <- tmp[-c(1:icol)]
##                }
##                i2 <- ni2
##            }
##            return(i2)
##        }
##        mdl$pentacleAdjustOutputs <- pentacleAdjustOutputs
##        save(mdl, file=Rdatafile )
##        return(T)
##        }
##        """ )
##        self.Rlearner = ro.globalenv['buildModel']

        self.trainList = []
        try:
            ifile = open (self.vpath+'/itrain.txt')
            for line in ifile:
                piece = line.split('\t')
                self.trainList.append(piece)
            ifile.close()
        except:
            pass
    
        
    def standardize (self, moli):
        """Applies a structure normalization protocol

           DUMMY method. At present it does nothing

           The name of the output molecules is built as a+'original name'

           Returns a tuple containing:
           1) True/False: depending on the result of the method
           2) (if True ) The name of the output molecule
              (if False) The error message
        """
        
        molo = 'a'+moli

        shutil.copy(moli,molo)
        
        return (True,molo)


    def protonate (self, moli, pH):
        """Adjusts the ionization state of the molecule "moli" 

           In this implementation, it ises blabber_sd from Molecular Discovery
           The result is a tuple containing:
           1) True/False: describes the success of the protonation for this compound
           2) (if True ) The name of the protonated molecules and its formal charge
              (if False) The error message
        """

        molo = 'b'+moli

        stderrf = open (os.devnull, 'w')
        stdoutf = open (os.devnull, 'w')     

        call = [opt+'blab_etox/blabber_sd', moli,
                '-p',  str(pH),
                '-o',  molo]

        try:
            retcode = subprocess.call(call,stdout=stdoutf, stderr=stderrf)
        except:
            return (False, 'Blabber execution error', 0.0)
        
        stdoutf.close()
        stderrf.close()
        
        if retcode != 0:
            return (False, 'Blabber execution error', 0.0)
        
        try:
            finp = open (molo)
        except:
            return (False, 'Blabber output not found', 0.0)

        charge = 0
        for line in finp:
            if line.startswith ('M  CHG'):
                items = line.split()
                if int(items[2]):
                    for c in range (4,len(items),2): charge+=int(items[c])
                break   
        finp.close()
        return (True, molo, charge)


    def convert3D (self, moli):
        """Converts the 2D structure of the molecule "moli" to 3D

           In this implementation, it uses CORINA from Molecular Networks
           The result is a tuple containing:
           1) suucTrue/False: describes the success of the 3D conversion for this compound
           2) (if True ) The name of the 3D molecule
              (if False) The error mesage
        """

        molo = 'c'+moli
        
        stderrf = open (os.devnull, 'w')
        stdoutf = open (os.devnull, 'w')
        
        call = [opt+'corina/corina',
                '-dwh','-dori',
                '-ttracefile=corina.trc',
                '-it=sdf', moli,
                '-ot=sdf', molo]

        try:
            retcode = subprocess.call(call, stdout=stdoutf, stderr=stderrf)
        except:
            return (False, 'Corina execution error')

        stdoutf.close()
        stderrf.close()
            
        if retcode != 0:
            return (False, 'Corina execution error')

        if not os.path.exists(molo):
            return (False, 'Corina output not found')

        return (True,molo)


    def checkIdentity (self, mol):
        """ Checks if the compound "mol" is part of the training set 

            We used InChiKeys without the last three chars (ionization state) to make the comparison

            This version uses OpenBabel (OB)

            TODO: Migrate to RCDKit
        """
##      DISABLED TEMPORARILY FOR DEVELOPMENT, UNCOMMENT IN THE PRODUCTION VERSION
##        conv = ob.OBConversion()
##        conv.SetInAndOutFormats("sdf", "inchi")
##        conv.SetOptions("Kw", conv.OUTOPTIONS)
##        moli = ob.OBMol()
##        conv.ReadFile(moli, mol)
##        ik = conv.WriteString(moli)
##
##        # alternative version using pybel (shorter but I was unable to supress warnings)
##        #   moli = pybel.readfile ('sdf',mol).next()
##        #   ik = moli.write('INCHIKEY',')
##         
##        ik = ik[:-3] # remove the right-most part expressing ionization
##        
##        for l in self.trainList:
##            if ik in l[0]:
##                #print 'the query compound is in the training set'
##                return (True, float(l[1]))
          
        return (False, mol)


    def computeMD (self, mol):
        """ Computes the Molecular Descriptors for compound "mol"

            In this implementation we run Pentacle with default settings
            It returms a tuple that contains
            1) True or False, indicating the success of the computation
            2) A vector of floats (if True) with the GRIND descriptors
               An Error message (if False)
        """
        
        molr = mol[:-4] # mol root; name without extension. I can safely asume the extension is ".sdf"

        removefile ( '/var/tmp/'+molr )
        removefile ( '/var/tmp/'+molr+'.ppf' )
        
        t = open ('template-md','w')
        t.write ('name '+molr+'\n')
        t.write ('input_file '+mol+' sdf\n')
        t.write ('mif_computation grid\n')
        t.write ('mif_discretization amanda\n')
        t.write ('mif_encoding  macc2\n')
        t.write ('probe DRY\nprobe O\nprobe N1\nprobe TIP\n')
        t.write ('dynamic yes\n')
        t.write ('export_data csv\n')
        t.close()
        
        call = [opt+'pentacle_etox/pentacle',
                '-c','template-md']  

        stdoutf = open ('stdout.txt','w')
        stderrf = open (os.devnull, 'w')

        try:
            retcode = subprocess.call(call,stdout=stdoutf,stderr=stderrf)
        except:
            return (False, 'Pentacle execution error' )

        stdoutf.close()
        stderrf.close()

        removefile ( '/var/tmp/'+molr )
        removefile ( '/var/tmp/'+molr+'.ppf' )

        try:
            stdoutf = open('stdout.txt')
        except:
            return (False, 'Pentacle std output not found')

        for line in stdoutf:
            if 'Error' in line :
                return (False, line)
        stdoutf.close()
            
        try:
            fpr = open (molr+'.csv')
        except:
            return (False, 'Pentacle results not found')
        
        line = fpr.readline()
        fpr.close()
        
        #temp = line.split(',')    
        #md = [float(x) for x in temp[1:]]

        # remove mol name using line.partition(',') and extracting the right piece [2]
        md = np.loadtxt(StringIO(line.partition(',')[2]),delimiter=',')
        
        return (True,md)


    def computePR (self, md, charge):
        """ Computes the prediction for compound "mol"

            This function makes use of the molecular descriptor vector (md) to project the compound in the model
            The model has been loaded previously as an R object
        """

##        try:
##            pr = self.Rmodel(ro.FloatVector(md))
##        except:
##            return(False, "Prediction not calculated")

        model = pls()
        model.loadModel(self.vpath+'/modelPLS.npz')
        success, result = model.project(self.adjustPentacle(md,4,model.nvarx),self.numLV)

        if success:
            yp = result[0]
            return (True, yp[0])
        else:
            return (False, result)


    def computeAD (self, md, pr, detail):
        """Carries out a protocol for determining how fat is the query compound from the model

           Provisionally, implements a temporary version of the ADAN method

           Returns a tuple that contains
           1) True or False, indicating the success of the computation
           2) (if True ) a number between 0 and 5 with the number of criteria broken
              (if False) an error message 
           

        """
        f = file (self.vpath+'/tscores.npz','rb')
        nlv = np.load(f)
        p95dcentx = np.load(f)
        p95dclosx = np.load(f)
        p95dmodx = np.load(f)
        p95dcenty = np.load(f)
        p95dclosy = np.load(f)
        centx = np.load(f)
        centy = np.load(f)
        T = np.load(f)
        Y = np.load(f)
        f.close()

        model = pls ()
        model.loadModel(self.vpath+'/modelPLS.npz')
        success, result = model.project(self.adjustPentacle(md,4,model.nvarx),nlv)

        y = pr
        t = result[1]
        d = result[2]

        AD=dict()
        # compute distance t to centroid x
        dcentx = np.sqrt(np.sum(np.square(centx-t)))
        AD['dcentx']=dcentx>p95dcentx

        # compute min distance to T    
        dclosx=1e20
        for i in range (model.nobj):
            dtemp = np.sqrt(np.sum(np.square(T[i,:]-t)))
            if dtemp < dclosx : dclosx = dtemp
        AD['dclosx']=dclosx>p95dclosx

        # compute distance to modx
        AD['dmodx']=d[-1]>p95dmodx

        # compute distance to centroid y
        dcenty = np.abs(y-centy)
        AD['dcenty']=dcenty>p95dcenty

        # compute min distance to Y
        dclosy=1e20
        for i in range (model.nobj):
            dtemp = np.abs(Y[i]-y)
            if dtemp < dclosy : dclosy = dtemp
        AD['dclosy']=dclosy>p95dclosy

        print "DCENTX %6.3f (%6.3f)\n" % (dcentx,p95dcentx),
        print "DCLOSX %6.3f (%6.3f)\n" % (dclosx,p95dclosx),
        print "DCMODX %6.3f (%6.3f)\n" % (d[-1],p95dmodx),
        print "DCENTY %6.3f (%6.3f)\n" % (dcenty,p95dcenty),
        print "DCLOSY %6.3f (%6.3f)\n" % (dclosy,p95dclosy)

        return (True,sum(AD.values()))


    def computeRI (self, ad):
        """Calculates a Reliability Index for the given prediction

           Provisionally it returns a tuple that contains
           1) True or False, indicating the success of the computation
           2) (if True ) - the model SDEP when the #broken AD criteria is 0 or 1
                         - twice the model SDEP when the #broken AD criteria is 2 or 3
                         - 0.0 otherwyse (this must be interpreted as a "model NA")
              (if False) An error message
        """

        model = pls ()
        model.loadModel(self.vpath+'/modelPLS.npz')

        ri=0.0
        if ad<4:
            model = pls ()
            model.loadModel(self.vpath+'/modelPLS.npz')
            sdep = model.SDEP[model.Av-1]
            if ad<2:
                ri = sdep    # 0 or 1 criteria broken
            else:
                ri = sdep*2  # 2 or 3 criteria broken
        
        return (True,ri)


    def normalize (self, mol):
        """Preprocesses the molecule "mol" by running a workflow that:

        - Normalizes the 2D structure (DUMMY)
        - Adjusts the ionization state 
        - Converts the structure to 3D

        The result is a tuple containing:
        1) True/False: describes the success of the normalization
        2) (if True ) The name of the normalized molecule and its formal charge
           (if False) The error mesage
        """

        success, resulta = self.standardize (mol)
        if not success: return (False, resulta)
        
        success, resultb, charge = self.protonate (resulta, self.pH)
        if not success: return (False, resultb)

        success, resultc = self.convert3D (resultb)
        if not success: return (False, resultc)

        return (True,resultc,charge)
    

    def predict (self, molN,detail):
        """Runs the prediction protocol

        """

        # alias
        mol, charge = molN[0], molN[1]

        # default return values
        pr=ri=ad=(False,0.0)

        ic = self.checkIdentity (mol)
        if ic[0]: return (ic, (True, 100.0), (True,100.0))
     
        md = self.computeMD (mol)
        if not md[0]: return (pr,ad,ri)
        
        pr = self.computePR (md[1],charge)
        if not pr[0]: return (pr,ad,ri)

        ad = self.computeAD (md[1], pr[1], detail)
        if not ad[0]: return (pr,ad,ri)

        ri = self.computeRI (ad[1])
        return (pr,ad,ri)


    def getInChi (self, mol):
        """ Computes the InChiKey for the compound "mol"

            We used InChiKeys without the last three chars (ionization state) to make the comparison
        """
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("sdf", "inchi")
        conv.SetOptions("Kw", conv.OUTOPTIONS)
        moli = ob.OBMol()
        conv.ReadFile(moli, mol)
        ik = conv.WriteString(moli)

        return (True,ik[:-3])


    def getBio (self, mol):
        """ Extracts the value of the experimental biological property from the compound "mol" 

            Such value must be identified by the tag <activity>
        """

        moli = pybel.readfile ('sdf',mol).next()
        data = moli.data
        if data.has_key('activity'):
            return (True, float(data['activity']))
        else:
            return (False, 'not found')


    def adjustPentacle (self, row, nprobes, Bcol):
        """Adjust the row of GRIND descriptors in "row" to Bcol size, asuming that we have nprobes
           GRID probes, applying a procrustean transform for each block (correlogram) 

           Both the row and the returning array are NumPy float64 arrays 
        """

        Acol = len(row)        # orignal num of columns
        blocks = [0,1,3,6,10]
        nblocks = blocks[nprobes]
        
        if Acol == Bcol:
            return row

        deltaA = Acol/nblocks  # num col in original
        deltaB = Bcol/nblocks  # num col in new

        nn = np.empty(0,dtype='float64')   
        if Acol > Bcol:
            for i in range (nblocks):
                start = i*deltaA
                nn=np.hstack((nn,row[start:(start+deltaB)]))
        else:
            zero = np.zeros(deltaB-deltaA,dtype='float64')
            for i in range (nblocks):
                start = i*deltaA
                nn=np.hstack((nn,row[start:(start+deltaA)]))
                nn=np.hstack((nn,zero))
                
        return nn


    def build (self, data):
        """Uses the data extracted from the training series to build a model, using the Rlearner object 

           This function also creates the "itrain.txt" file that describes the training series, including InChiKey of the compounds
        """

        ncol = 0
        xx = []
        yy = []

        # obtain X and Y
        for success, i in data:
            if len(i[1])>ncol: ncol = len(i[1])
            xx.append(i[1])
            yy.append(i[2])       
        nrow = len (xx)

        # new
        Y = np.array (yy)
        X = np.empty ((nrow,ncol),dtype=np.float64)
        i=0
        for row in xx:
            X[i,:]=self.adjustPentacle(row,4,ncol)
            i+=1
            
        nrows, ncols = np.shape(X)
        print nrows, ncols

        model = pls ()
        model.build (X,Y,self.numLV)
        model.validateLOO(self.numLV)
        print 'R2: %6.4f Q2: %6.4f SDEP: %6.4f' % \
              (model.SSYac[3],model.Q2[3],model.SDEP[3])
        model.saveModel (self.vpath+'/modelPLS.npz')
        
        # translate X to numpy format
        
##        xflat = np.empty(0,dtype='float64')
##        for row in xx:
##            xflat = np.hstack((xflat,self.adjustPentacle(row,4,ncol)))       
##        Y = ro.FloatVector(yy)
##        X = ro.r.matrix(ro.FloatVector(xflat),nrow)
##        try:
##            pr = self.Rlearner(X, Y, self.vpath+'/model.Rdata', self.numLV, 4)
##        except:
##            return(False, "R error building the model")

            
        # compute PP on X
        model = pls ()
        model.build (X,Y,targetSSX=0.4)

        nlv = model.Am

        # initialize arrays
        T = np.empty((nrows,nlv),dtype=np.float64)
        centx = np.empty(nlv,dtype=np.float64)
        dcentx = np.empty(nrows, dtype=np.float64)
        dclosx = np.empty(nrows, dtype=np.float64)
        dcenty = np.empty(nrows, dtype=np.float64)
        dclosy = np.empty(nrows, dtype=np.float64)
        
        # extract t 
        for a in range (nlv):
            ta = model.t[a]
            T[:,a]=ta
            centx[a]=np.mean(ta)

        i95 = np.ceil(nrows*0.95)-1
        # compute distances to X centroid and percentil 95
        for i in range(nrows):
            dcentx[i] = np.sqrt(np.sum(np.square(centx-T[i,:])))

        dcentx = np.sort(dcentx)
        p95dcentx = dcentx[i95]

        # compute closer distances in X and percentil 95
        for i in range (nrows):
            dclosx[i]=1e20
            for j in range (nrows):
                if j != i:
                    dtemp = np.sqrt(np.sum(np.square(T[j,:]-T[i,:])))
                    if dtemp < dclosx[i] : dclosx[i] = dtemp
                    
        dclosx = np.sort(dclosx)
        p95dclosx = dclosx[i95]

        # compute percentil 96 dmodx
        dmodx = np.array(model.dmodx[-1])
        p95dmodx=dmodx[i95]
        
        # compute distance to Y centroid and percentil 95
        centy = np.mean(Y)
        dcenty = np.abs(Y-centy)
        dcenty = np.sort(dcenty)
        p95dcenty = dcenty[i95]

        # compute mutual distance in Y and percentil 95
        for i in range (nrows):
            dclosy[i]=1e20
            for j in range (nrows):
                if j != i:
                    dtemp = np.abs(Y[i]-Y[j])
                    if dtemp < dclosy[i] : dclosy[i] = dtemp

        dclosy = np.sort(dclosy)
        p95dclosy = dclosy[i95]

        # write in a file, Am -> critical distances -> centroid -> t
        f = file (self.vpath+'/tscores.npz','wb')
        np.save(f,nlv)
        np.save(f,p95dcentx)
        np.save(f,p95dclosx)
        np.save(f,p95dmodx)
        np.save(f,p95dcenty)
        np.save(f,p95dclosy)
        np.save(f,centx)
        np.save(f,centy)
        np.save(f,T)
        np.save(f,Y)
        f.close()

        # write itraining data, we must include PP computed on X
        ftrain = open (self.vpath+'/itrain.txt','w')
        for success, i in data:
            ftrain.write (i[0]+'\t'+str(i[2])+'\n') # now it writes i0 (InChi) + i2 (activity)
        ftrain.close()


        return(True, "Model built")


    def extract (self, mol):
        """Process the compound "mol" for obtaining
           1) InChiKey (string)
           2) Molecular Descriptors (NumPy float64 array)
           3) Biological Activity (float)

           Returns a Tuple as (True,('InChi',array[0.0,0.0, 0.0],4.56))
        """
        
        moli = mol[0]
        i1=''
        i2=[]
        i3=0.0
        
        success, i1 = self.getInChi(moli)
        if not success:
            print 'error in InChi'
            return (False,(i1,i2,i3))
        
        success, i2 = self.computeMD (moli)
        if not success:
            print 'error in MD '
            return (False,(i1,i2,i3))

        success, i3 = self.getBio (moli)
        if not success:
            print 'error in Bio'
            return (False,(i1,i2,i3))
        
        return (True,(i1,i2,i3))
