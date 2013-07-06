# -*- coding: utf-8 -*-

##    Description    PLS toolkit using NIPALS algorithm 
##                   
##    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
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


import numpy as np
import sys
from scale import center, scale
from qualit import *

class pls:

    def __init__ (self):
        self.X = None
        self.Y = None
        
        self.Am = 0  # model dimensionality
        self.Av = 0  # number of LV validated. Notice that Av <= Am
        self.nobj = 0
        self.nvarx = 0
        
        self.mux = None
        self.muy = None

        self.wgx = None

        self.autoscale = False
        
        self.t = []  # scores
        self.p = []  # loadings
        self.w = []  # weights
        self.c = []  # inner relationship
        self.SSXex = []   # SSX explained
        self.SSXac = []   # SSX accumulated
        self.SSYex = []   # SSY explained
        self.SSYac = []   # SSY accumulated
        self.SDEC = []    # SD error of the calculations
        self.dmodx = []   # distance to model
        
        self.SSY  = []    # SSY explained
        self.SDEP = []    # SD error of the predictions
        self.Q2   = []    # cross-validated R2

        self.cutoff = None
        
        self.TP = 0
        self.TN = 0
        self.FP = 0
        self.FN = 0

    def saveModel(self,filename):
        """Saves the whole model to a binary file in numpy .npy format

        """

        f = file(filename,'wb')

        np.save(f,self.Am)
        np.save(f,self.Av)
        np.save(f,self.nobj)
        np.save(f,self.nvarx)
        
        np.save(f,self.mux)
        np.save(f,self.muy)
        np.save(f,self.wgx)

        np.save(f,self.autoscale)
        
        for a in range(self.Am):
            np.save(f,self.t[a])
            np.save(f,self.p[a])
            np.save(f,self.w[a])
            np.save(f,self.c[a])
            np.save(f,self.SSXex[a])
            np.save(f,self.SSXac[a])
            np.save(f,self.SSYex[a])
            np.save(f,self.SSYac[a])
            np.save(f,self.SDEC[a])
            np.save(f,self.dmodx[a])

        for a in range(self.Av):
            np.save(f,self.SSY[a])
            np.save(f,self.SDEP[a])
            np.save(f,self.Q2[a])

        np.save(f,self.cutoff)

        np.save(f,self.TP)
        np.save(f,self.TN)
        np.save(f,self.FP)
        np.save(f,self.FN)
        
        f.close()

            
    def loadModel(self,filename):
        """Loads the whole model from a binary file in numpy .npy format

        """

        f = file(filename,'rb')
        
        self.Am = np.load(f)
        self.Av = np.load(f)
        self.nobj = np.load(f)
        self.nvarx = np.load(f)
        
        self.mux = np.load(f)
        self.muy = np.load(f)
        self.wgx = np.load(f)

        self.autoscale = np.load(f)
        
        for a in range(self.Am):
            self.t.append (np.load(f))
            self.p.append (np.load(f))
            self.w.append (np.load(f))
            self.c.append (np.load(f))
            self.SSXex.append (np.load(f))
            self.SSXac.append (np.load(f))
            self.SSYex.append (np.load(f))
            self.SSYac.append (np.load(f))
            self.SDEC.append (np.load(f))
            self.dmodx.append (np.load(f))

        for a in range(self.Av): 
            self.SSY.append (np.load(f))
            self.SDEP.append (np.load(f))
            self.Q2.append (np.load(f))

        self.cutoff = np.load(f)

        self.TP = np.load(f)
        self.TN = np.load(f)
        self.FP = np.load(f)
        self.FN = np.load(f)

        f.close()                     

        
    def build (self, X, Y, targetA=0, targetSSX=0.0, autoscale=False):
        """Build a new PLS model with the X and Y numpy matrice provided using NIPALS algorithm

           The dimensionality of the model can be defined either providing
           1. directly the number of LV to extract (targetA)
           2. the fraction of SSX that the model will explain (targetSSX)

           The X and Y matrices are centered but no other scaling transform is applied

           Does not return anything, but updates internals vectors and variables
        """
        nobj, nvarx= np.shape(X)

        self.nobj = nobj
        self.nvarx = nvarx
        self.X = X
        self.Y = Y

        X, mux = center(X)
        Y, muy = center(Y)
        X, wgx = scale(X, autoscale)

        self.mux = mux
        self.muy = muy
        self.wgx = wgx

        self.autoscale = autoscale
        
        SSXac=0.0
        SSYac=0.0

        SSX0,SSY0, null = self.computeSS(X,Y)
        
        SSXold=SSX0
        SSYold=SSY0

        a=0
        while True:
            t, p, w, c = self.extractLV(X, Y)
                
            self.t.append(t) 
            self.p.append(p)
            self.w.append(w)
            self.c.append(c)
            
            X, Y = self.deflateLV(X, Y, t, p, c)
            
            SSXnew, SSYnew, dmodx = self.computeSS(X, Y)

            SSXex = (SSXold-SSXnew)/SSX0
            SSXac+=SSXex

            SSYex = (SSYold-SSYnew)/SSY0
            SSYac+=SSYex

            SDEC = np.sqrt(SSYnew/nobj)

            dof = nvarx-a
            if dof <= 0 : dof = 1
            dmodx = [np.sqrt(d/dof) for d in dmodx] 

            SSXold=SSXnew
            SSYold=SSYnew

            self.SSXex.append(SSXex)
            self.SSXac.append(SSXac)
            self.SSYex.append(SSYex)
            self.SSYac.append(SSYac)
            self.SDEC.append(SDEC)
            self.dmodx.append(dmodx)

            a+=1
                
            if targetA>0:
                if a==targetA : break

            if targetSSX>0.0:
                if SSXac>targetSSX: break

        self.Am=a

    def validateLOO (self, A):
        """ Validates A dimensions of an already built PLS model, using Leave-One-Out cross-validation

            Returns nothing. The results of the cv (SSY, SDEP and Q2) are stored internally
        """

        if self.X == None or self.Y == None:
            return 
        
        X = self.X
        Y = self.Y

        nobj,nvarx = np.shape (X)

        SSY0 = 0.0
        for i in range (nobj):
            SSY0+=np.square(Y[i]-np.mean(Y))

        SSY = np.zeros(A,dtype=np.float64)
        YP = np.zeros ((nobj,A+1),dtype=np.float64)
        
        for i in range (nobj):
            
            # build reduced X and Y matrices removing i object
            Xr = np.delete(X,i,axis=0)
            Yr = np.delete(Y,i)

            Xr,muxr = center(Xr)
            Xr,wgxr = scale (Xr, self.autoscale)
           
            Yr,muyr = center(Yr)

            xp = np.copy(X[i,:])
            xp -= muxr
            xp *= wgxr

            # predicts y for the i object, using A LV
            yp = self.getLOO(Xr,Yr,xp,A)      
            yp += muyr

            # updates SSY with the object i errors
            YP[i,0]=Y[i]
            
            for a in range(A):
                SSY[a]+= np.square(yp[a]-Y[i])
                YP[i,a+1]=yp[a]

        self.SSY  = SSY        
        self.SDEP = [np.sqrt(i/nobj) for i in SSY]
        self.Q2   = [1.00-(i/SSY0) for i in SSY]
        
        self.Av = A

        return (YP)


    def project (self, x, A):
        """projects query object x into current model using A LV

           Returns
           y:    vector of predicted Y values using growing number of LV
           t:    vector of scores
           d:    SSX for every dimension
        """

        if A > self.Am:
            return (False, 'Too many LV')
                
        x-=self.mux
        x*=self.wgx

        y=np.zeros(A,dtype=np.float64)
        t=np.zeros(A,dtype=np.float64)
        d=np.zeros(A,dtype=np.float64)

        for a in range (A):        
            t[a] = np.dot(x,self.w[a])
            y[a] += t[a]*self.c[a]
            x -= self.p[a]*t[a]
            dof = (self.nvarx-a)
            if dof <= 0 : dof = 1
            d[a] = np.sqrt(np.dot(x.T,x)/dof) 

        y+=self.muy

        return (True, (y, t, d))

        
    def extractLV (self, X, Y):
        """Extracts a single LV from the provided X and Y matrices using NIPALS algorithm

           This method assumes that both X and Y are centered. No deflation is applied
           
           Returns
           t:    vector of scores
           p:    vector of loadings
           w:    vector of weights
           c:    inner relationship
        """
        
        nobj,nvarx = np.shape (X)
        w = np.zeros(nvarx, dtype=np.float64)
        p = np.zeros(nvarx, dtype=np.float64)
        t = np.zeros(nobj , dtype=np.float64)

        uu = np.dot(Y.T,Y)
        for j in range(nvarx):
            w[j] = np.dot(Y.T,X[:,j])/uu
            
        ww = np.sqrt(np.dot(w.T,w))
        w/=ww

        for i in range(nobj):
            t[i] = np.dot(w.T,X[i,:])

        tt = np.dot(t.T,t)

        for j in range(nvarx):
            p[j] = np.dot(t.T,X[:,j])/tt

        c = np.dot(t.T,Y)/tt

        return t, p, w, c


    def computeSS (self, X, Y):
        """Computes the Sum-Of-Squares for provided X and Y matrices

           Returns
           SSX:    sum of squates of the X matrix
           SSY:    sum of squares of the Y matrix
           d:      vector with the SSX for every object 
        """
        
        nobj,nvarx = np.shape (X)
        
        SSX=SSY=0.0

        d = np.zeros(nobj,dtype=np.float64)
        
        for i in range (nobj):
            objx = X[i,:]
            SSX += np.dot(objx.T,objx)
            SSY += np.square(Y[i])
            d[i] = np.dot(objx.T,objx)
            
        return SSX, SSY, d


    def deflateLV (self, X, Y, t, p, c):
        """Deflates both the X and Y matrices, using the provided t, p and c vectors

           Returns deflated X and Y
        """
        
        nobj,nvarx = np.shape (X)
        
        for i in range (nobj):
            ti=t[i]
            X[i,:] -= (ti*p)
            Y[i] -= ti*c
        
        return X,Y


    def getLOO (self, X, Y, x, A):
        """Builds a model of A dimension with the provided X and Y matrices, yielding a prediction y for the query object x.
           Typically used as inner loop in LOO CV method.

           Notice that both X and Y must be centered, while x must have been centered with the model averages

           Returns the predicted y value for the query object x
        """
        
        y = np.zeros(A,dtype=np.float64)
        
        for a in range(A):
            t, p, w, c = self.extractLV(X, Y)
            if a>0 : y[a]=y[a-1]
           
            tt = np.dot(x,w)    
            y[a] += tt*c
            x -= p*tt
            
            X, Y = self.deflateLV(X, Y, t, p, c)
        return y

    def calcOptCutoff (self, ycutoff = 0.5, nsteps = 100):

        by = []
        yp = np.zeros (self.nobj,dtype=np.float64)
        for i in range (self.nobj):
            success, result = self.project(self.X[i,:],self.Am) # just for final #LV
            if success:
                yp[i] = result[0][-1]        
            by.append (self.Y[i] > ycutoff)   # we asume here 1=True, 0=False

        bestv  = 1.0e20
        bestc  = 0.0           
        cutoff = 0.0
        bTP=bTN=bFP=bFN=0
        
        for step in range (nsteps):
            cutoff+=1.0/nsteps
            TP=TN=FP=FN=0

            for i in range(len(yp)):                
                if by[i]:
                    if yp[i] > cutoff:
                        TP+=1
                    else:
                        FN+=1
                else:
                    if yp[i] > cutoff:
                        FP+=1
                    else:
                        TN+=1

            sens = sensitivity (TP, FN)
            spec = specificity (TN, FP)
            
            if abs(sens-spec) < bestv :
                bestv = abs(sens-spec)
                bestc = cutoff
                bTP = TP
                bTN = TN
                bFP = FP
                bFN = FN

        self.cutoff = bestc
        self.TP = bTP
        self.TN = bTN
        self.FP = bTP
        self.FN = bFN
        
        #return (bestc, (bTP,bTN,bFP,bFN))
          

################################################################################################


def readData (filename):
    """Reads numpy X and Y matrices from a file in GOLPE .dat format, asuming a single Y value at the end

       Returns X and Y as a numpy matrices
    """

    f = open (filename)
    line=f.readline()
    line=f.readline()
    nvar=int(line)
    line=f.readline()
    nobj=int(line)

    X = np.zeros((nobj,nvar-1),dtype=np.float64)
    Y = np.zeros(nobj,dtype=np.float64)
    for i in range(nobj):
        line = f.readline()
        line = f.readline()
        for j in range(nvar-1):
            line = f.readline()
            X[i,j]=float(line)
        line = f.readline()
        Y[i]=float(line)

    f.close()
    return X, Y     


if __name__ == "__main__":
    # this is only testing code that can be used as an example of use

    # loads data
    X, Y = readData ('Biopsycho_2A_activity.dat')

    # builds a PLS model
    mypls = pls ()
    mypls.build(X,Y,targetA=5,autoscale=False)     
    mypls.validateLOO(5)
    mypls.saveModel('modelPLS.npy')

    # everything complete, print the results
    for a in range (mypls.Am):
        print "SSXex %6.4f SSXac %6.4f " % \
              (mypls.SSXex[a], mypls.SSXac[a]),
        print "SSYex %6.4f SSYac %6.4f SDEC %6.4f" % \
              (mypls.SSYex[a], mypls.SSYac[a], mypls.SDEC[a])
     
    for a in range (mypls.Av):
        print 'A:%2d  SSY: %6.4f Q2: %6.4f SDEP: %6.4f' % \
              (a+1,mypls.SSY[a],mypls.Q2[a],mypls.SDEP[a])

##    # reloads the data
##    x, y = readData ('test01.dat')
##    nobj,nvarx= np.shape(x)
##
##    # creates a new PLS object, reading the model saved above
##    pls2 = pls ()
##    pls2.loadModel('modelPLS.npy')
##
##    # projects the data on the loaded model
##    for i in range(nobj):
##        success, result = pls2.project(x[i,:],3)
##        if success:
##            yp, tp, dmodx = result
##            #print yp, tp, dmodx[0], pls2.dmodx[0][i]
##            print pls2.dmodx[0][i]
##        else:
##            print result
  
