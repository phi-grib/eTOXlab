# -*- coding: utf-8 -*-
#
#    Description    PCA toolkit using NIPALS algorithm
#                   
#
#    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
#
#    (c) PhI 2013

import numpy as np

def center (X):
    """Centers the numpy matrix (X) provided as argument"""   

    mu = np.mean(X, axis=0)
    return X-mu, mu

    
def extractPC (X):
    """Computes a single PC from the numpy X matrix (X) provided as argument
       using the NIPALS algorithm

       NIPALS-PCA is iterative and runs until convergence. Criteria used here are:
       - less than 100 iterarios
       - changes in any p value <= 1.0E-9

       Returns two numpy vectors
       t:    scores
       p:    loadings
    """

    nobj,nvar = np.shape (X)
    p = np.zeros(nvar)
    pold = np.zeros(nvar)
    t = np.zeros(nobj)
    
    ttmax = 0.00
    for k in range(nvar):
        obj = X[:,k]
        tt = np.dot(obj.T,obj)
        if tt>ttmax:
            ttmax = tt
            tti = k

    t=np.copy(X[:,tti])
 
    for iter in range (100):  # max 100 iterations

        # (ii) p' = t'X/t't
        for k in range(nvar):
            p[k] = np.dot(t.T,X[:,k])/ np.dot(t.T,t)

        # (iii) normalice P to length 1
        p /= np.sqrt(np.dot(p.T,p))

        # (iv) t = Xp/p'p)
        for j in range(nobj):
            t[j] = np.dot(X[j,:],p) / np.dot(p.T,p)

        # check convergence
        if max(pold-p) > 1.0e-9:  # convergence criteria set to 1.0e-9
            pold = np.copy(p)
        else:
            print 'converges after '+str(iter+1)+' steps'
            break

    return t,p


def deflatePC (X, t, p):
    """Deflates the numpy X matrix (X) using the numpy scores (t) and loadings (p) 
       using the NIPALS algorithm

       Returns
       X:     the deflated X matrix
       SSX:   Sum-of-squares of the X matrix before the deflation
       SSXex: Sum-of-squares of the scores vector, hence explained by this PC
    """
    nobj,nvar = np.shape (X)
    SSX=0.0
    for i in range (nobj):
        obj = X[i,:]
        SSX += np.dot(obj.T,obj)
        obj -= (t[i]*p)
    SSXex = np.dot(t.T,t)
    return X, SSX, SSXex

    
def projectPC (X, mu, p, a):
    """The numpy X matrix (X) is projected into an existing PCA model to extract a single PC

       This call is repeated A times (one for each model dimension) passing the deflated X matrix in
       each call
      
       The value of a is only used to check if this is the first call. If true, the matrix is centered using
       the model mean vector (mu)
       
       Returns three numpy objects  
       X:    deflated X matrix
       t:    scores
       d:    distance to model   
    """
    
    nobj,nvar = np.shape (X)

    if a==0: X-=mu
    
    t = np.zeros(nobj)
    d = np.zeros(nobj)
    
    for i in range (nobj):
        obj = X[i,:]

        # obtain scores for object
        t[i] = np.dot(obj,p)

        # deflates X
        obj -=p*t[i]

        # DModX computed after deflating
        d[i] = np.dot(obj.T,obj)
        d[i] = np.sqrt(d[i]/(nvar-a)) # must be divided by SSX/DOF for the model!

    return X, t, d


def readData (filename):
    """Reads a numpy X matrix from a file in GOLPE .dat format

       Returns the X matrix as a numpy matrix
    """

    f = open (filename)
    line=f.readline()
    line=f.readline()
    nvar=int(line)
    line=f.readline()
    nobj=int(line)

    X = np.zeros((nobj,nvar),dtype=np.float64)
    for i in range(nobj):
        line = f.readline()
        line = f.readline()
        for j in range(nvar):
            line = f.readline()
            X[i,j]=float(line)

    f.close()
    return X   


if __name__ == "__main__":
    
    # this is only testing code that can be used as an example of use

    # loads data
    X = readData('test01.dat')
    Q = readData ('test01.dat')
    print X

    # center X
    X, mu = center(X)
    print X

    SSXpca=0.0
    SSXtot=0.0

    for a in range(4):
        # extracts LV
        t, p = extractPC(X)
        print t
        print p

        # deflates X
        X, SSX, SSXex = deflatePC(X,t,p)
        if a==0: SSXtot = SSX

        SSXpca += SSXex
        print SSXex/SSXtot, SSXpca/SSXtot

        Q, tq, dq = projectPC (Q, mu, p, a)
        print tq

