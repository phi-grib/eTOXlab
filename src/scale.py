# -*- coding: utf-8 -*-
#
#    Description    Scaling tools for PCA and PLS 
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

def scale (X, autoscale):
    """Scales the numpy matrix (X) provided as argument by
       dividing the values by the standard deviation

       The standard deviation is computed as SSX/(N-1)

       Return the scaled matrix and the inverse of the sd (weights)
    """
    nobj, nvar= np.shape(X)

    wg = np.ones (nvar,dtype=np.float64)
    if not autoscale:
        return X, wg
    
    st = np.std (X, axis=0, ddof=1)

    for i in range (nvar):
        sti=st[i]
        if sti<1.0e-7:
            wg[i]=0.0
        else:
            wg[i]/=sti

    #wg=np.select ([st<1.0e-7,st>=1.0e-7],[0.0,wg/st])
            
    #wg/=st
    #wg[st<1.0e-7]=0.0 # the weight of variables with small var is set to 0

    return X*wg, wg
