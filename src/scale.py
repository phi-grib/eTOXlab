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
    for s, w in zip(st,wg):
        if s>10e-6: w/=s

    return X*wg, wg
