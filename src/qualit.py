# -*- coding: utf-8 -*-
#
#    Description    Scaling tools for PCA and PLS 
#                   
#
#    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
#
#    (c) PhI 2013

from math import sqrt

def sensitivity (TP, FN):
    if (TP+FN) > 0 :
        return (float(TP) / float(TP + FN))
    else:
        return float(0)

def specificity (TN, FP):
    if (TN+FP) > 0 :
        return (float(TN) / float(TN + FP))
    else:
        return float(0)

def MCC (TP, TN, FP, FN):
    d = float(TP+FP)*float(TP+FN)*float(TN+FP)*float(TN+FN)
    if d > 0.0:
        return ((float(TP*TN)-float(FP*FN)) / sqrt(d))
    else:
        return float(0)
