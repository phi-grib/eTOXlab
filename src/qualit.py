# -*- coding: utf-8 -*-

##    Description    tools for qualitative endpoints
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
