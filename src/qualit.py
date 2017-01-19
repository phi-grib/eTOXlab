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

import numpy as np
import matplotlib.pyplot as plt
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

def FourfoldDisplay(TP, TN, FP, FN, label, name):
    """ Draws confusion matrix graphical representaion

    """

    width = np.pi / 2.0
    theta = np.radians([0,90,180,270])
    table = [FP,TP,FN,TN]
    plt.figure("RF-Qualitative_validation")    
    plt.clf()
    ax = plt.subplot(121, polar=True, adjustable='box', aspect=1)    
    bars = ax.bar(theta, table, width=width, color=["red", "green", "red", "green"])
    plt.title( label + ' Confusion Matrix')

    ax.set_xticklabels(["","FP (%s)" % str(FP), "",  "TP (%s)" % str(TP), "", "FN (%s)" % str(FN), 
                        "",  "TN (%s)" % str(TN)], fontsize=20)
    ax.set_yticks([])
    ax.grid(False)
    ax.axes.spines['polar'].set_visible(False)

    ax2 = plt.subplot(122, adjustable='box', aspect=3)
    plt.ylim([0,1])
    plt.title( 'Sensitivity and Specifity')
    bar_width = 0.5
    y = [0, sensitivity(TP,FN), specificity(TN,FP), 0]
    index = np.arange(4)
    ax2.bar(index, y, bar_width)
    #ax.offset(0.5)
    plt.xticks( index + bar_width / 2.0, ("", 'Sensitivity', 'Specifity', ""))

    plt.savefig(name)
