#!/usr/bin/env python

# -*- coding: utf-8 -*-
#
#    Description    eTAM component for running a predictive model
#                   
#
#    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
#
#    (c) PhI 2013

import sys
import os
import getopt

from utils import splitSDF
from utils import lastVersion
from utils import writeError

def predict (molecules, verID=-1, detail=False):
    """Top level prediction function

       molecules:  SDFile containing the collection of 2D structures to be predicted
       verID:      version of the model that will be used. Value -1 means the last one
       detail:     level of detail of the prediction. If True the structure of the
                   closest compond will be returned
    """
    
    # compute version (-1 means last) and point normalize and predict to the right version
    vpath = lastVersion (verID)
    if vpath:
        sys.path.append(vpath)
        from imodel import imodel
        model = imodel(vpath)
    else:
        return (False,"No versions directory found")
    
    # split SDFfile into individual molecules
    molList = splitSDF (molecules)
    if not molList:
        return (False,"No molecule found in %s; SDFile format not recognized" % molecules)

    # iterate molList and normalize + predict every molecule
    pred=[]
    for mol in molList:
        molN  = model.normalize (mol)
        if not molN[0]:
##            print 'error in normalization'
            pred.append(molN)
            continue
        
        predN = model.predict (molN[1:], detail)
        pred.append((True,predN))

    return (True, pred)

def writePrediction (pred):
    """Writes the result of the prediction into a log file and prints some of them in the screen
    """
    # print predicted value or 'NA'  
    if pred[0]:
        for x in pred[1]:
            if x[0]:
                for y in x[1]:
                    if y[0]:
                        if y[1] is float:
                            print "%8.3f" % y[1],
                        else:
                            print y[1]
                    else:
                        print y
                print
            else:
                print x
    else:
        print pred

##    f = open ('results.txt','w')
##    f.write (str(pred))
##    f.close()
##    
##    if pred[0]:
##        for x in pred[1]:
##            if x[0]:
##                y = x[1][0]
##                if y[0]:
##                    print "%.3f" % y[1]
##                else:
##                    print 'NA'
##            else:
##                print 'NA'
##    else:
##        print 'NA'

def predict_eTOXsys (tfile,version=-1):
    
    pred = predict (tfile,version)
    
##    f = open ('results.txt','w')
##    f.write (str(pred))
##    f.close()

    # results are stored in a list of tuples (stat,val,msg)
    results = []
    
    if pred[0]:
        #loop for compounds
        for x in pred[1]:
            val = ''
            msg = ''
            stat = 1
            if x[0]:
                y = x[1][0]
                if y[0]:
                    #val = float(y[1])
                    val = y[1]
                    stat = 0
                else:
                    msg = str(y[1])
            else:
                msg = str(x[1])
            results.append((val,stat,msg))
    else:
        stat = 1
        val = ''
        msg = str(pred[1])
        results.append((val,stat,msg))

    return results

    
def usage ():
    """Prints in the screen the command syntax and argument"""
    
    print 'predict [-f filename.sdf][-v 1]'

def main ():
    
    ver = -1
    mol = 'test1.sdf'  # remove!

    try:
       opts, args = getopt.getopt(sys.argv[1:], 'f:v:h')

    except getopt.GetoptError:
       writeError('Error. Arguments not recognized')
       usage()
       sys.exit(1)

    if args:
       writeError('Error. Arguments not recognized')
       usage()
       sys.exit(1)
        
    if len( opts ) > 0:
        for opt, arg in opts:
            if opt in '-f':
                mol = arg
            elif opt in '-v':
                ver = int(arg)
            elif opt in '-h':
                usage()
                sys.exit(0)
        
    result=predict (mol,ver)

    writePrediction (result)
        
if __name__ == '__main__':

##    results = predict_eTOXsys()
##    for val,stat, msg in results:
##        print 'val:', val, 'stat :', stat, 'msg :', msg
    
    main()
    sys.exit(0)
