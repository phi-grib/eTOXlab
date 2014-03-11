#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    eTOXlab component for runing a predictive model
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

import sys
import os
import getopt
import cPickle as pickle

from utils import lastVersion
from utils import writeError
from utils import removefile

def predict (endpoint, molecules, verID=-1, api=0, detail=False):
    """Top level prediction function

       molecules:  SDFile containing the collection of 2D structures to be predicted
       verID:      version of the model that will be used. Value -1 means the last one
       detail:     level of detail of the prediction. If True the structure of the
                   closest compond will be returned
    """
    
    # compute version (-1 means last) and point normalize and predict to the right version
    vpath = lastVersion (endpoint,verID)
    if not vpath:
        return (False,"No versions directory found")

    if api>0:
        head, tail = os.path.split (vpath)
        if tail == 'version0000':
            return (False, 'no published model found')
    
    sys.path.append(vpath)
    from imodel import imodel

    # load model
    model = imodel(vpath)

    datList = []
    datList = model.loadData ()
   
    i=0
    pred = []
    mol=''
    fout = None

    # open SDFfile and iterate for every molecule
    try:
        f = open (molecules,'r')
    except:
        return (False,"No molecule found in %s; SDFile format not recognized" % molecules)
    
    for line in f:
        if not fout or fout.closed:
            i += 1
            mol = 'm%0.10d.sdf' % i
            fout = open(mol, 'w')

        fout.write(line)
    
        if '$$$$' in line:
            fout.close()

            ## workflow for molecule i (mol) ###########
            success, result  = model.normalize (mol)
            if not success:
                pred.append((False, result))
                continue

            molFile   = result[0]
            molName   = result[1]
            molCharge = result[2]

            predN = model.predict (molFile, molName, molCharge, detail)
            
            pred.append((True, predN))
            ############################################

            removefile(mol)

    return (True, pred)

def presentPrediction (pred):
    
    """Writes the result of the prediction into a log file and prints some of them in the screen
    """
    
    if pred[0]:
        for x in pred[1]:
            if x[0]:
                for y in x[1]:
                    if y[0]:
                        if isinstance(y[1], float):
                            print "%8.5f" % y[1],
                        else:
                            print y[1],
                    else:
                        print y,
                print
            else:
                print x
    else:
        print pred

def presentPredictionWS2 (pred):
    
    """Writes the result of the prediction into a log file and prints some of them in the screen
    """

    with open('result.txt','w') as fp:
        if pred[0]:
            for compound in pred[1]:  # loop for compounds
                if compound[0]:
                    vaTuple = compound[1][0]
                    adTuple = compound[1][1]
                    riTuple = compound[1][2]

                    fp.write(str(vaTuple[1]))
                    fp.write('\t')
                    fp.write(str(adTuple[1]))
                    fp.write('\t')
                    fp.write(str(riTuple[1]))
                    fp.write('\n')
    

def presentPredictionWS (pred):
    
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

    pkl = open('results.pkl', 'wb')
    pickle.dump(results, pkl)
    pkl.close()

def testimodel():
    try:
        from imodel import imodel
    except:
        return

    print 'please remove file imodel.py or imodel.pyc from eTOXlab/src'
    sys.exit(1)
    
def usage ():
    """Prints in the screen the command syntax and argument"""
    
    print 'predict -e endpoint [-f filename.sdf][-v 1|last][-a|b]'

def main ():

    endpoint = None
    ver = -99
    api = 0
    mol = None

    try:
       opts, args = getopt.getopt(sys.argv[1:], 'abe:f:v:h')

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

            if opt in '-e':
                endpoint = arg               
            elif opt in '-f':
                mol = arg
            elif opt in '-v':
                if 'last' in arg:
                    ver = -1
                else:
                    try:
                        ver = int(arg)
                    except ValueError:
                        ver = -99
            elif opt in '-a':
                api = 1
                mol = './input_file.sdf'
                ver = -1
                # calls from web services might not have PYTHONPATH updated
                sys.path.append ('/opt/RDKit/')
                sys.path.append ('/opt/standardiser/standardise20140206/')
            elif opt in '-b':
                api = 2
                mol = './input_file.sdf'
                ver = -1
                # calls from web services might not have PYTHONPATH updated
                sys.path.append ('/opt/RDKit/')
                sys.path.append ('/opt/standardiser/standardise20140206/')    
            elif opt in '-h':
                usage()
                sys.exit(0)

    if ver == -99:
        usage()
        sys.exit (1)

    if not mol:
        usage()
        sys.exit (1)
        
    if not endpoint:
        usage()
        sys.exit (1)

    # make sure imodel has not been copied to eTOXlab/src. If this were true, this version will
    # be used, instead of those on the versions folder producing hard to track errors and severe
    # misfunction
    testimodel()
    
    result=predict (endpoint,mol,ver,api)
    
    if api==0:
        presentPrediction (result)
    elif api==1:
        presentPredictionWS (result)
    elif api==2:
        presentPredictionWS2 (result)

    sys.exit(0)
        
if __name__ == '__main__':
    
    main()
