#!/usr/bin/env python

# -*- coding: utf-8 -*-
#
#    Description    eTAM component for creating a new predictive model
#                   
#
#    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
#
#    (c) PhI 2013

import sys
import os
import getopt
import shutil

from utils import splitSDF
from utils import lastVersion
from utils import nextVersion
from utils import writeError


def build (molecules, verID=-1):
    """Top level buildind function

       molecules:  SDFile containing the collection of 2D structures to be predicted
       verID:      version of the model that will be used. Value -1 means the last one

    """
    # getMolecule

    # clone directory
    va = lastVersion (verID)
    vb = nextVersion ()
    
    if not va:
        return (False,"No versions directory found")

    shutil.copytree(va,vb)

    # compute version (-1 means last) and point normalize and predict to the right version
    #try:
    sys.path.append(vb)
    from imodel import imodel
    model = imodel (vb)
    #from inormalize import inormalize
    #from ibuild import iextract, ibuild
    #except:
    #    return (False,"Wrong version number (%d)" % verID)

    # split SDFfile into individual molecules
    molList = splitSDF (molecules)
    if not molList:
        return (False,"No molecule found in %s; SDFile format not recognized" % molecules)

    # iterate molList and normalize + predict every molecule
    datList=[]
    for mol in molList:
        molN = model.normalize (mol)
        if not molN[0]:
           datList.append((False,'normalization error'))
           print 'error in normalize'
           continue

        success, infN = model.extract (molN[1:])
        if not success:
           datList.append((False,'data extraction error'))
           print 'error in extract'
           continue

        datList.append((True,infN))

    for imol,idat in zip(molList,datList):
        if not idat[0]:
           molList.remove(imol)
           datList.remove(idat)

    results = model.build (datList)
    return (results)

def writeResults (result):
    """Writes the result of the model building
    """
    print result 

def usage ():
    """Prints in the screen the command syntax and argument"""
    
    print 'build [-f filename.sdf][-v 1]'

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
        
    result=build (mol,ver)

    writeResults (result)
        
if __name__ == '__main__':   
    main()
    sys.exit(0)
