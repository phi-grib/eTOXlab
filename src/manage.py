#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    eTAM component for managing predictive model
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

import sys
import os
import getopt
import shutil
import cPickle as pickle

from utils import sandVersion
from utils import nextVersion
from utils import lastVersion
from utils import writeError
from utils import wkd


def publishVersion (endpoint, tag):
    """Top level buildind function

       molecules:  SDFile containing the collection of 2D structures to be predicted
       verID:      version of the model that will be used. Value -1 means the last one

    """

    # clone directory
    va = sandVersion (endpoint)
    vb = nextVersion (endpoint)
    
    if not va:
        return (False,"No versions directory found")

    shutil.copytree(va,vb)

    modelInfo = open (vb+'/info.pkl','rb')
    infoID = pickle.load(modelInfo)
    infoSeries = pickle.load(modelInfo)
    infoMD = pickle.load(modelInfo)
    infoModel = pickle.load(modelInfo)
    infoResult = pickle.load(modelInfo)
    modelInfo.close()

    for i in range(len(infoID)):
        if infoID[i][0]=='version': 
            infoID.remove (infoID[i])
            infoID.insert (i,('version', int (vb[-4:])))

    if not tag: tag = 'none'

    infoID.append (('tag',tag))

    modelInfo = open (vb+'/info.pkl','wb')
    pickle.dump(infoID, modelInfo)
    pickle.dump(infoSeries, modelInfo)
    pickle.dump(infoMD, modelInfo)
    pickle.dump(infoModel, modelInfo)
    pickle.dump(infoResult, modelInfo)
    modelInfo.close() 

    return (True, vb)

def createVersion (endpoint):

    ndir = wkd +'/'+endpoint
    
    # check if there is already a tree for this endpoint
    if os.path.isdir (ndir):
        return (False, 'This endpoint already exists')
    try:
        os.mkdir (ndir)
    except:
        return (False,'unable to create directory '+ndir)

    ndir+='/version0000'
    try:
        os.mkdir (ndir)
    except:
        return (False,'unable to create directory '+ndir)   
    try:
        shutil.copy(wkd+'/tmpl-imodel.py',ndir+'/imodel.py')
    except:
        return (False,'unable to create imodel.py at '+ndir)

    return (True,'version created OK')

           
def removeVersion (endpoint):

    va = sandVersion (endpoint)
    vb = lastVersion (endpoint,-1)

    if va == vb:
        return (False, 'no more removable versions')

    try:
        shutil.rmtree (vb, ignore_errors=True)
    except:
        return (False, 'unable to remove '+vb)
    
    return (True,'version '+vb+' removed OK')

     
def infoVersion (endpoint,ver,style):

    vb = lastVersion (endpoint,ver)
    
    if not os.path.isfile (vb+'/info.pkl'):
        return (False,'model information file not found')
    
    modelInfo = open (vb+'/info.pkl','rb')
    infoID = pickle.load(modelInfo)
    infoSeries = pickle.load(modelInfo)
    infoMD = pickle.load(modelInfo)
    infoModel = pickle.load(modelInfo)
    infoResult = pickle.load(modelInfo)
    modelInfo.close()

    shortList = ['version','MD','model','nobj','R2','Q2','sens','spec','MCC']  # edit to build a reasonable list of items to show

    if style in 'long':
        for i in infoID:
            print i[0]+':'+str(i[1])+'  ',
        print
        for i in infoSeries:
            print i[0]+':'+str(i[1])+'  ',
        print
        for i in infoMD:
            print i[0]+':'+str(i[1])+'  ',
        print
        for i in infoModel:
            print i[0]+':'+str(i[1])+'  ',
        print
        for i in infoResult:
            print i[0]+':'+str(i[1])+'  ',
        print
            
    elif style in 'short':
        for i in infoID:
            if i[0] in shortList :  print i[0]+':'+str(i[1])+'  ',
        for i in infoSeries:
            if i[0] in shortList :  print i[0]+':'+str(i[1])+'  ',
        for i in infoMD:
            if i[0] in shortList :  print i[0]+':'+str(i[1])+'  ',
        for i in infoModel:
            if i[0] in shortList :  print i[0]+':'+str(i[1])+'  ',
        for i in infoResult:
            if i[0] in shortList :  print i[0]+':'+str(i[1])+'  ',
            
    print
    
    return (True,'OK')

def info (endpoint,ver,style):

    items = []

    print 'working: ', wkd
    
    itemswkd = os.listdir(wkd)
    itemswkd.sort()
    for iendpoint in itemswkd:
       
        if not os.path.isdir(wkd+'/'+iendpoint): continue
        if endpoint:
            if iendpoint != endpoint: continue

        print 'endpoint: '+iendpoint

        itemend = os.listdir(wkd+'/'+iendpoint)
        itemend.sort()
        for iversion in itemend:
            
            if not os.path.isdir(wkd+'/'+iendpoint+'/'+iversion): continue
            if not iversion.startswith('version'): continue
            vi = int (iversion[-4:])
            if ver>-99:
                if vi != ver: continue
            
            inform = infoVersion(iendpoint, vi, style)
            items.append (inform)

    correct = 0
    for i in items:
        if i[0] : correct+=1

    if correct == len(items):
        return (True, 'All requested models informed OK')
    else:
        return (False, '%d models found, %d reported correctly' % (len(items), correct))

def get (endpoint, ver, piece):

    if 'model' in piece:
        piece_name = '/imodel.py'
    elif 'series' in piece:
        piece_name = '/training.sdf'
        
    vb = lastVersion (endpoint,ver)
    try:
        shutil.copy(vb+piece_name,'./')
    except:
        return (False,'Unable to copy '+piece_name)

    return (True, 'File retrieved OK')

def usage ():
    """Prints in the screen the command syntax and argument"""
    
    print 'manage  --publish|new|remove|info=[short|long]|get=[model|series]] -e endpoint [-v 1|last] [-t tag]'

def main ():

    endpoint = None
    tag = None
    action = None
    infoList = None
    ver = -99
    
    try:
       opts, args = getopt.getopt(sys.argv[1:], 'e:v:t:h', ['publish','new','remove','info=', 'get='])

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
            elif opt in '-t':
                tag = arg
            elif opt in '-v':
                if 'last' in arg:
                    ver = -1
                else:
                    try:
                        ver = int(arg)
                    except ValueError:
                        ver = -99
            elif opt in '--publish':
                action = 'publish'
            elif opt in '--new':
                action = 'new'
            elif opt in '--remove':
                action = 'remove'
            elif opt in '--info':
                action = 'info'
                infoStyle = arg
            elif opt in '--get':
                action = 'get'
                getPiece = arg
            elif opt in '-h':
                usage()
                sys.exit(0)

    if not action:
        usage()
        sys.exit (1)
            
        
    ## publish
    if 'publish' in action:

        if not endpoint:
            usage()
            sys.exit (1)

        result = publishVersion (endpoint, tag)

    ## new
    if 'new' in action:

        if not endpoint:
            usage()
            sys.exit (1)
            
        result = createVersion (endpoint)

    ## remove
    if 'remove' in action:

        if not endpoint:
            usage()
            sys.exit (1)
            
        result = removeVersion (endpoint)

    ## info
    if 'info' in action:

        style='long'  # default
        for i in ['short','long']:
            if infoStyle in i:
                style = i
            
        result = info (endpoint,ver,style)

    ## get
    if 'get' in action:

        piece='model'  # default
        for i in ['model','series']:
            if getPiece in i:
                piece = i
            
        result = get (endpoint,ver,piece)
        
    if result[0]:
        print result[1]
        sys.exit(0)
    else:
        print "ERROR: "+result[1]
        sys.exit(1)

        
if __name__ == '__main__':
    
    main()
