#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    eTOXlab component for managing predictive model
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
import shutil
import cPickle as pickle

from utils import sandVersion
from utils import nextVersion
from utils import lastVersion
from utils import writeError
from utils import wkd
from utils import VERSION

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

def createVersion (endpoint, tag):

    ndir = wkd +'/'+endpoint
    
    # check if there is already a tree for this endpoint
    if os.path.isdir (ndir):
        return (False, 'This endpoint already exists')
    try:
        os.mkdir (ndir)
    except:
        return (False,'unable to create directory '+ndir)

    try:
        f = open (ndir+'/service-label.txt','w')
        f.write (tag+'\n')
        f.close()
    except:
        return (False, 'unable to create service label')

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

    if style in 'long':
        for i in infoID:
            if 'version' in i:
                print i[1]
        for i in infoID:
            if not 'version' in i: print '  %-10s'%i[0],' : '+str(i[1])
        for i in infoSeries:
            print '  %-10s'%i[0],' : '+str(i[1])
        for i in infoMD:
            print '  %-10s'%i[0],' : '+str(i[1])
        for i in infoModel:
            print '  %-10s'%i[0],' : '+str(i[1])
        for i in infoResult:
            print '  %-10s'%i[0],' : '+str(i[1])
        print
            
    elif style in 'short':
        iversion = 2*' ' 
        iMD      = 8*' ' 
        imod     = 16*' '
        imol = isen = ispe = iMCC = ir2 = iq2 = 4*' '
        
        for i in infoID:
            if 'version' in i: iversion = '%-2s'%(i[1])
        for i in infoMD:
            if 'MD' in i: iMD = '%-8s'%(i[1])
        for i in infoModel:
            if 'model' in i : imod =  '%-16s'%(i[1])
        for i in infoResult:
            if 'nobj' in  i: imol = '%4d'%int(i[1])
            if 'sens' in i : isen =  '%4.2f'%(float(i[1]))
            elif 'spec' in i : ispe =  '%4.2f'%(float(i[1]))
            elif 'MCC' in i : iMCC =  '%4.2f'%(float(i[1]))
            elif 'R2' in i : ir2 =  '%4.2f'%(float(i[1]))
            elif 'Q2' in i : iq2 =  '%4.2f'%(float(i[1]))

        if ir2 == '    ':
            print iversion+'  MD:'+iMD+'  mod:'+imod+'  mol:'+imol+'  sen:'+isen+'  spe:'+ispe+'  MCC:'+iMCC
        else:
            print iversion+'  MD:'+iMD+'  mod:'+imod+'  mol:'+imol+'  R2 :'+ir2+'  Q2 :'+iq2
    
    return (True,'OK')

def info (endpoint,ver,style):

    items = []
    
    itemswkd = os.listdir(wkd)
    itemswkd.sort()
    for iendpoint in itemswkd:
       
        if not os.path.isdir(wkd+'/'+iendpoint): continue
        
        if endpoint:
            if iendpoint != endpoint: continue

        tag = ''
        try:
            f = open (wkd+'/'+iendpoint+'/service-label.txt','r')
            tag = f.readline ()[:-1]
            f.close()
        except:
            pass
        
        print 78*'-'
        print iendpoint+' ['+tag+']'
        

        itemend = os.listdir(wkd+'/'+iendpoint)
        itemend.sort()

        vi = -99
        for iversion in itemend:
            
            if not os.path.isdir(wkd+'/'+iendpoint+'/'+iversion): continue
            if not iversion.startswith('version'): continue
            vi = int (iversion[-4:])
            if ver>-99:
                if vi != ver: continue
            
            inform = infoVersion(iendpoint, vi, style)
            items.append (inform)
            
        if ver == -1:
            if vi == -99 : break # in case no version was found exit
            inform = infoVersion(iendpoint, vi, style)
            items.append (inform)

        #print 78*'-'

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
    
    print 'manage  --publish|new|remove|version|info=[short|long]|get=[model|series]] -e endpoint [-v 1|last] [-t tag]'

def main ():

    endpoint = None
    tag = None
    action = None
    infoList = None
    ver = -99
    
    try:
       opts, args = getopt.getopt(sys.argv[1:], 'e:v:t:h', ['publish','new','remove','version','info=', 'get='])

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
            elif opt in '--version':
                print 'version: '+VERSION
                sys.exit(0)

    if not action:
        usage()
        sys.exit (1)
            
        
    ## publish
    if 'publish' in action:

        if not endpoint:
            print 'please provide the name of the endpoint'
            sys.exit (1)

        if ver != -99:
            print 'publish uses version 0 to create a new version. No version must be specified'
            sys.exit (1)

        result = publishVersion (endpoint, tag)

    ## new
    if 'new' in action:

        if not endpoint:
            print 'please provide the name of the endpoint'
            sys.exit (1)

        if not tag:
            print 'please provide the label of the eTOXsys service'
            sys.exit (1)
            
        result = createVersion (endpoint,tag)

    ## remove
    if 'remove' in action:

        if not endpoint:
            print 'please provide the name of the endpoint'
            sys.exit (1)

        if ver != -99:
            print 'remove always removed last published version. No version must be specified'
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

        if not endpoint:
            print 'please provide the name of the endpoint'
            sys.exit (1)

        if ver == -99 :
            print 'please provide the version number or "last" '
            sys.exit (1)
            
        piece=None  
        for i in ['model','series']:
            if getPiece in i:
                piece = i

        if not piece:
            print 'please enter what you want to obtain: "model" or "series" '
            sys.exit (1)
            
        result = get (endpoint,ver,piece)
        
    if result[0]:
        print result[1]
        sys.exit(0)
    else:
        print "ERROR: "+result[1]
        sys.exit(1)
        
if __name__ == '__main__':
    
    main()
