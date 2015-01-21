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
import tarfile

from utils import sandVersion
from utils import nextVersion
from utils import lastVersion
from utils import removefile
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

    if not os.path.isfile(va+'/info.pkl'):
        return (False,"No suitable model found")

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

    # add the type of endpoint (quantitative or qualitative) after the tag
    # this is needed by views2 to publish appropriately the model type
    ndir = wkd +'/'+endpoint

    if not os.path.isfile(ndir+'/service-label.txt'):
        return (False, 'unable to open service-label.txt file')
    
    f = open (ndir+'/service-label.txt','r')
    tag = f.readline()
    f.close()

    f = open (ndir+'/service-label.txt','w')
    f.write (tag)
    ytype='undefined'
    for i in infoID:
        if 'dependent' in i:
                ytype = i[1]
    f.write (ytype+'\n')
    f.close()

    return (True, vb)


def exposeVersion (endpoint, ver):
    edir = wkd +'/'+endpoint
    
    # check if there is already a tree for this endpoint
    if not os.path.isdir (edir):
        return (False, 'This endpoint do not exists')

    if os.path.isfile (edir+'/service-version.txt'):
        f = open (edir+'/service-version.txt','r')
        oldver = int(f.readline())
        f.close()
        if oldver == ver:
            removefile(edir+'/service-version.txt')
            return (True, 'version un-exposed OK')
        
    try:
        f = open (edir+'/service-version.txt','w')
        f.write (str(ver))
        f.close()
    except:
        return (False, 'unable to create version label')

    return (True, 'version exposed OK')

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

def createConfVersion (endpoint, tag):

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
        shutil.copy(wkd+'/tmpl-conf-imodel.py',ndir+'/imodel.py')
    except:
        return (False,'unable to create imodel.py at '+ndir)

    return (True,'version (for confidential model development) created OK')
           
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

def exportEndpoint (endpoint):

    currentPath = os.getcwd ()
    exportfile = currentPath+'/'+endpoint+'.tgz'
    
    os.chdir(wkd)
    if not os.path.isdir(endpoint):
        return (False, 'endpoint directory not found')

    tar = tarfile.open(exportfile, 'w:gz')
    tar.add(endpoint+'/service-label.txt')
    
    itemend = os.listdir(endpoint)
    itemend.sort()

    modconfY=False
    modconfN=False
    for iversion in itemend:
        if not os.path.isdir(wkd+'/'+endpoint+'/'+iversion): continue
        if not iversion.startswith('version'): continue

        # confidential models are recognized by file "distiledPLS.txt"
        if os.path.isfile (wkd+'/'+endpoint+'/'+iversion+'/distiledPLS.txt'):
            tar.add(endpoint+'/'+iversion+'/imodel.py')
            tar.add(endpoint+'/'+iversion+'/distiledPLS.txt')
            tar.add(endpoint+'/'+iversion+'/info.pkl')
            modconfY=True
        else:
            tar.add(endpoint+'/'+iversion)
            modconfN=True
            
    tar.close()

    # avoid exporting mixtures of confidential and non-confidential models. This could be a major security breach
    if modconfY and modconfN:
        os.remove (exportfile)
        return (False, 'endpoint '+endpoint+' contains a mixture of confidential and non-confidential models. EXPORT ABORTED')
   
    return (True,'endpoint '+endpoint+' exported as '+endpoint+'.tgz')

def importEndpoint (endpoint):

    if os.path.isdir (wkd+'/'+endpoint):
        return (False, 'endpoint already existing')

    importfile = endpoint+'.tgz'
    
    if not os.path.isfile (importfile):
        return (False, 'importing package '+importfile+' not found')
    
    shutil.copy (endpoint+'.tgz',wkd)   
    os.chdir(wkd)
   
    tar = tarfile.open(importfile, 'r:gz')
    tar.extractall()
    tar.close()

    os.remove (importfile)
    
    return (True,'endpoint '+endpoint+' imported OK')

def infoVersion (endpoint,ver,style,isWS):

    vb = lastVersion (endpoint,ver)
    unk = not os.path.isfile (vb+'/info.pkl')
##    
##    if :
##        if ver == 0:
##            print '*     no model info available'
##        else:
##            print '%-4s  no model info available'%ver
##        return (True, 'OK')
##    
##        #return (False,'model information file not found')

    if not unk:
        modelInfo = open (vb+'/info.pkl','rb')
        infoID = pickle.load(modelInfo)
        infoSeries = pickle.load(modelInfo)
        infoMD = pickle.load(modelInfo)
        infoModel = pickle.load(modelInfo)
        infoResult = pickle.load(modelInfo)
        modelInfo.close()

    if style in 'long':

        if unk :
            if ver == 0: print '*'
            else       : print '%-4s'%ver
            print '  web service :', isWS
            print '  no model info available'
            return (True, 'OK')
            
        for i in infoID:
            if 'version' in i: print i[1]
        
        for i in infoID:
            if not 'version' in i: print '  %-10s'%i[0],' : '+str(i[1])

        print "  web service :", isWS
        
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
        iversion = 4*' ' 
        iMD      = 8*' ' 
        imod     = 16*' '
        imol = isen = ispe = iMCC = ir2 = iq2 = 4*' '

        iconf = ''
        if isWS : ws = ' @ '
        else    : ws = '   '

        if unk :
            print '%-4s'%ver + ws + 'no model info available'
            return (True, 'OK')
        
        for i in infoID:
            if 'version' == i[0] : iversion = '%-4s'%(i[1])
            if 'confident' == i[0] and i[1]=='True' : iconf = '  confident'
        for i in infoMD:
            if 'MD' == i[0] : iMD = '%-8s'%(i[1])
        for i in infoModel:
            if 'model' == i[0] : imod =  '%-16s'%(i[1])
        for i in infoResult:
            if 'nobj' == i[0] : imol = '%4d'%int(i[1])
            if 'sens' == i[0] : isen =  '%4.2f'%(float(i[1]))
            elif 'spec' == i[0] : ispe =  '%4.2f'%(float(i[1]))
            elif 'MCC' == i[0] : iMCC =  '%4.2f'%(float(i[1]))
            elif 'R2' == i[0] : ir2 =  '%4.2f'%(float(i[1]))
            elif 'Q2' == i[0] : iq2 =  '%4.2f'%(float(i[1]))

        if ir2 == '    ':
            print iversion+ws+'MD:'+iMD+'  mod:'+imod+'  mol:'+imol+'  sen:'+isen+'  spe:'+ispe+'  MCC:'+iMCC+iconf
        else:
            print iversion+ws+'MD:'+iMD+'  mod:'+imod+'  mol:'+imol+'  R2 :'+ir2+'  Q2 :'+iq2+iconf

    
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

        wsID = -999
        try:
            f = open (wkd+'/'+iendpoint+'/service-version.txt','r')
            wsID = int(f.readline ())
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
            
            inform = infoVersion(iendpoint, vi, style, vi==wsID)
            items.append (inform)
            
        if ver == -1:
            if vi == -99 : break # in case no version was found exit
            inform = infoVersion(iendpoint, vi, style, vi==wsID)
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
    
    print 'manage  --publish|expose|new|conf|remove|import|export|version|info=[short|long]|get=[model|series]] -e endpoint [-v 1|last] [-t tag]'
    

def printResult (result):

    if result[0]:
        print result[1]
        sys.exit(0)
    else:
        print "\nERROR:", result[1]
        sys.exit(1)

def main ():

    endpoint = None
    tag = None
    action = None
    infoList = None
    ver = -99
    
    try:
       opts, args = getopt.getopt(sys.argv[1:], 'e:v:t:h', ['publish','expose','new','conf','remove','export','import','version','info=', 'get='])

    except getopt.GetoptError:
       usage()
       printResult ((False, "Arguments not recognized"))

    if args:
       usage()
       printResult ((False, "Arguments not recognized"))
        
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
            elif opt in '--expose':
                action = 'expose'
            elif opt in '--new':
                action = 'new'
            elif opt in '--conf':
                action = 'conf'
            elif opt in '--remove':
                action = 'remove'
            elif opt in '--import':
                action = 'import'
            elif opt in '--export':
                action = 'export'
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
        printResult ((False, 'wrong command syntax'))
                    
    ## publish
    if 'publish' in action:

        if not endpoint:
            printResult ((False, 'please provide the name of the endpoint'))

        if ver != -99:
            printResult ((False, 'publish uses version 0 to create a new version. No version must be specified'))

        result = publishVersion (endpoint, tag)
        printResult (result)
        
    ## expose
    if 'expose' in action:

        if not endpoint:
            printResult ((False, 'please provide the name of the endpoint'))

        if ver == -99:
            printResult ((False, 'please provide the version number or "last"'))

        if ver == 0:
            printResult ((False, 'version 0 (the sandbox) cannot be published'))
            
        result = exposeVersion (endpoint, ver)
        printResult (result)

    ## new
    if 'new' in action:

        if not endpoint:
            printResult ((False, 'please provide the name of the endpoint'))

        if not tag:
            printResult ((False, 'please provide the label of the eTOXsys service'))
            
        result = createVersion (endpoint,tag)
        printResult (result)
        
    ## conf
    if 'conf' in action:

        if not endpoint:
            printResult ((False, 'please provide the name of the endpoint'))

        if not tag:
            printResult ((False, 'please provide the label of the eTOXsys service'))
            
        result = createConfVersion (endpoint,tag)
        printResult (result)

    ## remove
    if 'remove' in action:

        if not endpoint:
            printResult ((False, 'please provide the name of the endpoint'))

        if ver != -99:
            printResult ((False, 'remove always removed last published version. No version must be specified'))
            
        result = removeVersion (endpoint)
        printResult (result)

    ## import
    if 'import' in action:

        if not endpoint:
            printResult ((False, 'please provide the name of the endpoint'))
            
        result = importEndpoint (endpoint)
        printResult (result)
        
    ## export
    if 'export' in action:

        if not endpoint:
            printResult ((False, 'please provide the name of the endpoint'))
            
        result = exportEndpoint (endpoint)
        printResult (result)
        
    ## info
    if 'info' in action:

        style='long'  # default
        for i in ['short','long']:
            if infoStyle in i:
                style = i
            
        result = info (endpoint,ver,style)
        printResult(result)

    ## get
    if 'get' in action:

        if not endpoint:
            printResult ((False, 'please provide the name of the endpoint'))

        if ver == -99 :
            printResult ((False, 'please provide the version number or "last" '))
            
        piece=None  
        for i in ['model','series']:
            if getPiece in i:
                piece = i

        if not piece:
            printResult ((False, 'please enter what you want to obtain: "model" or "series" '))
            
        result = get (endpoint,ver,piece)
        printResult(result)
        
        
if __name__ == '__main__':
    
    main()
