# -*- coding: utf-8 -*-
#
#    Description    eTAM utility functions 
#                   
#
#    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
#
#    (c) PhI 2013

import os
import shutil
import string
import random

#opt = '/opt/'   # full path to the external software 
#wkd = '/root/soft/eTAM/src'    # full path of the working directory
opt = os.environ['ETAM_OPT']
wkd = os.path.dirname(os.path.abspath(__file__))
print 'this is the ETAM_HOME: ', wkd
#wkd = os.environ['ETAM_HOME']

def removefile(file):
    """Removes silently files or whole paths.

       It checks the existence of the file in the argument and if it does not exists do not issue any error
    """
    if not os.path.exists(file):
        return
    if os.path.isfile(file):
        try:
            os.remove(file)
        except OSError:
            pass
    elif os.path.isdir(file):
        try:
            shutil.rmtree( file )
        except OSError :
            pass

def randomName (size=10, chars=string.ascii_uppercase + string.digits):
    name = ''
    for i in range (size):
        name+=random.choice(chars)
    return name

def splitSDF (molecules):
    """Splits a SDFile containing 1 to n molecules

       Every molecule is written to a ""m0000000001.sdf"" file and the full list of names is returned in a list
    """

    try:
        finp = open( molecules)
    except:
        return None

    molList = []
    count = 0
    fout = None

    for line in finp:
        if not fout or fout.closed:
            count += 1
            molList.append ("m%0.10d.sdf" % count)
            fout = open(molList[-1], "w")

        fout.write(line)
        
        if "$$$$" in line:
            fout.close()

    finp.close()
    fout.close()

    return molList

def lastVersion (verID):
    """Returns the path to the directory where the verID version of the model is located

       If a value of -1 is provided, analyzes the "version0000" directories and returns the last one
    """

    if verID==-1:       # we request the last version
        v = [int(x[7:]) for x in os.listdir (wkd) if x.startswith("version")]
        if not v:
            return None
        verID = max(v)

    return wkd+"/version%0.4d" % verID

def nextVersion ():
    """Returns the path to the directory where the verID version of the model is located

       If a value of -1 is provided, analyzes the "version0000" directories and returns the last one
    """
    v = [int(x[7:]) for x in os.listdir (wkd) if x.startswith("version")]
    if not v:
        return None
    verID = max(v)
    verID+=1

    return wkd+"/version%0.4d" % verID

def writeError (error):
    """Print an error message"""
    
    print error

