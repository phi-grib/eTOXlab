#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    eTOXlab component for creating a new predictive model
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

from utils import *
  
def licenseTest (endpoint, verID):
    
    if (verID!=None):
        vv = lastVersion (endpoint, verID)  # full path to endpoint+version or last if -1 is provided
    
    # load model
    try:
        sys.path.append(vv)
        from imodel import imodel
        model = imodel (vv)
    except:
        return (False, 'unable to load imodel')

    if not model:
        return (False, 'unable to load imodel')
    
    result = model.licenseTesting()

    return (result)


def usage ():
    """Prints in the screen the command syntax and argument"""
    
    print 'ERROR: licensetest -e endpoint [-v 1|last]'

def main ():

    endpoint = None
    ver = None

    try:
       opts, args = getopt.getopt(sys.argv[1:], 'e:v:h')

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
            elif opt in '-v':
                if 'last' in arg:
                    ver = -1
                else:
                    try:
                        ver = int(arg)
                    except ValueError:
                        ver = None

    if endpoint == None or ver==None:
        usage()
        sys.exit(1)
    
    result=licenseTest (endpoint,ver)

    print result

    sys.exit(0)
        
if __name__ == '__main__':
    
    main()
