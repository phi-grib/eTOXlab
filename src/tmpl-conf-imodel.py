# -*- coding: utf-8 -*-

##    Description    eTOXlab model template
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

from model import model

class imodel(model):
    def __init__ (self, vpath):
        model.__init__(self, vpath)

        ##
        ## General settings
        ##
        self.buildable = True
        self.quantitative = True
        self.confidential = True
        self.identity = False
        self.experimental = False
        self.SDFileName = 'name'
        self.SDFileActivity = 'activity'
        self.SDFileExperimental = ''

        
        ##
        ## Normalization settings
        ##
        self.norm = False
        self.normStand = True
        self.normNeutr = True
        self.normNeutrMethod = 'moka'
        self.normNeutr_pH = 7.4
        self.norm3D = True

        ##
        ## Molecular descriptor settings
        ##
        self.MD = 'padel'                         # 'padel'|'pentacle'|'adriana'
        self.padelMD = ['-2d']                       # '-2d'|'-3d'
        self.padelMaxRuntime = None
        self.padelDescriptor = None
        self.pentacleProbes = ['DRY','O','N1','TIP']
        self.pentacleOthers = []

        ##
        ## Modeling settings
        ##
        self.model = 'pls'
        self.modelLV = 2
        self.modelAutoscaling = True
        self.modelCutoff = 'auto'
        self.selVar = False
        #self.selVarMethod = GOLPE
        self.selVarLV = 2
        #self.selVarCV = 'LOO'
        self.selVarRun = 2
        self.selVarMask = None

        ##
        ## View settings
        ##
        self.viewType = None    # 'pca' | 'property' | 'project'
        self.viewBackground = False
        self.viewReferenceEndpoint = None
        self.viewReferenceVersion = 0

        self.plotPCAColor = 'red'
        self.plotPCAMarkerShape = 'D'
        self.plotPCAMarkerSize = 40
        self.plotPCAMarkerLine = 0
        
        self.plotPRPColor = 'red'
        self.plotPRPMarkerShape = 'D'
        self.plotPRPMarkerSize = 40
        self.plotPRPMarkerLine = 0
        
        self.plotPRJColor = 'DModX'    # DModX | [color] (e.g. red)
        self.plotPRJMarkerShape = 'o'
        self.plotPRJMarkerSize = 50
        self.plotPRJMarkerLine = 1
        
        self.plotBGColor = '#aaaaaa'
        self.plotBGMarkerShape = 'o'
        self.plotBGMarkerSize = 20
        self.plotBGMarkerLine = 0

        ##
        ## Path to external programs
        ##
        self.mokaPath = '/opt/blabber/blabber4eTOX/'
        self.padelPath = '/opt/padel/padel218ws/'
        self.padelURL = 'http://localhost:9000/computedescriptors?params=' 
        self.pentaclePath = '/opt/pentacle/pentacle106eTOX/'
        self.adrianaPath = '/opt/AdrianaCode/AdrianaCode226/'
        self.corinaPath = '/opt/corina/corina3494/'
        self.javaPath = '/usr/java/jdk1.7.0_51/'
        self.RPath = '/opt/R/R-3.0.2/'
        self.standardiserPath = '/opt/standardiser/standardise20140206/'
