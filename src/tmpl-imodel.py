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
        self.confidential = False
        
        ##
        ## Normalization settings
        ##
        self.norm = True
        self.normStand = True
        self.normNeutr = True
        self.normNeutrMethod = 'moka'
        self.normNeutr_pH = 7.4
        self.norm3D = True

        ##
        ## Molecular descriptor settings
        ##
        self.MD = 'pentacle'                         # 'padel'|'pentacle'
        self.padelMD = ['-3d']                       # '-2d'|'-3d'
        self.padelMaxRuntime = None
        self.padelDescriptor = None
        self.pentacleProbes = ['DRY','O','N1','TIP']
        self.pentacleOthers = []

        ##
        ## Modeling settings
        ##
        self.model = 'pls'
        self.modelLV = 2
        self.modelAutoscaling = False
        self.modelCutoff = 'auto'
