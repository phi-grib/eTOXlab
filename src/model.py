# -*- coding: utf-8 -*-

##    Description    eTOXlab model class
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

import os
import sys
import shutil
import subprocess
import cPickle as pickle
import time
import urllib2

import matplotlib
from pylab import *
import matplotlib.pyplot as plt

import numpy as np

from pls import pls
from pca import pca
from StringIO import StringIO
from utils import removefile
from utils import randomName
from utils import updateProgress
from utils import writeError
from utils import wkd

from rdkit import Chem
from rdkit import RDLogger
from standardise import standardise
from qualit import *
from rdkit.Chem import Descriptors

class model:
            
    def __init__ (self, vpath):

        self.vpath = vpath
        self.tdata = []

        ##
        ## General settings
        ##
        self.buildable = False
        self.quantitative = False
        self.confidential = False
        self.identity = False
        self.SDFileName = ''
        self.SDFileActivity = ''
        
        ##
        ## Normalization settings
        ##
        self.norm = False
        self.normStand = False
        self.normNeutr = False
        self.normNeutrMethod = None
        self.normNeutr_pH = None
        self.norm3D = False

        ##
        ## Molecular descriptor settings
        ##
        self.MD = None
        self.padelMD = [] 
        self.padelMaxRuntime = None
        self.padelDescriptor = None
        self.pentacleProbes = []
        self.pentacleOthers = []

        ##
        ## Modeling settings
        ##
        self.model = None
        self.modelLV = None
        self.modelAutoscaling = None
        self.modelCutoff = None
        self.selVar = None
        #self.selVarMethod = None
        self.selVarLV = None
        #self.selVarCV = None
        self.selVarRun = None
        self.selVarMask = None

        ##
        ## Path to external programs
        ##
        self.mokaPath = None
        self.padelPath = None
        self.padelURL = None 
        self.pentaclePath = None
        self.adrianaPath = None
        self.corinaPath = None
        self.javaPath = None
        self.RPath = None
        self.standardiserPath = None
        
        # Info lists serve only to store properties of new models
        # and inform the users. This list does not set model properties
        self.infoID = []
        self.infoSeries = []
        self.infoMD = []
        self.infoModel = []
        self.infoResult = []

        ##
        ## View settings
        ##
        self.viewType = 'pca'
        self.viewBackground = False
        self.viewReferenceEndpoint = None
        self.viewReferenceVersion = 0


    def licenseTesting (self):
        
        if self.norm and self.normNeutr and self.normNeutrMethod=='moka':
            if not os.path.isfile (self.mokaPath+'/license.txt'):
                print 'No suitable license found for Moka software. Aborting'
                return (False)
                
        if self.MD == 'pentacle':
            
            if not os.path.isfile (self.pentaclePath+'/license.txt'):
                print 'No suitable license found for Pentacle software. Aborting'
                return (False)

        return (True)

        


##################################################################
##    SHARED (PREDICTION & BUILD) METHODS
##################################################################    

    def loadData (self):
        """Gets all model settings stored for the last model and compares with current settings.
           In case of any disagreement, a False is returned.

           Please note that this method does not intend to set up the model settings (__init__)
        """

        if self.confidential:
            return False
        
        if not os.path.isfile (self.vpath+'/tdata.pkl'):
            return False

        try:
            f = open (self.vpath+'/tdata.pkl','rb')
        except:
            return False
        
        norm = pickle.load(f)
        if norm != self.norm:
            return False

        if norm:
            normStand = pickle.load(f)
            if normStand != self.normStand:
                return False
            
            normNeutr = pickle.load(f)
            if normNeutr != self.normNeutr:
                return False

            if normNeutr:
                normNeutrMethod = pickle.load(f)
                if normNeutrMethod != self.normNeutrMethod:
                    return False

                normNeutr_pH = pickle.load(f)
                if normNeutr_pH != self.normNeutr_pH:
                    return False
                
            norm3D = pickle.load(f)
            if norm3D != self.norm3D:
                return False
        
        MD = pickle.load(f)
        if MD != self.MD:
            return False
        
        if 'pentacle' in MD:
            pentacleProbes = pickle.load(f)
            if pentacleProbes != self.pentacleProbes:
                return False

            pentacleOthers = pickle.load(f)
            if pentacleOthers != self.pentacleOthers:
                return False

        elif 'padel' in MD:
            padelMD = pickle.load(f)
            if padelMD != self.padelMD:
                return False
        
        self.tdata = pickle.load(f)
        f.close()

        return True

    def saveData (self):
        """Saves all model settings in file tdata.pkl, so they can be compared with actual model settings in future runs

           Please note that this method does not intend to save information for setting up the model. This is carried out
           by the __init__ method
        """

        if self.confidential:
            return
        
        try:
            f = open (self.vpath+'/tdata.pkl','wb')
        except:
            return
                
        pickle.dump(self.norm, f)
        if self.norm:
            pickle.dump(self.normStand, f)
            pickle.dump(self.normNeutr, f)
            if self.normNeutr:
                pickle.dump(self.normNeutrMethod, f)
                pickle.dump(self.normNeutr_pH, f)
            pickle.dump(self.norm3D, f)
            
        pickle.dump(self.MD,f)
        if 'pentacle' in self.MD:
            pickle.dump(self.pentacleProbes,f)
            pickle.dump(self.pentacleOthers,f)
        elif 'padel' in self.MD:
            pickle.dump(self.padelMD,f)
            
        pickle.dump(self.tdata, f)

        # remove variables that might not be applicable any longer, like FFD excluded variables
        removefile (self.vpath+'/ffdexcluded.pkl')

        f.close()

    
    def savePropertyData (self):
        """Saves visualization matrix of property data in file pdata.pkl

           Please note that this method does not intend to save information for setting up the model. This is carried out
           by the __init__ method
        """

        if self.confidential:
            return
        
        try:
            f = open (self.vpath+'/pdata.pkl','wb')
        except:
            return
            
        pickle.dump(self.tdata, f)
        f.close ()

    def loadPropertyData (self):
        """Gets all model settings stored for the last model and compares with current settings.
           In case of any disagreement, a False is returned.

           Please note that this method does not intend to set up the model settings (__init__)
        """

        if self.confidential:
            return False
        
        if not os.path.isfile (self.vpath+'/pdata.pkl'):
            return False

        try:
            f = open (self.vpath+'/pdata.pkl','rb')
        except:
            return False
        
        self.tdata = pickle.load(f)
        f.close()

        return True

    def loadVisualData (self):
        """Gets all model settings stored for the last model and compares with current settings.
           In case of any disagreement, a False is returned.

           Please note that this method does not intend to set up the model settings (__init__)
        """

        if self.confidential:
            return False
        
        if not os.path.isfile (self.vpath+'/tdata.pkl'):
            return False

        try:
            f = open (self.vpath+'/tdata.pkl','rb')
        except:
            return False
        
        norm = pickle.load(f)

        if norm:
            normStand = pickle.load(f)
            normNeutr = pickle.load(f)
            if normNeutr:
                normNeutrMethod = pickle.load(f)
                normNeutr_pH = pickle.load(f)
            norm3D = pickle.load(f)
        
        MD = pickle.load(f)
        
        if 'pentacle' in MD:
            pentacleProbes = pickle.load(f)
            pentacleOthers = pickle.load(f)
        elif 'padel' in MD:
            padelMD = pickle.load(f)
        
        self.tdata = pickle.load(f)
        f.close()

        return True        

        
    def loadSeriesInfo (self):
        """Gets information about the series used to build the model (model.infoSeries) stored in file info.pkl
           This information is used only for informative purposes and does not change the model properties
        """
        try:
            modelInfo = open (self.vpath+'/info.pkl','rb')
        except:
            return (False)
        infoID = pickle.load(modelInfo)
        self.infoSeries = pickle.load(modelInfo)
        modelInfo.close()
        return (True)

        
    def adjustPentacle (self, row, nprobes, Bcol):
        """Adjust the row of GRIND descriptors in "row" to Bcol size, asuming that we have nprobes
           GRID probes, applying a procrustean transform for each block (correlogram) 

           Both the row and the returning array are NumPy float64 arrays 
        """
        
        Acol = len(row)        # orignal num of columns
        blocks = [0,1,3,6,10]
        nblocks = blocks[nprobes]
        
        if Acol == Bcol:
            return row

        deltaA = Acol/nblocks  # num col in original
        deltaB = Bcol/nblocks  # num col in new

        nn = np.empty(0,dtype='float64')   
        if Acol > Bcol:
            for i in range (nblocks):
                start = i*deltaA
                nn=np.hstack((nn,row[start:(start+deltaB)]))
        else:
            zero = np.zeros(deltaB-deltaA,dtype='float64')
            for i in range (nblocks):
                start = i*deltaA
                nn=np.hstack((nn,row[start:(start+deltaA)]))
                nn=np.hstack((nn,zero))
                
        return nn

    def computeMDAdriana (self, mol, clean=True):

        molr = randomName(20)+'.csv'
        
        call = [self.adrianaPath+'adriana',
                '-i', mol,
                '-o', molr, 
                self.adrianaPath+'adrianaFull.prj']   #TODO use diverse projects and select here  

        stdoutf = open (os.devnull, 'w')
        stderrf = open (os.devnull, 'w')

        try:
            retcode = subprocess.call(call,stdout=stdoutf,stderr=stderrf)
        except:
            stdoutf.close()
            stderrf.close()
            return (False, 'AdrianaCode execution error' )

        stdoutf.close()
        stderrf.close()
            
        try:
            fpr = open (molr)
        except:
            return (False, 'Adriana results not found')
        
        line = fpr.readline()
        line = fpr.readline()
        fpr.close()

        md = np.genfromtxt(StringIO(line.partition(',')[2]),delimiter=',')
        md = np.nan_to_num(md)
        
        if clean:
            removefile (molr)
        
        return (True,md)
                
    def computeMDPentacle (self, mol, clean=True):
        """ Computes Pentacle Molecular Descriptors for compound "mol"

            It returms a tuple that contains
            1) True or False, indicating the success of the computation
            2) A vector of floats (if True) with the GRIND descriptors
               An Error message (if False)
        """
        
        molr = randomName(20)
        
        t = open ('template-md','w')
        t.write ('name '+molr+'\n')
        t.write ('input_file '+mol+' sdf\n')
        t.write ('mif_computation grid\n')
        t.write ('mif_discretization amanda\n')
        t.write ('mif_encoding  macc2\n')
        for key in self.pentacleOthers:
            t.write (key+'\n')
        for probe in self.pentacleProbes:
            t.write ('probe '+probe+'\n')
        t.write ('dynamic yes\n')
        t.write ('export_data csv\n')
        
        t.close()       
        
        call = [self.pentaclePath+'pentacle',
                '-c','template-md']  

        stdoutf = open ('stdout.txt','w')
        stderrf = open (os.devnull, 'w')

        try:
            retcode = subprocess.call(call,stdout=stdoutf,stderr=stderrf)
        except:
            removefile ( '/var/tmp/'+molr )
            removefile ( '/var/tmp/'+molr+'.ppf' )
            stdoutf.close()
            stderrf.close()
            return (False, 'Pentacle execution error' )

        stdoutf.close()
        stderrf.close()

        removefile ( '/var/tmp/'+molr )
        removefile ( '/var/tmp/'+molr+'.ppf' )

        try:
            stdoutf = open('stdout.txt')
        except:
            return (False, 'Pentacle std output not found')

        for line in stdoutf:
            if 'Error' in line :
                stdoutf.close()
                return (False, line)
        stdoutf.close()
            
        try:
            fpr = open (molr+'.csv')
        except:
            return (False, 'Pentacle results not found')
        
        line = fpr.readline()
        fpr.close()
        
        line = line.partition(',')[2] # removes mol name

        if len(line):
            lfile = StringIO(line)
            md = np.loadtxt(lfile,delimiter=',')
            lfile.close()
        else:
            return (False,'error in Pentacle')

        if clean:
            removefile (molr+'.csv')
            removefile ('template-md')
            removefile ('stdout.txt')
        
        return (True,md)


    def computeMDPadelws (self, mol, clean=False):
        """ Computes PaDEL Molecular Descriptors for compound "mol"

            In this implementation we use PaDEL making calls to a web service
            in the background to avoid startig up the Java VM again and again
            
            It returms a tuple that contains
            1) True or False, indicating the success of the computation
            2) A vector of floats (if True) with the PaDEL descriptors
               An Error message (if False)
        """
        try:
            shutil.rmtree('padel')
        except:
            pass
        
        os.mkdir ('padel')
        shutil.copy (mol,'padel')

        homepath = os.getcwd()
        call = ['-dir', homepath+'/padel',
                '-file',homepath+'/padel.txt']

        for key in self.padelMD:
            call.append (key)

        if self.padelMaxRuntime:
            call.append ('-maxruntime')
            call.append (str(self.padelMaxRuntime))

        if self.padelDescriptor:
            dname,fname = os.path.split(self.padelDescriptor)
            dfile = self.vpath+'/'+fname
            if not os.path.isfile(dfile):
                shutil.copy(self.padelDescriptor,dfile)
            call.append ('-descriptortypes')
            call.append (dfile)
    
        params = "|".join(call)
        try:
            url = self.padelURL+params
            req  = urllib2.Request(url)
            resp = urllib2.urlopen(req)
            the_page = resp.read() 
            #print the_page
            retcode = 0
        except urllib2.HTTPError as e:
            return (False, 'PaDEL execution HTTPError' )
        except urllib2.URLError as e:
            return (False, 'PaDEL execution URLError' )
        except:
            return (True, 'PaDEL execution error' )

        try:
            fpr = open ('padel.txt','r')
        except:
            return (False, 'PaDEL results not found')
        
        line = fpr.readline()
        line = fpr.readline()
        fpr.close()
        
        md = np.genfromtxt(StringIO(line.partition(',')[2]),delimiter=',')
        md = np.nan_to_num(md)
        
        if clean:
            try:
                shutil.rmtree('padel')
                removefile ('padel.txt')
            except:
                pass
        
        return (True,md)

    def computeMDPadelcl (self, mol, clean=False):
        """ Computes PaDEL Molecular Descriptors for compound "mol"

            In this implementation we use PaDEL calling the Java VM. This is
            less efficient than computeMDPadelws for large collections of compounds
            
            It returms a tuple that contains
            1) True or False, indicating the success of the computation
            2) A vector of floats (if True) with the PaDEL descriptors
               An Error message (if False)
        """
        try:
            shutil.rmtree('padel')
        except:
            pass
        
        os.mkdir ('padel')
        shutil.copy (mol,'padel')

        call = [self.javaPath+'bin/java','-Djava.awt.headless=true','-jar',
                self.padelPath+'PaDEL-Descriptor.jar',
                '-dir','./padel',
                '-file','padel.txt']

        for key in self.padelMD:
            call.append (key)

        stdoutf = open (os.devnull, 'w')
        stderrf = open (os.devnull, 'w')
        try:
            retcode = subprocess.call(call,stdout=stdoutf,stderr=stderrf)
        except:
            return (True, 'PaDEL execution error' )
        finally:
            stdoutf.close()
            stderrf.close()

        try:
            fpr = open ('padel.txt','r')
        except:
            return (False, 'PaDEL results not found')
        
        line = fpr.readline()
        line = fpr.readline()
        fpr.close()
        
        md = np.genfromtxt(StringIO(line.partition(',')[2]),delimiter=',')
        md = np.nan_to_num(md)
        
        if clean:
            try:
                shutil.rmtree('padel')
                removefile ('padel.txt')
            except:
                pass
        
        return (True,md)

    def computeMDlogpmw (self, mol, clean=True):

        md = np.zeros (2, dtype='float64')
        
        try:
            suppl = Chem.SDMolSupplier(mol)
            mi = suppl.next()

            if mi is None:
                return (False, 'wrong input format')

            md[0]=Descriptors.MolLogP(mi)
            md[1]=Descriptors.MolWt(mi)

        except:
            return (False, 'error computing logP and MW')

        return (True, md)
    
    def computeMD (self, mol, clean=True):
        """ Computes Molecular Descriptors for compound "mol"

            This is just a wrapper method that calls the appropriate
            computeMDxxxxx metho, depending on the MD used by this model
        """
        
        if 'pentacle' in self.MD:
            success, md = self.computeMDPentacle (mol, clean)
        elif 'padel' in self.MD:
            success, md = self.computeMDPadelws (mol, clean)
        elif 'adriana' in self.MD:
            success, md = self.computeMDAdriana (mol, clean)

##        print md
        
        return (success, md)
    
##################################################################
##    NORMALIZE METHODS
##################################################################    

    def getMolName (self, molFile):
        """ Find the name in the SDFile

            We first used the field defined by SDFileName parameter, then the firt line and last a seq number
        """

        try:
            suppl=Chem.SDMolSupplier(molFile)
            mi = suppl.next()

            if not mi:
                return (False, 'wrong input format')

        except:
            return (False, 'wrong input format')
        
        name = ''
        if self.SDFileName:
            if mi.HasProp (self.SDFileName):
                name = mi.GetProp(self.SDFileName)
        if not name:
            name = mi.GetProp('_Name')
        if not name:
            name = molFile[:-4]

        return (True, name)
    
    def standardize (self, moli, clean=True):
        """Applies a structure normalization protocol provided by Francis Atkinson (EBI)

           The name of the output molecules is built as a+'original name'

           Returns a tuple containing:
           1) True/False: depending on the result of the method
           2) (if True ) The name of the output molecule
              (if False) The error message
        """
        molo = 'a'+moli

        suppl=Chem.SDMolSupplier(moli)
        m = suppl.next()
            
        try:
            parent = standardise.apply(Chem.MolToMolBlock(m))
        except standardise.StandardiseException as e:
            return (False, e.name)
            
        fo = open (molo,'w')
        fo.write(parent)

        if self.SDFileActivity:        
            if m.HasProp(self.SDFileActivity):
                activity = m.GetProp(self.SDFileActivity)
                fo.write('>  <'+self.SDFileActivity+'>\n'+activity+'\n\n$$$$')

        if m.HasProp('activity'):
            activity = m.GetProp('activity')
            fo.write('>  <activity>\n'+activity+'\n\n$$$$')
            
        fo.close()

        if clean:
            removefile (moli)

        return (True,molo)

    def protonate (self, moli, pH, clean=True):
        """Adjusts the ionization state of the molecule "moli" 

           In this implementation, it uses blabber_sd from Molecular Discovery
           The result is a tuple containing:
           1) True/False: describes the success of the protonation for this compound
           2) (if True ) The name of the protonated molecules and its formal charge
              (if False) The error message
        """

        molo = 'b'+moli

        stderrf = open (os.devnull, 'w')
        stdoutf = open (os.devnull, 'w')     

        call = [self.mokaPath+'blabber_sd', moli,
                '-p',  str(pH),
                '-o',  molo]

        try:
            retcode = subprocess.call(call,stdout=stdoutf, stderr=stderrf)
        except:
            return (False, 'Blabber execution error', 0.0)
        
        stdoutf.close()
        stderrf.close()

        # for the newer versions use !=0, for the old one use ==0
        if retcode == 0:
            return (False, 'Blabber execution error', 0.0)

        try:
            finp = open (molo)
        except:
            return (False, 'Blabber output not found', 0.0)

        charge = 0
        for line in finp:
            if line.startswith ('M  CHG'):
                items = line.split()
                if int(items[2]):
                    for c in range (4,len(items),2): charge+=int(items[c])
                break   
        finp.close()

        if clean:
            removefile (moli)
        
        return (True, molo, charge)


    def convert3D (self, moli, clean=True):
        """Converts the 2D structure of the molecule "moli" to 3D

           In this implementation, it uses CORINA from Molecular Networks
           The result is a tuple containing:
           1) suucTrue/False: describes the success of the 3D conversion for this compound
           2) (if True ) The name of the 3D molecule
              (if False) The error mesage
        """

        molo = 'c'+moli
        
        stderrf = open (os.devnull, 'w')
        stdoutf = open (os.devnull, 'w')
        
        call = [self.corinaPath+'corina',
                '-dwh','-dori',
                '-ttracefile=corina.trc',
                '-it=sdf', moli,
                '-ot=sdf', molo]

        try:
            retcode = subprocess.call(call, stdout=stdoutf, stderr=stderrf)
        except:
            return (False, 'Corina execution error')

        stdoutf.close()
        stderrf.close()
            
        if retcode != 0:
            return (False, 'Corina execution error')

        if not os.path.exists(molo):
            return (False, 'Corina output not found')

        if clean:
            removefile(moli)
            removefile('corina.trc')
            
        return (True,molo)


    def checkIdentity (self, mol, ypcutoff=0.5):
        """ Checks if the compound "mol" is part of the training set 

            We used InChiKeys without the last three chars (ionization state) to make the comparison

            This version uses RDKit
        """
        
        # this dissables warnings
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.ERROR)

        try:
            suppl = Chem.SDMolSupplier(mol)
            mi = suppl.next()
        except:
            return (False, 'wrong input format')

        if mi is None:
            return (False, 'wrong input format')

        try:
            ik = Chem.InchiToInchiKey(Chem.MolToInchi(mi))
        except:
            return (False, 'failed to obtain InChiKey')

        
        ik = ik[:-3] # remove the right-most part expressing ionization
        
        for l in self.tdata:   # the InChi is the element 1 of the tuple and the Activity is the element 4
            
            if ik in l[1]:
            
                yp = float (l[4])
                
                if self.quantitative:
                    return (True, yp)
                if (yp < ypcutoff):
                    return (True, 'negative')
                else:
                    return (True, 'positive')
        
        return (False, mol)

    def saveNormalizedMol (self, mol):
        """Appends normalized moleculed to a SDFile containing normalized structures for all the compounds in
           the training series

           return the possition, so every molecule can be easily extracted from the file
        """

        if self.confidential:
            return
        
        fi = open (mol)
        fo = open (self.vpath+'/tstruct.sdf','a')
        alpha = fo.tell()
        
        for line in fi:
            fo.write(line)
        fi.close()
        fo.close()

        return alpha

    def normalize (self, mol):
        """Preprocesses the molecule "mol" by running a workflow that:

        - Normalizes the 2D structure (DUMMY)
        - Adjusts the ionization state 
        - Converts the structure to 3D

        The result is a tuple containing:
        1) True/False: describes the success of the normalization
        2) (if True ) The name of the normalized molecule and its formal charge
           (if False) The error mesage
        """

        charge = 0.0 #fallback

        success, result = self.getMolName (mol)
        if not success: return (False, result)

        molName = result

        if not self.norm:
            return (True, (mol, molName, charge))
        
        if self.normStand:
            success, resulta = self.standardize (mol)
            if not success: return (False, resulta)
        else:
            resulta = mol

        if self.normNeutr:
            success, resultb, charge = self.protonate (resulta, self.normNeutr_pH)
            if not success: return (False, resultb)
        else:
            resultb = resulta

        if self.norm3D:
            success, resultc = self.convert3D (resultb)
            if not success: return (False, resultc)
        else:
            resultc = resultb
        
        return (True,(resultc, molName, charge))


##################################################################
##    PREDICT METHODS
##################################################################    

    def computePredictionOther (self, md, charge):
        # empty method to be overriden
        return (False, 'not implemented')

    
    def computePredictionPLS (self, md, charge):
        """ Computes the prediction for compound "mol"

            This function makes use of the molecular descriptor vector (md) to project the compound in the model
            The model has been loaded previously as an R object
        """
        
        model = pls()
        if not self.confidential:
            model.loadModel(self.vpath+'/modelPLS.npy')
        else:
            model.loadDistiled(self.vpath+'/distiledPLS.txt')
            
        if 'pentacle' in self.MD:
            md = self.adjustPentacle(md,len(self.pentacleProbes),model.nvarx)
            
        success, result = model.project(md,self.modelLV)

        if not success:
            return (success, result)

        yp = result[0][-1]

        if self.quantitative:
            return (True, yp)

        # qualitative
        if len (model.cutoff) == 0:
            return (False, 'cutoff not defined')
        
        if yp < model.cutoff[-1]: # use last cutoff
            return (True, 'negative')
        else:
            return (True, 'positive')


    def computePrediction (self, md, charge):
        """ Computes the prediction for compound "mol"

            This function makes use of the molecular descriptor vector (md) to project the compound in the model
            The model has been loaded previously as an R object
        """      
        success = False
        result = 0.0
        
        if self.model == 'pls':
            success, result = self.computePredictionPLS (md, charge)
        else :
            success, result = self.computePredictionOther (md, charge)

        return (success, result)



    def computeAD (self, md, pr, detail):
        """Carries out a protocol for determining how far is the query compound from the model

           Provisionally, implements a temporary version of the ADAN method. Note that this requires that
           ADAN was run on the model, building an ad hoc PLS model (adan.npy) and extracting information
           (tscores.npy) 

           Returns a tuple that contains
           1) True or False, indicating the success of the computation
           2) (if True ) a number between 0 and 6 with the number of criteria broken
              (if False) an error message 
           
        """
        f = file (self.vpath+'/tscores.npy','rb')
        nlv = np.load(f)
        p95dcentx = np.load(f)
        p95dclosx = np.load(f)
        p95dmodx = np.load(f)
        p95dcenty = np.load(f)
        p95dclosy = np.load(f)
        p95dpredy = np.load(f)
        centx = np.load(f)
        centy = np.load(f)
        T = np.load(f)
        Y = np.load(f)
        squareErr = np.load(f)
        f.close()

        model = pls ()
        model.loadModel(self.vpath+'/adan.npy')

        if 'pentacle' in self.MD:
            md = self.adjustPentacle(md,len(self.pentacleProbes),model.nvarx)
        
        success, result = model.project(md,model.Am)

        y = pr
        t = result[1]
        d = result[2]

        AD=dict()
        # compute distance t to centroid x (A)
        dcentx = np.sqrt(np.sum(np.square(centx-t)))
        AD['dcentx']=dcentx>p95dcentx

        # compute min distance to T (B)
        dclosx = 1e20
        dclosi = 0
        for i in range (model.nobj):
            dtemp = np.sqrt(np.sum(np.square(T[i,:]-t)))
            if dtemp < dclosx :
                dclosx = dtemp
                dclosi = i
        
        AD['dclosx']=dclosx>p95dclosx

        # compute distance to modx (C)
        AD['dmodx']=d[-1]>p95dmodx

        # compute distance to centroid y, only for quantitative (D)
        if self.quantitative:
            dcenty = np.abs(y-centy)
            AD['dcenty']=dcenty>p95dcenty
        else:
            AD['dcenty']=False

        # compute Y of the closer compound (E)
        if self.quantitative:
            dclosy = np.abs(y-Y[dclosi])
            AD['dclosy']=dclosy>p95dclosy
        else:
            AD['dclosy']=False

        # compute SDEP of 5% closer neighbours (F) and percentil 95
        if self.quantitative:
            p5 = int(np.rint(0.05*model.nobj))
            if p5 < 1 : p5 = int(1)
        
            closerDis = np.empty(p5,dtype=np.float64)
            closerErr = np.empty(p5,dtype=np.float64)
            
            for j in range (p5):
                closerDis[j]=1e20
                
            for i in range (model.nobj):
                dtemp = np.sum(np.square(T[i,:]-t))    
                if dtemp < np.amax(closerDis) :
                    imindis=np.argmax(closerDis)
                    closerDis[imindis] = dtemp
                    closerErr[imindis] = squareErr[i]
                        
            dpredy=np.sqrt(np.sum(closerErr)/p5)
            AD['dpredy']=dpredy>p95dpredy
        else:
            AD['dpredy']=False
        
##        print '[',
##        if AD['dcentx'] : print '1 ',
##        else : print '0 ',
##        if AD['dclosx'] : print '1 ', 
##        else : print '0 ',
##        if AD['dmodx' ] : print '1 ',
##        else : print '0 ',
##        if AD['dcenty'] : print '1 ',
##        else : print '0 ',
##        if AD['dclosy'] : print '1 ',
##        else : print '0 ',
##        if AD['dpredy'] : print '1 ',
##        else : print '0 ',
##        print '] %d' % sum (AD.values())
        
##        print "DCENTX %6.3f (%6.3f)\n" % (dcentx,p95dcentx),
##        print "DCLOSX %6.3f (%6.3f)\n" % (dclosx,p95dclosx),
##        print "DCMODX %6.3f (%6.3f)\n" % (d[-1],p95dmodx),
##        
##        if self.quantitative:
##            print "DCENTY %6.3f (%6.3f)\n" % (dcenty,p95dcenty),
##            print "DCLOSY %6.3f (%6.3f)\n" % (dclosy,p95dclosy),
##            print "DPREDY %6.3f (%6.3f)\n" % (dpredy,p95dpredy)

        return (True,sum(AD.values()))

    def computeCI (self, ad):
        """Calculates a Reliability Index for the given prediction

           Provisionally it returns a tuple that contains
           1) True or False, indicating the success of the computation
           2) (if True ) - the model SDEP when the #broken AD criteria is 0 or 1
                         - twice the model SDEP when the #broken AD criteria is 2 or 3
                         - 0.0 otherwyse (this must be interpreted as a "model NA")
              (if False) An error message
        """

        # this only works for quantitative models, the qualitative version
        # has not been implemented

        if not self.quantitative:
            return (False, 0)
        
        ci=0.0
        
        if ad<4:
            model = pls ()
            model.loadModel(self.vpath+'/modelPLS.npy')
            ci = 1.96*model.SDEP[model.Av-1]
            if ad<2:
                pass        # 0 or 1 criteria broken
            else:
                ci *= 2.0   # 2 or 3 criteria broken
        
        return (True,ci)
  

    def predict (self, molFile, molName, molCharge, detail, clean=True):
        """Runs the prediction protocol

        """
        # default return values
        molPR=molCI=molAD=(False,0.0)
        
        if self.identity:
            success, result = self.checkIdentity (molFile)
            
            if success:
                molPR = (success, result)    # the value of the training set
                molAD = (True, 0)            # no ADAN rules broken
                molCI = (True, 0.0)          # CI is 0.0 wide
                return (molPR,molAD,molCI)
     
        success, molMD = self.computeMD (molFile)
        if not success: return (molPR,molAD,molCI)

        # MD are passed as copy because the scaling changes them and they
        # need to be reused for computing the AD
        
        success, pr = self.computePrediction (molMD.copy(),molCharge)
        molPR = (success, pr)
        if not success: return (molPR,molAD,molCI)
        
        if not self.confidential:
            success, ad = self.computeAD (molMD.copy(), pr, detail)
            molAD = (success, ad)
            if not success: return (molPR, molAD ,molCI)

            success, ci = self.computeCI (ad)
            molCI = (success, ci)

        if clean: removefile(molFile)
        
        return (molPR,molAD,molCI)

##################################################################
##    EXTRACT METHODS
##################################################################          
 
    def getInChi (self, mi):
        """ Computes the InChiKey for the compound 

            We used InChiKeys without the last three chars (ionization state) to make the comparison
        """

        lg = RDLogger.logger()
        lg.setLevel(RDLogger.ERROR)
        
        try:
            ik = Chem.InchiToInchiKey(Chem.MolToInchi(mi))
        except:
            return (False, 'Failed to obtain InChiKey')
        return (True,ik[:-3])


    def getBio (self, mi):
        """ Extracts the value of the experimental biological property from the compound 

            Such value must be identified by the tag <activity>
        """

        bio = None
        if self.SDFileActivity:
            if mi.HasProp(self.SDFileActivity):
                bio = mi.GetProp(self.SDFileActivity)
                
        if not bio:
            if mi.HasProp('activity'):
                bio = mi.GetProp('activity')

        if not bio:
            return (False, 'Biological activity not found')
        
        try:
            nbio = float (bio)
        except:
            return (False, 'Activity cannot be converted to numerical format')
        return (True, nbio)



    def extract (self, molFile, molName, molCharge, molPos, clean=True):
        """Process the compound "mol" for obtaining
           1) InChiKey (string)
           2) Molecular Descriptors (NumPy float64 array)
           3) Biological Activity (float)

           Returns a Tuple as (True,('InChi',array[0.0,0.0, 0.0],4.56))
        """

        ##  data[0] molname
        ##  data[1] InChi
        ##  data[2] MD
        ##  data[3] charge            +++
        ##  data[4] activity
        ##  data[5] index in tstruct  +++

        if not self.buildable:
            return (False, 'this model cannot by built automatically')

        molInChi=''
        molMD=[]
        molActivity=0.0
        
        suppl = Chem.SDMolSupplier(molFile)
        mol = suppl.next()
        
        success, molInChi = self.getInChi(mol)
        if not success:
            return (False,(molName,molInChi,molMD,molCharge,molActivity,molPos))
        
        success, molActivity = self.getBio(mol)
        if not success:
            return (False,(molName,molInChi,molMD,molCharge,molActivity,molPos))

        success, molMD = self.computeMD(molFile)
        if not success:
            return (False,(molName,molInChi,molMD,molCharge,molActivity,molPos))

        self.tdata.append( (molName,molInChi,molMD,molCharge,molActivity,molPos) )
        
        if clean:
            removefile (molFile)
                            
        return (True,'extraction OK')


    def extractView (self, molFile, molName, molCharge, clean=True):
        """Process the compound "mol" for obtaining
           2) Molecular Descriptors (NumPy float64 array)

           Returns a Tuple as (True,('InChi',array[0.0,0.0, 0.0],4.56))
        """

        molInChi=''
        molMD=[]
        molActivity=0.0
        
        suppl = Chem.SDMolSupplier(molFile)

        if self.viewType == 'property':
            success, molMD = self.computeMDlogpmw(molFile)
        else :
            success, molMD = self.computeMD(molFile)
            
        if not success:
            return (False,(molName,molInChi,molMD,molCharge,molActivity, 0))   # last argument (molPos) replaced by 0

        self.tdata.append( (molName,molInChi,molMD,molCharge,molActivity, 0) ) # last argument (molPos) replaced by 0 
        
        if clean:
            removefile (molFile)
                            
        return (True,'extraction OK')
    

##################################################################
##    BUILD METHODS
##################################################################    

    def setSeries (self, molecules, numMol):
        self.infoSeries = []
        self.infoSeries.append ( ('series',molecules) )
        self.infoSeries.append ( ('nmol',numMol) )
        
    def getMatrices (self):
        """ Returns NumPy X and Y matrices extracted from tdata. In case of Pentacle MD, it also adjusts the X vectors
        """
        ncol = 0
        xx = []
        yy = []
        
        # obtain X and Y from tuple elements 2 (MD) and 4 (Activity)
        for i in self.tdata:
            if len(i[2])>ncol: ncol = len(i[2])
            xx.append(i[2])
            yy.append(i[4])

        nrow = len (xx)
        
        Y = np.array (yy,dtype=np.float64)
        X = np.empty ((nrow,ncol),dtype=np.float64)
      
        i=0
        for row in xx:
            if 'pentacle' in self.MD:
                row=self.adjustPentacle(row,len(self.pentacleProbes),ncol)
            X[i,:]=np.array(row,dtype=np.float64)
            i+=1

        return X, Y
        
    def buildPLS (self, X, Y):
        """ Builds a PLS model using the internal PLS implementation. The number of LV and the autoscaling
            are model settings defined in __init__

            returns the PLS model object
        """
        
        model = pls ()

        if self.selVar:
            iRuns = 0

            resSet = []
            
            if (os.path.isfile(self.vpath+'/ffdexcluded.pkl')):
                resfile = open (self.vpath+'/ffdexcluded.pkl', 'rb')

                oslv = pickle.load(resfile) # checking if the old FFD used the same number of LV
                resn = pickle.load(resfile)

                for i in range (resn):
                    resSet.append (pickle.load(resfile))
                resfile.close()
                
                if resn >= self.selVarRun:
                    res = resSet[self.selVarRun-1]
                    iRuns = self.selVarRun
                else:
                    res = resSet[-1]
                    iRuns = resn
                        
                nobj,nvarx = np.shape (X)

                # make sure that the old FFD mask has the same nvarx and was generated using the same
                # number of LV
                if (len(res) != nvarx) or (oslv != self.selVarLV):
                    
                    iRuns = 0
                    removefile (self.vpath+'/ffdexcluded.pkl')
                    
                else :

                    X = model.excludeVar (X, res)    
                    self.selVarMask = res
                    
                    print '\napplying previous FFD...'
                    
                    model.build (X,Y,self.modelLV,autoscale=self.modelAutoscaling)
                    
                    model.validateLOO(self.modelLV)
                    
                    for i in range (self.modelLV):
                       print 'LV%2d R2:%5.3f Q2:%5.3f SDEP:%7.3f' % \
                            (i+1,model.SSYac[i],model.Q2[i],model.SDEP[i])
                                   
            while True:
                if iRuns >= self.selVarRun :
                    break
                
                print 'FFD var selection... (please be patient)'

                res, nexcluded = model.varSelectionFFD (X,Y,self.selVarLV,self.modelAutoscaling)
                X = model.excludeVar (X, res)
                self.selVarMask = res
                
                iRuns += 1
                resSet.append(res)

                print '\n',nexcluded,' variables excluded'
                
                model.build (X,Y,self.modelLV,autoscale=self.modelAutoscaling)

                model.validateLOO(self.modelLV)

                for i in range (self.modelLV):
                   print 'LV%2d R2:%5.3f Q2:%5.3f SDEP:%7.3f' % \
                         (i+1,model.SSYac[i],model.Q2[i],model.SDEP[i])

                if not nexcluded :
                    break
           
            # refresh the FFD mask in any case
            
            resfile = open (self.vpath+'/ffdexcluded.pkl','wb')
            pickle.dump (self.selVarLV, resfile)
            pickle.dump (len(resSet),resfile)
            for ri in resSet:
                pickle.dump(ri,resfile)
            resfile.close()

        else:
            model.build (X,Y,self.modelLV,autoscale=self.modelAutoscaling)

        self.infoModel = []
        if self.quantitative:
            self.infoModel.append( ('model','PLS-R  (NIPALS)') )
        else:
            self.infoModel.append( ('model','PLS-DA (NIPALS)') )
        self.infoModel.append( ('LV', self.modelLV ))
        return model
    
    def diagnosePLS_R (self, model):
        """ Runs CV diagnostic on the PLS-R model provided as argument

            The results are:
            - printed (R2, Q2, SDEP)
            - stored in infoResult (R2, Q2, SDEP)
            - images with the recalculated and predicted results for all the objects and LV
            - text files with the recalculated and predicted results for all the objects and LV
        """

        print 'cross-validating...'

        if self.selVar :
            model.X = model.excludeVar (model.X, self.selVarMask)
            #model.build (X,Y,self.modelLV,autoscale=self.modelAutoscaling)
            
        yp = model.validateLOO(self.modelLV)
        for i in range (self.modelLV):
            print 'LV%2d R2:%5.3f Q2:%5.3f SDEP:%7.3f' % \
                  (i+1,model.SSYac[i],model.Q2[i],model.SDEP[i])
            
        self.infoResult = []    
        self.infoResult.append( ('nobj',model.nobj) )
        self.infoResult.append( ('R2','%5.3f' % model.SSYac[self.modelLV-1]) )
        self.infoResult.append( ('Q2','%5.3f' % model.Q2[self.modelLV-1]) )
        self.infoResult.append( ('SDEP','%5.3f' % model.SDEP[self.modelLV-1]) )
        
        yr = model.recalculate()

        # generate rec vs experimental and pred vs experimental for all model dimensions
        for i in range(self.modelLV):
            nvar = str(i+1)
            
            fig1=plt.figure()
            plt.xlabel('experimental y')
            plt.ylabel('predicted LV'+nvar)
            plt.title('Predicted')
            plt.plot(yp[:,0],yp[:,i+1],"ro")
            fig1.savefig("pls-predicted-LV"+nvar+".png", format='png')
            
            fig2=plt.figure()
            plt.xlabel('experimental y')
            plt.ylabel('recalculated LV'+nvar)
            plt.title('Recalculated')
            plt.plot(yr[:,0],yr[:,i+1],"ro")
            fig2.savefig("pls-recalculated-LV"+nvar+".png", format='png')
        #plt.show()

        # write a file with experimental Y (yp[0]) vs LOO predicted Y 
        fp=open ('pls-predicted.txt','w')
        fr=open ('pls-recalculated.txt','w')

        # write simple header
        # indicate that the first column corresponds with experimental values
        header = 'Name Yexp '
        for i in range (self.modelLV):
            header += 'Y-LV%d '%(i+1)
        header += '\n'

        fp.write (header)
        fr.write (header)
        
        for i in range (model.nobj):
            fp.write(self.tdata[i][0]+' ')
            fr.write(self.tdata[i][0]+' ')
            for j in range (self.modelLV+1):
                fp.write('%.3f '% yp[i][j])
                fr.write('%.3f '% yr[i][j])
            fp.write('\n')
            fr.write('\n')
        fp.close()
        fr.close()

        return (yp[:,-1])

    def diagnosePLS_DA (self, model):
        """ Runs CV diagnostic on the PLS-DA model provided as argument. In case no cutoff is provided, it also
            computes an appropriate cutoff that balances model sensibility and specificity

            The results are:
            - printed (TP, TN, FP, FN, spec, sens, MCC)
            - stored in infoResult (TP, TN, FP, FN, spec, sens, MCC)
        """

        if 'auto' == self.modelCutoff:
            model.calcOptCutoff ()
        else:
            model.calcConfussion(self.modelCutoff)

        for i in range (self.modelLV):
            sens = sensitivity(model.TP[i],model.FN[i])
            spec = specificity(model.TN[i],model.FP[i])
            mcc  = MCC(model.TP[i],model.TN[i],model.FP[i],model.FN[i])

            print "rec  LV:%-2d cutoff:%4.2f TP:%3d TN:%3d FP:%3d FN:%3d spec:%5.3f sens:%5.3f MCC:%5.3f" % (i+1,
                    model.cutoff[i], model.TP[i], model.TN[i], model.FP[i], model.FN[i], spec, sens, mcc)

        print 'cross-validating...'
        yp = model.predConfussion()
        
        for i in range (self.modelLV):
            sensp = sensitivity(model.TPpred[i],model.FNpred[i])
            specp = specificity(model.TNpred[i],model.FPpred[i])
            mccp  = MCC(model.TPpred[i],model.TNpred[i],model.FPpred[i],model.FNpred[i])

            print "pred LV:%-2d cutoff:%4.2f TP:%3d TN:%3d FP:%3d FN:%3d spec:%5.3f sens:%5.3f MCC:%5.3f" % (i+1,
                    model.cutoff[i], model.TPpred[i], model.TNpred[i], model.FPpred[i], model.FNpred[i], specp, sensp, mccp)
        
        self.infoResult = []    
        self.infoResult.append( ('nobj',model.nobj) )
        self.infoResult.append( ('cutoff',str(self.modelCutoff) ) )
        
        self.infoResult.append( ('sens','%5.3f' % sens ) )
        self.infoResult.append( ('spec','%5.3f' % spec ) )
        self.infoResult.append( ('MCC' ,'%5.3f' % mcc ) )
        
        self.infoResult.append( ('sens pred','%5.3f' % sensp ) )
        self.infoResult.append( ('spec pred','%5.3f' % specp ) )
        self.infoResult.append( ('MCC  pred' ,'%5.3f' % mccp ) )

        return (yp)
        
    def ADAN (self, X, Y, yp):
        """Runs ADAN method for assessing the reliability of the prediction 

        """

        nrows, ncols = np.shape(X)
        
        # compute PP on X
        model = pls ()
        model.build (X,Y,targetSSX=0.8, autoscale=self.modelAutoscaling)
        model.saveModel (self.vpath+'/adan.npy')
        
        nlv = model.Am

        # initialize arrays
        T = np.empty((nrows,nlv),dtype=np.float64)
        centx = np.empty(nlv,dtype=np.float64)
        dcentx = np.empty(nrows, dtype=np.float64)
        dclosx = np.empty(nrows, dtype=np.float64)
        dcenty = np.empty(nrows, dtype=np.float64)
        dclosy = np.empty(nrows, dtype=np.float64)
        dpredy = np.empty(nrows, dtype=np.float64)
        
        # extract t 
        for a in range (nlv):
            ta = model.t[a]
            T[:,a]=ta
            centx[a]=np.mean(ta)

        i95 = np.round(nrows*0.95)-1
        # compute distances to X centroid (A) and percentil 95 
        for i in range(nrows):
            dcentx[i] = np.sqrt(np.sum(np.square(centx-T[i,:])))

        dcentx = np.sort(dcentx)
        p95dcentx = dcentx[i95]
   
        # compute closer distances in X (B) and percentil 95 
        for i in range (nrows):
            dclosx[i]=1e20
            closj=0
            for j in range (nrows):
                if j != i:
                    dtemp = np.sqrt(np.sum(np.square(T[j,:]-T[i,:])))
                    if dtemp < dclosx[i] :
                        dclosx[i] = dtemp
                        closj = j
            dclosy[i]=np.abs(Y[i]-Y[closj])
                    
        dclosx = np.sort(dclosx)
        p95dclosx = dclosx[i95]

        # compute DModX (C) and percentil 95
        dmodx = np.array(model.dmodx[-1])
        dmodx = np.sort(dmodx)
        p95dmodx=dmodx[i95]
        
        # compute distance to Y centroid (D) and percentil 95
        if self.quantitative:
            centy = np.mean(Y)
            dcenty = np.abs(Y-centy)
            dcenty = np.sort(dcenty)
            p95dcenty = dcenty[i95]
        else:
            centy = 0.0
            p95dcenty = 0.0       

        # compute Y of the closer compound (E) and percentil 95
        if self.quantitative:
            dclosy = np.sort(dclosy)
            p95dclosy = dclosy[i95]
        else:
            p95dclosy = 0.0

        # compute SDEP of 5% closer neighbours (F) and percentil 95
        if self.quantitative:
            p5 = int(np.round(nrows*0.05))
            if p5 < 1 : p5 = 1
            
            closerDis = np.empty(p5,dtype=np.float64)
            closerErr = np.empty(p5,dtype=np.float64)
            squareErr = np.empty(nrows,dtype=np.float64)
            squareErr = np.square(Y-yp)
            
            for i in range (nrows):
         
                closerDis.fill(1e20)
                    
                for j in range (nrows):
                    if j != i:
                        dtemp = np.sum(np.square(T[j,:]-T[i,:]))      
                        if dtemp < np.amax(closerDis) :
                            imindis=np.argmax(closerDis)
                            closerDis[imindis] = dtemp
                            closerErr[imindis] = squareErr[j]
                            
                dpredy[i]=np.sqrt(np.sum(closerErr)/float(p5))

            dpredy = np.sort(dpredy)
            
            p95dpredy = dpredy[i95]
        else:
            # TODO: the values in yp can be used to compute criteria G (equivalent to F)
            squareErr = np.empty(nrows,dtype=np.float64)
            p95dpredy = 0.0


##        print "DCENTX %6.3f \n" % (p95dcentx),
##        print "DCLOSX %6.3f \n" % (p95dclosx),
##        print "DCMODX %6.3f \n" % (p95dmodx),
##        print "DCENTY %6.3f \n" % (p95dcenty),
##        print "DCLOSY %6.3f \n" % (p95dclosy),
##        print "DPREDY %6.3f \n" % (p95dpredy)

        # write in a file, Am -> critical distances -> centroid -> t
        f = file (self.vpath+'/tscores.npy','wb')
        np.save(f,nlv)
        np.save(f,p95dcentx)
        np.save(f,p95dclosx)
        np.save(f,p95dmodx)
        np.save(f,p95dcenty)
        np.save(f,p95dclosy)
        np.save(f,p95dpredy)
        np.save(f,centx)
        np.save(f,centy)
        np.save(f,T)
        np.save(f,Y)
        np.save(f,squareErr)
        f.close()

        return (True, "Model OK")


    def cleanConfidentialFiles (self):
        """Removes from the current model version all non-relevant files that can compromise
           structures of the training series, for models built using the Confidential option

           Note: if imodel.py is called again (e.g. by predict), the file imodel.pyc can appear again
           TODO: implement a export method in manage designed ad hoc for Confidential models
        """

        preserve = ['imodel.py',
                    'distiledPLS.txt',
                    'info.pkl']
        
        if self.MD == 'padel' and self.padelDescriptor:
            dname,fname = os.path.split(self.padelDescriptor)
            preserve.append(fname)
        
        for item in os.listdir (self.vpath):
            if item in preserve : continue
            removefile (self.vpath+'/'+item)
        

    def build (self):
        """Uses the data extracted from the training series to build a model, using the Rlearner object 

           This function also creates the "itrain.txt" file that describes the training series, including InChiKey of the compounds
        """
        if not self.buildable:
            return (False, 'this model cannot by built automatically')
 
        X,Y = self.getMatrices ()

        nobj, nvarx = np.shape(X)
        if (nobj==0) or (nvarx==0) : return (False, 'no MD generated')

        nobj = np.shape(Y)
        if (nobj==0) : return (False, 'no activity found')

        if 'pls' in self.model:
            model = self.buildPLS (X,Y)

            if self.quantitative:
                yp = self.diagnosePLS_R (model)
            else:
                yp = self.diagnosePLS_DA (model)

            if self.confidential:
                model.saveDistiled (self.vpath+'/distiledPLS.txt')
            else:
                model.saveModel (self.vpath+'/modelPLS.npy')
        else:
            return (False, 'modeling method not recognised')

        if self.confidential:
            self.cleanConfidentialFiles()
            return (True, 'Confidential Model OK')
        else:
            success, result = self.ADAN (X,Y,yp)
            return (success, result)
        
##################################################################
##    VIEW METHODS
##################################################################   


    def viewPlotBackground (self):

        backname = wkd + '/' + self.viewReferenceEndpoint + '/version%0.4d' % self.viewReferenceVersion
        if   self.viewType == 'property':
            backname += '/backproperty.txt'
        elif self.viewType == 'pca':
            backname += '/backpca.txt'

        if not os.path.exists (backname):
            return (False)

        try:
            f = open (backname, 'r')  
        except:
            return (False)
        
        for line in f:
            lxy= line.split()
            try:
                x = float (lxy[-1])
            except:
                continue
            plt.scatter(lxy[-2],lxy[-1], c='#aaaaaa', marker='o', s=30, linewidths=0)

        f.close()

        return (True)

        
    def viewPCA (self):
        """Uses the data extracted from the training series to build a model, using the Rlearner object 

           This function also creates the "itrain.txt" file that describes the training series, including InChiKey of the compounds
        """
        
        X,Y = self.getMatrices ()

        model = pca ()

        model.build (X, 2)
        model.saveModel (self.vpath+'/pcmodel.npy')

        fig1=plt.figure(figsize=(9,6))
        plt.xlabel('PC 1')
        plt.ylabel('PC 2')
        
        if self.viewBackground :
            self.viewPlotBackground()
        
        plt.scatter(model.t[0],model.t[1], c='red', marker='D', s=40, linewidths=0)

        if os.path.isfile ('pca-scores12.png'):
            removefile ('pca-scores12.png')
            
        fig1.savefig("pca-scores12.png", format='png')

        # write a file with experimental Y (yp[0]) vs LOO predicted Y 
        ft=open ('pca-scores12.txt','w')

        # write simple header
        ft.write ('Name PC1 PC2\n')
        
        for i in range (model.nobj):
            ft.write(self.tdata[i][0]+' ')
            ft.write('%.3f '% model.t[0][i])
            ft.write('%.3f '% model.t[1][i])
            ft.write('\n')
        ft.close()

        shutil.copy ('./pca-scores12.txt', self.vpath+'/backpca.txt')
        
        return (True, 'pca-scores12.png')


    def viewProperty (self):
        """Uses the data extracted from the training series to build a model, using the Rlearner object 

           This function also creates the "itrain.txt" file that describes the training series, including InChiKey of the compounds
        """
        X,Y = self.getMatrices ()
        
        fig1=plt.figure(figsize=(9,6))

        plt.xlabel('log P')
        plt.ylabel('MW')

        if self.viewBackground : self.viewPlotBackground()

        plt.scatter (X[:,0],X[:,1], c='red', marker='D', s=40, linewidths=0)

        if os.path.isfile ('generic.png'):
            removefile ('generic.png')
        
        fig1.savefig("generic.png", format='png')

        # write a file with experimental Y (yp[0]) vs LOO predicted Y 
        ft=open ('generic.txt','w')

        # write simple header
        ft.write ('Name logP MW\n')

        nrows, ncols = np.shape(X)
        
        for i in range (nrows):
            ft.write(self.tdata[i][0]+' ')
            ft.write('%.3f '% X[i,0])
            ft.write('%.3f '% X[i,1])
            ft.write('\n')
        ft.close()
        
        shutil.copy ('./generic.txt', self.vpath+'/backproperty.txt')
        
        return (True, 'generic.png')
        

    def viewProject (self):
        """Uses the data extracted from the training series to build a model, using the Rlearner object 

           This function also creates the "itrain.txt" file that describes the training series, including InChiKey of the compounds
        """

        X,Y = self.getMatrices ()

        XX = []
        tt = []
        dd = []
        
        model = pca ()

        refname = wkd+ '/' + self.viewReferenceEndpoint + '/version%0.4d/pcmodel.npy'    % self.viewReferenceVersion
        
        if not os.path.isfile (refname):
            return (False, 'no reference PCA model found')

        model.loadModel (refname)
        
        # projects the data on the loaded model
            
        for md in X:
            if self.MD == 'pentacle' :
                XX.append(self.adjustPentacle(md,len(self.pentacleProbes),model.nvar))
            else:
                XX.append(md)
         
        for i in range(2):
            success, result = model.projectPC(XX,i)
            if success:
                XX, t, dmodx = result
                tt.append(t)
                dd.append(dmodx)
            else:
                return (False, 'no projection')

        fig1=plt.figure(figsize=(9,6))
        plt.xlabel('PC 1 (projected)')
        plt.ylabel('PC 2 (projected)')

        if self.viewBackground : 
            plt.scatter(model.t[0], model.t[1], c='#aaaaaa', marker='o', s=20, linewidths=0)
            
        plt.scatter(tt[0],tt[1], c=dd[-1], marker='o', s=50)
        
        plt.colorbar()

        if os.path.isfile ('pca-scores12.png'):
            removefile ('pca-scores12.png')
            
        fig1.savefig("pca-scores12.png", format='png')

        # write a file with experimental Y (yp[0]) vs LOO predicted Y 
        ft=open ('pca-scores12.txt','w')

        # write simple header
        ft.write ('Name PC1 PC2 dmodx\n')
        
        for i in range (len(tt)):
            ft.write(self.tdata[i][0]+' ')
            ft.write('%.3f '% tt[0][i])
            ft.write('%.3f '% tt[1][i])
            ft.write('%.3f '% dd[-1][-1])
            ft.write('\n')
        ft.close()

        #shutil.copy ('./pca-scores12.txt', self.vpath+'/backpca.txt')
        
        return (True, 'pca-scores12.png')
  

    def view (self):

        success = False
        
        if self.viewType == 'pca':
            success = self.viewPCA () 
        elif self.viewType == 'property':
            success = self.viewProperty ()
        elif self.viewType == 'project':
            success = self.viewProject ()

        return (success)
        
    
        
##################################################################
##    LOG METHODS
##################################################################    


    def log (self):
        """Stores information about the Series, the MD, the Model and the Results to the info.pkl file

           This information us used only for inform the user and not to recover the status of the model 
        """

        self.infoID = []
        self.infoID.append (('version', '*'))
        self.infoID.append (('date', time.asctime(time.localtime(time.time()))))

        self.infoID.append (('buildable', str(self.buildable) ))
        if self.quantitative:
            self.infoID.append (('dependent', 'quantitative'))
        else:
            self.infoID.append (('dependent', 'qualitative'))
            
        self.infoMD = []

        if 'pentacle' in self.MD:
            self.infoMD.append( ('MD','Pentacle') )
            for probe in self.pentacleProbes:
                self.infoMD.append ( ('probe',probe) )
            for key in self.pentacleOthers:
                self.infoMD.append ( ('key',key) )
        elif 'padel' in self.MD:
            self.infoMD.append( ('MD','PaDEL') )
            if self.padelDescriptor:
                self.infoMD.append( ('descriptors', self.padelDescriptor) )
            if self.padelMaxRuntime:
                self.infoMD.append( ('max runtime', str(self.padelMaxRuntime)) )
        elif 'adriana' in self.MD:
            self.infoMD.append( ('MD','Adriana') )
            
        try:
            modelInfo = open (self.vpath+'/info.pkl','wb')
        except:
            return (False, 'Failed to write model log')
                
        pickle.dump(self.infoID, modelInfo)
        pickle.dump(self.infoSeries, modelInfo)
        pickle.dump(self.infoMD, modelInfo)
        pickle.dump(self.infoModel, modelInfo)
        pickle.dump(self.infoResult, modelInfo)
        
        modelInfo.close()
        
        return (True, "Model OK")

##################################################################
##    WORKFLOW METHODS
##################################################################

    def buildWorkflow(self, molecules):

        if not self.licenseTesting ():
            return (False, 'licenses not found')
        
        if not self.buildable:
            success, result = self.log ()
            if not success:
                return (False, result)
            return (result)

        # load data, if stored, or compute it from the provided SDFile

        dataReady = False

        if not molecules:
            dataReady = self.loadData ()

            if not self.loadSeriesInfo ():
                self.setSeries ('training.sdf', len(self.tdata))

        if not dataReady: # datList was not completed because load failed or new series was set

            # estimate number of molecules inside the SDFile
            nmol = 0
            try:
                f = open (self.vpath+'/training.sdf','r')
            except:
                return (False,"Unable to open file %s" % molecules)
            for line in f:
                if '$$$$' in line: nmol+=1
            f.close()

            if not nmol:
                return (False,"No molecule found in %s:  SDFile format not recognized" % molecules)

            self.setSeries (molecules, nmol)

            print (molecules, nmol)  # DEBUG ONLY. REMOVE!!!!

            i = 0
            fout = None
            mol = ''

            # open SDFfile and iterate for every molecule
            f = open (self.vpath+'/training.sdf','r')

            # clean normalized structures
            removefile (self.vpath+'/tstruct.sdf')
            
            updateProgress (0.0)

            for line in f:
                if not fout or fout.closed:
                    i += 1
                    mol = 'm%0.10d.sdf' % i
                    fout = open(mol, 'w')

                fout.write(line)

                if '$$$$' in line:
                    fout.close()

                    ## workflow for molecule i (mol) ############
                    success, result = self.normalize (mol)
                    if not success:
                       writeError('error in normalize: '+result)
                       continue

                    molFile   = result[0]
                    molName   = result[1]
                    molCharge = result[2]
                    molPos    = self.saveNormalizedMol(molFile)

                    success, infN = self.extract (molFile,molName,molCharge,molPos)
                    if not success:
                       writeError('error in extract: '+ str(infN))
                       continue

                    updateProgress (float(i)/float(nmol))
                    ##############################################

                    removefile (mol)

            f.close()
            if fout :
                fout.close()

            self.saveData ()

        # build the model with the datList stored data

        success, result = self.build ()
        if not success:
            return (False, result)

        success, result = self.log ()
        if not success:
            return (False, result)

        return (result)


    def predictWorkflow(self, molecules, detail, progress):

        if not self.licenseTesting ():
            return (False, 'licenses not found')
        
        datList = []
        datList = self.loadData ()
        
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
                success, result  = self.normalize (mol)
                if not success:
                    pred.append((False, result))
                    continue

                molFile   = result[0]
                molName   = result[1]
                molCharge = result[2]

                predN = self.predict (molFile, molName, molCharge, detail)

                pred.append((True, predN))
                ############################################

                if progress:
                    sys.stdout.write('completed: %d\n'%i)
                    sys.stdout.flush()

                removefile(mol)

        return (True, pred)


    def viewWorkflow(self, molecules):

        # load data, if stored, or compute it from the provided SDFile
        
        dataReady = False

        if not molecules:

            if self.viewType == 'pca' or self.viewType == 'project':
                dataReady = self.loadVisualData ()
            else:
                dataReady = self.loadPropertyData ()

            if not dataReady:
                molecules = self.vpath+'/training.sdf'

##            if not self.loadSeriesInfo ():
##                self.setSeries ('training.sdf', len(self.tdata))
            
        if not dataReady: # datList was not completed because load failed or new series was set
            
            # estimate number of molecules inside the SDFile
            nmol = 0
            try:
                f = open (molecules,'r')
            except:
                return (False,"Unable to open file %s" % molecules)
            for line in f:
                if '$$$$' in line: nmol+=1
            f.close()

            if not nmol:
                return (False,"No molecule found in %s:  SDFile format not recognized" % molecules)

            #self.setSeries (molecules, nmol)

            i = 0
            fout = None
            mol = ''

            # open SDFfile and iterate for every molecule
            f = open (molecules,'r')

            updateProgress (0.0)

            for line in f:
                if not fout or fout.closed:
                    i += 1
                    mol = 'm%0.10d.sdf' % i
                    fout = open(mol, 'w')

                fout.write(line)

                if '$$$$' in line:
                    fout.close()

                    ## workflow for molecule i (mol) ############
                    success, result = self.normalize (mol)
                    if not success:
                       writeError('error in normalize: '+result)
                       continue

                    molFile   = result[0]
                    molName   = result[1]
                    molCharge = result[2]

                    success, infN = self.extractView (molFile,molName,molCharge)
                    if not success:
                       writeError('error in extract: '+ str(infN))
                       continue

                    updateProgress (float(i)/float(nmol))
                    ##############################################

                    removefile (mol)

            f.close()
            if fout :
                fout.close()

            if self.viewType == 'property':
                self.savePropertyData ()

        success = self.view ()

        return (success)
