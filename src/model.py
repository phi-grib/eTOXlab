# -*- coding: utf-8 -*-

##    Description    eTAM model class
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
from StringIO import StringIO
from utils import removefile
from utils import opt
from utils import randomName

from rdkit import Chem
from rdkit import RDLogger
from standardise import standardise
from qualit import *

class model:
            
    def __init__ (self, vpath):

        self.vpath = vpath

        self.trainList = []
        try:
            ifile = open (self.vpath+'/itrain.txt')
            for line in ifile:
                piece = line.split('\t')
                self.trainList.append(piece)
            ifile.close()
        except:
            pass

        ##
        ## General settings
        ##
        self.buildable = False
        self.quantitative = False
        
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
        self.modelAutoscaling = False
        self.modelCutoff = None

        # Info lists serve only to store properties of new models
        # and inform the users. This list does not set model properties
        self.infoID = []
        self.infoSeries = []
        self.infoMD = []
        self.infoModel = []
        self.infoResult = []


##################################################################
##    SHARED (PREDICTION & BUILD) METHODS
##################################################################    

    def loadData (self):
        datList = []

        if not os.path.isfile (self.vpath+'/data.pkl'):
            return datList

        try:
            f = open (self.vpath+'/data.pkl','rb')
        except:
            return datList

        norm = pickle.load(f)
        if norm != self.norm:
            return datList

        if norm:
            normStand = pickle.load(f)
            if normStand != self.normStand:
                return datList
            
            normNeutr = pickle.load(f)
            if normNeutr != self.normNeutr:
                return datList

            if normNeutr:
                normNeutrMethod = pickle.load(f)
                if normNeutrMethod != self.normNeutrMethod:
                    return datList

                normNeutr_pH = pickle.load(f)
                if normNeutr_pH != self.normNeutr_pH:
                    return datList
                
            norm3D = pickle.load(f)
            if norm3D != self.norm3D:
                return datList

            MD = pickle.load(f)
            if MD != self.MD:
                return datList
            
            if 'pentacle' in MD:
                pentacleProbes = pickle.load(f)
                if pentacleProbes != self.pentacleProbes:
                    return datList

                pentacleOthers = pickle.load(f)
                if pentacleOthers != self.pentacleOthers:
                    return datList

            elif 'padel' in MD:
                padelMD = pickle.load(f)
                if padelMD != self.padelMD:
                    return datList
                
        datList = pickle.load(f)
        f.close()

        return datList

    def saveData (self,datList):
        try:
            f = open (self.vpath+'/data.pkl','wb')
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
        pickle.dump(datList, f)

        f.close()

        
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
    
    def computeMDPentacle (self, mol, clean=True):
        """ Computes the Molecular Descriptors for compound "mol"

            In this implementation we run Pentacle with default settings
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
        
        call = [opt+'pentacle_etox/pentacle',
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
            url = 'http://localhost:9000/computedescriptors?params='+params
            #print "call "+ url
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
        try:
            shutil.rmtree('padel')
        except:
            pass
        
        os.mkdir ('padel')
        shutil.copy (mol,'padel')

        call = [opt+'jdk/bin/java','-Djava.awt.headless=true','-jar',
                opt+'padel/PaDEL-Descriptor.jar',
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
    
    def computeMD (self, mol, clean=True):

        if 'pentacle' in self.MD:
            success, md = self.computeMDPentacle (mol, clean)
        elif 'padel' in self.MD:
            success, md = self.computeMDPadelws (mol, clean)

        return (success, md)
    
##################################################################
##    NORMALIZE METHODS
##################################################################    
        
    def standardize (self, moli, clean=True):
        """Applies a structure normalization protocol provided by Francis Atkinson (EBI)

           The name of the output molecules is built as a+'original name'

           Returns a tuple containing:
           1) True/False: depending on the result of the method
           2) (if True ) The name of the output molecule
              (if False) The error message
        """
        molo = 'a'+moli

        #print molo
        
        try:
            suppl=Chem.SDMolSupplier(moli)
            m = suppl.next()

            if m is None:
                return (False, 'wrong input format')
            
        except:
            return (False, 'wrong input format')
            
        try:
            parent = standardise.apply(Chem.MolToMolBlock(m))
        except standardise.StandardiseException as e:
            return (False, e.name)
            
        fo = open (molo,'w')
        fo.write(parent)
        
        if m.HasProp('activity'):
            activity = m.GetProp('activity')
            fo.write('>  <activity>\n'+activity+'\n\n$$$$')
            
        fo.close()

        if clean:
            removefile (moli)

        return (True,molo)

    def protonate (self, moli, pH, clean=True):
        """Adjusts the ionization state of the molecule "moli" 

           In this implementation, it ises blabber_sd from Molecular Discovery
           The result is a tuple containing:
           1) True/False: describes the success of the protonation for this compound
           2) (if True ) The name of the protonated molecules and its formal charge
              (if False) The error message
        """

        molo = 'b'+moli

        stderrf = open (os.devnull, 'w')
        stdoutf = open (os.devnull, 'w')     

        call = [opt+'blab_etox/blabber_sd', moli,
                '-p',  str(pH),
                '-o',  molo]

        try:
            retcode = subprocess.call(call,stdout=stdoutf, stderr=stderrf)
        except:
            return (False, 'Blabber execution error', 0.0)
        
        stdoutf.close()
        stderrf.close()

        # for the newer versions use !=0, for the old one use ==0
        if retcode != 0:
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
        
        call = [opt+'corina/corina',
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
        
        for l in self.trainList:
            if ik in l[0]:

                yp = float (l[1])
                
                if self.quantitative:
                    return (True, yp)
                if (yp < ypcutoff):
                    return (True, 'negative')
                else:
                    return (True, 'positive')
        
        return (False, mol)


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

        if not self.norm:
            return (True, (mol, charge))
        
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

        return (True,(resultc,charge))


##################################################################
##    PREDICT METHODS
##################################################################    
   

    def computePR (self, md, charge):
        """ Computes the prediction for compound "mol"

            This function makes use of the molecular descriptor vector (md) to project the compound in the model
            The model has been loaded previously as an R object
        """
        
        model = pls()
        model.loadModel(self.vpath+'/modelPLS.npy')

        if 'pentacle' in self.MD:
            md = self.adjustPentacle(md,len(self.pentacleProbes),model.nvarx)
            
        success, result = model.project(md,self.modelLV)

        if success:
            yp = result[0][-1] # yp is an array of modelLV elements, pick the last

            if not self.quantitative:
                if model.cutoff is None:
                    return (False, 'cutoff not defined')
                if yp < model.cutoff[-1]: # use last cutoff
                    return (True, 'negative')
                else:
                    return (True, 'positive')
            else:
                return (True, yp)
        else:
            return (False, result)


    def computeAD (self, md, pr, detail):
        """Carries out a protocol for determining how fat is the query compound from the model

           Provisionally, implements a temporary version of the ADAN method

           Returns a tuple that contains
           1) True or False, indicating the success of the computation
           2) (if True ) a number between 0 and 5 with the number of criteria broken
              (if False) an error message 
           

        """
        f = file (self.vpath+'/tscores.npy','rb')
        nlv = np.load(f)
        p95dcentx = np.load(f)
        p95dclosx = np.load(f)
        p95dmodx = np.load(f)
        p95dcenty = np.load(f)
        p95dclosy = np.load(f)
        centx = np.load(f)
        centy = np.load(f)
        T = np.load(f)
        Y = np.load(f)
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
        # compute distance t to centroid x
        dcentx = np.sqrt(np.sum(np.square(centx-t)))
        AD['dcentx']=dcentx>p95dcentx

        # compute min distance to T    
        dclosx=1e20
        for i in range (model.nobj):
            dtemp = np.sqrt(np.sum(np.square(T[i,:]-t)))
            if dtemp < dclosx : dclosx = dtemp
        AD['dclosx']=dclosx>p95dclosx

        # compute distance to modx
        AD['dmodx']=d[-1]>p95dmodx

        # compute distance to centroid y, only for quantitative
        if self.quantitative:
            dcenty = np.abs(y-centy)
            AD['dcenty']=dcenty>p95dcenty
        else:
            AD['dcenty']=False

        # compute min distance to Y
        if self.quantitative:
            dclosy=1e20
            for i in range (model.nobj):
                dtemp = np.abs(Y[i]-y)
                if dtemp < dclosy : dclosy = dtemp
            AD['dclosy']=dclosy>p95dclosy
        else:
            AD['dclosy']=False

##        print "DCENTX %6.3f (%6.3f)\n" % (dcentx,p95dcentx),
##        print "DCLOSX %6.3f (%6.3f)\n" % (dclosx,p95dclosx),
##        print "DCMODX %6.3f (%6.3f)\n" % (d[-1],p95dmodx),
##        print "DCENTY %6.3f (%6.3f)\n" % (dcenty,p95dcenty),
##        print "DCLOSY %6.3f (%6.3f)\n" % (dclosy,p95dclosy)

        return (True,sum(AD.values()))

    def computeRI (self, ad):
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
        
        ri=0.0
        
        if ad<4:
            model = pls ()
            model.loadModel(self.vpath+'/modelPLS.npy')
            sdep = model.SDEP[model.Av-1]
            if ad<2:
                ri = sdep    # 0 or 1 criteria broken
            else:
                ri = sdep*2  # 2 or 3 criteria broken
        
        return (True,ri)
  

    def predict (self, molN, detail, clean=True):
        """Runs the prediction protocol

        """

        # alias
        mol, charge = molN[0], molN[1]

        # default return values
        pr=ri=ad=(False,0.0)

        ic = self.checkIdentity (mol)
        if ic[0]: return (ic, (True, 0), (True, 0.0))
     
        md = self.computeMD (mol)
        if not md[0]: return (pr,ad,ri)
        
        pr = self.computePR (md[1],charge)
        if not pr[0]: return (pr,ad,ri)

        ad = self.computeAD (md[1], pr[1], detail)
        if not ad[0]: return (pr,ad,ri)

        ri = self.computeRI (ad[1])

        if clean:
            removefile (mol)
            
        return (pr,ad,ri)


##################################################################
##    EXTRACT METHODS
##################################################################          

    def getInChi (self, mol):
        """ Computes the InChiKey for the compound "mol"

            We used InChiKeys without the last three chars (ionization state) to make the comparison
        """

        lg = RDLogger.logger()
        lg.setLevel(RDLogger.ERROR)
        
        suppl = Chem.SDMolSupplier(mol)
        mi = suppl.next()
        try:
            ik = Chem.InchiToInchiKey(Chem.MolToInchi(mi))
        except:
            return (False, 'Failed to obtain InChiKey')
        return (True,ik[:-3])


    def getBio (self, mol):
        """ Extracts the value of the experimental biological property from the compound "mol" 

            Such value must be identified by the tag <activity>
        """
        
        suppl = Chem.SDMolSupplier(mol)
        mi = suppl.next()

        if mi.HasProp('activity'):
            bio = mi.GetProp('activity')
            try:
                nbio = float (bio)
            except:
                return (False, 'Activity cannot be converted to numerical format')
            return (True, nbio)
        else:
            return (False, 'Biological activity not found as <activity>')


    def extract (self, mol, clean=True):
        """Process the compound "mol" for obtaining
           1) InChiKey (string)
           2) Molecular Descriptors (NumPy float64 array)
           3) Biological Activity (float)

           Returns a Tuple as (True,('InChi',array[0.0,0.0, 0.0],4.56))
        """

        if not self.buildable:
            return (False, 'this model cannot by built automatically')

        moli = mol[0]
        i1=''
        i2=[]
        i3=0.0
        
        success, i1 = self.getInChi(moli)
        if not success:
            return (False,(i1,i2,i3))
        
        success, i2 = self.computeMD (moli)
        if not success:
            return (False,(i1,i2,i3))

        success, i3 = self.getBio (moli)
        if not success:
            return (False,(i1,i2,i3))

        if clean:
            removefile (moli)
            
        return (True,(i1,i2,i3))


##################################################################
##    BUILD METHODS
##################################################################    

    def setSeries (self, molecules, numMol):
        self.infoSeries = []
        self.infoSeries.append ( ('series',molecules) )
        self.infoSeries.append ( ('nmol',numMol) )

    def saveTraining (self, data):
        ftrain = open (self.vpath+'/itrain.txt','w')
        for success, i in data:          
            ftrain.write (i[0]+'\t'+str(i[2])+'\n') #  i0 (InChi) + i2 (activity)
        ftrain.close()
        
    def getMatrices (self, data):  
        ncol = 0
        xx = []
        yy = []
        
        # obtain X and Y
        for success, i in data:
            if len(i[1])>ncol: ncol = len(i[1])
            xx.append(i[1])
            yy.append(i[2])

        nrow = len (xx)
        
        Y = np.array (yy)
        X = np.empty ((nrow,ncol),dtype=np.float64)
      
        i=0
        for row in xx:
            if 'pentacle' in self.MD:
                row=self.adjustPentacle(row,len(self.pentacleProbes),ncol)
            X[i,:]=np.array(row)
            i+=1

        return X, Y
    
    def buildPLS (self, X, Y):
        model = pls ()
        model.build (X,Y,self.modelLV,autoscale=self.modelAutoscaling)

        self.infoModel = []
        if self.quantitative:
            self.infoModel.append( ('model','PLS-R  (NIPALS)') )
        else:
            self.infoModel.append( ('model','PLS-DA (NIPALS)') )
        self.infoModel.append( ('LV', self.modelLV ))
        return model
    
    def diagnosePLS_R (self, model):
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
            plt.ylabel('experimental y')
            plt.xlabel('predicted LV'+nvar)
            plt.title('Predicted')
            plt.plot(yp[:,i+1],yp[:,0],"ro")
            fig1.savefig("pls-predicted-LV"+nvar+".png", format='png')
            
            fig2=plt.figure()
            plt.ylabel('experimental y')
            plt.xlabel('recalculated LV'+nvar)
            plt.title('Recalculated')
            plt.plot(yr[:,i+1],yr[:,0],"ro")
            fig2.savefig("pls-recalculated-LV"+nvar+".png", format='png')
        #plt.show()

        # write a file with experimental Y (yp[0]) vs LOO predicted Y 
        fp=open ('pls-predicted.txt','w')
        fr=open ('pls-recalculated.txt','w')

        for i in range (model.nobj):
            for j in range (self.modelLV+1):
                fp.write('%.3f '% yp[i][j])
                fr.write('%.3f '% yr[i][j])
            fp.write('\n')
            fr.write('\n')
        fp.close()
        fr.close()      

    def diagnosePLS_DA (self, model, data):

        if 'auto' == self.modelCutoff:
            model.calcOptCutoff ()
        else:
            model.calcConfussion(self.modelCutoff)

        for i in range (self.modelLV):
            sens = sensitivity(model.TP[i],model.FN[i])
            spec = specificity(model.TN[i],model.FP[i])
            mcc  = MCC(model.TP[i],model.TN[i],model.FP[i],model.FN[i])

            print "LV:%-2d cutoff:%4.2f TP:%3d TN:%3d FP:%3d FN:%3d spec:%5.3f sens:%5.3f MCC:%5.3f" % (i+1,
                    model.cutoff[i], model.TP[i], model.TN[i], model.FP[i], model.FN[i], spec, sens, mcc)

        self.infoResult = []    
        self.infoResult.append( ('nobj',model.nobj) )
        self.infoResult.append( ('cutoff',str(self.modelCutoff) ) )
        self.infoResult.append( ('sens','%5.3f' % sens ) )
        self.infoResult.append( ('spec','%5.3f' % spec ) )
        self.infoResult.append( ('MCC' ,'%5.3f' % mcc ) )
        
    def ADRI (self, X, Y):

        nrows, ncols = np.shape(X)
        
        # compute PP on X
        model = pls ()
        model.build (X,Y,targetSSX=0.4)
        model.saveModel (self.vpath+'/adan.npy')
        
        nlv = model.Am

        # initialize arrays
        T = np.empty((nrows,nlv),dtype=np.float64)
        centx = np.empty(nlv,dtype=np.float64)
        dcentx = np.empty(nrows, dtype=np.float64)
        dclosx = np.empty(nrows, dtype=np.float64)
        dcenty = np.empty(nrows, dtype=np.float64)
        dclosy = np.empty(nrows, dtype=np.float64)
        
        # extract t 
        for a in range (nlv):
            ta = model.t[a]
            T[:,a]=ta
            centx[a]=np.mean(ta)

        i95 = np.ceil(nrows*0.95)-1
        # compute distances to X centroid and percentil 95
        for i in range(nrows):
            dcentx[i] = np.sqrt(np.sum(np.square(centx-T[i,:])))

        dcentx = np.sort(dcentx)
        p95dcentx = dcentx[i95]

        # compute closer distances in X and percentil 95
        for i in range (nrows):
            dclosx[i]=1e20
            for j in range (nrows):
                if j != i:
                    dtemp = np.sqrt(np.sum(np.square(T[j,:]-T[i,:])))
                    if dtemp < dclosx[i] : dclosx[i] = dtemp
                    
        dclosx = np.sort(dclosx)
        p95dclosx = dclosx[i95]

        # compute percentil 96 dmodx
        dmodx = np.array(model.dmodx[-1])
        p95dmodx=dmodx[i95]
        
        # compute distance to Y centroid and percentil 95
        centy = np.mean(Y)
        dcenty = np.abs(Y-centy)
        dcenty = np.sort(dcenty)
        p95dcenty = dcenty[i95]

        # compute mutual distance in Y and percentil 95
        for i in range (nrows):
            dclosy[i]=1e20
            for j in range (nrows):
                if j != i:
                    dtemp = np.abs(Y[i]-Y[j])
                    if dtemp < dclosy[i] : dclosy[i] = dtemp

        dclosy = np.sort(dclosy)
        p95dclosy = dclosy[i95]

        # write in a file, Am -> critical distances -> centroid -> t
        f = file (self.vpath+'/tscores.npy','wb')
        np.save(f,nlv)
        np.save(f,p95dcentx)
        np.save(f,p95dclosx)
        np.save(f,p95dmodx)
        np.save(f,p95dcenty)
        np.save(f,p95dclosy)
        np.save(f,centx)
        np.save(f,centy)
        np.save(f,T)
        np.save(f,Y)
        f.close()

        return (True, "Model OK")


    def build (self, data):
        """Uses the data extracted from the training series to build a model, using the Rlearner object 

           This function also creates the "itrain.txt" file that describes the training series, including InChiKey of the compounds
        """
        if not self.buildable:
            return (False, 'this model cannot by built automatically')
        
        self.saveTraining (data)
 
        X,Y = self.getMatrices (data)

        if 'pls' in self.model:
            model = self.buildPLS (X,Y)

            if self.quantitative:
                self.diagnosePLS_R (model)
            else:
                self.diagnosePLS_DA (model,data)

            model.saveModel (self.vpath+'/modelPLS.npy')
        else:
            return (False, 'modeling method not recognised')

        success, result = self.ADRI (X,Y)
             
        return (success, result)


##################################################################
##    LOG METHODS
##################################################################    


    def log (self):

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
