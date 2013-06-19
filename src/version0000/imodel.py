from model import model

import os
import numpy as np
from pls import pls

class imodel(model):
    def __init__ (self, vpath):
        model.__init__(self, vpath)

        # set common model settings
        self.numLV = 1
        self.pH = 4.8

        # set pentacle settings
        self.pentacleProbes = ('DRY','O','N1')
        self.pentacleOthers = ('macc2_window 1.6','step 1.3')                               

        if os.path.isfile(self.vpath+'/cutoff.npy'):
            f = file (self.vpath+'/cutoff.npy','rb')
            self.cutoff = np.load(f)
            f.close()
        else:
            self.cutoff = 0.0

##        self.norma = True
##        self.norma-standard = True
##        self.norma-neutraliz = True
##        self.norma-3D = True
##        self.MD = Pentacle
##
##        self.model = PLS
##        self.model-autoscaling = False


    def extract (self, mol, clean=True):

        charge = mol[1]
        
        base = model.extract (self, mol, clean)

        return (base[0], (base[1][0], base[1][1], charge, base[1][2]))


    def build (self, data):

        ncol = 0
        xx = []
        yy = []

        # obtain X and Y
        for success, i in data:
            if i[2]<1 :
                continue
            if len(i[1])>ncol: ncol = len(i[1])
            xx.append(i[1])
            yy.append(i[3])

        nrow = len (xx)

        # new
        Y = np.array (yy)
        X = np.empty ((nrow,ncol),dtype=np.float64)
        i=0
        for row in xx:
            X[i,:]=self.adjustPentacle(row,len(self.pentacleProbes),ncol)
            i+=1
            
        nrows, ncols = np.shape(X)

        model = pls ()
        model.build (X,Y,self.numLV)
        model.validateLOO(self.numLV)
        
## evaluate the quality of the model. TODO       
##        for i in range (self.numLV):
##            print 'LV%2d R2: %6.4f Q2: %6.4f SDEP: %6.4f' % \
##                  (i+1,model.SSYac[i],model.Q2[i],model.SDEP[i])
            
        model.saveModel (self.vpath+'/modelPLS.npy')

        # write itraining data
        ftrain = open (self.vpath+'/itrain.txt','w')
        for success, i in data:
            ftrain.write (i[0]+'\t'+str(i[2])+'\n') # now it writes i0 (InChi) + i2 (activity)
        ftrain.close()
        
        # obtain the best cutoff according to some criteria. TODO
        self.cutoff = 0.44
    
        f = file (self.vpath+'/cutoff.npy','wb')
        np.save(f,self.cutoff)
        f.close()

        return (True, (X,Y))

        
    def predict (self, molN, detail, clean=True):

        if molN[1] < 1:
            return ((True,'negative'), (True, 100.0), (True, 100.0))
            if clean:
                removefile(mol)
        
        pr, ad, ri = model.predict (self, molN, detail, clean)
        
        if pr[0]:
            if pr[1]>= self.cutoff:
                npr = (True, 'positive')
            else:
                npr = (True, 'negative')
        else:
            npr = pr

        return (npr, ad, ri)


        


    
