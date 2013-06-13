from model import model

class imodel(model):
    def __init__ (self, vpath):
        model.__init__(self, vpath)

        # set common model settings
        self.numLV = 2
        self.pH = 7.4

##        self.norma = True
##        self.norma-standard = True
##        self.norma-neutraliz = True
##        self.norma-3D = True
##
##        self.MD = Pentacle
##
##        self.model = PLS
##        self.model-autoscaling = False
