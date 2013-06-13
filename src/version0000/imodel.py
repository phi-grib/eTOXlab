from model import model

class imodel(model):
    def __init__ (self, vpath):
        model.__init__(self, vpath)

        # set common model settings
        self.numLV = 2
        self.pH = 7.4

