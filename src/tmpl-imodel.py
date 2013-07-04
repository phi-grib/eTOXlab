from model import model

class imodel(model):
    def __init__ (self, vpath):
        model.__init__(self, vpath)

        ##
        ## General settings
        ##
        self.buildable = True
        self.quantitative = True
        
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
        self.pentacleProbes = ['DRY','O','N1','TIP']
        self.pentacleOthers = []

        ##
        ## Modeling settings
        ##
        self.model = 'pls'
        self.modelLV = 2
        self.modelAutoscaling = False
