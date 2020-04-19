import OmicsIntegrator as oi
import leidenalg as la

class hypha:

    def __init__(self):
        '''Constructor for hypha'''
        self.forestParams={}
        self.communityParams={}
        self.networks=[]

    def runForest(self):
        print('running forest')
        oi.graph()

    def runCommunity(self):
        print('running community detection)')
