from igraph import *
import leidenalg as la
import pandas as pd

class hypha:


    def __init__(self,proteinWeights,interactome):
        '''Constructor for hypha'''
        self.proteins = proteinWeights
        self.interactome = interactome
        self.networks = {}
        for pat in proteinWeights.keys():
            self.networks[pat] =self._getSubnet(proteinWeights[pat].keys())

    def _getSubnet(self,nodenames):
        """
        Searches for nodes in the interactome that have nodes of a specific name,
        then returns subnetwork
        """

        nodelist=[]
        for n in nodenames:
            try:
                nodelist.append(self.interactome.vs.find(name = n))
            except ValueError:
                print('Node '+n+' not found')
       # print("Reducing",len(nodenames),"proteins to",len(nodelist)," nodes")
        return(self.interactome.subgraph(nodelist))

    def runCommunity(self):
        print('running community detection')
        optimizer = la.Optimiser()
        netlist = [self.networks[p] for p in self.networks.keys()]
        [membership,improv] =la.find_partition_multiplex([netlist],la.ModularityVertexPartition)
        return(membership,improv)

    def networkStats():
        print("We have "+len(self.networks)+"networks")
        #plot number of nodes in each network with number of edges

    def plotNetwork():
        """
        Uses some sort of graph visualization to plot network

        Returns
        -------
        None.

        """


    def mapHGNCtoNetwork():
        tab = pd.read_csv('data/human.name_2_string.tsv','\t',skiprows=[0],header=None)
        res = dict(zip(tab[1],tab[2]))
        return(res)

    def mapNetworkToHGNC():
        tab = pd.read_csv('data/human.name_2_string.tsv','\t',skiprows=[0],header=None)
        res = dict(zip(tab[2],tab[1]))
        return(res)
