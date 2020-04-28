from igraph import *
import leidenalg as la
import pandas as pd
import OmicsIntegrator as oi

class hypha:


    def __init__(self,proteinWeights,interactome):
        '''Constructor for hypha'''
        self.proteins = proteinWeights
        self.interactome = interactome
        self.networks = {}
        nodeSet = set()
        for pat in proteinWeights.keys():
            print('Running PCSF for patient'+pat)
            pdf = pd.DataFrame(list(proteinWeights[pat].items()))
            pdf.columns = ['name','prize']
            graph  = self.interactome
            graph._prepare_prizes(pdf)
            verts, edges = graph.pcsf()
            nodeSet.add(verts)
            print('found ',len(verts),' vertices and ',len(edges),' edges for patient ',pat)
            #     self.networks[pat] =self._getSubnet(proteinWeights[pat].keys())
            forest, augmented_forest = graph.output_forest_as_networkx(verts, edges)
            self.networks[pat] = forest

    def makeGraph(gfile):
        print('making single graph')
        graph = oi.Graph(gfile)
        return(graph)

    def _getSubnet(self,nodenames):
        """
        Searches for nodes in the interactome that have nodes of a specific name,
        then returns inferred subgraph, nothing fancy
        """

        nodelist=[]
        for n in nodenames:
            try:
                nodelist.append(self.interactome.vs.find(name = n))
            except ValueError:
                print('Node '+n+' not found')
       # print("Reducing",len(nodenames),"proteins to",len(nodelist)," nodes")
        sg = self.interactome.subgraph(nodelist)
        print("Subgraph has",len(sg.vs),"nodes")
        return(sg)

    def _getForest(self,nodeweights):
        """
        Uses the omics integrator package to build a weighted subgrpah, inferring
        additional nodes and returning a larger subgraph on the interatome
        """
        #map nodes to interactome
        #

    def runCommunity(self):
        #first we have to create a union/interaction
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
