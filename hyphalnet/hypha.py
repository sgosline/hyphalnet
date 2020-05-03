from igraph import *
import leidenalg as la
import pandas as pd
import OmicsIntegrator as oi
import networkx as nx
import goenrich
import numpy as np

class hypha:
    """hypha class defines multiple networks and communities """
    def __init__(self, proteinWeights, interactome):
        '''Constructor for hypha'''
        self.proteins = proteinWeights #dictionary of protein weights for each forest
        self.interactome = interactome #interactome
        self.forests = {} #forests
        self.node_counts = dict() #nodes found in forests and their counts
        self.communities = {} #communities from hypha
        for pat in proteinWeights.keys():
            print('Running PCSF for patient'+pat)
            forest = self._getForest(list(proteinWeights[pat].items()))
            for node in list(forest.nodes()):
                if node in list(self.node_counts.keys()):
                    self.node_counts[node] += 1
                else:
                    self.node_counts[node] = 1
            self.forests[pat] = forest
        print('Ready to create hypha across ', len(self.node_counts), 'nodes and',\
              len(proteinWeights), 'forests')

    def make_graph(gfile):
        print('making networkx graph for OmicsIntegrator')
        graph = oi.Graph(gfile)
        return graph

    def _get_subnet(self, nodenames):
        """
        Searches for nodes in the interactome that have nodes of a specific name,
        then returns inferred subgraph, nothing fancy
        """
        nodelist = []
        for n in nodenames:
            try:
                nodelist.append(self.interactome.vs.find(name = n))
            except ValueError:
                print('Node '+n+' not found')
        sg = self.interactome.subgraph(nodelist)
        print("Subgraph has", len(sg.vs), "nodes")
        return sg

    def _nx2igraph(self, nx_graph):
        """Helper function that maps networkx object to igraph"""
        g = Graph(directed=False)
        alln = self.node_counts.keys()
        for n in alln:
            g.add_vertex(n)
        for e in nx_graph.edges():
            g.add_edge(e[0], e[1])
        return g

    def _getForest(self, nodeweights):
        """
        Uses the omics integrator package to build a weighted subgrpah, inferring
        additional nodes and returning a larger subgraph on the interatome
        """
        #map nodes interactome
        pdf = pd.DataFrame(nodeweights)
        pdf.columns = ['name', 'prize']
        graph = self.interactome
        graph._prepare_prizes(pdf)
        verts, edges = graph.pcsf()
        forest, augmented_forest = graph.output_forest_as_networkx(verts, edges)
        return forest #TODO: add in parameter shift for 0 size networks

    def runCommunityWithMultiplex(self):
        """After forests are computed can run community detection"""
        print('running community detection')
        optimizer = la.Optimiser()
        netlist = []
        for nx in self.forests.values():
            netlist.append(self._nx2igraph(nx))
        [membership, improv] = la.find_partition_multiplex(netlist,\
                                                           la.ModularityVertexPartition)
        comm_df = pd.DataFrame({'Node': list(self.node_counts.keys()),\
                           'NumForests': list(self.node_counts.values()),\
                           'Community': membership})
        self.communities = dict(comm_df.groupby('Community')['Node'].apply(list))
        return comm_df

    def forest_stats(self):
        """Compute statistics for node membership and forests"""
        print("We have "+str(len(self.forests))+" forests")
        #plot number of nodes in each network with number of edges
        dflist = []
        for pat, forest in self.forests.items():
            protweights = pd.DataFrame(self.proteins[pat].items(), \
                                       columns=['Gene', 'DisWeight'])
           # print(nx.get_node_attributes(forest,'prize'))
            nodevals = pd.DataFrame(nx.get_node_attributes(forest, 'prize').items(),\
                                    columns=['Gene', 'ForestPrize'])
            #print(nodevals)
            pat_df = protweights.merge(nodevals, on='Gene', how='outer').fillna(0, downcast='infer')
            pat_df['Patient'] = pat
            dflist.append(pat_df)
        full_df = pd.concat(dflist)
        return full_df

    def saveCommunityToFile(self, prefix='',doAllGraphs=False):
        """
        Currently writes everything to cytoscape files

        Returns
        -------
        None.
        """
        print('Reducing network to node set')
        gred = nx.from_pandas_edgelist(self.interactome.interactome_dataframe,\
                                       'protein1', 'protein2', 'cost')
        #comS = set()
        for com, nodes in self.communities.items():
            cred = gred.subgraph(list(nodes))
            nx.set_node_attributes(cred, self.node_counts)
            print('Adding membership to network option')
            #nx.set_node_attributes(gred,membership,'Community')
            print('saving to html')
            oi.output_networkx_graph_as_graphml_for_cytoscape(cred, \
                                                              output_dir=".", \
                                                              filename=prefix+"_"+\
                                                              str(com)+\
                                                              "_pcsf_results.graphml.gz")
            oi.output_networkx_graph_as_ginteractive_html(cred, \
                                                         output_dir=".", \
                                                         filename=prefix+"_"+\
                                                         str(com)+\
                                                         "_pcsf_results.html")
        return None

    def map_hgnc_to_ensPep():
        tab = pd.read_csv('data/human.name_2_string.tsv', '\t', skiprows=[0], header=None)
        res = dict(zip(tab[1], tab[2]))
        return res

    def map_ensPep_to_hgnc():
        tab = pd.read_csv('data/human.name_2_string.tsv', '\t', skiprows=[0], header=None)
        res = dict(zip(tab[2], tab[1]))
        return res

    def graph2graph_distance(g1, g2):
        """Computes edit distance between two networkx objects"""
        print('Computing distance between graphs by xxx')

    def distance_to_communities(self, g_query):
        """computes distance between networkx graph and all communities in hypha"""
        print('Computing difference from graph to all communities')
