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

    def jaccard_distance(self, ns_1, ns_2):
        """Computes jaccard distance between two networkx objects"""
        #print('Computing distance between graphs by jaccard')
        u_size = len(ns_1.union(ns_2))
        if u_size == 0:
            return 1.0
        return 1-len(ns_1.intersection(ns_2))/u_size

    def distance_to_networks(self, g_query):
        """ Compute distance to all networks assuming this networks is net2"""
        g_nodes = set(g_query.nodes())
        net_dist = {}
        for pat, forest in self.forests.items():
            net_dist[pat] = self.jaccard_distance(set(forest.nodes()), g_nodes)
        net_df = pd.DataFrame(net_dist.items(), columns=['net2', 'distance'])
        net_df['net1_type'] = 'forest'
        net_df['net2_type'] = 'forest'
        print(net_df)
        return net_df

    def distance_to_communities(self, g_query):
        """computes distance between entry graph and communities in hypha"""
        print('Computing difference from graph to all communities')
        g_nodes = set(g_query.nodes())
        comm_dist = {}
        for comm, nodes in self.communities.items():
            comm_dist[comm] = self.jaccard_distance(set(nodes), g_nodes)
        comm_df = pd.DataFrame(comm_dist.items(), columns=['net2', 'distance'])
        comm_df['net2_type'] = 'community'
        comm_df['net1_type'] = 'forest'
        print(comm_df)
        return comm_df

    def within_distances(self):
        distvals = []
        for pat, forest in self.forests.items():
            net_net_dist = self.distance_to_networks(forest)
            net_net_dist['net1'] = pat
            #now do distance to hypha
            net_hyph_dist = self.distance_to_communities(forest)
            net_hyph_dist['net1'] = pat
            distvals.append(net_net_dist)
            distvals.append(net_hyph_dist)
        return pd.concat(distvals)

    def inter_distance(self, hyp2):
        """ Computes distance to other hypha"""
        #first compute networks to other networks
        df_list = []
        for pat, forest in hyp2.forests.items():
            net_dist = self.distance_to_networks(forest)
            comm_dist = self.distance_to_communities(forest)
            comm_dist['net1'] = pat
            net_dist['net1'] = pat
            df_list.append(comm_dist)
            df_list.append(net_dist)
        return pd.concat(df_list)


        #the compute networks to communities

        #the compute communities to communities
