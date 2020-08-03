"""
Docstring goes here
"""
import pickle
import os
from igraph import *
import leidenalg as la
import pandas as pd
from pcst_fast import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plot
from sklearn.metrics import confusion_matrix, normalized_mutual_info_score

matplotlib.rcParams['pdf.fonttype'] = 42

class hyphalNetwork:
    """
    The hypha class represents a set of individual sample networks/graphs and
    creates the communities shared between each
    """
    def __init__(self, proteinWeights, interactome, beta=1, g=3, do_forest=False):
        """
        Hypha class
        Based a dictionary of proteins with weights and an interactome,
        derives individual protein networks for each dictionary and then
        finds the communities shared between them
        """
        self.proteins = proteinWeights #dictionary of protein weights for each forest
        ##do some interactome processing
        self.interactome = self._weightByDegree(interactome, g)
        self.idict = self._getPPIdict()
        # Add costs to each edge, proportional to the degrees of the nodes it connects, modulated by parameter g.

        self.forests = dict()#forests
        self.community_enrichment = dict() # enriched terms
        self.forest_enrichment = dict() #enriched terms
        self.orig_weights = proteinWeights
        self.node_counts = dict() #nodes found in forests and their counts
        for pat in list(proteinWeights.keys()):
            print('Building tree for sample '+pat+' with', len(proteinWeights[pat]), 'proteins')
            #            forest = self._getForest(list(proteinWeights[pat].items()))
            nbeta = beta
            fast_forest = self._getFastForest(proteinWeights[pat], nbeta, do_forest)
            while (len(fast_forest['vert']) < 2) and (nbeta < 1000):
                print("Tree only has 1 node, trying to increase beta by factor of 10")
                nbeta = nbeta*10.0
                fast_forest = self._getFastForest(proteinWeights[pat], nbeta, do_forest)
            if len(fast_forest['vert']) == 1:
                print("Could not create tree for", pat)
                next
            print("Built tree of", len(fast_forest['vert']), 'proteins at beta', nbeta)
            self.forests[pat] = fast_forest
            #print(fast_forest['vert'])
            for ni in fast_forest['vert']:
                #fast_forest.vs['name']:
                #list(forest.nodes()):
                node = self.interactome.vs['name'][ni]
                if node in list(self.node_counts.keys()):
                    self.node_counts[node] += 1
                else:
                    self.node_counts[node] = 1
        self.communities, self.comm_graphs = self.runCommunityWithMultiplex()
        ##now we create various statistics to compare communities
        #what is the distance (jaccard) between trees and communities?
        self.distVals = self.within_distances() #compute distances between forests
        #what is the score of communities for each sample?
        self.assignedCommunities = self.community_to_samples()
        print('Created hypha across ', len(self.node_counts), 'nodes and',\
              len(self.forests), 'forests')

    def _to_file(self, fname):
        """ enables saving to file"""
        pickle.dump(self, open(fname, 'wb'))

    def _getPPIdict(self):
        """Converts igraph to dictionary that can be used by `pcst_fast`
        Args:
           None
        Returns:
          A dictionary with `edges`, `nodes`, `weight`, `cost`, and `degree` items
        """
        ig = self.interactome
        print("Saving igraph to dict")
        edges = np.array([e.tuple for e in ig.es])
        pcst_dict = {'edges':edges, 'nodes':ig.vs['name'], 'weight':ig.es['weight'],\
                     'cost':ig.es['cost'], 'degree':ig.degree}
        return pcst_dict

    def _weightByDegree(self, ig, g):
        """ Re-weights edges based on degree and gamma (g) parameter
        Args:
          ig: iGraph graph
          g: Gamma paraameter
        Returns:
          An iGraph object wiht updated 'cost' elements

        """
        print('Adding degree updates to cost')
        N = len(ig.vs)
        node_degrees = ig.degree()
        edge_penalties = np.array([])
        #for e in ig.es:
        #    a,b = e.tuple
        edge_penalties = np.array([(10**g) * node_degrees[e.tuple[0]] * node_degrees[e.tuple[1]] /
                                             ((N - node_degrees[e.tuple[0]] - 1) * \
                                              (N - node_degrees[e.tuple[1]] - 1) + \
                                               node_degrees[e.tuple[0]] * node_degrees[e.tuple[1]]) \
                                              for e in ig.es])

        ig.es['cost'] = 1.0 - np.array(ig.es['weight']) + edge_penalties
        return ig

    def _getFastForest(self, nodeweights, beta, do_forest=False, w=4):
        """
        uses the pcst_fast package to build a weighted subgraph
        inferring nodes and returning a larger subgraph
        """
        #map nodes to indices in some base file
        edges = self.idict['edges']
        nodes = self.idict['nodes']
        cost = self.idict['cost']

        orig_node_count = len(nodes)
        orig_edge_count = len(edges)

        dummy = -1
        weights = np.zeros(orig_node_count)
        if do_forest:
            weights = np.append(weights, 0)
            dummy = len(nodes)

        for n in set(nodes).intersection(nodeweights.keys()):
            ni = nodes.index(n)
            weights[ni] = nodeweights[n]*beta
            if do_forest: ##TODO: move this to igraph!!!
                edges = np.append(edges, [[dummy, ni]], axis=0)
                cost = np.append(cost, w)
        if do_forest:
            print('Running pcst with', len(edges), 'edges compared to original', \
                  orig_edge_count, 'for', len(nodeweights), 'terminals')
        vert, edge = pcst_fast(edges, weights, cost, dummy, 1, 'strong', 0)
        #return as igraph
        return {'vert':vert, 'edge':edge} ##IS IT WORTH RETURNING A GRAPH?

    def runCommunityWithMultiplex(self):
        """
        Primary method for now
        TODO: update with moer features
        Currently takes all forests in hypha (defined by nodes on shared edge set)
        and finds communities
        """
        print('running community detection')
        optimizer = la.Optimiser()
        netlist = []
        all_nodes = set()
        for pat, vals in self.forests.items():
            all_nodes.update([self.interactome.vs['name'][i] for i in vals['vert']])
        print("Have", len(all_nodes), 'total nodes')
        for_graph = self.interactome#.subgraph(all_nodes)
        for nx_g in self.forests.values(): ##i'm not convince the forests have the same node/edge indices
            tmp_g = for_graph.subgraph_edges(nx_g['edge'],\
                                             delete_vertices=False)#Graph()#nx_g.copy() ###this is apointer, make sur
            netlist.append(tmp_g)
        [membership, improv] = la.find_partition_multiplex(netlist,\
                                                           #la.RBERVertexPartition)
                                                           la.ModularityVertexPartition)
        comm_df = pd.DataFrame({'Node': for_graph.vs['name'],\
                                'Community': membership})
        comm_counts = comm_df.groupby("Community")['Node'].count()
        comm_dict = dict(comm_df.groupby('Community')['Node'].apply(list))
        red_list = [comm_dict[c] for c in comm_counts.index[comm_counts > 5]]
        red_dict = dict(zip(comm_counts.index[comm_counts > 5], red_list))
        red_graph = {}
        for comm, vals in red_dict.items():
            rgraph = self.interactome.subgraph(vals)
            red_graph[comm] = rgraph
            print("Community", comm, " graph has", len(vals),'proteins and',\
                  len(rgraph.components()), 'component')
        return red_dict, red_graph

    def node_stats(self):
        """
        For each node computes how many trees it is in, and which community it is in
        """
        comms = []
        all_prots = set()
        for samp, pweights in self.orig_weights.items():
            all_prots.update(pweights.keys())
            for node, val in pweights.items():
                ndict = {'Sample': samp, 'Node': node, 'OrigWeight': val}
                if node in self.node_counts.keys():
                    ndict['NumForests'] = self.node_counts[node]
                    if node in [self.interactome.vs['name'][i] for i in self.forests[samp]['vert']]:
                        ndict['inForest'] = True
                    else:
                        ndict['inForest'] = False
                    for comm, nodelist in self.communities.items():
                        if node in nodelist:
                            ndict['Community'] = comm
                    if 'Community' not in ndict.keys():
                        ndict['Community'] = None
                else:
                    ndict['inForest'] = False
                    ndict['Community'] = None
                    ndict['NumForests'] = None
                comms.append(ndict)
        other_nodes = set(self.node_counts.keys()).difference(all_prots)
        print('Including', len(other_nodes), 'proteins that are not originally weighted')
        for node in other_nodes:
            for samp, forest in self.forests.items():
                if node in [self.interactome.vs['name'][i] for i in forest['vert']]:
                    ndict = {'Sample': samp, 'Node': node, 'OrigWeight': 0.0,\
                             'inForest': True}
                    ndict['NumForests'] = self.node_counts[node]
                else:
                    next
                for comm, nodelist in self.communities.items():
                    if node in nodelist:
                        ndict['Community'] = comm
                if 'Community' not in ndict.keys():
                    ndict['Community'] = None
            comms.append(ndict)
        return pd.DataFrame(comms)

    def community_stats(self, prefix=''):
        """Computes stats about community and returns data frame
        These statistics include:
        - number of nodes in each community
        - number of forests that are 'closest' to each community
        - number of terms enriched for each community
        - will create graphML file to be loaded into cytoscape if prefix is not empty
        - NMI between communities and original graph
        Args:
           prefix: String to use to write hte graph and plot to file
        Returns:
           Data frame with columns describing statistics for each community
        """

        ##get the assigned communities
        for_d = pd.DataFrame(zip(self.assignedCommunities.keys(),\
                                 self.assignedCommunities.values()),\
                             columns=['forest', 'community'])

        num_forests = dict(for_d.groupby('community')['forest'].apply(len))
        stat_list = []
        for comm, nodes in self.communities.items():
            res = {'Community': comm,\
                              'Nodes': len(nodes)}
            if comm in num_forests.keys():
                res['Forests'] = num_forests[comm]
            else:
                res['Forests'] = 0
            if str(comm) in self.community_enrichment.keys():
                res['Enriched terms'] = len(self.community_enrichment[str(comm)])
            else:
                res['Enriched terms'] = 1.0

            stat_list.append(res)
        if prefix != '':
            self.to_graph(prefix=prefix)
        ##plot the bar plot of communities and networks
        pdf = pd.DataFrame(stat_list)
        if prefix != '':
            tp = pdf[['Forests', 'Nodes', 'Enriched terms']]
            tp.set_index(pdf['Community'])
            tp.plot.scatter(x='Nodes', y='Forests', c='Enriched terms', \
                            colormap='viridis', fontsize=12)
            plot.savefig(prefix+'_communityStats.pdf')
        return pdf

    def assign_enrichment(self, enrich_df, type='community'):
        """
        Reads in a data frame containing the go/kegg enrichment and assigns each
        to a community or patient
        """
        if type == 'community':
            e_dict = enrich_df.groupby('Community')['name'].apply(list)
            self.community_enrichment = e_dict
        else:
            e_dict = enrich_df.groupby('Patient')['name'].apply(list)
            self.forest_enrichment = e_dict
       # print(e_dict)

    def community_to_samples(self):
        """
        maps the communities to the closest samples and returns a list of
        forests for each community
        """
        d_val = self.distVals
        closest = {}
        #print(self.forests)
        for samp in self.forests.keys():
            #self.forests.keys():
            #print(samp)
            res = d_val[(d_val.net1 == samp) & (d_val.net2_type == 'community')].nsmallest(1, 'distance').index[0]
            co = list(d_val.at[res, 'net2'])
            #print(samp,co)
            closest.update({samp:co[1]})
        return closest

    def forest_stats(self, forest_names=[], to_file=False):
        """
        Compute statistics for node membership and forests
        """
        print("We have "+str(len(self.forests))+" forests")
        #plot number of nodes in each network with number of edges
        dflist = []
        #create dictionary of all the forests we want
        if len(forest_names) > 0:
            f_dict = {pat:forest for pat, forest in self.forests.items() if pat in forest_names}
        else:
            f_dict = self.forests

        for pat, forest in f_dict.items():
            protweights = pd.DataFrame(self.proteins[pat].items(), \
                                       columns=['Gene', 'DisWeight'])
           # print(nx.get_node_attributes(forest,'pri          ze'))
        #    nodevals = pd.DataFrame(nx.get_node_attributes(forest, 'prize').items(),\
        #                            columns=['Gene', 'ForestPrize'])
            #print(nodevals)
        #    pat_df = protweights.merge(nodevals, on='Gene', how='outer').fillna(0, downcast='infer')
            pat_df = protweights
            pat_df['Patient'] = pat
            if pat in self.forest_enrichment.keys():
                pat_df['Enriched Terms'] = len(self.forest_enrichment[pat])
            else:
                pat_df['Enriched Terms'] = 0.0
            dflist.append(pat_df)
        full_df = pd.concat(dflist)
        return full_df

    def to_graph(self, prefix=''):
        """
        This function interprets the community membership in terms of the underlying
        network and trees that built the network

        Returns
        -------
        igraph Graph object
        """
        print('Creating graph with community annotations')
        try:
            os.mkdir('./'+prefix+'_graphs/')
        except OSError:
            print("Directory already exists")

        graph_list = []
        full_graph = Graph(directed=False)
        for comm, comm_graph in self.comm_graphs.items():
            comm_graph.vs['Community'] = comm
            comm_graph.vs['numForests'] = [self.node_counts[d] for d in comm_graph.vs['name']]
            graph_list.append(comm_graph.to_undirected())
#            new_graph = gred.subgraph(nodes).simplify(combine_edges=dict(weight='max', name='first', numForests='max'))
            print('Saving community', comm, 'with', len(comm_graph.vs), 'nodes and', len(comm_graph.components()), 'components')
            comm_graph.write_graphml(prefix+'_graphs/community'+str(comm)+'_'+prefix+'_graph.graphml')
        try:
            full_graph = full_graph.disjoint_union(graph_list)
            print('Full graph has',len(full_graph.vs),'nodes and',len(full_graph.es),'edges')
        except:
            print("Couldn't create graph")
        return full_graph

    def distance_to_networks(self, g_query):
        """
        Network distance helper function
        Takes a query network and computes distance to all
        networks within hypha
        """
        g_nodes = self.interactome.vs['name']#set(g_query.vs['name'])
        q_nodes = [g_nodes[i] for i in g_query['vert']]
        net_dist = {}
        for pat, forest in self.forests.items():
            net_dist[pat] = jaccard_distance(set([g_nodes[i] for i in forest['vert']]),\
                                             set(q_nodes))
        net_df = pd.DataFrame(net_dist.items(), columns=['net2', 'distance'])
        net_df['net1_type'] = 'forest'
        net_df['net2_type'] = 'forest'
        return net_df

    def distance_to_communities(self, g_query):
        """
        computes distance between entry graph and communities in hypha
        """
        #print('Computing difference from graph to all communities')
        g_nodes = set([self.interactome.vs['name'][i] for i in g_query['vert']])#g_query.vs['name'])
        comm_dist = {}
        for comm, nodes in self.communities.items():
            comm_dist[comm] = jaccard_distance(set(nodes), g_nodes)
        comm_df = pd.DataFrame(comm_dist.items(), columns=['net2', 'distance'])
        comm_df['net2_type'] = 'community'
        comm_df['net1_type'] = 'forest'
        return comm_df

    def within_distances(self):
        """Compute distance between all forests within the hypha"""
        print("Computing differences between forests and communities")
        distvals = []
        for pat, forest in self.forests.items():
            net_net_dist = self.distance_to_networks(forest)
            net_net_dist['net1'] = pat
            #now do distance to hypha
            net_hyph_dist = self.distance_to_communities(forest)
            net_hyph_dist['net1'] = pat
            distvals.append(net_net_dist)
            distvals.append(net_hyph_dist)
       # print(distvals)
        return pd.concat(distvals)

    def intra_distance(self, hyp2):
        """ Computes distance to other hypha"""
        #first compute networks to other networks
        df_list = []
        for pat, forest in hyp2.forests.items():
            #iterate through each element in the other hypha
            net_dist = self.distance_to_networks(forest)
            net_dist['net1'] = pat
            comm_dist = self.distance_to_communities(forest)
            comm_dist['net1'] = pat
            df_list.append(comm_dist)
            df_list.append(net_dist)
#        for comm, nodelist in hyp2.communities.items(): TODO
        return pd.concat(df_list)

def getIgraph(full_graph):
    """
    Helper function that takes the edge list interactome and returns igraph
    """

    #this takes mapping the edges back to the original node values which can
    #be compared across graphs
    edges = full_graph['edges']
    nodes = full_graph['nodes']
    weight = full_graph['edgeWeights']

    gr = Graph(directed=False)
    gr.add_vertices(nodes)
    gr.add_edges(edges)
    for edge in gr.es:
        edge['name'] = '_'.join(gr.vs.select(edge.tuple)['name'])
    gr.es['weight'] = weight
   # print(gr.es.attributes())
   # print("Simplifying")
    gr = gr.simplify(combine_edges=dict(weight='max', name='first'))
   # print(gr.es.attributes())
    return gr


def make_graph_from_dict(gfile):
    """
   Makes graph from dictionary of nodes and edges
    """
    print("Loading edge and node lists from dictionary")
    gd = pickle.load(open(gfile, 'rb'))
    return gd

def load_from_file(filename):
    """Loads hypha object"""
    print("Loading hypha "+filename)
    return pickle.load(open(filename, 'rb'))

def jaccard_distance(ns_1, ns_2):
    """Computes jaccard distance between two networkx objects"""
    #print('Computing distance between graphs by jaccard')
    u_size = len(ns_1.union(ns_2))
    if u_size == 0:
        return 1.0
    return 1-len(ns_1.intersection(ns_2))/u_size

def map_hgnc_to_ensPep():
    tab = pd.read_csv('data/human.name_2_string.tsv', '\t', skiprows=[0], header=None)
    res = dict(zip(tab[1], tab[2]))
    return res

def map_ensPep_to_hgnc():
    tab = pd.read_csv('data/human.name_2_string.tsv', '\t', skiprows=[0], header=None)
    res = dict(zip(tab[2], tab[1]))
    return res

def computeCommunityNMI(comm_dict1, comm_dict2):
    '''
    Computes the normalized mutual information metric for two community assignments
    Used to compare two communities
    '''
    c1_vals = dict()
    c2_vals = dict()
    for comm, gene in comm_dict1.items():
        for g in gene:
            c1_vals[g] = comm
    for comm, gene in comm_dict2.items():
        for g in gene:
            c2_vals[g] = comm

    all_genes = set(list(c1_vals.keys())).intersection(set(list(c2_vals.keys())))

    c1_vec = [c1_vals[g] for g in all_genes]
    c2_vec = [c2_vals[g] for g in all_genes]

    nmi = normalized_mutual_info_score(c1_vec, c2_vec)
    print(nmi)
    return 1-nmi

def communityFromGraph(ig):
    """
    Loads graph from edge list pKL file and calculates communities
    """
    #dfile = make_graph_from_dict(gfile)
    #ig = pickle.load(open(dfile, 'rb'))#getIgraph(dfile)
    partition = la.find_partition(ig, la.ModularityVertexPartition)
    comm_dict = {}
    for p in range(len(partition)):
        comm_dict[p] = ig.vs.select(partition[p])['name']
    print("Found", len(comm_dict), 'communities for the primary interactome')
    return comm_dict
