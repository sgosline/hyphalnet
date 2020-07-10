from igraph import *
import leidenalg as la
import pandas as pd

from pcst_fast import *
import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plot
from sklearn.metrics import confusion_matrix, normalized_mutual_info_score

matplotlib.rcParams['pdf.fonttype'] = 42

#import plotly.express as px

class hyphalNetwork:
    """
    The hypha class represents a set of individual sample networks/graphs and
    creates the communities shared between each
    """
    def __init__(self, proteinWeights, interactome):
        """
        Hypha class
        Based a dictionary of proteins with weights and an interactome,
        derives individual protein networks for each dictionary and then
        finds the communities shared between them
        """

        self.proteins = proteinWeights #dictionary of protein weights for each forest
        self.interactome = interactome #interactome
        self.forests = dict()#forests
        self.community_enrichment = {} # enriched terms
        self.forest_enrichment = dict() #enriched terms
        self.node_counts = dict() #nodes found in forests and their counts
        for pat in list(proteinWeights.keys()):
            print('Running PCSF for sample '+pat)
            #            forest = self._getForest(list(proteinWeights[pat].items()))
            fast_forest = self._getFastForest(proteinWeights[pat])
            self.forests[pat] = fast_forest
            for node in fast_forest.vs['name']:
                #list(forest.nodes()):
                if node in list(self.node_counts.keys()):
                    self.node_counts[node] += 1
                else:
                    self.node_counts[node] = 1
        self.communities = self.runCommunityWithMultiplex()
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

    def _get_subnet(self, nodenames):
        """
        Searches for nodes in the interactome that have nodes of a specific name,
        then returns inferred subgraph, nothing fancy
        TODO: deprecate this if we no longer need
        """
        nodelist = []
        for nn in nodenames:
            try:
                nodelist.append(self.interactome.vs.find(name = nn)) #TODO: change from OI
            except ValueError:
                print('Node '+nn+' not found')
        sub_g = self.interactome.subgraph(nodelist)
        print("Subgraph has", len(sub_g.vs), "nodes")
        return sub_g

    def _nx2igraph(self, nx_graph):
        """
        Helper function that maps networkx object to igraph
        TODO: deprecrate
        """
        new_g = Graph(directed=False)
        alln = self.node_counts.keys()
        for n in alln:
            new_g.add_vertex(n)
        for e in nx_graph.edges():
            new_g.add_edge(e[0], e[1])
        return new_g

    def _getForest(self, nodeweights):
        """
        Uses the omics integrator package to build a weighted subgrpah, inferring
        additional nodes and returning a larger subgraph on the interatome
        DEPRECATED
        """
        #map nodes interactome
        pdf = pd.DataFrame(nodeweights)
        pdf.columns = ['name', 'prize']
        graph = self.interactome
        graph._prepare_prizes(pdf)
        verts, edges = graph.pcsf()
        forest, augmented_forest = graph.output_forest_as_networkx(verts, edges)
        return forest #TODO: add in parameter shift for 0 size networks

    def _getFastForest(self, nodeweights):
        """
        uses the pcst_fast package to build a weighted subgraph
        inferring nodes and returning a larger subgraph
        """
        #map nodes to indices in some base file
        edges = self.interactome['edges']
        nodes = self.interactome['nodes']
        cost = self.interactome['cost']
       # print(nodeweights)
        weights = []
        for n in nodes:
            #print(n)
            if n in nodeweights.keys():
                weights.append(nodeweights[n])
            else:
                weights.append(0.0)
        vert, edge = pcst_fast(edges, weights, cost, -1,1,'gw',0)
        #return as igraph
        gr = Graph()
        #this takes mapping the edges back to the original node values which can
        #be compared across graphs
        enodes = [nodes[v] for v in vert]
        edge_i = [edges[e] for e in edge]
        all_e = [[nodes[e[0]],nodes[e[1]]] for e in edge_i]

        gr.add_vertices(enodes) ##ADD in all nodes!
        gr.add_edges(all_e)
        ##now remove the zero-vertex nodes
        term=[e for e in enodes if e in nodeweights.keys()]
        print("Created tree from", len(nodeweights), 'proteins with',\
              len(gr.vs), 'nodes (',len(term),'terminals) and', len(gr.es), 'edges')
        return gr

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
        [all_nodes.update(ig.vs['name']) for ig in self.forests.values()]
        print("Have", len(all_nodes), 'total nodes')
        for nx_g in self.forests.values():
            tmp_g = nx_g.copy() ###this is apointer, make sure to copy!!
            tmp_g.add_vertices([a for a in all_nodes.difference(set(nx_g.vs['name']))]) ##need to a in missing vertsxs
            print("Graph now has", len(tmp_g.es), 'edges and', len(tmp_g.vs), 'nodes')
            print("Orig graph still has",len(nx_g.vs),'nodes')
            netlist.append(tmp_g)#self._nx2igraph(nx_g)) ##forest is already igraph
        [membership, improv] = la.find_partition_multiplex(netlist,\
                                                           la.ModularityVertexPartition)
        comm_df = pd.DataFrame({'Node': list(self.node_counts.keys()),\
                           'Community': membership})
        comm_counts = comm_df.groupby("Community")['Node'].count()
        comm_dict = dict(comm_df.groupby('Community')['Node'].apply(list))
        red_list = [comm_dict[c] for c in comm_counts.index[comm_counts > 1]]
        red_dict = dict(zip(comm_counts.index[comm_counts > 1], red_list))
        for comm, vals in red_dict.items():
            print("Community", comm, 'has', len(vals), 'nodes')
        return red_dict

    def node_stats(self):
        """
        For each node computes how many trees it is in, and which community it is in
        """
        comms = []
        for no, vals in self.node_counts.items():
            ndict = {'Node': no,'NumForests': vals}
            for comm, nodelist in self.communities.items():
                if no in nodelist:
                    ndict['Community'] = comm
            comms.append(ndict)
        return pd.DataFrame(comms)

    def community_stats(self, prefix=''):
        """
        computes stats about community and returns data frame
        These statistics include:
        - number of nodes in each community
        - number of forests that are 'closest' to each community
        - number of terms enriched for each community
        - will create graphML file to be loaded into cytoscape if prefix is not empty
        - NMI between communities and original graph
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
        if type=='community':
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
        for samp in self.forests.keys():
            res = d_val[(d_val.net1 == samp) & (d_val.net2_type == 'community')].nsmallest(1, 'distance').index[0]
            co = d_val.at[res, 'net2']
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
        Currently plots community to graph file

        Returns
        -------
        None.
        """
        print('Creating graph with community annotations')
        #first build entire graph
        gr = getIgraph(self.interactome)
        #now reduce graph to only those nodes in the community
        rednodes = set()
        [rednodes.update(cn) for cn in self.communities.values()]
        gred = gr.subgraph(rednodes)
        print("Reducing full interactome of",len(gr.es),\
              "edges and",len(gr.vs),"nodes to one with",len(gred.es),\
              "edges and",len(gred.vs),"nodes")
        node_comm={}
        for comm,nodelist in self.communities.items():
            for n in nodelist:
                node_comm[n]=comm
        #nodes = self.communities[commName]
        #cred = gred.subgraph(list(nodes))
        nc = {node:self.node_counts[node] for node in list(rednodes)}
        #nx.set_node_attributes(cred, nc, name='Number of forests')
        comms=[]
        ntrees=[]
        for n in gr.vs['name']:
            if n in node_comm.keys():
                comms.append(node_comm[n])
            else:
                comms.append('None')
            if n in nc.keys():
                ntrees.append(nc[n])
            else:
                ntrees.append(0)

        gred.vs['Community'] = comms
        gred.vs['NumTrees'] = ntrees
        print('Adding membership to network option')
        #nx.set_node_attributes(gred,membership,'Community')
        gred.write_graphmlz(prefix+'_communityGraph.graphml.gz')
        return gred

    def distance_to_networks(self, g_query):
        """
        Network distance helper function
        Takes a query network and computes distance to all
        networks within hypha
        """
        g_nodes = set(g_query.vs['name'])
        net_dist = {}
        for pat, forest in self.forests.items():
            net_dist[pat] = jaccard_distance(set(forest.vs['name']), g_nodes)
        net_df = pd.DataFrame(net_dist.items(), columns=['net2', 'distance'])
        net_df['net1_type'] = 'forest'
        net_df['net2_type'] = 'forest'
        return net_df

    def distance_to_communities(self, g_query):
        """
        computes distance between entry graph and communities in hypha
        """
        #print('Computing difference from graph to all communities')
        g_nodes = set(g_query.vs['name'])
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
        return pd.concat(distvals)

    def intra_distance(self, hyp2):
        """ Computes distance to other hypha"""
        #first compute networks to other networks
        df_list = []
        for pat, forest in hyp2.forests.items():
            #iterate through each element in the other hypha
            net_dist = self.distance_to_networks(forest)
            comm_dist = self.distance_to_communities(forest)
            comm_dist['net1'] = pat
            net_dist['net1'] = pat
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

    gr = Graph()
    gr.add_vertices(nodes)
    gr.add_edges(edges)
    gr.es['weight'] = weight
    return gr


def make_graph_from_dict(gfile):
    """
   Makes graph from dictionary of nodes and edges
    """
    print("Loading edge and node lists from dictionary")
    gd = pickle.load(open(gfile,'rb'))
#    print(gd.keys())
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

def communityFromGraphFile(gfile):
    """
    Loads graph from edge list pKL file and calculates communities
    """
    dfile = make_graph_from_dict(gfile)
    ig = getIgraph(dfile)
    partition = la.find_partition(ig, la.ModularityVertexPartition)
    comm_dict={}
    for p in range(len(partition)):
        comm_dict[p]=ig.vs.select(partition[p])['name']
    print("Found",len(comm_dict),'communities for the primary interactome')
    return comm_dict
