from igraph import *
import leidenalg as la
import pandas as pd
import OmicsIntegrator as oi
import networkx as nx
import goenrich

class hypha:
    """hypha class defines multiple networks and communities """
    def __init__(self, proteinWeights, interactome):
        '''Constructor for hypha'''
        self.proteins = proteinWeights #dictionary of protein weights for each forest
        self.interactome = interactome #interactome
        self.forests = {} #forests
        self.node_set = set() #nodeset across all forests
        self.communities = {} #communities from hypha
        for pat in proteinWeights.keys():
            print('Running PCSF for patient'+pat)
            forest = self._getForest(list(proteinWeights[pat].items()))
            self.forests[pat] = forest
        for net in self.forests.values():
            self.node_set = self.node_set.union(set(net.nodes()))
        print('Ready to create hypha across ', len(self.node_set), 'nodes and',\
              len(proteinWeights), 'forests')

    def make_graph(gfile):
        print('making graph for OmicsIntegrator')
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
       # print("Reducing",len(nodenames),"proteins to",len(nodelist)," nodes")
        sg = self.interactome.subgraph(nodelist)
        print("Subgraph has", len(sg.vs), "nodes")
        return sg

    def _nx2igraph(self, nx_graph):
        g = Graph(directed=False)
        alln = self.node_set
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
        #first we have to create a union/interaction
        print('running community detection')
        optimizer = la.Optimiser()
        netlist = []
        for nx in self.forests.values():
            netlist.append(self._nx2igraph(nx))
        [membership, improv] = la.find_partition_multiplex(netlist, la.ModularityVertexPartition)
        df =pd.DataFrame(list(zip(self.node_set, membership)), columns=['Node','Partition'])
        self.communities = dict(df.groupby('Partition')['Node'].apply(list))
        return df

    def forestStats():
        print("We have "+len(self.forests)+"forests")
        #plot number of nodes in each network with number of edges

    def saveCommunityToFile(self, prefix=''):
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
        for com,nodes in self.communities.items():
            cred = gred.subgraph(list(nodes))
            print('Adding membership to network option')
            #nx.set_node_attributes(gred,membership,'Community')
            print('saving to html')
            oi.output_networkx_graph_as_graphml_for_cytoscape(cred, \
                                                              output_dir=".", \
                                                              filename=prefix+"_"+\
                                                              str(com)+\
                                                              "_pcsf_results.graphml.gz")

            oi.output_networkx_graph_as_interactive_html(cred, \
                                                         output_dir=".", \
                                                         filename=prefix+"_"+\
                                                         str(com)+\
                                                         "_pcsf_results.html")
        return(None)

    def map_hgnc_to_ensPep():
        tab = pd.read_csv('data/human.name_2_string.tsv', '\t', skiprows=[0], header=None)
        res = dict(zip(tab[1], tab[2]))
        return res

    def map_ensPep_to_hgnc():
        tab = pd.read_csv('data/human.name_2_string.tsv', '\t', skiprows=[0], header=None)
        res = dict(zip(tab[2], tab[1]))
        return res

    def get_go_enrichment(genelist, background):
        O = goenrich.obo.ontology('db/go-basic.obo')
        gene2go = goenrich.read.gene2go('db/gene2go.gz')
        gene2go = gene2go.loc[gene2go['GeneID'].isin(background)]
        values = {k: set(v) for k,v in gene2go.groupby('GO_ID')['GeneID']}
        background_attribute = 'gene2go'
        goenrich.enrich.propagate(O, values, background_attribute)
        df = goenrich.enrich.analyze(O, genelist, background_attribute)
        return(df)


    def go_enrich_forests(self, ncbi_mapping):
        """        Iterates through forests and gets enrichment """
        enrich={}
        for pat, forest in self.forests.items():
            print(pat)
            nodenames = []
            background = []
            for fn in forest.nodes():
                try:
                    nodenames.append(ncbi_mapping[fn])
                except KeyError:
                    print("No key", fn)
            for ns in self.node_set:
                try:
                    background.append(ncbi_mapping[ns])
                except KeyError:
                    print("No key", ns)

            try:
                enrich[pat] = get_go_enrichment(nodenames, background)
            except:
                print('error?')
        return(enrich)

    def go_enrich_communities(self, ncbi_mapping):
        """ GEts enrichment for individual communities"""
        print ('Doing community thing')
