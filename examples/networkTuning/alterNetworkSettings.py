# -*- coding: utf-8 -*-
"""
The goal of this script is to alter the network algorithm and assess how to
obtain the communities that best span the forests.

@author: gosl241
"""

import argparse
import sys
sys.path.insert(0, "../../")

import hyphalnet.hypha as hyp
from hyphalnet.hypha import hyphalNetwork
import hyphalnet.proteomics as prot
import hyphalnet.hyphEnrich as hyEnrich
import hyphalnet.hyphaeStats as hyStats
import pandas as pd
from pcst_fast import *
from igraph import *
import leidenalg as la

class kvdictAppendAction(argparse.Action):
    """
    argparse action to split an argument into KEY=VALUE form
    on the first = and append to a dictionary.
    """
    def __call__(self, parser, args, values, option_string=None):
        assert(len(values) == 1)
        values = values[0].split(',')
        for val in values:
            try:
                (k, v) = val.split("=", 2)
            except ValueError as ex:
                raise argparse.ArgumentError(self, \
                                             "could not parse argument \"{values[0]}\" as k=v format")
            d = getattr(args, self.dest) or {}
            d[k] = v
            setattr(args, self.dest, d)


#Parser information for command line
parser = argparse.ArgumentParser(description="""Get data from the proteomic \
                                 data commons and build community networks""")
parser.add_argument('--saveGraphs', dest='toFile', action='store_true',\
                    default=False, help='Flag to save networks to file')
parser.add_argument('--getDistances', dest='getDist', action='store_true',\
                    default=False, help='Get and save distances')
parser.add_argument('--fromFile',dest='fromFile', nargs=1,\
                    action=kvdictAppendAction,metavar='KEY=VALUE',\
                    help='Key/value params for extra files')


def build_hyphae_from_data():
    """ Temp function to load data from local directory"""
    ##this is the framework for the PDC data parser.
    norms = prot.normals_from_manifest('data/PDC_biospecimen_manifest_05112020_184928.csv')

#    bcData = prot.parsePDCfile('data/TCGA_Breast_BI_Proteome.itraq.tsv')
    bcData = prot.parsePDCfile('data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv')
    lungData = prot.parsePDCfile('data/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv')
    colData = prot.parsePDCfile('data/CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv')
    gbmData = prot.parsePDCfile('data/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv')

    normPats = {'brca': set([a for a in bcData['Patient'] if a in norms['Breast Invasive Carcinoma']]),\
                'coad': set([a for a in colData['Patient'] if a in norms['Colon Adenocarcinoma']]),\
                'luad': set([a for a in lungData['Patient'] if a in norms['Lung Adenocarcinoma']]),\
                'gbm': set([a for a in gbmData['Patient'] if a in norms['Other']])}

    g = hyp.make_graph_from_dict(gfile)
    namemapper = None #hyp.mapHGNCtoNetwork()

    ##here we get the top values for each patient
    patVals = {'brca':prot.getProtsByPatient(bcData, namemapper),\
               'luad':prot.getProtsByPatient(lungData, namemapper),\
             'coad':prot.getProtsByPatient(colData, namemapper),\
             'gbm':prot.getProtsByPatient(gbmData, namemapper)}

    #here we get the top most distinguished from normals
    patDiffs = {'brca': prot.getTumorNorm(bcData, normPats['brca'], namemapper),
                'luad': prot.getTumorNorm(lungData, normPats['luad'], namemapper),
                'coad': prot.getTumorNorm(colData, normPats['coad'], namemapper),
                'gbm': prot.getTumorNorm(gbmData, normPats['gbm'], namemapper)}
    #now we want to build network communities for each
    hyphae = dict()

    for key in patDiffs:
        this_hyp = hyphalNetwork(patDiffs[key], g)
        hyphae[key] = this_hyp
        this_hyp._to_file(key+'_hypha.pkl')
    return hyphae

def runCommunityWithMultiplex(forests,node_counts,vpm,ga):
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
    [all_nodes.update(ig.vs['name']) for ig in forests.values()]
    node_weights = [node_counts[n] for n in all_nodes]
    print("Have", len(all_nodes), 'total nodes')
    g = ga.copy()
    g = g.subgraph(g.select(all_nodes))
    for nx_g in forests.values(): ##i'm not convince the forests have the same node/edge indices
        tmp_g = Graph()#nx_g.copy() ###this is apointer, make sure to copy!!
        tmp_g.add_vertices([a for a in all_nodes])#[a for a in all_nodes.difference(set(nx_g.vs['name']))]) ##need to a in missing vertsxs
        tmp_g.vs['node_size']= node_weights
        for ed in nx_g.es:
            eps = nx_g.vs.select(ed.tuple)['name']
            tmp_g.add_edge(eps[0], eps[1],weights=ed['weights'])

        print("Graph now has", len(tmp_g.es), 'edges and', len(tmp_g.vs), 'nodes')
        ##compare to copy approach
        other_tmp = nx_g.copy()
        other_tmp.add_vertices([a for a in all_nodes.difference(set(nx_g.vs['name']))])
        print("Copied graph has", len(other_tmp.es), 'edges and', len(other_tmp.vs), 'nodes')
        #print("Orig graph still has", len(nx_g.vs), 'nodes')
        netlist.append(tmp_g)#self._nx2igraph(nx_g)) ##forest is already igraph
    [membership, improv] = la.find_partition_multiplex(netlist,vpm)
                                                          # la.CPMVertexPartition)
    #                                                           la.ModularityVertexPartition)
    print('Found',max(membership),'communities, now refining')
    membership = vpm(g,initial_membership=membership, weights='weight')
    print("NOw have",max(membership),'communiies')
    comm_df = pd.DataFrame({'Node': list(node_counts.keys()),\
                            'Community': membership})
    comm_counts = comm_df.groupby("Community")['Node'].count()
    comm_dict = dict(comm_df.groupby('Community')['Node'].apply(list))
    red_list = [comm_dict[c] for c in comm_counts.index[comm_counts > 5]]
    red_dict = dict(zip(comm_counts.index[comm_counts > 5], red_list))
    for comm, vals in comm_dict.items():
        print("Community", comm, 'has', len(vals), 'nodes')
    return comm_dict

def getFastForest(nodeweights, interactome):
    """
        uses the pcst_fast package to build a weighted subgraph
        inferring nodes and returning a larger subgraph
        """
    #map nodes to indices in some base file
    edges = interactome['edges']
    nodes = interactome['nodes']
    cost = interactome['cost']
    eweights = interactome['edgeWeights']
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
    enodes = []
    node_vals = []
    for v in vert:
        enodes.append(nodes[v])
        if nodes[v] in nodeweights.keys():
            node_vals.append(nodeweights[nodes[v]])
        else:
            node_vals.append(0.0)
    edge_i = [edges[e] for e in edge]
    all_e = [[nodes[e[0]], nodes[e[1]]] for e in edge_i]

    gr.add_vertices(enodes) ##ADD in all nodes! ## add in node_sizes
    gr.vs['node_size'] = node_vals
    gr.add_edges(all_e) ##add in weights
    gr.es['weights'] = [eweights[e] for e in edge]
    ##now add in edge names
    for edge in gr.es:
        edge['name'] = '_'.join(gr.vs.select(edge.tuple)['name'])

    ##now remove the zero-vertex nodes
    term = [e for e in enodes if e in nodeweights.keys()]
    print("Created tree from", len(nodeweights), 'proteins with',\
          len(gr.vs), 'nodes (', len(term), 'terminals) and', len(gr.es), 'edges')
    return gr


def loadFromFile(file_name_dict):
    hyphae = dict()
    for key, fname in file_name_dict.items():
        hyphae[key] = hyp.load_from_file(fname)
    return hyphae



norms = prot.normals_from_manifest('../cancerProtFeatures/data/PDC_biospecimen_manifest_05112020_184928.csv')

#    bcData = prot.parsePDCfile('data/TCGA_Breast_BI_Proteome.itraq.tsv')
bcData = prot.parsePDCfile('../cancerProtFeatures/data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv')
lungData = prot.parsePDCfile('../cancerProtFeatures/data/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv')
colData = prot.parsePDCfile('../cancerProtFeatures/data/CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv')
gbmData = prot.parsePDCfile('../cancerProtFeatures/data/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv')

namemapper = None

normPats = {'brca': set([a for a in bcData['Patient'] if a in norms['Breast Invasive Carcinoma']]),\
                'coad': set([a for a in colData['Patient'] if a in norms['Colon Adenocarcinoma']]),\
                'luad': set([a for a in lungData['Patient'] if a in norms['Lung Adenocarcinoma']]),\
                'gbm': set([a for a in gbmData['Patient'] if a in norms['Other']])}

   #here we get the top most distinguished from normals
patDiffs = {'brca': prot.getTumorNorm(bcData, normPats['brca'], namemapper),
            'luad': prot.getTumorNorm(lungData, normPats['luad'], namemapper),
            'coad': prot.getTumorNorm(colData, normPats['coad'], namemapper),
            'gbm': prot.getTumorNorm(gbmData, normPats['gbm'], namemapper)}

g = hyp.make_graph_from_dict(gfile)


forests = {}
node_counts={}
for pat in list(patDiffs['brca'].keys())[0:10]:
    print('Running PCSF for sample '+pat)
    #            forest = self._getForest(list(proteinWeights[pat].items()))
    fast_forest = getFastForest(patDiffs['brca'][pat],g)
    forests[pat] = fast_forest
    for node in fast_forest.vs['name']:
        if node in list(node_counts.keys()):
            node_counts[node] += 1
        else:
            node_counts[node] = 1


mod = runCommunityWithMultiplex(forests,node_counts,la.ModularityVertexPartition,g)
#cpm = runCommunityWithMultiplex(forests,node_counts,la.CPMVertexPartition,g)
rand = runCommunityWithMultiplex(forests,node_counts,la.RBERVertexPartition,g)
sig = runCommunityWithMultiplex(forests,node_counts,la.SignificanceVertexPartition,g)


##plot community size

def analysis():

    for key, this_hyp in hyphae.items():
        this_hyp.node_stats().to_csv(key+'_nodelist.csv')
        this_hyp.forest_stats().to_csv(key+'_TreeStats.csv')
        this_hyp.community_stats(prefix=key).to_csv(key+'_communityStats.csv')

    #now compute graph distances to ascertain fidelity
    if args.getDist:
        res = hyStats.compute_all_distances(hyphae)
        res.to_csv('panCancerDistances.csv')
        nmi = hyStats.compute_all_nmi(hyphae, gfile)
        nmi.to_csv('panCancerNMI.csv')

if __name__ == '__main__':
    main()
