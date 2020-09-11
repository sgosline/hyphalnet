# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:05:29 2020

@author: gosl241
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:05:29 2020
@author: gosl241
"""

import argparse
import sys

import hyphalnet.hypha as hyp
from hyphalnet.hypha import hyphalNetwork
import hyphalnet.proteomics as prot
import hyphalnet.hyphEnrich as hyEnrich
import hyphalnet.hyphaeStats as hyStats
import pandas as pd
import synapseclient
from synapseclient import File, Table, table
import pickle
from scipy import stats

syn = synapseclient.Synapse()
syn.login()
data = syn.tableQuery("SELECT * FROM syn22172602").asDataFrame()


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
parser.add_argument('--enrich', dest='doEnrich', action='store_true',\
                    default=False, help='Flag to do GO enrichment')
parser.add_argument('--saveGraphs', dest='toFile', action='store_true',\
                    default=False, help='Flag to save networks to file')
parser.add_argument('--getDistances', dest='getDist', action='store_true',\
                    default=False, help='Get and save distances')
parser.add_argument('--fromFile', dest='fromFile', nargs=1,\
                    action=kvdictAppendAction,metavar='KEY=VALUE',\
                    help='Key/value params for extra files')


def tumor_genes(data_frame, group, subgroup, value):
    significant_mut = data_frame[abs(data_frame[value]) > 0]
    mutation_dictionary = (significant_mut.groupby(group).apply(lambda x: dict(zip(x[subgroup], x[value]))).to_dict())
    return mutation_dictionary

# alpha / 2 = 0.01 /2 = 0.005 the zalpha/2 = 2.58
# Reject Ho when Z <= -2.58 or Z >= 2.58
 #Calculating z score

def significant_prots(data_frame, group, subgroup, value,quantThresh=0.01):
    pquants = pd.DataFrame({'thresh':data_frame.groupby(group)[value].quantile(1.0-quantThresh)})
    tdat = data_frame.merge(pquants, on=group)
    tdat = tdat.assign(topProt=tdat[value] > tdat['thresh'])
    significant = tdat[tdat['topProt']]

    #data_frame['zscore'] = stats.zscore(data_frame[value])
    #significant = data_frame[data_frame['zscore'] >= 2.58]
    gene_dictionary = (significant.groupby(group).apply(lambda x: dict(zip(x[subgroup], x[value]))).to_dict())
    return gene_dictionary

def store_network_files(parent_id='syn22269875'):
    '''
    store_network_files takes the files and uploads them to synapse
    '''


def loadFromFile(file_name_dict):
    hyphae = dict()
    for key, fname in file_name_dict.items():
        hyphae[key] = hyp.load_from_file(fname)
    return hyphae

def main():

    gfile = '../../data/igraphPPI.pkl'
    g = pickle.load(open(gfile, 'rb'))#hyp.make_graph_from_dict(gfile)

    args = parser.parse_args()
    beta = 0.5
    proteomics_dictionary = significant_prots(data, 'AML sample', 'Gene', 'LogFoldChange')
    gene_dictionary = tumor_genes(data, 'AML sample', 'Gene', 'Tumor VAF')
    if args.fromFile is None:
        hyphae = dict()
        hyphae['mutations'] = hyphalNetwork(gene_dictionary, g.copy(), beta)
        hyphae['proteomics'] = hyphalNetwork(proteomics_dictionary, g.copy(), beta)
        for key, this_hyp in hyphae.items():
            this_hyp._to_file(key+'_amlPatientData_hypha.pkl')
    else:
        hyphae = loadFromFile(args.fromFile)

        #now compute graph distances to ascertain fidelity
    if args.getDist:
        res = hyStats.compute_all_distances(hyphae)
        res.to_csv('amlNetworkdistances.csv')
        tab = table.build_table("AML Network Distances", 'syn22128879',res)
        syn.store(tab)
        nmi = hyStats.compute_all_nmi(hyphae, g)
        nmi.to_csv('amlNMI.csv')
        syn.store(File('amlNMI.csv',parent='syn22269875'))
        #store distances
    for key, this_hyp in hyphae.items():
        node_stats = this_hyp.node_stats()
        node_stats.to_csv(key+'_nodelist.csv')
        tab = table.build_table("AML Network Nodes", 'syn22128879', node_stats)
        syn.store(tab)
        if args.doEnrich:
            if len(this_hyp.forest_enrichment) == 0:
                for_e = hyEnrich.go_enrich_forests(this_hyp)#SG, ncbi)
                this_hyp.assign_enrichment(for_e, type='forest')
                for_e.to_csv(key+'enrichedForestGoTerms.csv')
                syn.store(File(key+'enrichedForestGoTerms.csv',parent='syn22269875'))
                this_hyp._to_file(key+'_amlPatientData_hypha.pkl')
            if len(this_hyp.community_enrichment) == 0:
                com_e = hyEnrich.go_enrich_communities(this_hyp)
                this_hyp.assign_enrichment(com_e, type='community')
                com_e.to_csv(key+'enrichedCommunityGOterms.csv')
                syn.store(File(key+'enrichedCommunityGOterms.csv',parent='syn22269875'))
                this_hyp._to_file(key+'_amlPatientData_hypha.pkl')
            ##next: compare enrichment between patients mapped to communities
        this_hyp.community_stats(prefix=key).to_csv(key+'_communityStats.csv')
        this_hyp.forest_stats().to_csv(key+'_TreeStats.csv')
        for files in [key+'_amlPatientData_hypha.pkl',key+'_communityStats.csv',key+'_TreeStats.csv']:
            syn.store(File(files,parent='syn22269875'))

if __name__ == '__main__':
    main()
