# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:05:29 2020

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
import synapseclient
from scipy import stats

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
parser.add_argument('--fromFile',dest='fromFile', nargs=1,\
                    action=kvdictAppendAction,metavar='KEY=VALUE',\
                    help='Key/value params for extra files')

gfile='../../data/pcstDictPPI.pkl'

syn = synapseclient.Synapse()
syn.login()
data = syn.tableQuery("SELECT * FROM syn22172602").asDataFrame()


#pull CPTAC data to dictionary of dictionaries 
def dictionary_creator(data_frame, group,subgroup,value):
    genes_dict = {k: f.groupby(subgroup)[value].apply(list).to_dict() for k, f in data_frame.groupby(group)}
    return genes_dict


# alpha / 2 = 0.01 /2 = 0.005 the zalpha/2 = 2.58
# Reject Ho when Z <= -2.58 or Z >= 2.58
 #Calculating z score threshold for significant genes
def significant_genes(data_frame, group, subgroup, value):
    data_frame['zscore'] = stats.zscore(data_frame[value])
    significant = data[abs(data['zscore']) >= 2.58]
    gene_dictionary = dictionary_creator(significant, group, subgroup, value)
 #   return gene_dictionary
    
    g = hyp.make_graph_from_dict(gfile)
    
    hyphae = dict()
    #for key in significant_genes:
    this_hyp = hyphalNetwork(significant_genes, g)
    hyphae['protein'] = this_hyp
    this_hyp._to_file(key+'_hypha.pkl')
    return hyphae


def loadFromFile(file_name_dict):
    hyphae = dict()
    for key, fname in file_name_dict.items():
        hyphae[key] = hyp.load_from_file(fname)
    return hyphae

def main():
    args = parser.parse_args()

    if args.fromFile is None:
        hyphae = significant_genes(data, 'AML sample', 'Gene', 'LogFoldChange')
    else:
        hyphae = loadFromFile(args.fromFile)

    for key, this_hyp in hyphae.items():
        this_hyp.node_stats().to_csv(key+'_nodelist.csv')
        if args.doEnrich:
            if len(this_hyp.forest_enrichment)==0:
                for_e = hyEnrich.go_enrich_forests(this_hyp)
                this_hyp.assign_enrichment(for_e, type='forest')
                for_e.to_csv(key+'enrichedForestGoTerms.csv')
                this_hyp._to_file(key+'_hypha.pkl')
            if len(this_hyp.community_enrichment)==0:
                com_e = hyEnrich.go_enrich_communities(this_hyp)
                this_hyp.assign_enrichment(com_e, type='community')
                this_hyp._to_file(key+'_hypha.pkl')
                com_e.to_csv(key+'enrichedCommunityGOterms.csv')
            ##next: compare enrichment between patients mapped to communities
        this_hyp.forest_stats().to_csv(key+'_TreeStats.csv')
        this_hyp.community_stats(prefix=key).to_csv(key+'_communityStats.csv')

    #now compute graph distances to ascertain fidelity
    if args.getDist:
        res = hyStats.compute_all_distances(hyphae)
        res.to_csv('proteomicDistances.csv')
        nmi = hyStats.compute_all_nmi(hyphae, gfile)
        nmi.to_csv('proteomicNMI.csv')

if __name__ == '__main__':
    main()
