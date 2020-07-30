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


gfile = '../../data/igraphPPI.pkl'
g = pickle.load(open(gfile, 'rb'))#hyp.make_graph_from_dict(gfile)



#Convert any dataframe to nested dict
def nested_dict(df):
    #Check for number of columns
    if len(df.columns) == 1:
        if df.values.size == 1: 
            return [df.values[0][0]]
        return df.values.squeeze().tolist()
    grouped = df.groupby(df.columns[0])
    d = {k: nested_dict(g.iloc[:,1:]) 
         for k,g in grouped}
    return d

# alpha / 2 = 0.01 /2 = 0.005 the zalpha/2 = 2.58
# Reject Ho when Z <= -2.58 or Z >= 2.58
 #Calculating z score

def significant_genes(data_frame, group, subgroup, value):
    data_frame['zscore'] = stats.zscore(data_frame[value])
    significant = data_frame[data_frame['zscore'] >= 2.58]
    gene_dictionary = (significant.groupby(group).apply(lambda x: dict(zip(x[subgroup], x[value]))).to_dict())


    hyphae = dict()
    beta=0.5
    key = 'proteomics'
    this_hyp = hyphalNetwork(gene_dictionary, g.copy())
    hyphae[key] = this_hyp
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

        #now compute graph distances to ascertain fidelity
    if args.getDist:
        res = hyStats.compute_all_distances(hyphae)
        res.to_csv('proteomicdistances.csv')
        nmi = hyStats.compute_all_nmi(hyphae, g)
        nmi.to_csv('proteomic.csv')
    for key, this_hyp in hyphae.items():
        this_hyp.node_stats().to_csv(key+'_nodelist.csv')
        if args.doEnrich:
            if len(this_hyp.forest_enrichment) == 0:
                for_e = hyEnrich.go_enrich_forests(this_hyp)#SG, ncbi)
                this_hyp.assign_enrichment(for_e, type='forest')
                for_e.to_csv(key+'enrichedForestGoTerms.csv')
                this_hyp._to_file(key+'_hypha.pkl')
            if len(this_hyp.community_enrichment) == 0:
                com_e = hyEnrich.go_enrich_communities(this_hyp)
                this_hyp.assign_enrichment(com_e, type='community')
                this_hyp._to_file(key+'_hypha.pkl')
                com_e.to_csv(key+'enrichedCommunityGOterms.csv')
            ##next: compare enrichment between patients mapped to communities
        this_hyp.community_stats(prefix=key).to_csv(key+'_communityStats.csv')
        this_hyp.forest_stats().to_csv(key+'_TreeStats.csv')

if __name__ == '__main__':
    main()