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
import pickle
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
parser.add_argument('--fromFile',dest='fromFile', nargs=1,\
                    action=kvdictAppendAction,metavar='KEY=VALUE',\
                    help='Key/value params for extra files')


gfile='../../data/pcstDictPPI.pkl'


syn = synapseclient.Synapse()
syn.login()
data = syn.tableQuery("SELECT * FROM syn22172602").asDataFrame()


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
    significant = data[abs(data['zscore']) >= 2.58]
    #replaced with new function
    gene_dictionary = nested_dict(significant[[group, subgroup, value]])
    
    g = hyp.make_graph_from_dict(gfile)
    
    hyphae = dict()
    beta=0.5
    for key, val in gene_dictionary.items():
        this_hyp = hyphalNetwork(val, g.copy())
        hyphae[key] = this_hyp
        this_hyp._to_file(key+'_hypha.pkl')
    print (hyphae)
    return hyphae


def loadFromFile(file_name_dict):
    hyphae = dict()
    for key, fname in file_name_dict.items():
        hyphae[key] = hyp.load_from_file(fname)
    return hyphae


def main():
    args = parser.parse_args()


    # First read in NCBI mapping for GO enrichmnet
    ncbi = pd.read_csv('../../data/gene_to_ncbi.txt', sep='\t', dtype={'NCBI Gene ID':str}).dropna()
    ncbi = dict(zip(ncbi['Approved symbol'], ncbi['NCBI Gene ID']))

    g = pickle.load(open(gfile, 'rb'))#hyp.make_graph_from_dict(gfile)
    if args.fromFile is None:
        hyphae = significant_genes(data, 'AML sample', 'Gene', 'LogFoldChange')
    else:
        hyphae = loadFromFile(args.fromFile)


    for key, this_hyp in hyphae.items():
        this_hyp.node_stats().to_csv(key+'_nodelist.csv')
        if args.doEnrich:
            if len(this_hyp.forest_enrichment) == 0:
                for_e = hyEnrich.go_enrich_forests(this_hyp, ncbi)
                this_hyp.assign_enrichment(for_e, type='forest')
                for_e.to_csv(key+'enrichedForestGoTerms.csv')
                this_hyp._to_file(key+'_hypha.pkl')
            this_hyp.forest_stats().to_csv(key+'_TreeStats.csv')

if __name__ == '__main__':
    main()