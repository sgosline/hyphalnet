"""
EValuate the ability of hyphal networks to correlate with drug response
"""


import argparse
import sys
#sys.path.insert(0, '../../')

import hyphalnet.hypha as hyp
from hyphalnet.hypha import hyphalNetwork
import hyphalnet.hyphEnrich as hyEnrich
import hyphalnet.hyphaeStats as hyStats
import pickle
import pandas as pd
import synapseclient as synapse

parser = argparse.ArgumentParser(description="""Get data from mutations""")


def getMutationalData():
    syn = synapse.login()
    tab = syn.tableQuery("SELECT * FROM syn22266278").asDataFrame()
    genes = tab.groupby('specimenID')['Gene'].apply(list)
    ad = tab.groupby('specimenID')['AD'].apply(list)
    pdict = {}
    for samp, genes in genes.items():
        pdict[samp] = dict(zip(genes, ad[samp]))
    return pdict

def main():
    args = parser.parse_args()
    beta =0.5
    #get mutational data
    mvals = getMutationalData()

    ##load up interactome
    gfile='../../data/igraphPPI.pkl'
    g = pickle.load(open(gfile, 'rb'))#hyp.make_graph_from_dict(gfile)

      # First read in NCBI mapping for GO enrichmnet
    #ncbi = pd.read_csv('../../data/gene_to_ncbi.txt', sep='\t', dtype={'NCBI Gene ID':str}).dropna()
    #ncbi = dict(zip(ncbi['Approved symbol'], ncbi['NCBI Gene ID']))

    this_hyp = hyphalNetwork(mvals, g)
    key = 'mpnstPDXmuts'
    this_hyp.node_stats().to_csv(key+'_nodelist.csv')
    for_e = hyEnrich.go_enrich_forests(this_hyp)
    this_hyp.assign_enrichment(for_e, type='forest')
    for_e.to_csv(key+'enrichedForestGoTerms.csv')
    com_e = hyEnrich.go_enrich_communities(this_hyp)
    this_hyp.assign_enrichment(com_e, type='community')
    this_hyp._to_file(key+'_hypha.pkl')
    com_e.to_csv(key+'enrichedCommunityGOterms.csv')
    this_hyp.community_stats(prefix=key).to_csv(key+'_communityStats.csv')

if __name__=='__main__':
    main()
