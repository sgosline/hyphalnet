"""
Evaluate the ability of hyphal networks to correlate with drug response
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
parser.add_argument('--graph', dest='graph', default='../../data/igraphPPI.pkl',\
                help='Path to pickled igraph interactome')

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
    beta = 0.5
    #get mutational data
    mvals = getMutationalData()

    ##load up interactome
    gfile = args.graph
    ##TODO: replace this with Docker image call
    g = pickle.load(open(gfile, 'rb'))
    key = 'mpnstPDXmuts'
    this_hyp = hyphalNetwork(mvals, g)
    this_hyp._to_file(key+'_hypha.pkl')

    ##read from file
    ###this is all we need to do in a single eval, then we can do tests later
    this_hyp.node_stats().to_csv(key+'_nodelist.csv')
    for_e = hyEnrich.go_enrich_forests(this_hyp)
    this_hyp.assign_enrichment(for_e, type='forest')
    for_e.to_csv(key+'enrichedForestGoTerms.csv')
    com_e = hyEnrich.go_enrich_communities(this_hyp)
    this_hyp.assign_enrichment(com_e, type='community')
    this_hyp._to_file(key+'_hypha.pkl')
    com_e.to_csv(key+'enrichedCommunityGOterms.csv')
    this_hyp.community_stats(prefix=key).to_csv(key+'_communityStats.csv')
    res = hyStats.compute_all_distances({'mutations':this_hyp})
    res.to_csv('panPDXDistances.csv')
    nmi = hyStats.compute_all_nmi({'mutations': this_hyp}, g)
    nmi.to_csv('panPDXNMI.csv')

if __name__=='__main__':
    main()
