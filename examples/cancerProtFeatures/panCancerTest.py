# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:05:29 2020

@author: gosl241
"""

import argparse
import sys
sys.path.insert(0, "../../")

from hyphalnet.hypha import hypha as hyp
import hyphalnet.proteomics as prot
import hyphalnet.hyphEnrich as hype

#Parser information for command line
parser = argparse.ArgumentParser(description="""Get data from the proteomic \
                                 data commons and build community networks""")
parser.add_argument('--enrich', dest='doEnrich', action='store_true',\
                    default=False, help='Flag to do GO enrichment')
parser.add_argument('--saveGraphs', dest='toFile', action='store_true',\
                    default=False, help='Flag to save networks to file')
parser.add_argument('--getDistances', dest='getDist', action='store_true',\
                    default=False, help='Get and save distances')

def main():
    args = parser.parse_args()

    ##this is the framework for the PDC data parser.
    bcData = prot.parsePDCfile('data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv')
    lungData = prot.parsePDCfile('data/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv')
    colData = prot.parsePDCfile('data/CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv')
    gbmData = prot.parsePDCfile('data/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv')

    gfile='../../../OmicsIntegrator2/interactomes/inbiomap.9.12.2016.full.oi2'
    g = hyp.make_graph(gfile)

    namemapper = None #hyp.mapHGNCtoNetwork()

    ##here we get the top values for each patient
    patVals = {'brca':prot.getProtsByPatient(bcData, namemapper),\
             'luad':prot.getProtsByPatient(lungData, namemapper),\
             'coad':prot.getProtsByPatient(colData, namemapper),\
             'gbm':prot.getProtsByPatient(gbmData, namemapper)}

    #now we want to build network communities for each
    hyphae = dict()

    for key in patVals:
        h = hyp(patVals[key], g)
        h.forest_stats().to_csv(key+'_forestStats.csv')
        members = h.runCommunityWithMultiplex()
        members.to_csv(key+'communities.csv')
        hyphae[key] = h
        if args.saveGraphs:
            h.saveCommunityToFile(prefix=key)
        if args.doEnrich:
            hype.go_enrich_forests(h, ncbi).to_csv(key+'enrichedForestGoTerms.csv')
            hype.go_enrich_communities(h, ncbi).to_csv(key+'enrichedCommunityGOterms.csv')

    #now compute graph distances to ascertain fidelity
    if args.getDist:
        res = compute_pairwise_distance(hyphae)

def compute_pairwise_distance(hyphae_dict):
    """Computes all vs all distances of forests and communities"""

if __name__ == '__main__':
    main()
