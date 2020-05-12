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
import pandas as pd

#Parser information for command line
parser = argparse.ArgumentParser(description="""Get data from the proteomic \
                                 data commons and build community networks""")
parser.add_argument('--enrich', dest='doEnrich', action='store_true',\
                    default=False, help='Flag to do GO enrichment')
parser.add_argument('--saveGraphs', dest='toFile', action='store_true',\
                    default=False, help='Flag to save networks to file')
parser.add_argument('--getDistances', dest='getDist', action='store_true',\
                    default=False, help='Get and save distances')

def compute_all_distances(hyp_dict):
    """
    Returns data frame with distances distances with columns
    col 1: dis 1
    col 2: dis 2
    col 3: net 1 name
    col 4: net 1 type
    col 5: net 2 name
    col 6: net 3 type
    col 7: distance
    """
    df_list = []
    for key1, hyp1 in hyp_dict.items():
        #compute forest distances
        within_dist = hyp1.within_distances()
        within_dist['hyp1'] = key1
        within_dist['hyp2'] = key1
        df_list.append(within_dist)
        for key2, hyp2 in hyp_dict.items():
            if key1 != key2:
                comm_net = hyp1.inter_distance(hyp2)
                comm_net['hyp1'] = key1
                comm_net['hyp2'] = key2
                df_list.append(pd.DataFrame(comm_net))
                #inter_distance is NOT symmetric
                comm_net = hyp2.inter_distance(hyp1)
                comm_net['hyp1'] = key2
                comm_net['hyp2'] = key1
                df_list.append(pd.DataFrame(comm_net))
    dist_df = pd.concat(df_list)
    return dist_df


def main():
    args = parser.parse_args()

    ##this is the framework for the PDC data parser.
    norms = prot.normals_from_manifest('data/PDC_biospecimen_manifest_05112020_184928.csv')

#    bcData = prot.parsePDCfile('data/TCGA_Breast_BI_Proteome.itraq.tsv')
    bcData = prot.parsePDCfile('data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv')
    lungData = prot.parsePDCfile('data/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv')
    colData = prot.parsePDCfile('data/CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv')
    gbmData = prot.parsePDCfile('data/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv')

    ##HACK: manually hard-coded in normal patients
    normPats = {'brca': set([a for a in bcData['Patient'] if a in norms['Breast Invasive Carcinoma']]),\
                'coad': set([a for a in colData['Patient'] if a in norms['Colon Adenocarcinoma']]),\
                'luad': set([a for a in lungData['Patient'] if a in norms['Lung Adenocarcinoma']]),\
                'gbm': set([a for a in gbmData['Patient'] if a in norms['Other']])}



    gfile='../../../OmicsIntegrator2/interactomes/inbiomap.9.12.2016.full.oi2'
    g = hyp.make_graph(gfile)

    namemapper = None #hyp.mapHGNCtoNetwork()
    ncbi = prot.map_ncbi_to_gene(bcData)

    ##here we get the top values for each patient
    patVals = {'brca':prot.getProtsByPatient(bcData, namemapper),\
               'luad':prot.getProtsByPatient(lungData, namemapper),\
             'coad':prot.getProtsByPatient(colData, namemapper),\
             'gbm':prot.getProtsByPatient(gbmData, namemapper)}

    patDiffs = {'brca': prot.getTumorNorm(bcData,normPats['brca'],namemapper),
                'luad': prot.getTumorNorm(lungData,normPats['luad'],namemapper),
                'coad': prot.getTumorNorm(colData,normPats['coad'],namemapper),
                'gbm': prot.getTumorNorm(gbmData,normPats['gbm'],namemapper)}
    #now we want to build network communities for each
    hyphae = dict()

    for key in patDiffs:
        h = hyp(patDiffs[key], g)
        h.forest_stats().to_csv(key+'_forestStats.csv')
        members = h.runCommunityWithMultiplex()
        members.to_csv(key+'communities.csv')
        hyphae[key] = h
        if args.toFile:
            h.saveCommunityToFile(prefix=key)
        if args.doEnrich:
            hype.go_enrich_forests(h, ncbi).to_csv(key+'enrichedForestGoTerms.csv')
            hype.go_enrich_communities(h, ncbi).to_csv(key+'enrichedCommunityGOterms.csv')

    #now compute graph distances to ascertain fidelity
    if args.getDist:
        res = compute_all_distances(hyphae)
        res.to_csv('panCancerDistances.csv')


if __name__ == '__main__':
    main()
