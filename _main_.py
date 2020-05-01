# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:05:29 2020

@author: gosl241
"""

import argparse
import hyphalnet as hypnet
from hyphalnet.hypha import hypha as hyp
from hyphalnet.pdcworkflow import pdcworkflow as pdc

parser = argparse.ArgumentParser(description="""Get data from the proteomic /
                                 data commons and build community networks""")

def main():
    args = parser.parse_args()

    ##this is the framework for the PDC data parser.
    bcData = pdc.parsePDCfile('data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv')
    lungData = pdc.parsePDCfile('data/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv')
    colData = pdc.parsePDCfile('data/CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv')
    gbmData = pdc.parsePDCfile('data/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv')

    namemapper = None #hyp.mapHGNCtoNetwork()


    gfile='../OmicsIntegrator2/interactomes/inbiomap.9.12.2016.full.oi2'
    g = hyp.make_graph(gfile)

    ##here we get the top values for each patient
    patVals={'brca':pdc.getProtsByPatient(bcData, namemapper),\
             'luad':pdc.getProtsByPatient(lungData, namemapper),\
             'coad':pdc.getProtsByPatient(colData, namemapper),\
             'gbm':pdc.getProtsByPatient(gbmData, namemapper)}

    #now we want to build network communities for each
    hyphae = dict()
    parts = dict()
    ncbi = pdc.map_ncbi_to_gene(bcData)
    for key in patVals:
        h = hyp(patVals[key], g)
        h.forest_stats().to_csv(key+'_forestStats.csv')
        members = h.runCommunityWithMultiplex()
        members.to_csv(key+'communities.csv')
#        h.saveCommunityToFile(prefix=key)
        h.go_enrich_forests(ncbi).to_csv(key+'enrichedForestGoTerms.csv')
        h.go_enrich_communities(ncbi).to_csv(key+'enrichedCommunityGOterms.csv')
        hyphae[key] = h

if __name__ == '__main__':
    main()
