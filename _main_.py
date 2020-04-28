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
    bcData = pdc.parsePDCfile('data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv')
    lungData = pdc.parsePDCfile('data/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv')
    colData = pdc.parsePDCfile('data/CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv')
    gbmData = pdc.parsePDCfile('data/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv')

    namemapper = hyp.mapHGNCtoNetwork()

    gfile='data/9606.protein.links.v11.0.txt'
    g = hyp.makeGraph(gfile)

    patVals={'brca':pdc.getProtsByPatient(bcData,namemapper),'luad':pdc.getProtsByPatient(lungData,namemapper),'coad':pdc.getProtsByPatient(colData,namemapper),'gbm':pdc.getProtsByPatient(gbmData,namemapper)}

    #now we want to build network communities for each
    hyphae=dict()
    parts = dict()
    for key in patVals:
        h =  hyp(patVals[key],g)
        hyphae[key]= h
#        parts[key]= h.runCommunity()

    #now we have hyphae for each disease



if __name__ == '__main__':
    main()
