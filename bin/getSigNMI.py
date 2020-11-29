'''
This script builds individual hyphalnetworks and computes the NMI
'''

import pandas as pd
import os
import proteomicsData as prot
import hyphalnet.hypha as hyp
from hyphalnet.hypha import hyphalNetwork
import hyphalnet.hyphaeStats as hyStats
import pickle
import argparse
import gzip
import proteomicsData as pdat
import re
import numpy as np

parser = argparse.ArgumentParser(description="""Collect pan can NMI values""")
parser.add_argument('--refName', dest='refName', default='CPTACpancan')
parser.add_argument('--quantile', dest='quant', default=0.01,  help='Threshold to use for top-expressed proteins')
parser.add_argument('--hypha', dest='hyph', help='Original pan can hypha')

def main():
    args = parser.parse_args()
    qval = args.quant

    ##first we run a helper function to make sure we have all cptac data
    fdict = pdat.cptacData()

    ##first get proteomics measurements
    allDat = pdat.getCancerData(fdict, qval, byType=True)

    ##get hyphal network and graph
    phyph = pickle.load(open(args.hyph, 'rb'))
    g = phyph.interactome

    beta = .5

    #build hyphal network of network communities

    hDict= {'panCan': phyph}
    for ct, dat in allDat.items():
        hDict[ct] = hyphalNetwork(dat, g, beta)

    nmi = hyStats.compute_all_nmi(hDict, g)
    nmi.to_csv(args.refName+'_nmi.csv')

if __name__=='__main__':
    main()
