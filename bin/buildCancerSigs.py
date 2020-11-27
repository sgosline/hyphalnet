'''
This script uses the cptac data package and to build a reference set of signatures
We also compute the distances between patients and the communities from the reference signatures
and compute the GO enrichment
'''

import pandas as pd
import os
import hyphalnet.hypha as hyp
from hyphalnet.hypha import hyphalNetwork
import hyphalnet.hyphaeStats as hyStats
import pickle
import argparse
import gzip
import cptac
import re
import numpy as np
import proteomicsData as pdata

parser = argparse.ArgumentParser(description="""Build hyphal network signatures for pan-can data""")
parser.add_argument('--ref', dest='refid', default='syn22392951', help='Synapse id of reference pkl')
parser.add_argument('--refName', dest='refName', default='CPTACpancan')
parser.add_argument('--quantile',dest='quant', default=0.01,  help='Threshold to use for top-expressed proteins')

def main():
    args = parser.parse_args()
    qval = args.quant

    ##first we run a helper function to make sure we have all cptac data
    fdict = pdata.cptacData()

    ##first get proteomics measurements
    allDat = pdata.getCancerData(fdict, qval, byType=False)
    patientData = pdata.getCombinedClinicalData(fdict)
    patientData.to_csv('clinicalData.csv')
    mutationData = pdata.getCombinedMutationData(fdict)
    mutationData.to_csv('mutationData.csv')

    ##make srue this file is built!
    g = pickle.load(open('../odata/igraphPPI.pkl','rb'))
    beta = .5

    #build hyphal network of network communities
    phyph = hyphalNetwork(allDat, g, beta)
    phpyh._to_file(args.refName+'_hypNet.pkl')

    #write out distances within communities
    res = phyph.distVals
    fname = args.refName+'_DistanceVals.csv'
    res.to_csv(fname)



if __name__=='__main__':
    main()
