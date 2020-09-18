'''
Gets data from synapse and processes it into a single dictionary for each sample
'''

import synapseclient as sc
import pandas as pd
import os
import hyphalnet.hypha as hyp
from hyphalnet.hypha import hyphalNetwork
import hyphalnet.hyphEnrich as hyEnrich
import hyphalnet.hyphaeStats as hyStats
import pickle
import argparse
import gzip

parser = argparse.ArgumentParser(description="""Get data from scProts""")
parser.add_argument('--ref', dest='refid', help='Synapse id of reference pkl')
parser.add_argument('--eval', dest='evalid',help='Synapse id of hypha to be evaluated')
parser.add_argument('--refName', dest='refName')
parser.add_argument('--evalName', dest='evalName')


def main():
    args = parser.parse_args()
    syn = sc.login()
    h1 = pickle.load(open(syn.get(args.refid).path,'rb'))
    h2 = pickle.load(open(syn.get(args.evalid).path,'rb'))

    newd= { args.refName: h1, args.evalName: h2}
    res = hyStats.compute_all_distances(newd)

    fname= 'hyp2hypDistances_'+args.refName+'_to_'+args.evalname+'.csv'
    res.to_csv(fname)


if __name__=='__main__':
    main()
