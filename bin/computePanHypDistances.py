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

parser = argparse.ArgumentParser(description="""Get data from scProts""")
parser.add_argument('--ref', dest='refid', default='syn22392951', help='Synapse id of reference pkl')
parser.add_argument('--refName', dest='refName', default='CPTACpancan')
parser.add_argument('--synProj',dest='synProj',help='id of synapse project to store result')
parser.add_argument('--eval',dest='evalFiles', nargs=1,\
                    action=kvdictAppendAction,metavar='KEY=VALUE',\
                    help='Key/value params for synapse ids')
def main():
    args = parser.parse_args()
    syn = sc.login()
    h1 = pickle.load(open(syn.get(args.refid).path, 'rb'))
    #h2 = pickle.load(open(syn.get(args.evalid).path, 'rb'))

    newd = {args.refName: h1}
    for key, val in args.evalFiles.items():
        newd[key] = pickle.load(open(syn.get(val).path,'rb'))

    res = hyStats.compute_all_distances(newd)

    fname = 'hyp2hypDistances_'+args.refName+'_to_data.csv'
    res.to_csv(fname)
    if args.synProj is not None:
        tab = sc.table.build_table("Original data to "+args.refName+' Distances',
                                   args.synProj, res)
        syn.store(tab)

def build_hyphae_from_data(qt, g, sample=False):
    """ Temp function to load data from local directory"""
    ##this is the framework for the PDC data parser.

    #now we want to build network communities for each
    hyphae = dict()
    patDiffs = loadCancerData(qt)
    beta = 0.5
    for key, vals in patDiffs.items():
        if sample:
            new_vals = {}
            for v in random.sample(list(vals), 5):
                new_vals[v] = vals[v]
            vals = new_vals
            print(len(vals))
        this_hyp = hyphalNetwork(vals, g.copy())
        hyphae[key+str(qt)] = this_hyp
#        this_hyp._to_file(key+str(qt)+'_hypha.pkl')
    return hyphae

if __name__=='__main__':
    main()
