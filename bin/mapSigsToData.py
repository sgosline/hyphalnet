'''
Basic script that reads in a data frame of samples, proteins, and values
And outputs a compressed hypha object to be analyzed by downstream methods
'''
import argparse
import sys

import hyphalnet.hypha as hyp
from hyphalnet.hypha import hyphalNetwork
import pickle
import pandas as pd

parser = argparse.ArgumentParser(description="""Get data from mutations""")
parser.add_argument('--hypha', dest='hyph', help='HyphalNetwork Signatures')
parser.add_argument('--inputData', dest='df', \
                    help='Data frame with three columns: `Sample`, `Gene`, `Value`')
parser.add_argument('--output', dest='output', default='hyphaOutput',\
                    help='Prefix to use for output file')


def df2dict(df):
    '''
    Converts data frame to dictionary
    '''
    genes = df.groupby('Sample')['Gene'].apply(list)
    ad = df.groupby('Sample')['Value'].apply(list)
    pdict = {}
    for samp, genes in genes.items():
        pdict[samp] = dict(zip(genes, ad[samp]))
    return pdict

def main():
    args = parser.parse_args()
    beta = 0.5
    #get mutational data
    mdf = pd.read_csv(args.df)
    key = args.output
    mvals = df2dict(mdf)
    ##load up interactome
    sigs = pickle.load(open(args.hyph, 'rb'))
    g = sigs.interactome

    this_hyp = hyphalNetwork(mvals, g)
    this_hyp._to_file(key+'_hypha.pkl')

    comm_dist = this_hyp.hypha_community_distances(sigs)
    ##now compute distances between panCan sigs and these networks
    comm_dist.to_csv()

if __name__=='__main__':
    main()
