'''
Sub-module designed to rid us of the omics integrator fake graph
'''
import pandas as pd
import pickle
import argparse
import sys
import igraph as ig
import networkx as nx

def buildMappingFromStringFile(fname):
    '''
    Builds a mapping dictionary from string file
    '''
    df = pd.read_csv(fname, sep='\t')
    name_dict = dict(zip(df['protein_external_id'],df['preferred_name']))
    return name_dict

def getEdgeWeightsFromString(fname,colname='experimental'):
    '''
    reads in edge weights from string file
    '''
    full_tab = pd.read_csv(fname,sep='\t')
    red_tab = full.tab['protein1','protein2',colname]
    return red_tab

def buildIgraphFromString(tab, mapping dest_dir):
    '''
    Builds igraph network from string file
    '''

def buildNxFromSring(tab, mapping dest_dir):
    colname='experimental'


def buildPCSTDictFromString(fname, mapping_file, dest_dir):
    '''
    creates a dictionary with the required input from pcst_fast package
    '''
    pcst_dict={'edges':[],'nodes':[],'edgeWeights':[]}
    #first read in mapping file
    node_map=buildMappingFromStringFile(mapping_file)

    open( dest_dir+"/string_pcst.pkl", "rb" )


def main():
    parser = argparse.ArgumentParser(description="""Process various graph formats and store appropriately""")
    parser.add_argument('--graphFile',dest='name',required=True,help='Path to file')
    parser.add_argument('--graphSource',default='string',dest='graphType',help='Source of graph file. Currently only supports `string`')
    parser.add_argument('--nodeMapping',dest='mapping',help='mapping file if needed')
    parser.add_argument('--outputFormat',dest='Desired output format')
    parser.add_argument('--dest',dest='destDir',help='Destination directory to store pkl file')
    args = parser.parse_args()

    if args.graphType == 'string':
        tab = getEdgeWeightsFromString(args.name)
        mapping = buildMappingFromStringfile(args.mapping)x




if __name__ == '__main__':
    main()
