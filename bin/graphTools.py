'''
Sub-module designed to rid us of the omics integrator fake graph
'''
import pandas as pd
import pickle
import argparse
import sys
import igraph as ig
#import networkx as nx

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
    full_tab = pd.read_csv(fname, sep=' ')
    red_tab = full_tab[['protein1', 'protein2', colname]]
    red_tab = red_tab[red_tab[colname] != 0]
    red_tab[colname] = [x/1000 for x in red_tab[colname]]
    return red_tab

def buildIgraphFromFile(fname, dest_dir):
    '''
    Builds igraph network from string file
    '''
    igg = ig.Graph.Read_Ncol(fname)
    igg.es['cost'] = [1-e for e in tab.es['weight']]
    pickle.dump(igg, open(dest_dir+'/igraphPPI.pkl', "wb"))
    return igg

#def buildNxFromFile(fname, dest_dir):
#    tnet = nx.readwrite.edgelist.read_weighted_edgelist(fname)
#    pickle.dump(tnet, open(dest_dir+'/networkPPI.pkl', "wb"))


def buildPCSTDictFromFile(fname, dest_dir,igg):
    '''
    creates a dictionary with the required input from pcst_fast package
    '''
    pcst_dict = {'edges':[], 'nodes':[], 'edgeWeights':[]}
    #first read in mapping file
    tab = pd.read_csv(fname, sep='\t')
    nodes = igg.vs['name']#list(set(tab.iloc[:0]).union(set(tab.iloc[:,1])))
    weights = tab.iloc[:,2]
    edges = []
    deg = igg.degree()
    for index, row in tab.iterrows():
        i1 = nodes.index(row[0])
        i2 = nodes.index(row[1])
        edges.append([i1,i2])
    pcst_dict = {'edges':edges, 'nodes':nodes, 'edgeWeights':weights, 'cost':1-weights, 'degree':deg}

    pickle.dump(pcst_dict, open(dest_dir+"/pcstDictPPI.pkl", "wb" ))

def write_to_file(tab, mapping):
    '''
    Write the table to a three column file - p1, p2, weight
    '''
    tab['protein1'] = [mapping[p] for p in tab['protein1']]
    tab['protein2'] = [mapping[p] for p in tab['protein2']]
    fname = 'tmpNet.tsv'
    tab.to_csv(fname, sep='\t', index=False, header=False)
    return fname

def main():
    parser = argparse.ArgumentParser(description="""Process various graph formats and store appropriately""")
    parser.add_argument('--graphFile', dest='name',required=True,help='Path to file')
    parser.add_argument('--graphSource', default='string',dest='graphType',help='Source of graph file. Currently only supports `string`')
    parser.add_argument('--nodeMapping', dest='mapping',help='mapping file if needed')
    parser.add_argument('--dest', dest='destDir',help='Destination directory to store pkl file')
    args = parser.parse_args()

    if args.graphType == 'string':
        tab = getEdgeWeightsFromString(args.name)
        mapping = buildMappingFromStringFile(args.mapping)

    fname=write_to_file(tab, mapping)
  #  buildNxFromFile(fname, args.destDir)
    ig=buildIgraphFromFile(fname, args.destDir)
    buildPCSTDictFromFile(fname, args.destDir, ig)


if __name__ == '__main__':
    main()
