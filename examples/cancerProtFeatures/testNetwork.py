'''
Create test suite that examines networks properties
'''

import argparse
import sys
sys.path.insert(0, "../../")


import hyphalnet.hypha as hyp
from hyphalnet.hypha import hyphalNetwork
import hyphalnet.proteomics as prot
import hyphalnet.hyphEnrich as hype
import pandas as pd
import leidenalg as la
import networkx as nx
from igraph import *


netDict={'brca':'brca_hypha.pkl','luad':'luad_hypha.pkl',\
         'coad':'coad_hypha.pkl','gbm':'gbm_hypha.pkl'}

hyphae = dict()
for key, fname in netDict.items():
    hyphae[key] = hyp.load_from_file(fname)

##now load in basic interactome
gfile = '../../../OmicsIntegrator2/interactomes/inbiomap.9.12.2016.full.oi2'
g = hyp.make_graph(gfile)

print("Maing igraph graph")
ig = Graph(directed=False)
for gi in g.nodes:
    ig.add_vertex(gi)

ig.add_edges(g.edges)

print("finding communities")
partition = la.find_partition(ig, la.ModularityVertexPartition)
comm_dict={}
for p in range(len(partition)):
    comm_dict[p]=ig.vs.select(partition[p])['name']

print("computing NMI")
df = pd.DataFrame()
##how do we compare eachx hypha community to interactome?
for n1, h1 in hyphae.items():
    for n2, h2 in hyphae.items():
        nmi = hyp.computeCommunityNMI(h1.communities, h2.communities)
        df = df.append({'hyp1':n1, 'hyp2':n2, 'NMIdist':nmi},ignore_index=True)
    nmi = hyp.computeCommunityNMI(h1.communities, comm_dict)
    df = df.append({'hyp1':n1, 'hyp2':'fullInteractome', 'NMIdist':nmi},ignore_index=True)
df.to_csv("allNMI.csv")
