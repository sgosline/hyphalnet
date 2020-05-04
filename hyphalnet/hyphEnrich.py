import goenrich
import numpy as np
import pandas as pd

def get_go_enrichment(genelist, background):
    O = goenrich.obo.ontology('db/go-basic.obo')
    gene2go = goenrich.read.gene2go('db/gene2go.gz')
    #  gene2go = gene2go.loc[gene2go['GeneID'].isin(background)]
    values = {k: set(v) for k, v in gene2go.groupby('GO_ID')['GeneID']}
    background_attribute = 'gene2go'
    goenrich.enrich.propagate(O, values, background_attribute)
    df = goenrich.enrich.analyze(O, np.array(genelist), background_attribute)
    df = df.dropna().loc[df['q'] < 0.05]
    return df


def go_enrich_forests(hyp, ncbi_mapping):
    """        Iterates through forests and gets enrichment """
    enrich = []
    background = []
    for ns in hyp.node_counts.keys():
        try:
            background.append(ncbi_mapping[ns])
        except KeyError:
            print("No key for background gene", ns)
    for pat, forest in hyp.forests.items():
        nodenames = []
        for fn in forest.nodes():
            try:
                nodenames.append(int(ncbi_mapping[fn]))
            except KeyError:
                print("No key for gene", fn)
        try:
            evals = get_go_enrichment(nodenames, background)
            evals['Patient'] = pat
            enrich.append(evals)
        except Exception as e:
            print(e)
    return pd.concat(enrich)

def go_enrich_communities(hyp, ncbi_mapping):
    """ Gets enrichment for individual communities"""
    print('Doing community thing')
    enrich = []
    background = []
    for ns in hyp.node_counts.keys():
        try:
            background.append(ncbi_mapping[ns])
        except KeyError:
            print("No key for background gene", ns)
    for comm, nodes in hyp.communities.items():
        print(comm)
        nodenames = []
        for fn in nodes:
            try:
                nodenames.append(int(ncbi_mapping[fn]))
            except KeyError:
                print("No key for gene", fn)
        try:
            evals = get_go_enrichment(nodenames, background)
            evals['Community'] = comm
            enrich.append(evals)
        except Exception as e:
            print(e)
    return pd.concat(enrich)
