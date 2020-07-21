import goenrich
import numpy as np
import pandas as pd


def get_kegg_enrichment(genelist, background):
    """Get KEGG based pathway enrichment"""
    O = goenrich.obo.ontology('db/go-basic.obo')
    gene2go = goenrich.read.gene2go('db/gene2go.gz')
    #  gene2go = gene2go.loc[gene2go['GeneID'].isin(background)]
    values = {k: set(v) for k, v in gene2go.groupby('GO_ID')['GeneID']}
    background_attribute = 'gene2go'
    goenrich.enrich.propagate(O, values, background_attribute)
    df = goenrich.enrich.analyze(O, np.array(genelist), background_attribute)
    df = df.dropna().loc[df['q'] < 0.05]
    return df

def get_go_enrichment(genelist, background):
    O = goenrich.obo.ontology('db/go-basic.obo')
    gene2go = goenrich.read.gene2go('db/gene2go.gz')
    #  gene2go = gene2go.loc[gene2go['GeneID'].isin(background)]
    values = {k: set(v) for k, v in gene2go.groupby('GO_ID')['GeneID']}
    background_attribute = 'gene2go'
    goenrich.enrich.propagate(O, values, background_attribute)
    df = goenrich.enrich.analyze(O, np.array(genelist), background_attribute)
    df = df.dropna().loc[df['q'] < 0.05]
 #   df = df.dropna().loc[df['namespace']=='biological_process']
    return df


def go_enrich_forests(hyp, ncbi_mapping):
    """
    Iterates through forests and gets enrichment
    """
    enrich = []
    background = []
    missed = 0
    for ns in hyp.interactome.vs['name']:
        #node_counts.keys():
        try:
            background.append(ncbi_mapping[ns])
        except KeyError:
            missed = missed+1
    print("Have background of", len(background), 'missing', missed)
    for pat, forest in hyp.forests.items():
        missed = 0
    #hyp.forests.items():
        nodenames = []
        for fn in list([hyp.interactome.vs['name'][i] for i in forest['vert']]):
            try:
                nodenames.append(int(ncbi_mapping[fn]))
            except KeyError:
                #print("No key for gene", fn)
                missed = missed+1
            except Exception as e:
                print('patient', pat, e, len(fn))
        print("Found matches for", len(nodenames), 'nodes for patient', pat, 'missing', missed)
        try:
            evals = get_go_enrichment(nodenames, background)
            evals['Patient'] = pat
            enrich.append(evals)
        except Exception as e:
            print('patient', pat, e)
    return pd.concat(enrich)

def go_enrich_communities(hyp, ncbi_mapping):
    """
    Gets enrichment for individual communities
    """
    print('Doing community thing')
    enrich = []
    background = []
    missed = 0
    for ns in hyp.interactome.vs['name']:
        #['nodes']:
        #hyp.node_counts.keys():
        try:
            background.append(ncbi_mapping[ns])
        except KeyError:
            #print("No key for background gene", ns)
            missed = missed+1
    print("Have background of", len(background), 'missing', missed)
    missed = 0
    for comm, nodes in hyp.communities.items():
      #  print(comm)
        nodenames = []
        for fn in nodes:
            try:
                nodenames.append(int(ncbi_mapping[fn]))
            except KeyError:
                #print("No key for gene", fn)
                missed = missed +1
        print("Found matches for", len(nodenames), 'nodes for community', comm, 'missing', missed)
        missed=0
        try:
            evals = get_go_enrichment(nodenames, background)
            evals['Community'] = str(comm)
            enrich.append(evals)
        except Exception as e:
            print('community', comm, e)
    return pd.concat(enrich)
