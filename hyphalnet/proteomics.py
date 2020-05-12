import pandas as pd
import re
import numpy as np


def downloadPDCfile():
    print("Not sure how this will work just yet")


def normals_from_manifest(fname):
    """ Parses PDC sample manifest"""
    dat = pd.read_csv(fname, sep=',')
    return dat.groupby("Disease Type")['Aliquot Submitter ID'].apply(list).to_dict()

def map_ncbi_to_gene(tdat):
    """ takes a  parsed file and returns dictionary of gene maps"""
    tdat = tdat.loc[~tdat['Gene'].isin(list(['Mean', 'Median', 'StdDev']))]
    return dict(zip(tdat['Gene'], [str(int(a)) for a in tdat['NCBIGeneID']]))

def parsePDCfile(fpath='data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv'):
    """
    Takes a PDC file ending in .tmt10.tsv or .itraq.tsv and creates
    tidied data frame with Gene, Patient, logratio and diffFromMean values

    Parameters
        ----------
        fpath : TYPE, optional
        DESCRIPTION. The default is 'data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv'.

        Return
        -------
        None.
        """

    dat = pd.read_csv(fpath, sep='\t')
    newdat = dat[['Gene', 'NCBIGeneID']]

    #retrieve log ratios
    pat = re.compile('.*[0-9]+\ Log Ratio')
    pats = list(filter(pat.match, dat.keys()))

    for pat in pats:
        up_pat = pat.replace(' Log Ratio', '')
        newdat[up_pat] = dat[pat]

        #now tidy data by log ratio by patient
    tdat = pd.melt(newdat, id_vars=['Gene', 'NCBIGeneID'],\
                   var_name='Patient', value_name='logratio')
    return tdat

def getProtsByPatient(tdf, namemapper=None, column='logratio', quantThresh=0.01):
    """
        Gets proteins from tidied data frame into dictionary for OI

        Parameters
        ----------
        tdf : Data frame
        Data frame indexed by Gene with 'Patient' column representing \
        patient name, and 'diffFromMean' column

        Returns
        -------
        Dictionary of dictionaries

        """
            #compute patient level quantiles
       # gene_means=tdat.groupby('Gene')['logratio'].mean()
    pquants = pd.DataFrame({'thresh':tdf.groupby("Patient")[column].quantile(1.0-quantThresh)})
    tdat = tdf.merge(pquants, on='Patient')

    if namemapper is not None:
        nm = pd.DataFrame.from_dict(namemapper, orient='index', columns=['Prot'])
        nm.loc[:, 'Gene'] = nm.index
        tdat = tdat.merge(nm, on='Gene')
    else:
        tdat.rename(columns={'Gene':'Prot'}, inplace=True)

    tdat = tdat.assign(topProt=tdat[column] > tdat['thresh'])
    selvals = tdat[tdat['topProt']]
    dprots = selvals.groupby('Patient')['Prot'].apply(list).to_dict()
    dvals = selvals.groupby("Patient")[column].apply(list).to_dict()
    res = {}
    for k in dprots.keys():
        res[k] = dict(zip(dprots[k], dvals[k]))

    return res


def getTumorNorm(tdf, normSamps, namemapper=None, column='logratio', quantThresh=0.01):
    """
    Gets per-patient tumor values compared to pooled normal
    TODO: update to do matched normal instead

    """

    tumSamps = set([a for a in tdf['Patient'] if a not in normSamps])
    normVals = tdf[tdf.Patient.isin(normSamps)]
    tumVals = tdf[tdf.Patient.isin(tumSamps)]

    meanVals = normVals.groupby('Gene')[column].apply(np.mean)

    tumMat = tumVals.pivot(index='Gene', columns='Patient', values=column)

    diffs = tumMat.subtract(meanVals, axis=0)
    diffs['Gene'] = diffs.index
    fd = diffs.melt(id_vars='Gene',value_vars=tumMat.columns, value_name='diffsToNormal',var_name='Patient')
    fd['absVal'] = np.abs(fd['diffsToNormal'])
    dquants = pd.DataFrame({'thresh':fd.groupby("Patient")['absVal'].quantile(1.0-quantThresh)})
    #which genes/patientss are abve that threshold
    fd = fd.merge(dquants, on='Patient')
    fd = fd.assign(topProt=fd['absVal'] > fd['thresh'])
    selvals = fd[fd['topProt']]
    dprots = selvals.groupby('Patient')['Gene'].apply(list).to_dict()
    dvals = selvals.groupby("Patient")['absVal'].apply(list).to_dict() #can't do neg prizes
    #return those values
    res = {}
    for k in dprots.keys():
        res[k] = dict(zip(dprots[k],dvals[k]))
    return res
