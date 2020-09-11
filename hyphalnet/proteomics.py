import pandas as pd
import re
import numpy as np


def downloadPDCfile():
    print("Not sure how this will work just yet")


def normals_from_manifest(fname):
    """
    Parses PDC sample manifest

    Parameters
    --------
    fname: chr, name of file
    Return
    --------
    Dictionary of normal samples by disease type
    """
    dat = pd.read_csv(fname, sep=',')
    nd = dat.groupby("Disease Type")['Aliquot Submitter ID'].apply(list).to_dict()
    print(nd.keys())
    return nd

def parsePDCfile(fpath='data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv'):
    """
    Takes a PDC file ending in .tmt10.tsv or .itraq.tsv and creates
    tidied data frame with Gene, Patient, logratio and diffFromMean values

    Parameters
        ----------
        fpath : chr, optional
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
        Gets proteins from tidied data frame into dictionary for HyphalNetwork

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




def getTumorNorm(tdf, normSamps, namemapper=None, column='logratio', quantThresh=0.01, doAbs=False):
    """
    Gets per-patient tumor values compared to pooled normal
    TODO: update to do matched normal instead

    """

    tumSamps = set([a for a in tdf['Patient'] if a not in normSamps])

    #separate data frame by tumor vs normal
    normVals = tdf[tdf.Patient.isin(normSamps)]
    tumVals = tdf[tdf.Patient.isin(tumSamps)]

    #TODO: match tumor/normal samples, for now just get mean for each gene
    meanVals = normVals.groupby('Gene')[column].apply(np.mean)

    #subtract mean log ratio vales
    tumMat = tumVals.pivot(index='Gene', columns='Patient', values=column)
    diffs = tumMat.subtract(meanVals, axis=0)
    diffs['Gene'] = diffs.index
    fd = diffs.melt(id_vars='Gene',value_vars=tumMat.columns,\
                    value_name='diffsToNormal', var_name='Patient')

    #now calculate absolute value to get top quantile
    fd['absVal'] = np.abs(fd['diffsToNormal'])
    dquants = pd.DataFrame({'absThresh':fd.groupby("Patient")['absVal'].quantile(1.0-quantThresh)})
    topquants = pd.DataFrame({'thresh':fd.groupby("Patient")['diffsToNormal'].quantile(1.0-quantThresh)})

    #which genes/patientss are abve that threshold
    fd = fd.merge(dquants, on='Patient').merge(topquants, on='Patient')

    if doAbs:
        fd = fd.assign(topProt=fd['absVal'] > fd['absThresh'])
    else:
        fd = fd.assign(topProt=fd['diffsToNormal'] > fd['thresh'])

    selvals = fd[fd['topProt']]
    dprots = selvals.groupby('Patient')['Gene'].apply(list).to_dict()
    dvals = selvals.groupby("Patient")['absVal'].apply(list).to_dict() #can't do neg prizes
    #return those values
    res = {}
    for k in dprots.keys():
        res[k] = dict(zip(dprots[k], dvals[k]))
    return res
