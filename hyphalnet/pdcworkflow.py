import pandas as pd
import re

class pdcworkflow:
    def __init__(self):
        '''Constructor for workflow'''
        self.datasets={}
        self.metadata={}

    def downloadPDFfile():
        print("Not sure how this will work just yet")

    def parsePDCfile(fpath='data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv'):
        """
        Takes a PDC file ending in .tmt10.tsv or .itraq.tsv and creates
        tidied data frame with Gene, Patient, logratio and diffFromMean values

        Parameters
        ----------
        fpath : TYPE, optional
            DESCRIPTION. The default is 'data/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv'.

        Returns
        -------
        None.

        """
        dat = pd.read_csv(fpath,sep='\t')
        newdat = dat[['Gene','NCBIGeneID']]

        #retrieve log ratios
        p=re.compile('.*[0-9]+\ Log Ratio')
        pats = list(filter(p.match,dat.keys()))
        for pat in pats:
            up_pat = pat.replace(' Log Ratio','')
            newdat[up_pat] = dat[pat]

        #now tidy data by log ratio by patient
        tdat = pd.melt(newdat,id_vars=['Gene','NCBIGeneID'],var_name='Patient',value_name='logratio')

        return(tdat)


    def getProtsByPatient(tdf,namemapper,column='logratio',quantThresh=0.01):
        """
        Gets proteins from tidied data frame into dictionary for OI

        Parameters
        ----------
        tdf : Data frame
            Data frame indexed by Gene with 'Patient' column representing patient name, and 'diffFromMean' column

        Returns
        -------
        Dictionary of dictionaries

        """
            #compute patient level quantiles
       # gene_means=tdat.groupby('Gene')['logratio'].mean()
        pquants = pd.DataFrame({'thresh':tdf.groupby("Patient")[column].quantile(1.0-quantThresh)})

        nm = pd.DataFrame.from_dict(namemapper,orient='index',columns=['Prot'])
        nm.loc[:,'Gene']=nm.index

        tdat = tdf.merge(pquants,on='Patient')
        tdat = tdat.merge(nm,on='Gene')
        tdat = tdat.assign(topProt=tdat[column]>tdat['thresh'])

        selvals=tdat[tdat['topProt']]

        dprots = selvals.groupby('Patient')['Prot'].apply(list).to_dict()
        dvals = selvals.groupby("Patient")[column].apply(list).to_dict()

        res = {}
        for k in dprots.keys():
            res[k] = dict(zip(dprots[k],dvals[k]))

        return(res)
