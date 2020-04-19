import pandas as pd
import re
import numpy as np

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
            newdat[pat.replace(' Log Ratio','')] = dat[pat]
        
        #now tidy data by log ratio by patient
        tdat = pd.melt(newdat,id_vars=['Gene','NCBIGeneID'],var_name='Patient',value_name='logratio')
        
        #compute gene level means and get values
        gene_means=tdat.groupby('Gene')['logratio'].mean()

        #divergence from mean
        newdat.index=newdat['Gene']
        gdiffs = newdat.drop(['NCBIGeneID','Gene'],1).subtract(gene_means,axis=0)
        gdiffs['Gene']=gdiffs.index
        rdat = tdat.merge(pd.melt(gdiffs,id_vars='Gene',var_name='Patient',value_name='diffFromMean'),on=['Gene','Patient'])
        

        #return tidied data frame
        return(rdat)
 
    def getProtsByPatient(tdf,column='diffFromMean',absthresh=1):
        """
        Gets proteins from tidied data frame into dictionary for OI

        Parameters
        ----------
        tdf : Data frame
            Data frame indexed by Gene with 'Patient' coluomn representing patient name, and 'diffFromMean' column

        Returns
        -------
        Dictionary of dictionaries 

        """
        selvals=tdf[abs(tdf[column])>absthresh]
        
        dvals = selvals.groupby('Patient')['Gene','logratio'].apply(dict)
         
        res = dict()
        for k,vals in dvals.iteritems():
            res[k] = dict(zip(vals['Gene'],vals['logratio']))
            
        return(res)