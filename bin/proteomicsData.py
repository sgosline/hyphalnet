'''
This script uses the cptac data package and to build a reference set of signatures
We also compute the distances between patients and the communities from the reference signatures
and compute the GO enrichment
'''

import pandas as pd
import os
import cptac
import numpy as np
import re

def getCombinedClinicalData(datdict):
    '''
    for a list of cancer types, collect all age, gender, stage (if possible), and race(if possible)
    not sure what else we can harmonize at this point
    TODO; add survival
    '''

    ##this dictionary describes the clinical data mapping
    '''
    mvals = {'brca':{'age':'Age.in.Month', 'gender':'Gender',\
                     'race':'Race', 'stage':'Stage'},\
                 'colon':{'age':'Age', 'gender':'Gender', 'stage':'Stage'},\
                 'ccrcc':{'gender':'gender', 'age':'age', 'race':'race', 'stage':''},\
                 'endometrial':{'gender':'Gender', 'race':'Race',\
                                'age':'Age', 'stage':'FIGO_stage'},\
                 'ovarian':{'race':'Participant_Race', 'stage':'Tumor_Stage_Ovary_FIGO',\
                            'age':'Participant_Procurement_Age', 'gender':'Participant_Gender'},\
                 'luad':{'age':'Age', 'gender':'Gender', 'race':'Ethnicity', 'stage':'Stage'}}
    '''
    mvals = {'brca':{'Age.in.Month':'age', 'Gender':'gender',\
                      'Race':'race', 'Stage':'stage'},\
             'colon':{'Age':'age', 'Gender':'gender', 'Stage':'stage'},\
             'ccrcc':{'gender':'gender', 'age':'age', 'race':'race', '':'stage'},\
             'endometrial':{'Gender':'gender', 'Race':'race', 'Age':'age', 'FIGO_stage':'stage'},\
             'ovarian':{'Participant_Race':'race', 'Tumor_Stage_Ovary_FIGO':'stage', \
                        'Participant_Procurement_Age':'age', 'Participant_Gender':'gender'},\
             'luad': {'Age':'age', 'Gender':'gender', 'Ethnicity':'race', 'Stage':'stage'}}

    fulldf = []
    for ct, dat in datdict.items():
        cols = [a for a in mvals[ct].keys() if a!=""]
        print(cols)
        clindat = dat.get_clinical()[cols]
        clindat = clindat.reset_index()
        clindat['CancerType'] = ct
        clindat = clindat.rename(columns=mvals[ct])
        fulldf.append(clindat)

    patData = pd.concat(fulldf)
    return patData


def getCombinedMutationData(fdict):
    '''
    We will also need to collect mutation data for analysis this needs to be in
    a data frame that we can write out to CSV to join with other available data
    '''
    return pd.DataFrame()

def getCancerData(fdict, qThresh=0.01, byType=False):

    '''
    Get cancer data from cptac, computes differentially expressed proteins between
    tumor and pooled normal samples
    '''
    ##these are the primary dictionaries to return
    tumNormDiffs = {} #dictionary of top `qThresh`
    for ct, dat in fdict.items():
        if byType:
            tumNormDiffs[ct] = {}
        print("Collecting "+ct+' tumor normal')
        dat = fdict[ct]
        cdf = dat.join_metadata_to_omics(metadata_df_name="clinical", \
                                         omics_df_name="proteomics", \
                                         metadata_cols=['Sample_Tumor_Normal'])
        tn = cdf['Sample_Tumor_Normal']

        prots = [a for a in cdf.columns if 'proteomics' in a[0]]
        if len(prots) == 0:
            prots = [a for a in cdf.columns if 'proteomics' in a]

        norm_rows = cdf.loc[cdf['Sample_Tumor_Normal'] == 'Normal'][prots]
        norm_means = np.zeros(len(norm_rows.columns))

        if len(norm_rows.index) == 0:
            print("No normal data, using top expressed proteins instead")
        else:
            norm_means = norm_rows.mean(axis=0)

        tum_rows = cdf.loc[cdf['Sample_Tumor_Normal'] == 'Tumor'][prots]-norm_means
        tumcols = tum_rows.columns
        tum_rows['Patient'] = tum_rows.index
        #reshape
        long_tab = pd.melt(tum_rows,id_vars='Patient').assign(Gene = lambda dataframe: dataframe['Name'].map(lambda Name: re.sub('_proteomics','',Name)))

        #groupby, get quantile
        pquants = pd.DataFrame({'thresh':long_tab.groupby("Patient")['value'].quantile(1-qThresh)})
        full_dat = long_tab.merge(pquants, on='Patient')
        #select over quantile
        full_dat = full_dat.assign(topProt=full_dat['value']>full_dat['thresh'])
        sel_dat = full_dat[full_dat['topProt']]
        pat_prots = sel_dat.groupby('Patient')['Gene'].apply(list).to_dict()
        pat_vals = sel_dat.groupby("Patient")['value'].apply(list).to_dict()

        test_count = 0 ##add for debugging
        for pval in pat_prots.keys():
            if(test_count > 10):
                continue
            test_count = test_count+1 #remove this once we have it working
            if byType:
                tumNormDiffs[ct][pval] = dict(zip(pat_prots[pval], pat_vals[pval]))
            else:
                tumNormDiffs[pval] = dict(zip(pat_prots[pval], pat_vals[pval]))
            #patientData[pval] = {'disease':ct}
        print('Added '+str(len(pat_prots))+' values to dictionary')
    return tumNormDiffs

def cptacData():
    '''
    We need to collect and load CPTAC data
    '''
    print("Loading cptac datasets")
        #we need to make sure all datasets are downloaded
    ##here are the cancers that are available without login information
    allcans = ['brca', 'ccrcc', 'colon', 'ovarian', 'luad',\
             #'hnscc','gbm','lscc',\
             'endometrial']
    print("Downloading cptac data")
    for ct in allcans:
        cptac.download(dataset=ct)
    #then we load them into a dictionary
    fdict = {'brca':cptac.Brca(), 'ccrcc':cptac.Ccrcc(),\
           'colon':cptac.Colon(), 'ovarian':cptac.Ovarian(),\
             #'hnscc':cptac.Hnscc(),'gbm':cptac.Gbm(), 'lscc':cptac.Lscc(),\
           'endometrial':cptac.Endometrial(), 'luad':cptac.Luad()}
    return fdict
