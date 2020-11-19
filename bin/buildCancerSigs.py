'''
This script uses the cptac data package and to build a reference set of signatures
We also compute the distances between patients and the communities from the reference signatures
and compute the GO enrichment
'''

import pandas as pd
import os
import hyphalnet.hypha as hyp
from hyphalnet.hypha import hyphalNetwork
import hyphalnet.hyphEnrich as hyEnrich
import hyphalnet.hyphaeStats as hyStats
import pickle
import argparse
import gzip
import cptac
import re
import numpy as np

parser = argparse.ArgumentParser(description="""Build hyphal network signatures for pan-can data""")
parser.add_argument('--ref', dest='refid', default='syn22392951', help='Synapse id of reference pkl')
parser.add_argument('--refName', dest='refName', default='CPTACpancan')
parser.add_argument('--synProj',dest='synProj',help='id of synapse project to store result')
parser.add_argument('--quantile',dest='quant', default=0.01,  help='Threshold to use for top-expressed proteins')

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

def getCancerData(qThresh=0.01):

    '''
    Get cancer data from cptac, computes differentially expressed proteins between
    tumor and pooled normal samples
    '''
    ##these are the primary dictionaries to return
    tumNormDiffs = {} #dictionary of top `qThresh`
    patientData = {} #dictionary of metadata

    ##make sur we have all the data
    #for each tumor type, colect data
    print("Downloading cptac data")
    allcans=['brca', 'ccrcc', 'colon', 'ovarian', 'luad',\
               #'hnscc','gbm','lscc',\
               'endometrial']
    #we need to make sure all datasets are downloaded
    for ct in allcans:
        cptac.download(dataset=ct)

    print("Loading cptac datasets")
    #then we load them into a dictionary
    fdict = {'brca':cptac.Brca(), 'ccrcc':cptac.Ccrcc(),\
           'colon':cptac.Colon(), 'ovarian':cptac.Ovarian(),\
             #'hnscc':cptac.Hnscc(),'gbm':cptac.Gbm(), 'lscc':cptac.Lscc(),\
           'endometrial':cptac.Endometrial(), 'luad':cptac.Luad()}

    patientData = getCombinedClinicalData(fdict)
    patientData.to_csv('clinicalData.csv')
    for ct in allcans:
        print("Collecting "+ct+' tumor normal')
        dat = fdict[ct]
        cdf = dat.join_metadata_to_omics(metadata_df_name="clinical", \
                                         omics_df_name="proteomics", \
                                         metadata_cols=['Sample_Tumor_Normal'])
        tn = cdf['Sample_Tumor_Normal']

        prots = [a for a in cdf.columns if 'proteomics' in a[0]]
        if len(prots)==0:
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
            tumNormDiffs[pval] = dict(zip(pat_prots[pval], pat_vals[pval]))
            #patientData[pval] = {'disease':ct}
        print('Added '+str(len(pat_prots))+' values to dictionary')
    return tumNormDiffs, patientData

def main():
    args = parser.parse_args()
    qval = args.quant
    allDat, patientData = getCancerData(qval)

    g = pickle.load(open('../odata/igraphPPI.pkl','rb'))
    beta = .5

    hDict = {'panCan' : hyphalNetwork(allDat, g, beta)}
    res = hyStats.compute_all_distances(hDict)
    ##eventually do GO enrichment

    fname = 'hyp2hypDistances_'+args.refName+'_to_data.csv'
    res.to_csv(fname)
#    if args.synProj is not None:
#        tab = sc.table.build_table("Original data to "+args.refName+' Distances',
#                                   args.synProj, res)
#        syn.store(tab)


if __name__=='__main__':
    main()
