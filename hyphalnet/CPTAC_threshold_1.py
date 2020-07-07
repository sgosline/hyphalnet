# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 10:37:37 2020

@author: linth
"""
import pandas as pd
from scipy import stats
data = pd.read_csv('CPTAC.csv')
data.info()
data.isnull().sum()
data.head()
# Creating protein dictionary
aml_sample = data['AML sample']
gene = data['Gene']
log_fold = data['LogFoldChange']
protein_dict = dict()
n = len(aml_sample)
for i in range(n):
    protein_dict[aml_sample[i]]={gene[i]:log_fold[i]}
protein_dict
 #Calculating z score
data['zscore'] = stats.zscore(data['LogFoldChange'])
data.head()
data['pvalue'] = stats.norm.cdf(data['zscore'])
significant_genes = data[data['pvalue'] < 0.01]
non_significant = data[data['pvalue'] > 0.01]

non_significant.head()
significant_genes.head()