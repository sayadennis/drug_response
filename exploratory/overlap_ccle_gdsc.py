import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

dn='/projects/b1131/saya/drug_response'

#####################################
#### Visualize cell line overlap ####
#####################################

ccle = pd.read_csv(f'{dn}/ccle/Cell_lines_annotations_20181226.txt', sep='\t')
gdsc = pd.read_csv(f'{dn}/gdsc/cell_lines_genomic/mutations_20191101.csv')
# pdx = pd.read_csv(f'{dn}/pdx/pdxe_mut_and_cn2.tsv', sep='\t')

# remove all spaces, hyphens, and underscores 
ccle_list = [re.sub(r'[-_ ]', '', x) for x in list(ccle['Name'].dropna().unique())]
gdsc_list = [re.sub(r'[-_ ]', '', x) for x in list(gdsc['model_name'].dropna().unique())]

ct=0
for item in ccle_list:
    if item in gdsc_list:
        ct+=1

venn2(subsets = (len(ccle_list)-ct, len(gdsc_list)-ct, ct), set_labels = ('CCLE', 'GDSC'))
plt.savefig(f'{dn}/data_summary/ccle_gdsc_cell_line_overlap_venn.png')

#########################################################
#### Align the RNA-seq data between the two datasets ####
#########################################################

ccle = pd.read_csv(f'{dn}/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct', sep='\t', skiprows=2) # 56202 genes 
gdsc = pd.read_csv(f'{dn}/gdsc/cell_lines_genomic/rnaseq_fpkm_20191101.csv') # 37282 genes 

ccle_genes = list(gdsc['Unnamed: 1'].dropna().unique())
gdsc_genes = list(ccle['Description'].dropna().unique())

## Assess the overlap in genes 
ol=[]
for item in gdsc_genes:
    if item in ccle_genes:
        ol.append(item)

print(f'{len(ol)} genes overlap between CCLE and GDSC. Keeping only these overlapping genes in following analysis.')

gdsc.iloc[0,1] = 'gene'; gdsc.columns = gdsc.iloc[0,:].values; gdsc.drop([0,1,2], axis=0, inplace=True)

gdsc = gdsc.iloc[[x in ol for x in gdsc['gene'].values],:]
ccle = ccle.iloc[[x in ol for x in ccle['Description'].values],:]

## Reformat row indices (gene names) to make the two matrices consistent 
# GDSC 
gdsc.set_index('gene', inplace=True, verify_integrity=True)
gdsc.drop('model_name', axis=1, inplace=True); gdsc.index.name = None

# CCLE
ccle.drop('Name', axis=1, inplace=True); ccle.set_index('Description', inplace=True); ccle.index.name = None # , verify_integrity=True 
ccle_singles = ccle.iloc[~ccle.index.duplicated(False),:]
dups = ccle.iloc[ccle.index.duplicated(False),:]
for gene in dups.index.unique():
    gene_expression = dups.loc[gene,:]
    ccle_singles = ccle_singles.append(pd.Series(gene_expression.median(axis=0), name=gene))

ccle = ccle_singles

## Reformat column names (cell line names) to make the two matrices consistent 
# GDSC 
gdsc.columns = [re.sub(r'[-_ ]', '', x.upper()) for x in gdsc.columns]
# CCLE
ccle.columns = [re.sub(r'[-_ ]', '', x.split('_')[0].upper()) for x in ccle.columns] # split because each column name is like "22RV1_PROSTATE"

# Assess cell line overlap again 
ol = []
gdsc_cl = gdsc.columns
ccle_cl = ccle.columns
for item in gdsc_cl:
    if item in ccle_cl:
        ol.append(item)
#

def unique(orig_list):
    unique_list = []
    for x in orig_list:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    return unique_list

ol = unique(ol)

gdsc = gdsc.loc[:,ol]
ccle = ccle.loc[:,ol]
#
gdsc.to_csv(f'{dn}/intermediary_data/gdsc_rna_aligned_to_ccle.csv')
ccle.to_csv(f'{dn}/intermediary_data/ccle_rna_aligned_to_ccle.csv')
