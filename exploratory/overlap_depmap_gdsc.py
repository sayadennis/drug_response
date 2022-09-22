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

dpmp = pd.read_csv('/projects/b1131/saya/drug_response/depmap/00_raw/sample_info.csv')
gdsc = pd.read_csv(f'{dn}/gdsc/cell_lines_genomic/mutations_20191101.csv')
# pdx = pd.read_csv(f'{dn}/pdx/pdxe_mut_and_cn2.tsv', sep='\t')

# remove all spaces, hyphens, and underscores 
dpmp_list = [re.sub(r'[-_ ]', '', x) for x in list(dpmp['stripped_cell_line_name'].dropna().unique())]
gdsc_list = [re.sub(r'[-_ ]', '', x) for x in list(gdsc['model_name'].dropna().unique())]

ct=0
for item in dpmp_list:
    if item in gdsc_list:
        ct+=1

venn2(subsets = (len(dpmp_list)-ct, len(gdsc_list)-ct, ct), set_labels = ('DepMap', 'GDSC'))
plt.savefig(f'{dn}/data_summary/depmap_gdsc_cell_line_overlap_venn.png')

#########################################################
#### Align the RNA-seq data between the two datasets ####
#########################################################

## Load expression data
dpmp = pd.read_csv(f'{dn}/depmap/00_raw/CCLE_expression.csv', index_col=0)
gdsc = pd.read_csv(f'{dn}/gdsc/cell_lines_genomic/rnaseq_fpkm_20191101.csv') # 37282 genes 

## Rename DepMap's expression matrix's cell lines (rows) from the "ACH-001113" IDs to the cell line names 
dp_meta = pd.read_csv('/projects/b1131/saya/drug_response/depmap/00_raw/sample_info.csv')

cell_names = dict(zip(dp_meta.DepMap_ID, dp_meta.stripped_cell_line_name))
dpmp = dpmp.rename(cell_names, axis=0)

## Rename DepMap's expression matrix's gene names (columns) from something like "TSPAN6 (7105)" to "TSPAN6"
dpmp.columns = [x.split()[0] for x in dpmp.columns]

## Find gene list for each matrix 
gdsc_genes = list(gdsc['Unnamed: 1'].dropna().unique())
dpmp_genes = list(dpmp.columns)

## Assess the overlap in genes 
ol=[]
for item in gdsc_genes:
    if item in dpmp_genes:
        ol.append(item)

print(f'{len(ol)} genes overlap between CCLE and GDSC. Keeping only these overlapping genes in following analysis.')

gdsc.iloc[0,1] = 'gene'; gdsc.columns = gdsc.iloc[0,:].values; gdsc.drop([0,1,2], axis=0, inplace=True)

gdsc = gdsc.iloc[[x in ol for x in gdsc['gene'].values],:]
dpmp = dpmp.iloc[:,[x in ol for x in dpmp.columns]]

## Reformat row indices (gene names) to make the two matrices consistent 
# GDSC 
gdsc.set_index('gene', inplace=True, verify_integrity=True)
gdsc.drop('model_name', axis=1, inplace=True); gdsc.index.name = None

# DepMap
# ccle.drop('Name', axis=1, inplace=True); ccle.set_index('Description', inplace=True); ccle.index.name = None # , verify_integrity=True 
dpmp = dpmp.T

## Reformat column names (cell line names) to make the two matrices consistent 
# GDSC 
gdsc.columns = [re.sub(r'[-_ ]', '', x.upper()) for x in gdsc.columns]
# DepMap
dpmp.columns = [re.sub(r'[-_ ]', '', x.upper()) for x in dpmp.columns]

# Assess cell line overlap again 
ol = []
gdsc_cl = gdsc.columns
dpmp_cl = dpmp.columns
for item in gdsc_cl:
    if item in dpmp_cl:
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
dpmp = dpmp.loc[:,ol]
#

# dpmp_singles = dpmp.iloc[:,~dpmp.columns.duplicated(False)]
# dups = dpmp.iloc[dpmp.index.duplicated(False),:]
# for gene in dups.index.unique():
#     gene_expression = dups.loc[gene,:]
#     dpmp_singles = dpmp_singles.append(pd.Series(gene_expression.median(axis=0), name=gene))

# dpmp = dpmp_singles


gdsc.to_csv(f'{dn}/intermediary_data/gdsc_rna_aligned_to_depmap.csv')
dpmp.to_csv(f'{dn}/intermediary_data/depmap_rna_aligned_to_gdsc.csv')
