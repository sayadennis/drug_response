import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

dn='/projects/b1131/saya/drug_response'

gdsc = pd.read_csv(f'{dn}/intermediary_data/gdsc_rna_aligned_to_ccle.csv', index_col=0)
ccle = pd.read_csv(f'{dn}/intermediary_data/ccle_rna_aligned_to_ccle.csv', index_col=0)

#############################################################################
#### Create narrowed down matrices using COSMIC genes and variable genes ####
#############################################################################

# Narrow down gene list by COSMIC genes 
cosmic1 = list(pd.read_csv(f'{dn}/Census_all_2021-11-22_tier1.csv')['Gene Symbol'])
cosmic2 = list(pd.read_csv(f'{dn}/Census_all_2021-11-22_tier2.csv')['Gene Symbol'])
variable = list(gdsc.std(axis=1).iloc[np.argsort(-1*gdsc.std(axis=1))][:2000].index) + list(ccle.std(axis=1).iloc[np.argsort(-1*ccle.std(axis=1))][:2000].index)

def unique(orig_list):
    unique_list = []
    for x in orig_list:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    return unique_list

variable = unique(variable)
print('Number of genes that are in either GDSC or CCLE\'s top 2000 most variable genes:', len(variable))

gdsc_cosmic1 = gdsc.loc[[x in cosmic1 for x in gdsc.index],:]
ccle_cosmic1 = ccle.loc[[x in cosmic1 for x in ccle.index],:]

gdsc_cosmic2 = pd.concat((gdsc.loc[[x in cosmic2 for x in gdsc.index],:], gdsc.loc[[x in cosmic1 for x in gdsc.index],:]), axis=0)
ccle_cosmic2 = pd.concat((ccle.loc[[x in cosmic2 for x in ccle.index],:], ccle.loc[[x in cosmic1 for x in ccle.index],:]), axis=0) 

gdsc_var = gdsc.loc[[x in variable for x in gdsc.index],:]
ccle_var = ccle.loc[[x in variable for x in ccle.index],:]

###########################################################
#### Loop through cell lines, record Pearson corr coef ####
###########################################################

## Files and lists to record coefficients in 
f = open(f'{dn}/data_summary/pearson_corr_ccle_gdsc.csv', 'w')
f1 = open(f'{dn}/data_summary/pearson_corr_ccle_gdsc_cosmic1.csv', 'w')
f2 = open(f'{dn}/data_summary/pearson_corr_ccle_gdsc_cosmic1-2.csv', 'w')
v = open(f'{dn}/data_summary/pearson_corr_ccle_gdsc_variablegenes.csv', 'w')

pearson = []
pearson1 = []
pearson2 = []
pearsonvar = []

## Loop 
for cell_line in ccle.columns:
    # pearson.append(gdsc[cell_line].corr(ccle[cell_line], method='pearson'))
    gdsc_exp = gdsc[cell_line].to_numpy(dtype=float)
    ccle_exp = ccle[cell_line].to_numpy(dtype=float)
    #
    gdsc_exp1 = gdsc_cosmic1[cell_line].to_numpy(dtype=float)
    ccle_exp1 = ccle_cosmic1[cell_line].to_numpy(dtype=float)
    #
    gdsc_exp2 = gdsc_cosmic2[cell_line].to_numpy(dtype=float)
    ccle_exp2 = ccle_cosmic2[cell_line].to_numpy(dtype=float)
    #
    gdsc_expvar = gdsc_var[cell_line].to_numpy(dtype=float)
    ccle_expvar = ccle_var[cell_line].to_numpy(dtype=float)
    try:
        # all genes 
        pearson_corr = pearsonr(gdsc_exp[~(np.isnan(gdsc_exp)|np.isnan(ccle_exp))], ccle_exp[~(np.isnan(gdsc_exp)|np.isnan(ccle_exp))])[0]
        pearson.append(pearson_corr)
        f.write(f'{cell_line},{pearson_corr}\n')
        # COSMIC tier 1 genes 
        pearson_corr1 = pearsonr(gdsc_exp1[~(np.isnan(gdsc_exp1)|np.isnan(ccle_exp1))], ccle_exp1[~(np.isnan(gdsc_exp1)|np.isnan(ccle_exp1))])[0]
        pearson1.append(pearson_corr1)
        f1.write(f'{cell_line},{pearson_corr1}\n')
        # COSMIC tier 2 genes 
        pearson_corr2 = pearsonr(gdsc_exp2[~(np.isnan(gdsc_exp2)|np.isnan(ccle_exp2))], ccle_exp2[~(np.isnan(gdsc_exp2)|np.isnan(ccle_exp2))])[0]
        pearson2.append(pearson_corr2)
        f.write(f'{cell_line},{pearson_corr2}\n')
        # top variable genes 
        pearson_corrvar = pearsonr(gdsc_expvar[~(np.isnan(gdsc_expvar)|np.isnan(ccle_expvar))], ccle_expvar[~(np.isnan(gdsc_expvar)|np.isnan(ccle_expvar))])[0]
        pearsonvar.append(pearson_corrvar)
        f.write(f'{cell_line},{pearson_corrvar}\n')
    except:
        continue

f.close(); f1.close(); f2.close()
## Plot the expression correlation ##

# All genes 
fig, ax = plt.subplots()
ax.violinplot(np.array(pearson))
ax.set_title('Correlation of CCLE/GDSC expression with all genes')
ax.set_ylabel('Pearson coefficients')
plt.tight_layout()
fig.savefig(f'{dn}/data_summary/pearson_corr_ccle_gdsc_violinplot_all.png')
plt.close()

# COSMIC genes 
fig, ax = plt.subplots()
ax.violinplot((np.array(pearson1), np.array(pearson2)), positions=[1,2])
ax.set_title('Correlation of CCLE/GDSC expression with COSMIC genes')
ax.set_xlabel('COSMIC genes used')
ax.set_xticks([1,2])
ax.set_xticklabels(['Tier 1', 'Tier 1&2'], rotation = 45, ha="right")
ax.set_ylabel('Pearson coefficients')
plt.tight_layout()
fig.savefig(f'{dn}/data_summary/pearson_corr_ccle_gdsc_violinplot_cosmic.png')
plt.close()

# Top variable genes
fig, ax = plt.subplots()
ax.violinplot(np.array(pearsonvar))
ax.set_title('Correlation of CCLE/GDSC expression with top 2000 variable genes')
ax.set_ylabel('Pearson coefficients')
plt.tight_layout()
fig.savefig(f'{dn}/data_summary/pearson_corr_ccle_gdsc_violinplot_topvariable.png')
plt.close()

