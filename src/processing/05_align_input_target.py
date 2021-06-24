import os
import sys
import numpy as np
import pandas as pd

## Define directories
ctrp_dn = '/projects/b1131/saya/drug_response/ctrp'
ccle_dn = '/projects/b1131/saya/drug_response/ccle'
inter_dn = '/projects/b1131/saya/drug_response/intermediary_data'
processed_dn = '/projects/b1131/saya/drug_response/processed_data'

ctrp_meta = pd.read_csv(os.path.join(ctrp_dn, 'v20.meta.per_cell_line.txt'), sep='\t')
cpd_meta = pd.read_csv(os.path.join(ctrp_dn, 'v20.meta.per_compound.txt'), sep='\t')
ccle_anno = pd.read_csv(os.path.join(ccle_dn, 'Cell_lines_annotations_20181226.txt'), sep='\t')

target_bin = pd.read_csv(os.path.join(inter_dn, 'ctrp_response_binary.csv'), index_col=0)
target_auc = pd.read_csv(os.path.join(inter_dn, 'ctrp_response_numeric.csv'), index_col=0)

vc = pd.read_csv(os.path.join(inter_dn, 'vc_mx_excludesilent.csv'), index_col=0)
exp = pd.read_csv(os.path.join(inter_dn, 'expression_mx.csv'), index_col=0)

###########################################################
#### Create lookup table for cell line identifications ####
###########################################################

ctrp_ids = ctrp_meta[['master_ccl_id', 'ccl_name']]
ccle_ids = ccle_anno[['CCLE_ID', 'depMapID']]

ctrp_ids = pd.concat(
    (ctrp_ids, 
    pd.DataFrame([ctrp_ids.loc[i,'ccl_name'].lower() for i in ctrp_ids.index], index=ctrp_ids.index, columns=['ccl_name_lower'])), 
    axis=1
)
ccle_ids = pd.concat(
    (ccle_ids,
    pd.DataFrame([ccle_ids.loc[i,'CCLE_ID'].split("_")[0].lower() for i in ccle_ids.index], index=ccle_ids.index, columns=['ccl_name_lower'])),
    axis=1
)

lookup = ctrp_ids.set_index('ccl_name_lower').join(ccle_ids.set_index('ccl_name_lower'), how='outer').reset_index()

## Show cases where relationship is not one-to-one
lookup.iloc[[np.sum(lookup['ccl_name']==x)>1 for x in lookup['ccl_name'].values],:]
lookup.iloc[[np.sum(lookup['master_ccl_id']==x)>1 for x in lookup['master_ccl_id'].values],:]

###############################
#### Get the intersect IDs ####
###############################

intersect = pd.DataFrame(columns=lookup.columns)

for i in lookup.index:
    has_target = (lookup.loc[i,'master_ccl_id'] in target_bin.index)
    has_exp = (lookup.loc[i,'CCLE_ID'] in exp.index)
    has_vc = (lookup.loc[i,'depMapID'] in vc.index)
    if (has_target & has_exp & has_vc):
        intersect = intersect.append(lookup.loc[i,:])

intersect.reset_index(drop=True, inplace=True)

###############################################################
#### Re-index feature and response matrices using depMapID ####
###############################################################

## For VC, just narrow down rows based on the intersect
vc_reindexed = vc.iloc[[item in intersect['depMapID'].values for item in vc.index],:]
vc_reindexed.sort_index(inplace=True)
vc_reindexed.to_csv(os.path.join(processed_dn, 'vc_mx_aligned.csv'))

## For RNA expression, reindex using the lookup table
# first narrow down based on intersect 
exp_reindexed = exp.iloc[[item in intersect['CCLE_ID'].values for item in exp.index],:]
# next, create a new index column and reindex with that
exp_reindexed['depMapID'] = None
for item in exp_reindexed.index:
    depmapid = lookup.iloc[item==lookup['CCLE_ID'].values]['depMapID'].values[0]
    exp_reindexed.loc[item,'depMapID'] = depmapid

exp_reindexed.set_index(['depMapID'], drop=True, inplace=True)
exp_reindexed.sort_index(inplace=True)

exp_reindexed.to_csv(os.path.join(processed_dn, 'expression_mx_aligned.csv'))

## For response variable 

# first narrow down based on intersect
target_bin = target_bin.iloc[[item in intersect['master_ccl_id'].values for item in target_bin.index],:]
target_auc = target_auc.iloc[[item in intersect['master_ccl_id'].values for item in target_auc.index],:]

# next, create a new index column and reindex with that
target_bin['depMapID'] = None
target_auc['depMapID'] = None
for item in target_bin.index: # they both have the same indices 
    depmapid = lookup.iloc[item==lookup['master_ccl_id'].values]['depMapID'].values[0]
    target_bin.loc[item,'depMapID'] = depmapid
    target_auc.loc[item,'depMapID'] = depmapid

target_bin.set_index(['depMapID'], drop=True, inplace=True)
target_auc.set_index(['depMapID'], drop=True, inplace=True)

target_bin.sort_index(inplace=True)
target_auc.sort_index(inplace=True)

target_bin.to_csv(os.path.join(processed_dn, 'response_binary_aligned.csv.csv'))
target_auc.to_csv(os.path.join(processed_dn, 'response_continuous_aligned.csv'))
