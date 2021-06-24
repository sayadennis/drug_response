import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ctrp_dn = '/projects/b1131/saya/drug_response/ctrp'
summary_dn = '/projects/b1131/saya/drug_response/data_summary'

#########################################################
#### Count how many cell lines come from what tissue ####
#########################################################

meta_cpd = pd.read_csv(os.path.join(ctrp_dn, 'v20.meta.per_compound.txt'), sep='\t', header=0)
meta_ccl = pd.read_csv(os.path.join(ctrp_dn, 'v20.meta.per_cell_line.txt'), sep='\t', header=0)
meta_exp = pd.read_csv(os.path.join(ctrp_dn, 'v20.meta.per_experiment.txt'), sep='\t', header=0)
meta_asy = pd.read_csv(os.path.join(ctrp_dn, 'v20.meta.per_assay_plate.txt'), sep='\t', header=0)
meta_med = pd.read_csv(os.path.join(ctrp_dn, 'v20.meta.media_comp.txt'), sep='\t', header=0)

data = pd.read_csv(os.path.join(ctrp_dn, 'v20.data.curves_post_qc.txt'), sep='\t', header=0)

summarytab = pd.DataFrame(
    None, 
    columns=["Primary site", "# cell lines", "# unique hist", "# unique cpd tested"]
)

# fill in the values for primary sites and number of cell lines 
summarytab["Primary site"] = meta_ccl["ccle_primary_site"].value_counts().index
summarytab["# cell lines"] = meta_ccl["ccle_primary_site"].value_counts().values
summarytab["# unique hist"] = meta_ccl["ccle_primary_site"].value_counts().values

# fill in the numbers of unique histology 
for site in meta_ccl["ccle_primary_site"].value_counts().index:
    subtab = meta_ccl[meta_ccl["ccle_primary_site"] == site]
    ct = len(subtab["ccle_primary_hist"].unique())
    summarytab.iloc[np.where(summarytab["Primary site"]==site)[0][0], np.where(summarytab.columns=="# unique hist")[0][0]] = ct

# fill in the numbers of unique compounds tested 
for site in meta_ccl["ccle_primary_site"].value_counts().index:
    ccl_ids = meta_ccl["master_ccl_id"].loc[meta_ccl["ccle_primary_site"]==site].values # get all CCL IDs of that site 
    exp_ids = meta_exp["experiment_id"].loc[[meta_exp["master_ccl_id"][i] in ccl_ids for i in range(meta_exp.shape[0])]].values # get all experiment IDs of that site
    cpd_ids = np.unique(data["master_cpd_id"].loc[[data["experiment_id"][i] in exp_ids for i in range(data.shape[0])]].values)
    ct = len(cpd_ids)
    summarytab.iloc[np.where(summarytab["Primary site"]==site)[0][0], np.where(summarytab.columns=="# unique cpd tested")[0][0]] = ct

summarytab.to_csv(os.path.join(summary_dn, 'cell_line_meta_summary.csv'), index=False, header=True)

pd.DataFrame(meta_ccl.groupby(["ccle_primary_site", "ccle_primary_hist"]).size()).to_csv(
    os.path.join(summary_dn, 'cell_line_meta_summary_grouped.csv') #, index=False, header=True
)

#####################################################
#### Plot histogram of cell viability AUC values ####
#####################################################

# Identify "experiment_id" where the cell line comes from breast cancer 
def find_experiment_ids(primary_site):
    ccl_ids = meta_ccl["master_ccl_id"].iloc[[x == primary_site for x in meta_ccl["ccle_primary_site"]]].values
    experiment_ids = meta_exp["experiment_id"].iloc[[x in ccl_ids for x in meta_exp["master_ccl_id"]]].unique()
    return experiment_ids

primary_site="breast"
breast_exp_ids = find_experiment_ids(primary_site)

curves = pd.read_csv(os.path.join(ctrp_dn, 'v20.data.curves_post_qc.txt'), sep='\t', header=0)
# select experiments that used breast cancer cell lines 
filt_curves = curves[[x in breast_exp_ids for x in curves["experiment_id"]]]

# integrated area under the sigmoid-fit concentration-response curve
plt.figure(figsize=(8,6))
plt.hist(filt_curves["area_under_curve"], bins=20)
plt.xlabel("Area under curve", size=16)
plt.ylabel("Counts", size=16)
plt.title("Histogram of integrated area under sigmoid-fit concentration-response curve", size=18)

plt.savefig(os.path.join(summary_dn, 'ctrp_response_auc_histogram.png'))
