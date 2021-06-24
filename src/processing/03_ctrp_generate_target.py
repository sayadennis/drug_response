import os
import sys
import numpy as np
import pandas as pd

din = "/projects/b1131/saya/drug_response/ctrp"
dout = "/projects/b1131/saya/drug_response/processed_data"

meta_exp = pd.read_csv(os.path.join(din, "v20.meta.per_experiment.txt"), sep="\t", header=0)
curves = pd.read_csv(os.path.join(din, "v20.data.curves_post_qc.txt"), sep="\t")

# Join the experiment metadata and response curves
merged = meta_exp.loc[:][["experiment_id", "run_id", "master_ccl_id"]].merge(
    curves.loc[:][["experiment_id", "area_under_curve", "master_cpd_id"]], on="experiment_id", how="outer"
)

# Save joined table to CSV
merged.loc[:][["master_ccl_id", "area_under_curve", "master_cpd_id"]].to_csv(
    os.path.join(dout, "ctrp_response_all_joined.csv"), header=True, index=False
)

# Next, make a matrix with row = cell lines, columns = compound
response_mx_num = pd.DataFrame(None, index=merged["master_ccl_id"].unique(), columns=merged["master_cpd_id"].unique())
response_mx_bin = pd.DataFrame(None, index=merged["master_ccl_id"].unique(), columns=merged["master_cpd_id"].unique())

for ccl_id in merged["master_ccl_id"].unique():
    for cpd_id in merged["master_cpd_id"].unique():
        response = merged["area_under_curve"].iloc[
            (np.array([x == ccl_id for x in merged["master_ccl_id"]]) & np.array([x == cpd_id for x in merged["master_cpd_id"]]))
        ].values
        if len(response) == 0:
            print("can't find response value at CCL {} and CPD {}".format(ccl_id, cpd_id))
            response = 0
        else:
            response = response[0]
            response_mx_num.loc[ccl_id][cpd_id] = response
            response_mx_bin.loc[ccl_id][cpd_id] = int(response > 9)

response_mx_num.to_csv(os.path.join(dout, "ctrp_response_numeric.csv"), header=True, index=True)
response_mx_bin.to_csv(os.path.join(dout, "ctrp_response_binary.csv"), header=True, index=True)
