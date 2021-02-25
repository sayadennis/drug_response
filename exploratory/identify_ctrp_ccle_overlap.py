import os
import sys
import numpy as np
import pandas as pd

ctrp_meta = pd.read_csv("/projects/b1131/saya/ctrp/v20.meta.per_cell_line.txt", sep="\t")
ccle_anno = pd.read_csv("/projects/b1131/saya/ccle/Cell_lines_annotations_20181226.txt", sep="\t")

## First, count overlapping cell lines in CTRP and CCLE 

ctrp_cclname = [name.lower() for name in list(ctrp_meta["ccl_name"])]
ccle_id = [name.lower() for name in list(ccle_anno["CCLE_ID"])]

ccle_id_split = {}

for ccle_name in ccle_id:
    identifier = ccle_name.split("_")[0]
    site = "".join(str(elem)+"_" for elem in ccle_name.split("_")[1:])[:-1]
    ccle_id_split[identifier] = site

ct_presence = 0
bool_presence = [] # boolean list indicating which CCLE cell lines are present in CTRP 
ct_sitematch = 0
only_ccle = 0

for key in ccle_id_split.keys(): # loop through CCLE
    if key in ctrp_cclname: # if CCLE ID is in CTRP
        ct_presence += 1 # count presence match 
        bool_presence.append(True)
        ctrp_site = ctrp_meta["ccle_primary_site"].iloc[[x == key for x in ctrp_cclname]] # record site name in CTRP 
        if ccle_id_split[key] == ctrp_site.values[0]: # if CTRP site name matches CCLE's site name
            ct_sitematch += 1 # count site match 
        else:
            print("Site doesn't match in cell line {}: \"{}\" vs. \"{}\"".format(key, ccle_id_split[key], ctrp_site.values[0]))
    else:
        only_ccle += 1
        bool_presence.append(False)

print("")
print("Number of cell lines annotated in CTRP: {}".format(len(ctrp_meta["ccl_name"].unique())))
print("Number of cell lines annotated in CCLE: {}".format(len(ccle_anno["CCLE_ID"].unique())))
print("")
print("Number of cell lines in both CTRP and CCLE: {}".format(ct_presence))
print("Number of cell lines with mismatch site information (including NA): {}".format(ct_presence-ct_sitematch))
print("")


## Do CTRP cell lines have CCLE mutational information? 
ccle_mut = pd.read_csv("/projects/b1131/saya/ccle/CCLE_DepMap_18q3_maf_20180718.txt", sep="\t")

dep_ids = ccle_anno["depMapID"].iloc[bool_presence]
bool_mut = [] # boolean indicating whether the cell line in the mutation table is present in CTRP
mut_ccl = ccle_mut["Broad_ID"].unique()
for dep_id in dep_ids:
    if dep_id in mut_ccl:
        bool_mut.append(True)
    else:
        bool_mut.append(False)

print("Number of CTRP cell lines that appear in mutation data: {}".format(sum(bool_mut)))
print("")

## how many rows per cell line are there in the mutational information table? 
num_rows = [] # number of mutation rows/observations that I found in 
for dep_id in dep_ids[bool_mut]:
    num_rows.append(len(ccle_mut[ccle_mut["Broad_ID"] == dep_id]))

quantile_names = ["min", "25%", "median", "75%", "max"]
quantile_nums = [0, 25, 50, 75, 100]

print("Distribution of the number of mutations per cell line:")
for i in range(len(quantile_names)):
    print(quantile_names[i], ":", np.percentile(num_rows, quantile_nums[i]))
print("")
