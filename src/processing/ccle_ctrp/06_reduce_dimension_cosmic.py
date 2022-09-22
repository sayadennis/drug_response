import numpy as np
import pandas as pd

## Define directories
din='/projects/b1131/saya/drug_response/processed_data'
dcos='/projects/b1131/saya/drug_response/intermediary_data'
dout=din

vc=pd.read_csv(f'{din}/vc_mx_aligned.csv', index_col=0)
exp=pd.read_csv(f'{din}/expression_mx_aligned.csv', index_col=0)

tier1=pd.read_csv(f'{dcos}/Census_all_2021-11-22_tier1.csv')
tier2=pd.read_csv(f'{dcos}/Census_all_2021-11-22_tier2.csv')

##############################
#### Reduce to only Tier1 ####
##############################

vc1=vc.iloc[:,[x in list(tier1['Gene Symbol'].values) for x in vc.columns]]
exp1=exp.iloc[:,[x in list(tier1['Gene Symbol'].values) for x in exp.columns]]

vc1.to_csv(f'{dout}/vc_mx_aligned_cosmic_tier1.csv')
exp1.to_csv(f'{dout}/exp_mx_aligned_cosmic_tier1.csv')

###################################
#### Reduce to Tier 1 + Tier 2 ####
###################################

vc2=vc.iloc[:,[x in list(tier2['Gene Symbol'].values) for x in vc.columns]]
exp2=exp.iloc[:,[x in list(tier2['Gene Symbol'].values) for x in exp.columns]]

pd.concat((vc1, vc2), axis=1).fillna(0).to_csv(f'{dout}/vc_mx_aligned_cosmic_tier1and2.csv')
pd.concat((exp1, exp2), axis=1).fillna(0).to_csv(f'{dout}/exp_mx_aligned_cosmic_tier1and2.csv')
