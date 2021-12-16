import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt

dn_ctrp='/projects/b1131/saya/drug_response/ctrp'
dn_rt='/projects/b1122/saya/reactome_pathways'
dn_out='/projects/b1131/saya/drug_response/processed_data'

meta_cpd = pd.read_csv(f'{dn_ctrp}/v20.meta.per_compound.txt', sep='\t')

#### Load the gene-to-pathway dictionary ####
with open(f'{dn_rt}/pathway_gene_dict.pkl', 'rb') as f:
    pathways = pickle.load(f)

##########################################
#### Features based on target protein ####
##########################################

# meta_cpd['gene_symbol_of_protein_target'] # 294 unique values, each contain multiple genes separated by semicolons 
prot_targets = {}
for i in meta_cpd.index:
    cpd_id = meta_cpd.loc[i,'master_cpd_id']
    targets = meta_cpd.loc[i,'gene_symbol_of_protein_target']
    if type(targets)==str:
        prot_targets[cpd_id] = targets.split(';')
    elif targets:
        prot_targets[cpd_id] = targets
    else:
        prot_targets[cpd_id] = None

cpd_feature = pd.DataFrame(0, index=meta_cpd['master_cpd_id'].values, columns=pathways.keys())

for cpd_id in prot_targets.keys():
    if type(prot_targets[cpd_id])==list:
        for target_gene in prot_targets[cpd_id]:
            for pathway in pathways.keys():
                if target_gene in pathways[pathway]:
                    cpd_feature.loc[cpd_id, pathway] = 1

cpd_feature = cpd_feature.iloc[:,cpd_feature.sum(axis=0).values!=0]

cpd_feature.to_csv(f'{dn_out}/cpd_features_reactome_pathways.csv', header=True, index=True)

fig, ax = plt.subplots()
ax.hist(cpd_feature.sum(axis=1), bins=20)
ax.set_ylabel('Number of compounds')
ax.set_xlabel('Number of pathways')
ax.set_title('Number of pathways associated with a given CPD')
fig.savefig('/projects/b1131/saya/drug_response/data_summary/histogram_cpd_feature_reactome_pathway.png')

###########################################
#### Features based on target activity ####
###########################################

meta_cpd['target_or_activity_of_compound'] # also sometimes multiple entries separated by semicolon 
no_target = [
    'screening hit',
    'natural product',
    'product of diversity oriented synthesis',
]
target_wordmatch = [
    'activator of ',
    'agonist of ',
    'analog of ',
    'antagonist of ',
    'inhibitor of ',
    'inducer of ',
    'modulator of ',
    'partial agonist of ',
    'promoter of ',
    're-activator of ',
]

for i in meta_cpd.index:
    match = []
    for phrase in target_wordmatch:
        if phrase in meta_cpd.loc[i,'target_or_activity_of_compound']:
            match.append(True)
        else:
            match.append(False)
    if ~np.any(match):
        print(meta_cpd.loc[i,'target_or_activity_of_compound'])



target_list = list(meta_cpd['target_or_activity_of_compound'].unique())
uniq_tar = []
for item in target_list:
    if ';' in item:
        indv_elements = item.split(';')
        for indv_element in indv_elements:
            if indv_element not in uniq_tar:
                uniq_tar.append(indv_element)
    else:
        uniq_tar.append(item)

sorted(uniq_tar)