import os
import numpy as np
import pandas as pd
import glob
from datetime import datetime

############################
#### Define directories ####
############################

din='/projects/b1131/saya/drug_response/depmap/00_raw'
dout='/projects/b1131/saya/drug_response/depmap/01_cleaned'

###########################################
#### Generate the variant count matrix ####
###########################################

# vc_fn='CCLE_mutations.csv'

# ## Define function 
# def write_vc_mx(fin, fout, exclude_class=['Silent'], exclude_type=[]):
#     '''
#     This function allows you to create a variant count matrix based on the CCLE DepMap mutaion call file.
#     exclude_class: specifies the Variant_Classification to be excluded from the count as a list of strings.
#         Possible classifications: 'Silent', 'Missense_Mutation', 'Splice_Site', 'Nonsense_Mutation', etc. (there's more)
#     exclude_type: specifies the Variant_Type to be excluded from the count as a list of strings.
#         Possible types: 'SNP', 'DNP', 'TNP', 'ONP', 'DEL', 'INS' (this is all)
#     '''
#     depmap = pd.read_csv(fin)
#     mx = pd.DataFrame(None)
#     # Loop through rows and record variant counts
#     for i in range(depmap.shape[0]):
#         varclass = depmap['Variant_Classification'].iloc[i]
#         vartype = depmap['Variant_Type'].iloc[i]
#         if ((varclass in exclude_class) | (vartype in exclude_type)):
#             continue
#         else:
#             gene = depmap['Hugo_Symbol'].iloc[i]
#             ccl = depmap['DepMap_ID'].iloc[i]
#             if gene in list(mx.columns):
#                 if ccl in list(mx.index):
#                     mx.loc[ccl][gene] += 1
#                 else:
#                     mx = mx.append(pd.DataFrame(0, index=[ccl], columns=mx.columns))
#                     mx.loc[ccl][gene] += 1
#             else:
#                 mx[gene] = 0 # add new column with default value zero 
#                 if ccl in list(mx.index):
#                     mx.loc[ccl][gene] += 1
#                 else:
#                     mx = mx.append(pd.DataFrame(0, index=[ccl], columns=mx.columns))
#                     mx.loc[ccl][gene] += 1
#     mx.to_csv(fout, header=True, index=True)
#     return 


# write_vc_mx(
#     fin = f'{din}/{vc_fn}', 
#     fout = f'{dout}/cleaned_vc_mx.csv'
# )

########################################
#### Generate the expression matrix ####
########################################

exp_fn='CCLE_expression.csv'

def write_exp_mx(fin, fout):
    '''
    This function allows you to create a RNA expression matrix based on the CCLE expression file.
    '''
    exp = pd.read_csv(fin, index_col=0) # rows = genes // columns = tissue or sample
    # exp = exp.T.drop('Name', axis=0) # dropping 'ENSG<num*>' gene names
    exp.to_csv(fout, header=True, index=True)
    return

write_exp_mx(
    fin = f'{din}/{exp_fn}', 
    fout = f'{dout}/cleaned_exp_mx.csv'
)
