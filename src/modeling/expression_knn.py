import os
import sys
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.neighbors import NearestNeighbors
# from sklearn.metrics import mean_absolute_error #, mean_absolute_percentage_error 
# from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, balanced_accuracy_score
# from sklearn.metrics import confusion_matrix

sys.path.append('drug_response/src/evaluation')
import evaluation_func

## Load data 
dn = '/projects/b1131/saya/drug_response/processed_data'

exp = pd.read_csv(os.path.join(dn, 'expression_mx_aligned.csv'), index_col=0)
auc = pd.read_csv(os.path.join(dn, 'response_continuous_aligned.csv'), index_col=0)
cpd_meta = pd.read_csv('/projects/b1131/saya/drug_response/ctrp/v20.meta.per_compound.txt', sep='\t')

## Train-test split 
exp_train, exp_test, auc_train, auc_test = train_test_split(exp, auc, test_size=0.2, random_state=42)

## Undersample the training data by randomly removing negative response values? 
keep_ix = np.random.randint(low=0, high=len(np.where(auc_train>9)[0]), size=len(np.where(auc_train<=9)[0])*4)  # indices of responses to keep (1:3 ratio)
neg_ixs = np.where(auc_train>9) # this is necessary in order to avoid error 
for ix in keep_ix:
    auc_train.iloc[neg_ixs[0][ix], neg_ixs[1][ix]] = None

## k-NN class 
knn = NearestNeighbors(n_neighbors=5, metric='cosine', algorithm='brute')
knn.fit(exp_train.values)
neigh_dist, neigh_ind = knn.kneighbors(exp_test)

## Record predictions on the test set 
auc_pred = pd.DataFrame(index=exp_test.index, columns=auc_test.columns)
for i in range(exp_test.shape[0]):
    auc_pred.iloc[i,:] = auc_train.iloc[neigh_ind[i,:]].mean(axis=0)

evaluation_func.evaluate(auc_test, auc_pred, outdir='drug_response/model_performance/exp_knn_all_cpd_balanced')

##################################################
#### Do the same with only FDA approved drugs ####
##################################################

fda_ids = np.array(cpd_meta['master_cpd_id'].iloc[cpd_meta['cpd_status'].values=='FDA'], dtype=str) # dtype has to be str

## Train-test split 
exp_train, exp_test, auc_train, auc_test = train_test_split(exp, auc.loc[:,fda_ids], test_size=0.2, random_state=42)

## Undersample the training data by randomly removing negative response values? 
keep_ix = np.random.randint(low=0, high=len(np.where(auc_train>9)[0]), size=len(np.where(auc_train<=9)[0])*4)  # indices of responses to keep (1:3 ratio)
neg_ixs = np.where(auc_train>9) # this is necessary in order to avoid error 
for ix in keep_ix:
    auc_train.iloc[neg_ixs[0][ix], neg_ixs[1][ix]] = None

## k-NN class 
knn = NearestNeighbors(n_neighbors=5, metric='cosine', algorithm='brute')
knn.fit(exp_train.values)
neigh_dist, neigh_ind = knn.kneighbors(exp_test)

## Record predictions on the test set 
auc_pred = pd.DataFrame(index=exp_test.index, columns=auc_test.columns)
for i in range(exp_test.shape[0]):
    auc_pred.iloc[i,:] = auc_train.iloc[neigh_ind[i,:]].mean(axis=0)

evaluation_func.evaluate(auc_test, auc_pred, outdir='drug_response/model_performance/exp_knn_fda_cpd_balanced')

#######################################
#### Try the COSMIC-reduced Tier 1 ####
#######################################

## Load data 
dn = '/projects/b1131/saya/drug_response/processed_data'
exp = pd.read_csv(os.path.join(dn, 'exp_mx_aligned_cosmic_tier1.csv'), index_col=0)

## Train-test split 
exp_train, exp_test, auc_train, auc_test = train_test_split(exp, auc.loc[:,fda_ids], test_size=0.2, random_state=42)

## Undersample the training data by randomly removing negative response values? 
keep_ix = np.random.randint(low=0, high=len(np.where(auc_train>9)[0]), size=len(np.where(auc_train<=9)[0])*4)  # indices of responses to keep (1:3 ratio)
neg_ixs = np.where(auc_train>9) # this is necessary in order to avoid error 
for ix in keep_ix:
    auc_train.iloc[neg_ixs[0][ix], neg_ixs[1][ix]] = None

## k-NN class 
knn = NearestNeighbors(n_neighbors=5, metric='cosine', algorithm='brute')
knn.fit(exp_train.values)
neigh_dist, neigh_ind = knn.kneighbors(exp_test)

## Record predictions on the test set 
auc_pred = pd.DataFrame(index=exp_test.index, columns=auc_test.columns)
for i in range(exp_test.shape[0]):
    auc_pred.iloc[i,:] = auc_train.iloc[neigh_ind[i,:]].mean(axis=0)

evaluation_func.evaluate(auc_test, auc_pred, outdir='drug_response/model_performance/exp_knn_fda_cpd_balanced_tier1')

###########################################
#### Try the COSMIC-reduced Tier 1 + 2 ####
###########################################

## Load data 
dn = '/projects/b1131/saya/drug_response/processed_data'
exp = pd.read_csv(os.path.join(dn, 'exp_mx_aligned_cosmic_tier1and2.csv'), index_col=0)

## Train-test split 
exp_train, exp_test, auc_train, auc_test = train_test_split(exp, auc.loc[:,fda_ids], test_size=0.2, random_state=42)

## Undersample the training data by randomly removing negative response values? 
keep_ix = np.random.randint(low=0, high=len(np.where(auc_train>9)[0]), size=len(np.where(auc_train<=9)[0])*4)  # indices of responses to keep (1:3 ratio)
neg_ixs = np.where(auc_train>9) # this is necessary in order to avoid error 
for ix in keep_ix:
    auc_train.iloc[neg_ixs[0][ix], neg_ixs[1][ix]] = None

## k-NN class 
knn = NearestNeighbors(n_neighbors=5, metric='cosine', algorithm='brute')
knn.fit(exp_train.values)
neigh_dist, neigh_ind = knn.kneighbors(exp_test)

## Record predictions on the test set 
auc_pred = pd.DataFrame(index=exp_test.index, columns=auc_test.columns)
for i in range(exp_test.shape[0]):
    auc_pred.iloc[i,:] = auc_train.iloc[neigh_ind[i,:]].mean(axis=0)

evaluation_func.evaluate(auc_test, auc_pred, outdir='drug_response/model_performance/exp_knn_fda_cpd_balanced_tier1and2')
