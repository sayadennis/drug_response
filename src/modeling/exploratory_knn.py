import os
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import mean_absolute_error, accuracy_score, balanced_accuracy_score, precision_score, recall_score, f1_score

## Load data 
dn = '/projects/b1131/saya/drug_response/processed_data'

ccl = pd.read_csv(os.path.join(dn, 'vc_mx_aligned.csv'), index_col=0)
auc = pd.read_csv(os.path.join(dn, 'response_continuous_aligned.csv'), index_col=0)
cpd_meta = pd.read_csv('/projects/b1131/saya/drug_response/ctrp/v20.meta.per_compound.txt', sep='\t')

## Train-test split 
ccl_train, ccl_test, auc_train, auc_test = train_test_split(ccl, auc, test_size=0.2, random_state=42)

## k-NN class 
knn = NearestNeighbors(n_neighbors=5, metric='cosine', algorithm='brute')
knn.fit(ccl_train.values)
neigh_dist, neigh_ind = knn.kneighbors(ccl_test)

## Record predictions on the test set 
auc_pred = pd.DataFrame(index=ccl_test.index, columns=auc_test.columns)
for i in range(ccl_test.shape[0]):
    auc_pred.iloc[i,:] = auc_train.iloc[neigh_ind[i,:]].mean(axis=0)

## Evaluate predictions based on MAE 
# convert to 1D array and create boolean array to remove NaNs 
y_test = np.array(auc_test.values.ravel(), dtype=float)
y_pred = np.array(auc_pred.values.ravel(), dtype=float)
not_nan_bool = (~np.isnan(y_test) & ~np.isnan(y_pred))

mean_absolute_error(y_test[not_nan_bool], y_pred[not_nan_bool])
# 1.27

## Evaluate predictions as a binary classification task (sensitive vs not)
# first load the binary target and create the same split to get the test set target
auc_bin = pd.read_csv(os.path.join(dn, 'response_binary_aligned.csv'), index_col=0)
_, _, auc_bin_train, auc_bin_test = train_test_split(ccl, auc_bin, test_size=0.2, random_state=42)

# convert to 1D array and create boolean array to remove NaNs 
y_pred_bin = np.array(y_pred>9, dtype=int)
y_test_bin = np.array(auc_bin_test.values.ravel(), dtype=float)
not_nan_bool = (~np.isnan(y_test_bin) & ~np.isnan(y_pred_bin))

print('accuracy\tbalanced accuracy\tprecision\trecall\tf1')
print('%s\t%s\t%s\t%s\t%s\n' % (
    accuracy_score(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool]),
    balanced_accuracy_score(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool]),
    precision_score(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool]),
    recall_score(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool]),
    f1_score(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool])
))
# accuracy	balanced accuracy	precision	recall	f1
# 0.91      0.70            	0.94	    0.96	0.95

# evaluate predictions based on ranks 

##################################################
#### Do the same with only FDA approved drugs ####
##################################################

fda_ids = np.array(cpd_meta['master_cpd_id'].iloc[cpd_meta['cpd_status'].values=='FDA'], dtype=str) # dtype has to be str

## Train-test split 
ccl_train, ccl_test, auc_train, auc_test = train_test_split(ccl, auc.loc[:,fda_ids], test_size=0.2, random_state=42)

## k-NN class 
knn = NearestNeighbors(n_neighbors=5, metric='cosine', algorithm='brute')
knn.fit(ccl_train.values)
neigh_dist, neigh_ind = knn.kneighbors(ccl_test)

## Record predictions on the test set 
auc_pred = pd.DataFrame(index=ccl_test.index, columns=auc_test.columns)
for i in range(ccl_test.shape[0]):
    auc_pred.iloc[i,:] = auc_train.iloc[neigh_ind[i,:]].mean(axis=0)

## Evaluate predictions based on MAE 
# convert to 1D array and create boolean array to remove NaNs 
y_test = np.array(auc_test.values.ravel(), dtype=float)
y_pred = np.array(auc_pred.values.ravel(), dtype=float)
not_nan_bool = (~np.isnan(y_test) & ~np.isnan(y_pred))

mean_absolute_error(y_test[not_nan_bool], y_pred[not_nan_bool])
# 1.36

## Evaluate predictions as a binary classification task (sensitive vs not)
# first load the binary target and create the same split to get the test set target
auc_bin = pd.read_csv(os.path.join(dn, 'response_binary_aligned.csv'), index_col=0)
_, _, auc_bin_train, auc_bin_test = train_test_split(ccl, auc_bin.loc[:,fda_ids], test_size=0.2, random_state=42)

# convert to 1D array and create boolean array to remove NaNs 
y_pred_bin = np.array(y_pred>9, dtype=int)
y_test_bin = np.array(auc_bin_test.values.ravel(), dtype=float)
not_nan_bool = (~np.isnan(y_test_bin) & ~np.isnan(y_pred_bin))

print('accuracy\tbalanced accuracy\tprecision\trecall\tf1')
print('%s\t%s\t%s\t%s\t%s\n' % (
    accuracy_score(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool]),
    balanced_accuracy_score(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool]),
    precision_score(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool]),
    recall_score(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool]),
    f1_score(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool])
))
# accuracy	balanced accuracy	precision	recall	f1
# 0.91	    0.71            	0.94    	0.96	0.95
