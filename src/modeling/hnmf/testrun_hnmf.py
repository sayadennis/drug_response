import os
import sys
import numpy as np
import pandas as pd
import pickle
import getopt
import datetime
import matplotlib.pyplot as plt
# from math import floor, ceil
# from scipy.sparse import coo_matrix

# from sklearn.linear_model import LogisticRegression
# from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, f1_score, balanced_accuracy_score
sys.path.append('drug_response/src/modeling/hnmf')
from hNMF import hNMF
import torch
import torch.nn as nn

import importlib
import warnings
warnings.filterwarnings(action='ignore')

print("Starting run ... {}".format(datetime.datetime.now()))

devstr = 'cpu'
config = 'please specify your own configuration string describing, e.g., germline mutations, filtering thresholds in pre-processing steps'
niter = 25
# lr = 0.01
seed = 1

## load data 
X1 = pd.read_csv('/projects/b1131/saya/drug_response/depmap/01_cleaned/cleaned_vc_mx.csv', header=0, index_col=0)
X2 = pd.read_csv('/projects/b1131/saya/drug_response/depmap/01_cleaned/cleaned_exp_mx.csv', header=0, index_col=0)

# reduce features 
X1 = X1.iloc[:,(X1==0).sum().values/X1.shape[0]<0.99] # select genes where at least 1% of cells have some mutation 
X2 = X2.iloc[:,np.argsort(X2.std()).values<5000] # select top 5000 most variable genes 

# align X1 and X2 
ol=[]
for item in X1.index:
    if item in X2.index:
        ol.append(item)

X1 = X1.loc[ol,:]
X2 = X2.loc[ol,:]

X1 = X1.astype(float).to_numpy()
X2 = X2.astype(float).to_numpy()
device = torch.device(devstr)

for lr in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    m = hNMF(X1, X2, floss1='kl', floss2='kl', verbose=True, n_iter=1000, weight_decay=1e-6)
    m.fit()

    plt.plot(m.show_report().loss)
    plt.xlabel('epochs')
    plt.ylabel('loss')
    plt.title('Loss curve of hNMF')
    plt.tight_layout()
    plt.savefig(f'20220411_test_hnmf_loss_curve_lr{lr}.png')
    plt.close()

# ncs = range(100,501,100) # range(50,501,50)

# print('nc,wcls,C,lr,tr bal acc,val bal acc,te acc,te bal acc,te precis,te recall,te f1,w2,b2,celoss') # ,best iter,mse,mse tr,mse val,mse te
# for lr in [0.001, 0.01]: # 0.00001, 0.0001, 
#     for nc in ncs:
#         for C in [0.01, 0.1, 1, 10, 100]: # , 0.1, 1, 10, 100, 1000, 0.001, 
#             for wcls in [0.01, 0.1, 1]: # , 2, 10, 50, 0.5, 
#                 #### to implement k-fold cross-validation, perform the splits here #### 
#                 #### run a for-loop here and record performance for every fold to average it later #### 
#                 fn = '%s/scanmap_k%d_wcls%s_C%s_lr%s.p' % (outdir, nc, wcls, C, lr) # scanmap%d/s%d/ # niter, seed, 
#                 m = ScanMap(
#                     np.vstack((X_train, X_val, X_test)), 
#                     cf = pts_tr, cfval = pts_val, y = y_train, yval = y_val, 
#                     k=nc, n_iter=niter, weight_decay=0, lr=lr, wcls=wcls, C=C, seed=2*722019+seed, fn=fn, device=device
#                 ) # 2*722019+seed is just to have a large odd number for seeding that is recommended for generating random numbers, fixed for reproducibility
#                 [X_tr_nmf, X_val_nmf, X_te_nmf, H] = m.fit_transform()

#                 chkpt = torch.load(fn)
#                 m.load_state_dict(chkpt['state_dict'])
#                 # best_iter = chkpt['epoch']
#                 accval = chkpt['best_val_acc']
#                 m.eval()

#                 y_tr_pred = m.predict(X_tr_nmf, pts_tr)
#                 y_te_pred = m.predict(X_te_nmf, pts_te)
#                 acctr = balanced_accuracy_score(y_train, y_tr_pred)
#                 accte = accuracy_score(y_test, y_te_pred)
#                 balaccte = balanced_accuracy_score(y_test, y_te_pred)
#                 preciste = precision_score(y_test, y_te_pred)
#                 recallte = recall_score(y_test, y_te_pred)
#                 f1te = f1_score(y_test, y_te_pred)
                
#                 w2 = np.square(m.state_dict()['fc.weight'].cpu().numpy()).sum(axis=None)
#                 b2 = np.square(m.state_dict()['fc.bias'].cpu().numpy()).sum(axis=None)
#                 w = 1 / pd.Series(y_train).value_counts(normalize=True).sort_index().to_numpy()
#                 vce = chkpt['celoss']

#                 # err = np.vstack((X_train, X_val, X_test)) - np.vstack((X_tr_nmf, X_val_nmf, X_te_nmf)) @ H
#                 # err_tr = X_train - X_tr_nmf @ H
#                 # err_val = X_val - X_val_nmf @ H
#                 # err_te = X_test - X_te_nmf @ H
#                 # mse = np.square(err).mean(axis=None)
#                 # mse_tr = np.square(err_tr).mean(axis=None)
#                 # mse_val = np.square(err_val).mean(axis=None)
#                 # mse_te = np.square(err_te).mean(axis=None)
                
#                 print('%d,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f' % 
#                     (nc, wcls, C, lr, acctr, accval, accte, balaccte, preciste, recallte, f1te, w2, b2, vce)) # best_iter, mse, mse_tr, mse_val, mse_te

# # print("Ending run ... {}".format(datetime.datetime.now()))
