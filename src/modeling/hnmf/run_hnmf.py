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

# print("Starting run ... {}".format(datetime.datetime.now()))

opts, extraparams = getopt.getopt(sys.argv[1:], 'X1:X2:l1:l2:k:', 
                                  ['X1=', 'X2=', 'X1loss=', 'X2loss=', 'max_k='])

print(sys.argv)
devstr = 'cuda'
config = 'please specify your own configuration string describing, e.g., germline mutations, filtering thresholds in pre-processing steps'
niter = 1000
# lr = 0.01
seed = 1

for o,p in opts:
    if o in ['-X1', '--X1']:
        X1_fn = p
    if o in ['-X2', '--X2']:
        X2_fn = p
    if o in ['-l1', '--X1loss']:
        X1_loss = p
    if o in ['-l2', '--X2loss']:
        X2_loss = p
    if o in ['-k', '--max_k']:
        max_k = int(p)

## load data 
X1 = pd.read_csv(X1_fn, header=0, index_col=0) # variant count (BELOW LINES NEED TO BE EDITED WHEN CHANGING DATA)
X2 = pd.read_csv(X2_fn, header=0, index_col=0) # expression (BELOW LINES NEED TO BE EDITED WHEN CHANGING DATA)

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

print('nc,lr,wd,loss')
for lr in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    for nc in np.arange(50, max_k, step=50):
        for wd in [1e-3, 1e-4, 1e-5, 1e-6]:
            m = hNMF(X1, X2, k=nc, lr=lr, floss1=X1_loss, floss2=X2_loss, verbose=False, n_iter=niter, weight_decay=wd)
            m.fit()
            loss = m.show_report().loss.iloc[len(m.show_report().loss)-1]
            print(f'{nc},{lr},{wd},{loss:.4f}')
