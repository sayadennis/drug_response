import os
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.metrics import mean_absolute_error #, mean_absolute_percentage_error
from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import coverage_error, label_ranking_average_precision_score

def confirm_numpy(x):
    # Function that will convert inputs into compatible format 
    if type(x) != np.ndarray:
        x = x.to_numpy(dtype=float)
    else:
        x = np.array(x, dtype=float)
    return x

def regres_error(y_test, y_pred):
    y_test, y_pred = confirm_numpy(y_test), confirm_numpy(y_pred)
    y_test = np.array(y_test.ravel(), dtype=float)
    y_pred = np.array(y_pred.ravel(), dtype=float)
    not_nan_bool = (~np.isnan(y_test) & ~np.isnan(y_pred))
    mae = mean_absolute_error(y_test, y_pred)
    # pmae = mean_absolute_percentage_error(y_test, y_pred)
    print('MAE: %s' % mae)
    # print('PMAE: %s' % pmae)
    return mae #, pmae

def performance_table(y_test_bin, y_pred_bin):
    pf = pd.DataFrame(0, index=['acc', 'prec', 'recall', 'f1'], columns=['binary', 'macro-averaged'])
    # fill in binary metrics
    pf.loc['acc','binary'] = accuracy_score(y_test_bin, y_pred_bin)
    pf.loc['prec', 'binary'] = precision_score(y_test_bin, y_pred_bin)
    pf.loc['recall', 'binary'] = recall_score(y_test_bin, y_pred_bin)
    pf.loc['f1', 'binary'] = f1_score(y_test_bin, y_pred_bin)
    # fill in macro averaged metrics 
    # pf.loc['acc','macro-averaged'] = accuracy_score(y_test_bin, y_pred_bin)
    pf.loc['prec', 'macro-averaged'] = precision_score(y_test_bin, y_pred_bin, average='macro')
    pf.loc['recall', 'macro-averaged'] = recall_score(y_test_bin, y_pred_bin, average='macro')
    pf.loc['f1', 'macro-averaged'] = f1_score(y_test_bin, y_pred_bin, average='macro')
    return pf

def evaluate_rank(y_test, y_pred):
    return

def evaluate(y_test, y_pred, outdir):
    # confirm numpy 
    y_test, y_pred = confirm_numpy(y_test), confirm_numpy(y_pred)
    # ravel so they are 1D arrays 
    y_test = np.array(y_test.ravel(), dtype=float)
    y_pred = np.array(y_pred.ravel(), dtype=float)
    # create binary predictions and targets 
    y_pred_bin = np.array(y_pred<9, dtype=int)
    y_test_bin = np.array(y_test<9, dtype=int)
    # find where pred or true is NaN 
    not_nan_bool = (~np.isnan(y_test) & ~np.isnan(y_pred))
    # save metrics
    mae = regres_error(y_test[not_nan_bool], y_pred[not_nan_bool]) # , pmae
    rank_met = evaluate_rank(y_test[not_nan_bool], y_pred[not_nan_bool])
    try:
        os.mkdir(outdir)
    except:
        FileExistsError(f'Directory already exists: {outdir}. \nPlease delete existing directory or specify a new directory.')
    with open(f'{outdir}/regression_metrics.txt', 'w') as f:
        f.write('MAE: %s' % mae)
        # f.write('PMAE: %s' % pmae)
        f.write('Rank based metric: %s\n' % rank_met)
    #
    performance_table(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool]).to_csv(f'{outdir}/classification_performance.csv')
    cfm = confusion_matrix(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool])
    pd.DataFrame(cfm, index=['true neg', 'true pos'], columns=['pred neg', 'pred pos']).to_csv(f'{outdir}/classification_confusion_matrix.csv')
    ncfm = confusion_matrix(y_test_bin[not_nan_bool], y_pred_bin[not_nan_bool], normalize='true')
    pd.DataFrame(ncfm, index=['true neg', 'true pos'], columns=['pred neg', 'pred pos']).to_csv(f'{outdir}/classification_confusion_matrix_norm.csv')
    return
