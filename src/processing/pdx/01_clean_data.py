"""
Data downloaded as an Excel file from https://www.nature.com/articles/nm.3954#Sec28
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dn='/projects/b1131/saya/drug_response/pdx'

rna = pd.read_csv(f'{dn}/RNAseq_fpkm.tsv', sep='\t', nrows=200)
cnv = pd.read_csv(f'{dn}/copy_number.tsv', sep='\t', nrows=200)
mutcn2 = pd.read_csv(f'{dn}/pdxe_mut_and_cn2.tsv', sep='\t')
pct_raw = pd.read_csv(f'{dn}/PCT_raw_data.tsv', sep='\t')
pct_met = pd.read_csv(f'{dn}/PCT_curve_metrics.tsv', sep='\t')


outcomes = ['BestResponse', 'Day_BestResponse', 'BestAvgResponse', 'Day_BestAvgResponse', 'TimeToDouble', 'Day_Last']
other = 'ResponseCategory'

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(14,7))
for i in range(6):
    axs[(i//3),(i%3)].hist(pct_met[outcomes[i]])
    axs[(i//3),(i%3)].set_title(outcomes[i])
    axs[0,0].set_ylabel('counts')
    axs[1,0].set_ylabel('counts')

fig.savefig(f'{dn}/metrics_summary.png')
plt.close()
