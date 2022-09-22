import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.decomposition import NMF

din = '/projects/b1131/saya/drug_response/processed_data'
dout = '/projects/b1131/saya/drug_response/data_summary'

mor64 = pd.read_csv(f'{din}/morgan_fingerprint_64.csv', index_col=0)
mor128 =  pd.read_csv(f'{din}/morgan_fingerprint_128.csv', index_col=0)

datadict = {
    64 : mor64,
    128 : mor128
}

#### load compound metadata to visualize #### 

cpd = pd.read_csv(f'/projects/b1131/saya/drug_response/ctrp/v20.meta.per_compound.txt', sep='\t')

#### visualize PCA ####
for bits in [64, 128]:
    pca = PCA()
    Xpca = pca.fit_transform(datadict[bits])
    # plot 
    plt.scatter(Xpca[:,0], Xpca[:,1])
    plt.title(f'PCA of chemicals Morgan Fingerprints ({bits} bits)')
    plt.savefig(f'{dout}/pca_compound_morgan_fingerprints{bits}.png')
    plt.close()

