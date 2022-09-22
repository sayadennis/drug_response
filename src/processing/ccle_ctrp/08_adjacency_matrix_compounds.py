import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

din='/projects/b1131/saya/drug_response/processed_data'
dout='/projects/b1131/saya/drug_response/data_summary'

## Load compound features 
reactome = pd.read_csv(f'{din}/cpd_features_reactome_pathways.csv', index_col=0)
morgan = pd.read_csv(f'{din}/morgan_fingerprint_128.csv', index_col=0)

Am = morgan.T.corr()
Rm = reactome.T.corr()

Am.replace(to_replace=1., value=pd.NA, inplace=True)

#### Try plotting a graph using NetworkX ####
for thres in [0.3, 0.4, 0.5]:
    Adj_m = (Am>thres).astype(int)
    rows, cols = np.where(Adj_m==1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    nx.draw(gr, node_size=20, alpha=0.5) # , node_size=500, with_labels=True
    plt.savefig(f'{dout}/networkx_visualization_thres{thres}.png')
    plt.close()

#### Plot histogram of correlation values ####

plt.hist(Am.values.ravel(), bins=30)
plt.title('Distribution of correlations of Morgan Fingerprints')
plt.xlabel('Correlation')
plt.ylabel('Counts')
plt.savefig(f'{dout}/histogram_compound_morgan128_adjacency_distribution.png')
plt.close()

