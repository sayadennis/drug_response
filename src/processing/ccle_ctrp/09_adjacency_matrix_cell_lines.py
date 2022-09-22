import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

din='/projects/b1131/saya/drug_response/processed_data'
dout='/projects/b1131/saya/drug_response/data_summary'

## Load compound features 
exp = pd.read_csv(f'{din}/expression_mx_aligned.csv', index_col=0)

Am = exp.T.corr()

Am.replace(to_replace=1., value=pd.NA, inplace=True)

#### Try plotting a graph using NetworkX ####
for thres in [0.8, 0.9, 0.95]:
    Adj_m = (Am>thres).astype(int)
    rows, cols = np.where(Adj_m==1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    nx.draw(gr, node_size=20, alpha=0.5) # , node_size=500, with_labels=True
    plt.savefig(f'{dout}/cell_line_expression_networkx_visualization_thres{thres}.png')
    plt.close()

#### Plot histogram of correlation values ####

plt.hist(Am.values.ravel(), bins=30)
plt.title('Distribution of correlations of RNA Expression')
plt.xlabel('Correlation')
plt.ylabel('Counts')
plt.tight_layout()
plt.savefig(f'{dout}/histogram_cell_line_expression_adjacency_distribution.png')
plt.close()

