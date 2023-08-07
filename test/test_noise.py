# %%
import pandas as pd
from anndata import AnnData

from test.SPFinderTester import SPFinderTester

spft = SPFinderTester()
count = pd.read_csv('./data/test2/simulate_exp.csv', sep=',', index_col=0)
position = pd.read_csv('./data/test2/simulate_position.csv', sep=',', index_col=0)
adata = AnnData(X=count)
adata.obs = position.reindex(adata.obs.index)
spft.set_adata(adata)
spft.normalize()
spft.fit_pattern(n_top_genes=200, n_comp=20)
spft.build_distance_array()
spft.cluster_gene(n_clusters=3, mds_components=30)
# %%
label_df = spft.genes_labels
spft.plot_pattern()
