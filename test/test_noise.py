# %%
import pandas as pd
from test.SPFinderTester import SPFinderTester
from anndata import AnnData

spft = SPFinderTester()

# %%
count = pd.read_csv('./data/test1/count.csv', sep=',', index_col=0)
position = pd.read_csv('./data/test1/position.csv', sep=',', index_col=0)
adata = AnnData(X=count.T)
adata.obs = position.reindex(adata.obs.index)

# %%
spft.set_adata(adata)
spft.normalize()
spft.fit_pattern(n_top_genes=50, n_comp=10)
spft.build_distance_array()
spft.cluster_gene(n_clusters=2, mds_components=30)
#%%

adata.X
