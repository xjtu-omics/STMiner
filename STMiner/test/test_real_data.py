import pandas as pd
from anndata import AnnData

from STMiner.SPFinder import SPFinder

count = pd.read_csv('E://SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', sep=',', index_col=0)
position = pd.read_csv('E://SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', sep=',', index_col=0)

adata = AnnData(X=count)
adata.obs = position.reindex(adata.obs.index)

spf = SPFinder()
spf.set_adata(adata)
spf.merge_bin(bin_width=50, amplification=1000)
spf.normalize()
spf.fit_pattern(n_top_genes=100, n_comp=20)
spf.cluster_gene(n_clusters=3, mds_components=30)
spf.genes_patterns
