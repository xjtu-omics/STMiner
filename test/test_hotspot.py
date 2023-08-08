import hotspot
import pandas as pd
import scanpy as sc
from anndata import AnnData

from IO.IOUtil import merge_bin_coordinate
from SPFinder import SPFinder

# %%
raw_df = pd.read_csv("E:/data/brain10x.csv",
                     sep=',',
                     index_col=0)
APP23_1 = raw_df.filter(regex=r'^APP23_I_')
APP23_1 = APP23_1.rename(columns=lambda x: x.replace('APP23_I_', ''))
adata = AnnData(X=APP23_1.T)
position_df = pd.read_csv('E:/data/RegionC/RegionC_10x_immunoloc.csv', sep=',', index_col=0)
# replace '-' to '.', example: 'AAACAAGTATCTCCCA-1' to 'AAACAAGTATCTCCCA.1'
position_df.index = position_df.index.str.replace('-', '.')
position_df.columns = position_df.columns.str.replace('pxl_row_in_immuno', 'x')
position_df.columns = position_df.columns.str.replace('pxl_col_in_immuno', 'y')
adata.obs = position_df.reindex(adata.obs.index)

# %%
adata.obs['x'] = merge_bin_coordinate(adata.obs['x'], adata.obs['x'].min(), bin_size=80)
adata.obs['y'] = merge_bin_coordinate(adata.obs['y'], adata.obs['y'].min(), bin_size=80)
sc.pp.filter_genes(adata, min_cells=200)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.calculate_qc_metrics(adata, inplace=True)

# %%
sc.tl.pca(adata, svd_solver='arpack')

# %%
sp = SPFinder(adata)

# %%
hs = hotspot.Hotspot(
    adata,
    model='danb',
    latent_obsm_key="X_pca",
    umi_counts_obs_key="total_counts"
)
# %%
hs.create_knn_graph(weighted_graph=False, n_neighbors=30)

# %%
hs_results = hs.compute_autocorrelations()
# %%
hs_genes = hs_results.loc[hs_results.FDR < 0.05].index  # Select genes

local_correlations = hs.compute_local_correlations(hs_genes, jobs=1)  # jobs for parallelization

# %%
modules = hs.create_modules(
    min_gene_threshold=30, core_only=True, fdr_threshold=0.05
)
module_scores = hs.calculate_module_scores()

# %%
module_scores
# %%
hs.plot_local_correlations()
