import hotspot
import pandas as pd
from anndata import AnnData
from SPFinder import SPFinder



#%%
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

#%%
hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model='danb',
    latent_obsm_key="X_pca",
    umi_counts_obs_key="total_counts"
)
