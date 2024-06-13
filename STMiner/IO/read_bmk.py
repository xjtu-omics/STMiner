import pandas as pd
import scanpy as sc
import os

from STMiner.IO.IOUtil import *


def read_bmk(file_path, bin_size=1):
    if not os.path.isdir(file_path):
        raise Exception("The path is not a directory")
    path = file_path
    adata = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)
    barcode_pos = pd.read_csv(f'{path}/barcodes_pos.tsv.gz', compression='gzip', sep='\t',
                              names=["barcode", "array_col", "array_row"], header=None)
    barcode_pos['in_tissue'] = 1
    barcode_pos = barcode_pos[['barcode', 'in_tissue', 'array_row', 'array_col']]
    adata.obs = barcode_pos
    adata.obs.index = adata.obs['barcode'].to_list()
    adata.obs = adata.obs[['in_tissue', 'array_row', 'array_col']]

    obsm = barcode_pos[['array_col', 'array_row']]
    obsm = obsm.to_numpy()
    adata.obsm["spatial"] = obsm
    if (('x' not in adata.obs.keys()) | ('y' not in adata.obs.keys())) and 'spatial' in adata.obsm.keys():
        position = adata.obsm['spatial']
        x_min = position[:, 0].min()
        y_min = position[:, 1].min()
        adata.obs['x'] = merge_bin_coordinate(position[:, 0], x_min, bin_size=bin_size)
        adata.obs['y'] = merge_bin_coordinate(position[:, 1], y_min, bin_size=bin_size)
        adata.misc = {'up': position[:, 1].min(),
                      'down': position[:, 1].max(),
                      'left': position[:, 0].min(),
                      'right': position[:, 0].max()}
    return adata
