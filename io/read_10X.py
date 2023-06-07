import scanpy as sc

from IO.IOUtil import *


def read_10X_h5ad(file, amplification=1, bin_size=1):
    """
    Read h5ad file to anndata and add more information for SEP pipline.
    :param bin_size:
    :type bin_size:
    :param file: path to h5ad file
    :type file: str
    :param amplification:
    :type amplification: int
    :return:
    :rtype: Anndata
    """
    amplification = np.int(amplification)
    adata = sc.read_h5ad(file)
    if (('x' not in adata.obs.keys()) | ('y' not in adata.obs.keys())) and 'spatial' in adata.obsm.keys():
        position = adata.obsm['spatial']
        x_min = position[:, 0].min() * amplification
        y_min = position[:, 1].min() * amplification
        adata.obs['x'] = merge_bin_coordinate(position[:, 0] * amplification, x_min, bin_size=bin_size)
        adata.obs['y'] = merge_bin_coordinate(position[:, 1] * amplification, y_min, bin_size=bin_size)
    return adata
