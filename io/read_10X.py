import numpy as np
import scanpy as sc


def read_10X_h5ad(file, amplification=1):
    """
    Read h5ad file to anndata and add more information for SEP pipline.
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
        adata.obs['x'] = _correct_coordinate(position[:, 0] * amplification)
        adata.obs['y'] = _correct_coordinate(position[:, 1] * amplification)
    return adata


def _correct_coordinate(array):
    return (array - array.min()).astype(np.int32)
