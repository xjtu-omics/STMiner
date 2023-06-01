import numpy as np
import scanpy as sc


def read_10X_h5ad(file, amplification=1):
    """

    :param file:
    :type file:
    :param amplification:
    :type amplification:
    :return:
    :rtype:
    """
    amplification = np.int(amplification)
    adata = sc.read_h5ad(file)
    if ('x' not in adata.obs.keys() | 'y' not in adata.obs.keys()) and 'spatial' in adata.obsm.keys():
        position = adata.obsm['spatial']
        adata.obs['x'] = _correct_coordinate(position[:, 0] * amplification)
        adata.obs['y'] = _correct_coordinate(position[:, 1] * amplification)
    return adata


def _correct_coordinate(array):
    return (array - array.min()).astype(np.int32)
