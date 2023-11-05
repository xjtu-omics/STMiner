import numpy as np
from scipy import sparse


def get_exp_array(adata, gene_name, remove_low_exp_spots=False):
    exp_array = adata[:, adata.var_names == gene_name].X
    if sparse.issparse(exp_array):
        data = np.array(exp_array.todense())
    else:
        data = np.array(exp_array)
    sparse_matrix = sparse.coo_matrix((data[:, 0], (np.array(adata.obs['x']), np.array(adata.obs['y']))))
    dense_array = sparse_matrix.todense()
    if remove_low_exp_spots:
        dense_array = np.maximum(dense_array - np.mean(dense_array[dense_array != 0]), 0)
    dense_array = np.array(np.round(dense_array), dtype=np.int32)
    return dense_array
