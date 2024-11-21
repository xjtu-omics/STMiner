import numpy as np
from scipy import sparse


def get_exp_array(adata, gene_name, remove_low_exp_spots=False):
    sparse_matrix = _preprocess(adata, gene_name)
    return _to_dense(remove_low_exp_spots, sparse_matrix)


def _to_dense(remove_low_exp_spots, sparse_matrix):
    if sparse_matrix is not None:
        dense_array = sparse_matrix.todense()
        if remove_low_exp_spots:
            dense_array = np.maximum(dense_array - np.mean(dense_array[dense_array != 0]), 0)
        dense_array = np.array(np.round(dense_array), dtype=np.int32)
        return dense_array
    else:
        print("matrix is None")


def _preprocess(adata, gene_name):
    if gene_name not in adata.var_names:
        print(f"Warning: {gene_name} is not in adata gene list.")
    else:
        exp_array = adata[:, adata.var_names == gene_name].X
        if sparse.issparse(exp_array):
            data = np.array(exp_array.todense())
        else:
            data = np.array(exp_array)
        sparse_matrix = sparse.coo_matrix((data[:, 0], (np.array(adata.obs['x']), np.array(adata.obs['y']))))
        return sparse_matrix
