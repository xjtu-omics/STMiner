import numpy as np
import pandas as pd
from tqdm import tqdm

from Algorithm.distance import get_exp_array


def get_high_var_genes(adata, n_top_genes=200):
    gene_list = list(adata.var.index)
    var_dict = {}
    for i in tqdm(gene_list, desc='Processing...'):
        matrix = get_exp_array(adata, i)
        var = np.var(matrix, dtype=np.float32)
        var_dict[i] = var
    mse_df = pd.DataFrame(var_dict.items(), columns=['gene', 'mse'])
    sorted_df = mse_df.sort_values(by='mse', ascending=False)
    gene_list = list(sorted_df['gene'][:n_top_genes])
    return gene_list
