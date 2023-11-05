import numpy as np
import pandas as pd
from tqdm import tqdm

from STMiner.Algorithm.AlgUtils import get_exp_array


def get_high_var_genes(adata, n_top_genes=200):
    """
    The variance is the average of the squared deviations from the mean,
    i.e., var = mean(x), where x = abs(a - a.mean())**2.
    :param adata: Anndata of spatial data
    :type adata: Anndata
    :param n_top_genes: number of top-genes selected by variances.
    :type n_top_genes: int
    :return: list with genes name
    :rtype: list
    """
    gene_list = list(adata.var.index)
    var_dict = {}
    for gene in tqdm(gene_list, desc='Processing...'):
        matrix = get_exp_array(adata, gene)
        var = np.var(matrix, dtype=np.float32)
        var_dict[gene] = var
    mse_df = pd.DataFrame(var_dict.items(), columns=['gene', 'mse'])
    sorted_df = mse_df.sort_values(by='mse', ascending=False)
    gene_list = list(sorted_df['gene'].iloc[:n_top_genes])
    return gene_list
