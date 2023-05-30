import anndata
import multiprocessing

import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm
from numba import njit
from scipy import sparse
from sklearn import mixture
from util import array_to_list
from Algorithm.Algorithm import *


def distribution_distance(gmm1, gmm2, method='hellinger'):
    """
    Calculates the distance between gmm1 and gmm2
    :param method:
    :type method:
    :param gmm1: first GMM model
    :type gmm1:
    :param gmm2: second GMM model
    :type gmm2:
    :return: Distance between gmm1 and gmm2
    :rtype:
    """
    gmm1_weights = gmm1.weights_
    gmm1_means = gmm1.means_
    gmm1_covs = gmm1.covariances_
    gmm2_weights = gmm2.weights_
    gmm2_means = gmm2.means_
    gmm2_covs = gmm2.covariances_
    n_components = gmm1_weights.size
    # calculate the distance
    distance_array = np.zeros((n_components, n_components))
    # TODO: other distance metrics
    if method == 'hellinger':
        for i in range(n_components):
            for j in range(n_components):
                distance_array[i, j] = get_hellinger_distance(gmm1_covs[i], gmm1_means[i], gmm2_covs[j], gmm2_means[j])
        distance = linear_sum(distance_array)
    return distance


def get_bh_distance(gmm1_covs, gmm1_means, gmm2_covs, gmm2_means):
    mean_cov = (gmm1_covs + gmm2_covs) / 2
    mean_cov_det = np.linalg.det(mean_cov)
    mean_cov_inv = np.linalg.inv(mean_cov)
    means_diff = gmm1_means - gmm2_means
    first_term = (means_diff.T @ mean_cov_inv @ means_diff) / 8
    second_term = np.log(mean_cov_det / (np.sqrt(np.linalg.det(gmm1_covs) * np.linalg.det(gmm2_covs)))) / 2
    result = first_term + second_term
    return result


@njit
def get_hellinger_distance(gmm1_covs, gmm1_means, gmm2_covs, gmm2_means):
    """
    Calculates the distance between two GMM models by hellinger distance.

    **Hellinger distance** (closely related to, although different from, the Bhattacharyya distance) is used to quantify
    the similarity between two probability distributions.

    Ref:
     - https://en.wikipedia.org/wiki/Hellinger_distance
    :param gmm1_covs: first gmm covariances
    :type gmm1_covs: np.Array_
    :param gmm1_means: first gmm means
    :type gmm1_means: Array
    :param gmm2_covs: second gmm covariances
    :type gmm2_covs: Array
    :param gmm2_means: second gmm means
    :type gmm2_means: Array
    :return: distance between two GMM models
    :rtype: np.float
    """
    mean_cov = (gmm1_covs + gmm2_covs) / 2
    mean_cov_det = np.linalg.det(mean_cov)
    mean_cov_inv = np.linalg.inv(mean_cov)
    gmm1_cov_det = np.linalg.det(gmm1_covs)
    gmm2_cov_det = np.linalg.det(gmm2_covs)
    means_diff = gmm1_means - gmm2_means
    first_term = np.exp(-(means_diff.T @ mean_cov_inv @ means_diff) / 8)
    second_term = ((np.power(gmm1_cov_det, 0.25)) * np.power(gmm2_cov_det, 0.25)) / np.sqrt(mean_cov_det)
    hellinger_distance = np.sqrt(np.abs(1 - first_term * second_term))
    return hellinger_distance


def fit_gmm(adata: anndata,
            gene_name: str,
            n_comp: int = 10,
            max_iter: int = 200):
    """
    Representation of a Gaussian mixture model probability distribution.
    Estimate the parameters of a Gaussian mixture distribution.

    Estimate model parameters with the EM algorithm.
    :param adata: Anndata of spatial data
    :type adata: Anndata
    :param gene_name: The gene name to fit
    :type gene_name: str
    :param n_comp: The number of mixture components.
    :type n_comp: int
    :param max_iter: The number of EM iterations to perform.
    :type max_iter: int
    :return: The fitted mixture.
    :rtype: GaussianMixture
    """
    data = np.array(adata[:, adata.var_names == gene_name].X.todense())
    sparse_matrix = sparse.coo_matrix((data[:, 0], (np.array(adata.obs['x']), np.array(adata.obs['y']))))
    dense_array = np.array(sparse_matrix.todense(), dtype=np.int32)
    result = array_to_list(dense_array)
    # Number of unique center must be larger than the number of components.
    if len(set(map(tuple, result))) > n_comp:
        gmm = mixture.GaussianMixture(n_components=n_comp, max_iter=max_iter)
        gmm.fit(result)
        return gmm
    else:
        return None


def fit_gmms(adata,
             gene_name_list: list,
             n_comp: int = 5,
             max_iter: int = 1000):
    """
    Same as fit_gmm, accepts a list of gene name.
    Representation of a Gaussian mixture model probability distribution.
    Estimate the parameters of a Gaussian mixture distribution.
    :param adata: Anndata of spatial data
    :type adata: Anndata
    :param gene_name_list: The python list, each element is a gene name in adata.
    :type gene_name_list: list
    :param n_comp: The number of mixture components.
    :type n_comp: int
    :param max_iter: The number of EM iterations to perform.
    :type max_iter: int
    :return: A Python dict of given genes list, key is gene name, value is GMM object.
    :return:
    :rtype: dict
    """
    gmm_dict = {}
    dropped_genes = []
    for gene_id in tqdm(gene_name_list, desc='Processing ...'):
        try:
            fit_result = fit_gmm(adata, gene_id, n_comp=n_comp, max_iter=max_iter)
            if fit_result is not None:
                gmm_dict[gene_id] = fit_result
            else:
                dropped_genes.append(gene_id)
        except Exception as e:
            error_msg = str(e)
            raise ValueError("Gene id: " + gene_id + "\nError: " + error_msg)
    print('Dropped genes:')
    print(dropped_genes)
    return gmm_dict


def fit_gmms_multiprocessing(adata: anndata,
                             gene_name_list: list,
                             n_comp: int = 5,
                             max_iter: int = 1000,
                             thread: int = 4) -> dict:
    """
    Same as fit_gmm, use multiple threads.
    Representation of a Gaussian mixture model probability distribution.
    Estimate the parameters of a Gaussian mixture distribution.

    Estimate model parameters with the EM algorithm.
    :param adata: Anndata of spatial data
    :type adata: Anndata
    :param gene_name_list: The genes list to fit
    :type gene_name_list: list
    :param n_comp: The number of mixture components.
    :type n_comp: int
    :param max_iter: The number of EM iterations to perform.
    :type max_iter: int
    :param thread: The number of threads to use, default:4
    :type thread: int
    :return: A Python dict of given genes list, key is gene name, value is GMM object.
    :rtype: dict
    """
    manager = multiprocessing.Manager()
    shared_dict = manager.dict()
    pool = multiprocessing.Pool(processes=thread)
    for i in gene_name_list:
        pool.apply_async(_fit_worker, args=(shared_dict, adata, i, n_comp, max_iter))
    pool.close()
    pool.join()
    normal_dict = dict(shared_dict)
    return normal_dict


def _fit_worker(shared_dict, adata, gene_name, n_comp, max_iter):
    shared_dict[gene_name] = fit_gmm(adata, gene_name, n_comp, max_iter)


def view_gmm(gmm, scope, bin_count=None):
    """
    View the fitted GMM model.
    :param bin_count:
    :type bin_count:
    :param scope:
    :type scope:
    :param gmm: fitted GMM model by sklearn.mixture.GaussianMixture
    :type gmm: sklearn.mixture.GaussianMixture
    """
    start = scope[0]
    end = scope[1]
    if bin_count is None:
        num = end - start
    else:
        num = bin_count
    x = np.linspace(start, end, num)
    y = np.linspace(start, end, num)
    x_range, y_range = np.meshgrid(x, y)
    x_y = np.column_stack([x_range.flat, y_range.flat])
    # calculate the density
    density = gmm.score_samples(x_y)
    density = density.reshape(x_range.shape)
    sns.heatmap(np.exp(density))
    plt.show()


def view_pattern(adata, gene_list, size=6):
    """
    Plot expression pattern of given genes.
    :param adata:
    :type adata:
    :param gene_list:
    :type gene_list:
    :param size:
    :type size:
    """
    tmp = adata[:, gene_list]
    sc.pl.scatter(tmp,
                  x='x',
                  y='y',
                  color='log1p_n_genes_by_counts',
                  size=size)
