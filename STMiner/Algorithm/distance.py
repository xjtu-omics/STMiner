import numpy as np
import ot
import pandas as pd
from numba import njit
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cosine
from tqdm import tqdm

from STMiner.Algorithm.AlgUtils import get_exp_array
from STMiner.Algorithm.algorithm import linear_sum
from STMiner.Utils.utils import issparse


def distribution_distance(first_gmm, second_gmm):
    """
    Calculates the distance between gmm1 and gmm2
    :param first_gmm: The first GMM model
    :type first_gmm:
    :param second_gmm: The second GMM model
    :type second_gmm:
    :return: Distance between gmm1 and gmm2
    :rtype: np.float64
    """
    gmm1_weights = first_gmm.weights_
    gmm1_means = first_gmm.means_
    gmm1_covs = first_gmm.covariances_
    gmm2_weights = second_gmm.weights_
    gmm2_means = second_gmm.means_
    gmm2_covs = second_gmm.covariances_
    n_components = gmm1_weights.size
    distance_array = np.zeros((n_components, n_components))
    for i in range(n_components):
        for j in range(n_components):
            hd = get_hellinger_distance(gmm1_covs[i], gmm1_means[i], gmm2_covs[j], gmm2_means[j])
            distance_array[i, j] = (gmm1_weights[i] + gmm2_weights[j]) * hd
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
    :rtype: np.float64
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


def build_gmm_distance_array(gmm_dict):
    """
    Generate a distance matrix by the given gmm dictionary.
    :param gmm_dict: gmm dictionary, key is gene name, value is GMM model.
    :type gmm_dict: dict
    :return: distance array
    :rtype: np.Array
    """
    gene_list = list(gmm_dict.keys())
    gene_counts = len(gene_list)
    distance_array = pd.DataFrame(0, index=gene_list, columns=gene_list, dtype=np.float64)
    for i in tqdm(range(gene_counts), desc='Building distance array...'):
        for j in range(i + 1, gene_counts):
            distance = distribution_distance(gmm_dict[gene_list[i]], gmm_dict[gene_list[j]])
            distance_array.loc[gene_list[i], gene_list[j]] = distance
            distance_array.loc[gene_list[j], gene_list[i]] = distance
    return distance_array


def compare_gmm_distance(gmm, gmm_dict):
    gene_list = list(gmm_dict.keys())
    gene_counts = len(gene_list)
    distance_dict = {}
    for gene in tqdm(range(gene_counts), desc='Comparing distance...'):
        distance = distribution_distance(gmm, gmm_dict[gene_list[gene]])
        distance_dict[gene_list[gene]] = distance
    df = pd.DataFrame.from_dict(distance_dict, orient='index')
    sorted_df = df.sort_values(0)
    return sorted_df


def build_cosine_similarity_array(adata, gene_list):
    gene_arr_dict = {}
    gene_list = list(gene_list)
    for gene in gene_list:
        gene_arr_dict[gene] = get_exp_array(adata, gene)
    gene_counts = len(gene_list)
    distance_array = pd.DataFrame(0, index=gene_list, columns=gene_list, dtype=np.float64)
    for i in tqdm(range(gene_counts), desc='Building distance array...'):
        for j in range(gene_counts):
            if i != j:
                distance = cosine(gene_arr_dict[gene_list[i]].flatten(),
                                  gene_arr_dict[gene_list[j]].flatten())
                distance_array.loc[gene_list[i], gene_list[j]] = distance
    return distance_array


def build_mse_distance_array(adata, gene_list):
    gene_arr_dict = {}
    gene_list = list(gene_list)
    for gene in gene_list:
        gene_arr_dict[gene] = get_exp_array(adata, gene)
    gene_counts = len(gene_list)
    distance_array = pd.DataFrame(0, index=gene_list, columns=gene_list, dtype=np.float64)
    for i in tqdm(range(gene_counts), desc='Building distance array...'):
        for j in range(gene_counts):
            if i != j:
                arr_i = gene_arr_dict[gene_list[i]]
                arr_j = gene_arr_dict[gene_list[j]]
                distance = mse(arr_i, arr_j)
                distance_array.loc[gene_list[i], gene_list[j]] = distance
    return distance_array


@njit()
def mse(arr_i, arr_j):
    return np.mean(np.square(arr_i - arr_j))


def build_mix_distance_array(adata, gmm_dict):
    """

    :param adata:
    :type adata:
    :param gmm_dict:
    :type gmm_dict:
    :return:
    :rtype:
    """
    gene_list = list(gmm_dict.keys())
    gene_counts = len(gene_list)
    distance_array = pd.DataFrame(0, index=gene_list, columns=gene_list, dtype=np.float64)
    for i in tqdm(range(gene_counts), desc='Building distance array...'):
        for j in range(gene_counts):
            if i != j:
                diff = adata[:, [gene_list[i]]].X - adata[:, [gene_list[j]]].X
                if issparse(diff):
                    mse_distance = np.mean(np.square(diff.todense()))
                else:
                    mse_distance = np.mean(np.square(diff))
                gmm_distance = distribution_distance(gmm_dict[gene_list[i]],
                                                     gmm_dict[gene_list[j]])
                distance_array.loc[gene_list[i], gene_list[j]] = mse_distance + gmm_distance
    return distance_array


def build_ot_distance_array(csr_dict, gene_list=None):
    genes_to_build = gene_list if gene_list is not None else list(csr_dict.keys())
    gene_counts = len(genes_to_build)
    distance_array = pd.DataFrame(0, index=genes_to_build, columns=genes_to_build, dtype=np.float64)
    for i in tqdm(range(gene_counts), desc='Building distance array...'):
        for j in range(i + 1, gene_counts):
            distance = calculate_ot_distance(csr_dict[genes_to_build[i]], csr_dict[genes_to_build[j]])
            distance_array.loc[genes_to_build[i], genes_to_build[j]] = distance
            distance_array.loc[genes_to_build[j], genes_to_build[i]] = distance
    return distance_array


def calculate_ot_distance(source, target) -> float:
    xsd = source.data
    xtd = target.data
    xs = np.array(source.nonzero()).T
    xt = np.array(target.nonzero()).T
    # uniform distribution on samples
    a = source.data[source.data > 0] / xsd.sum()
    b = target.data[target.data > 0] / xtd.sum()
    # loss matrix
    loss_matrix = ot.dist(xs, xt)
    distance = ot.emd2(a, b, loss_matrix, numItermax=200000)
    return distance


def domain_alignment(source_domain_dict, target_domain_dict):
    """
    Returns the distance between two domains, index is source domain and column is target domain.
    :param source_domain_dict: dict
    :param target_domain_dict: dict
    :return: pd.DataFrame
    """
    source_keys = list(source_domain_dict.keys())
    target_keys = list(target_domain_dict.keys())
    distance_array = pd.DataFrame(0, index=source_keys, columns=target_keys, dtype=np.float64)
    for i in range(len(source_domain_dict)):
        for j in range(len(target_domain_dict)):
            csr_source = csr_matrix(source_domain_dict[i])
            csr_target = csr_matrix(target_domain_dict[j])
            ot_dist = calculate_ot_distance(csr_source, csr_target)
            distance_array.loc[source_keys[i], target_keys[j]] = ot_dist
            distance_array.loc[source_keys[j], target_keys[i]] = ot_dist
    return distance_array
