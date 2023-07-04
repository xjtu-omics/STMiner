import networkx as nx
from scipy.spatial.distance import cdist
from scipy.spatial.distance import cosine

from Algorithm.distribution import *
from Utils.utils import issparse


class Graph:
    def __init__(self):
        pass


def build_gmm_distance_array(gmm_dict, method='weight_match'):
    """
    Generate a distance matrix by the given gmm dictionary.
    :param method: 'weight_match' or 'optimized_match' or 'emd, default: 'weight_match'
    :type method:
    :param gmm_dict: gmm dictionary, key is gene name, value is GMM model.
    :type gmm_dict: dict
    :return: distance array
    :rtype: np.Array
    """
    gene_list = list(gmm_dict.keys())
    gene_counts = len(gene_list)
    distance_array = pd.DataFrame(0, index=gene_list, columns=gene_list, dtype=np.float64)
    for i in tqdm(range(gene_counts), desc='Building distance array...'):
        for j in range(gene_counts):
            if i != j:
                distance = distribution_distance(gmm_dict[gene_list[i]], gmm_dict[gene_list[j]], method=method)
                distance_array.loc[gene_list[i], gene_list[j]] = distance
    return distance_array


def compare_gmm_distance(gmm, gmm_dict, method='weight_match'):
    gene_list = list(gmm_dict.keys())
    gene_counts = len(gene_list)
    distance_dict = {}
    for gene in tqdm(range(gene_counts), desc='Building distance array...'):
        distance = distribution_distance(gmm, gmm_dict[gene_list[gene]], method=method)
        distance_dict[gene_list[gene]] = distance
    df = pd.DataFrame.from_dict(distance_dict, orient='index')
    sorted_df = df.sort_values(0)
    return sorted_df


def build_cosine_similarity_array(adata, gene_list):
    gene_list = list(gene_list)
    gene_counts = len(gene_list)
    distance_array = pd.DataFrame(0, index=gene_list, columns=gene_list, dtype=np.float64)
    for i in tqdm(range(gene_counts), desc='Building distance array...'):
        for j in range(gene_counts):
            if i != j:
                distance = cosine(get_exp_array(adata, gene_list[i]).flatten(),
                                  get_exp_array(adata, gene_list[j]).flatten())
                distance_array.loc[gene_list[i], gene_list[j]] = distance
    return distance_array


def build_mse_distance_array(adata, gene_list):
    gene_list = list(gene_list)
    gene_counts = len(gene_list)
    distance_array = pd.DataFrame(0, index=gene_list, columns=gene_list, dtype=np.float64)
    for i in tqdm(range(gene_counts), desc='Building distance array...'):
        for j in range(gene_counts):
            if i != j:
                distance = np.mean(np.square(adata[:, [gene_list[i]]].X - adata[:, [gene_list[j]]].X))
                distance_array.loc[gene_list[i], gene_list[j]] = distance
    return distance_array


def build_mix_distance_array(adata, gmm_dict, method='optimized_match'):
    """

    :param adata:
    :type adata:
    :param gmm_dict:
    :type gmm_dict:
    :param method: optimized_match or weight_match
    :type method:
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
                                                     gmm_dict[gene_list[j]],
                                                     method=method)
                # TODO:
                distance_array.loc[gene_list[i], gene_list[j]] = mse_distance + gmm_distance
    return distance_array


def build_ot_distance_array(adata, gene_list):
    gene_list = list(gene_list)
    gene_counts = len(gene_list)
    distance_array = pd.DataFrame(0, index=gene_list, columns=gene_list, dtype=np.float64)
    for i in tqdm(range(gene_counts), desc='Building distance array...'):
        for j in range(gene_counts):
            if i != j:
                d = cdist(get_exp_array(adata, gene_list[i]), get_exp_array(adata, gene_list[j]))
                assignment = linear_sum_assignment(d)
                distance = d[assignment].sum()
                distance_array.loc[gene_list[i], gene_list[j]] = distance
    return distance_array


def build_graph(gmm_dict: dict,
                method='optimized_match',
                distance_threshold: int = 1):
    """
    Build graph by distance matrix
    :param method: weight_match or optimized_match
    :type method:
    :param gmm_dict:
    :type gmm_dict: dict
    :param distance_threshold:
    :type distance_threshold:
    :return:
    :rtype:
    """
    # build graph
    graph = nx.Graph()
    gene_list = list(gmm_dict.keys())
    gene_counts = len(gene_list)
    # add nodes
    for gene_id in gene_list:
        graph.add_node(gene_id)
    # calculate the weight and add edges
    for i in tqdm(range(gene_counts)):
        for j in range(gene_counts):
            if i != j and not graph.has_edge(gene_list[i], gene_list[j]):
                distance = distribution_distance(gmm_dict[gene_list[i]], gmm_dict[gene_list[j]], method=method)
                if distance < distance_threshold:
                    weight = 1 / distance
                    graph.add_edge(gene_list[i], gene_list[j], weight=weight)
    return graph


def cut_graph(graph):
    node = []
    for i in graph.nodes:
        if len(list(graph.neighbors(i))) == 0:
            node.append(i)
    for i in node:
        graph.remove_node(i)
    return graph


def cluster_graph(graph, k=2):
    similarity_matrix = nx.to_numpy_array(graph)
    clustering_model = SpectralClustering(n_clusters=k, affinity='precomputed')
    clustering_model.fit(similarity_matrix)
    return clustering_model.labels_
