import networkx as nx
from utils import issparse
from Algorithm.distribution import *


def build_gmm_distance_array(gmm_dict, method='weight_match'):
    """
    Generate a distance matrix by the given gmm dictionary.
    :param method:
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


def build_mix_distance_array(adata, gmm_dict):
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
                                                     method='weight_match')
                # TODO:
                distance_array.loc[gene_list[i], gene_list[j]] = mse_distance + gmm_distance
    return distance_array


def build_graph(gmm_dict: dict, distance_threshold: int = 1):
    """
    Build graph by distance matrix
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
                distance = distribution_distance(gmm_dict[gene_list[i]], gmm_dict[gene_list[j]])
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
