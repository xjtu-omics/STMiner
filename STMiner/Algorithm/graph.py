import networkx as nx
from sklearn.cluster import SpectralClustering

from STMiner.Algorithm.distance import distribution_distance
from STMiner.Algorithm.distribution import *


class Graph:
    def __init__(self):
        pass


def build_graph(gmm_dict: dict,
                distance_threshold: int = 1):
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
