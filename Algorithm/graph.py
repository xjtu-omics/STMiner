import networkx as nx
from distribution import *


def build_graph(gmm_dict):
    gmm_dict: dict
    # build graph
    graph = nx.Graph()
    gene_list = list(gmm_dict.keys())
    gene_counts = len(gene_list)
    for gene_id in gene_list:
        graph.add_node(gene_id)
    # calculate the weight
    for i in range(gene_counts):
        for j in range(gene_counts):
            if i != j and not graph.has_edge(gene_list[i], gene_list[j]):
                distance = distribution_distance(gmm_dict[gene_list[i]], gmm_dict[gene_list[j]])
                graph.add_edge(gene_list[i], gene_list[j], weight=distance)
