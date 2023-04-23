import networkx as nx


def build_graph(gmm_dict):
    gmm_dict: dict
    # build graph
    graph = nx.Graph()
    for gene_id in gmm_dict:
        graph.add_node(gene_id)
    # calculate the weight

