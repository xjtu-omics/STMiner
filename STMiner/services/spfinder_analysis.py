from typing import TYPE_CHECKING

import pandas as pd

from STMiner.Algorithm.algorithm import cluster
from STMiner.Algorithm.distance import (
    build_cosine_similarity_array,
    build_gmm_distance_array,
    build_mse_distance_array,
    build_ot_distance_array,
)


if TYPE_CHECKING:
    from STMiner.SPFinder import SPFinder


def build_distance_array_for_spfinder(spfinder: "SPFinder", method="gmm", gene_list=None):
    if gene_list is None:
        gene_list = list(spfinder.adata.var.index)

    if method == "gmm":
        spfinder.genes_distance_array = build_gmm_distance_array(spfinder.patterns)
    elif method == "mse":
        spfinder.genes_distance_array = build_mse_distance_array(spfinder.adata, gene_list)
    elif method == "cs":
        spfinder.genes_distance_array = build_cosine_similarity_array(
            spfinder.adata, gene_list
        )
    elif method == "ot":
        spfinder.genes_distance_array = build_ot_distance_array(
            spfinder.csr_dict, gene_list
        )
    else:
        raise ValueError("Unknown method, method should be one of gmm, mse, cs, ot")


def cluster_gene_for_spfinder(
    spfinder: "SPFinder",
    n_clusters: int,
    mds_components=20,
    use_highly_variable_gene=False,
    n_top_genes=500,
):
    if use_highly_variable_gene:
        df = pd.DataFrame(spfinder.genes_distance_array.mean(axis=1), columns=["mean"])
        df = df.sort_values(by="mean", ascending=False)
        hv_gene_list = list(df[:n_top_genes].index)
        spfinder.filtered_distance_array = spfinder.genes_distance_array.loc[
            hv_gene_list, hv_gene_list
        ]
        spfinder.genes_labels, spfinder.kmeans_fit_result, spfinder.mds_features = cluster(
            spfinder.filtered_distance_array,
            n_clusters=n_clusters,
            mds_components=mds_components,
        )
    else:
        spfinder.genes_labels, spfinder.kmeans_fit_result, spfinder.mds_features = cluster(
            spfinder.genes_distance_array,
            n_clusters=n_clusters,
            mds_components=mds_components,
        )
