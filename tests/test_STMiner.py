from STMiner.SPFinder import SPFinder
import anndata as ad
import numpy as np


def test_spfinder():
    spfinder = SPFinder()
    assert spfinder.adata is None
    print("[PASS] SPFinder initialized, adata is None")

    # Test set_adata
    adata = ad.read_h5ad("./tests/data/test.h5ad")
    spfinder.set_adata(adata)
    assert isinstance(spfinder.adata, ad.AnnData)
    assert (spfinder.adata.obs['x'].max() == 63)
    print(f"[PASS] set_adata succeeded, data shape: {spfinder.adata.shape}, "
          f"max x: {spfinder.adata.obs['x'].max()}")

    # Test merge_bin function
    spfinder.merge_bin(bin_width=2)
    assert (spfinder.adata.obs['x'].max() == 31)
    print(f"[PASS] merge_bin(bin_width=2) succeeded, merged max x: "
          f"{spfinder.adata.obs['x'].max()}")

    # Test IO & merge_bin
    spfinder.read_h5ad("tests/data/test.h5ad", bin_size=5, merge_bin=True)
    assert isinstance(spfinder.adata, ad.AnnData)
    print("[PASS] read_h5ad(bin_size=5, merge_bin=True) loaded and merged")

    # Check bin_size & merge_bin result
    assert spfinder.adata.shape == (282, 32268)
    print(f"[PASS] bin merge result verified, data shape: {spfinder.adata.shape}")

    # Test data preprocess
    spfinder.get_genes_csr_array(min_cells=200)
    assert len(spfinder.csr_dict) > 0
    assert np.mean(spfinder.csr_dict["rpl18a"].todense()) > 0
    print(f"[PASS] get_genes_csr_array preprocessed, {len(spfinder.csr_dict)} genes, "
          f"rpl18a mean: {np.mean(spfinder.csr_dict['rpl18a'].todense()):.4f}")

    # Test SVG detection
    spfinder.spatial_high_variable_genes()
    assert len(spfinder.global_distance) > 0
    print(f"[PASS] spatial variable gene (SVG) detection found "
          f"{len(spfinder.global_distance)} genes")

    # Check SVG detection result
    assert spfinder.global_distance.iloc[0, 0] == 'vmhcl'
    assert spfinder.global_distance.iloc[1, 0] == 'pvalb1'
    print(f"[PASS] SVG result verified, Top1: {spfinder.global_distance.iloc[0, 0]}, "
          f"Top2: {spfinder.global_distance.iloc[1, 0]}")

    # Test fitting
    spfinder.fit_pattern(
        n_comp=10, gene_list=list(spfinder.global_distance[:200]["Gene"])
    )
    assert len(spfinder.patterns) > 0
    print(f"[PASS] fit_pattern succeeded, {len(spfinder.patterns)} patterns fitted")

    # Test GMM OT
    spfinder.build_distance_array()
    assert len(spfinder.genes_distance_array) > 0
    print(f"[PASS] build_distance_array succeeded, distance matrix shape: "
          f"{spfinder.genes_distance_array.shape}")

    # Test gene cluster
    spfinder.cluster_gene(n_clusters=2, mds_components=20)
    assert len(spfinder.genes_labels) > 0
    print(f"[PASS] cluster_gene(n_clusters=2, mds_components=20) succeeded, "
          f"labels: {len(spfinder.genes_labels)}")
    spfinder.cluster_gene(n_clusters=2, mds_components=10)
    assert len(spfinder.genes_labels) > 0
    print(f"[PASS] cluster_gene(n_clusters=2, mds_components=10) succeeded, "
          f"labels: {len(spfinder.genes_labels)}")
    spfinder.cluster_gene(n_clusters=3, mds_components=20)
    assert len(spfinder.genes_labels) > 0
    print(f"[PASS] cluster_gene(n_clusters=3, mds_components=20) succeeded, "
          f"labels: {len(spfinder.genes_labels)}")

    # Test plot genes
    spfinder.plot.plot_gene('vmhcl', s=100, vmax=99)
    spfinder.plot.plot_gene('vmhcl', s=1, vmax=90)
    print("[PASS] plot_gene succeeded (s=100/vmax=99, s=1/vmax=90)")

    # Test plot pattern
    spfinder.get_pattern_array(vote_rate=0.2)
    spfinder.plot.plot_pattern(heatmap=False)
    print("[PASS] get_pattern_array + plot_pattern succeeded")

    # Test plot_gmm
    spfinder.plot_gmm("vmhcl")
    print("[PASS] plot_gmm('vmhcl') succeeded")

    # Test plot_intersection
    spfinder.get_pattern_array()
    spfinder.plot.plot_intersection(pattern_list=[1,2])
    print("[PASS] plot_intersection(pattern_list=[1, 2]) succeeded")

    # Test get_custom_pattern
    spfinder.get_custom_pattern(["pvalb1", "myhc4", "vmhcl"], n_components=10, vote_rate=0, mode="vote")
    assert (spfinder.custom_pattern is not None)
    print("[PASS] get_custom_pattern succeeded, custom_pattern is not None")

    # Test compare_gene_to_genes
    df = spfinder.compare_gene_to_genes("pvalb1")
    assert (df.loc['pvalb1'] is not None)
    print(f"[PASS] compare_gene_to_genes('pvalb1') succeeded, result shape: {df.shape}")

    print("=" * 50)
    print("All tests passed successfully!")

    
