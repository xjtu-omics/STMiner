from STMiner.SPFinder import SPFinder
import anndata as ad
import numpy as np


def test_spfinder():
    spfinder = SPFinder()
    assert spfinder.adata is None

    # Test set_adata
    adata = ad.read_h5ad("./tests/data/test.h5ad")
    spfinder.set_adata(adata)
    assert isinstance(spfinder.adata, ad.AnnData)
    assert (spfinder.adata.obs['x'].max() == 63)
    
    # Test merge_bin function
    spfinder.merge_bin(bin_width=2)
    assert (spfinder.adata.obs['x'].max() == 31)
    
    # Test IO & merge_bin
    spfinder.read_h5ad("tests/data/test.h5ad", bin_size=5, merge_bin=True)
    assert isinstance(spfinder.adata, ad.AnnData)    
    
    # Check bin_size & merge_bin result
    assert spfinder.adata.shape == (282, 32268)

    # Test data preprocess
    spfinder.get_genes_csr_array(min_cells=200)
    assert len(spfinder.csr_dict) > 0
    assert np.mean(spfinder.csr_dict["rpl18a"].todense()) > 0

    # Test SVG detection
    spfinder.spatial_high_variable_genes()
    assert len(spfinder.global_distance) > 0
    
    # Check SVG detection result
    assert spfinder.global_distance.iloc[0, 0] == 'vmhcl'
    assert spfinder.global_distance.iloc[1, 0] == 'pvalb1'

    # Test fitting
    spfinder.fit_pattern(
        n_comp=10, gene_list=list(spfinder.global_distance[:200]["Gene"])
    )
    assert len(spfinder.patterns) > 0

    # Test GMM OT
    spfinder.build_distance_array()
    assert len(spfinder.genes_distance_array) > 0

    # Test gene cluster
    spfinder.cluster_gene(n_clusters=2, mds_components=20)
    assert len(spfinder.genes_labels) > 0
    spfinder.cluster_gene(n_clusters=2, mds_components=10)
    assert len(spfinder.genes_labels) > 0
    spfinder.cluster_gene(n_clusters=3, mds_components=20)
    assert len(spfinder.genes_labels) > 0
    
    # Test plot genes
    spfinder.plot.plot_gene('vmhcl', s=100, vmax=99)
    spfinder.plot.plot_gene('vmhcl', s=1, vmax=90)

    # Test plot pattern
    spfinder.get_pattern_array(vote_rate=0.2)
    spfinder.plot.plot_pattern(heatmap=False)

    # Test plot_gmm
    spfinder.plot_gmm("vmhcl")

    # Test plot_intersection
    spfinder.get_pattern_array()
    spfinder.plot.plot_intersection(pattern_list=[1,2])

    # Test get_custom_pattern
    spfinder.get_custom_pattern(["pvalb1", "myhc4", "vmhcl"], n_components=10, vote_rate=0, mode="vote")
    assert (spfinder.custom_pattern is not None)
    
    # Test compare_gene_to_genes
    df = spfinder.compare_gene_to_genes("pvalb1")
    assert (df.loc['pvalb1'] is not None)

    
