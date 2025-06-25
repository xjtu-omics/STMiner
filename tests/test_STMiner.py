from STMiner.SPFinder import SPFinder
import anndata as ad
import numpy as np


def test_spfinder():
    spfinder = SPFinder()
    assert spfinder.adata is None

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
    
    # Test plot
    spfinder.plot.plot_gene('vmhcl', s=100, vmax=99)

