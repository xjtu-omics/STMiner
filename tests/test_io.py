from STMiner.SPFinder import SPFinder
import anndata as ad


def test_spfinder():
    spfinder = SPFinder()
    spfinder.read_h5ad("tests/data/test.h5ad")
    assert isinstance(spfinder.adata, ad.AnnData)
