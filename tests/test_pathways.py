import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from STMiner import SPFinder, KEGGFinder


@pytest.fixture
def sp():
    matrix = np.array([[0.125, 0.375, 9], [0.5, 1.5, 6], [2, 4, 3]], dtype=np.float64)
    obj = SPFinder(
        ad.AnnData(
            X=csr_matrix(matrix),
            obs=pd.DataFrame(
                {"x": [0, 1, 1], "y": [0, 0, 1]}, index=["s1", "s2", "s3"]
            ),
            var=pd.DataFrame(index=["A", "B", "C"]),
        )
    )
    obj.genes_labels = pd.DataFrame({"gene_id": ["A", "B", "C"], "labels": [0, 0, 1]})
    return obj


@pytest.fixture
def enrichment():
    return pd.DataFrame(
        [
            dict(
                pattern=0,
                query="pattern_0",
                source="KEGG",
                native="KEGG:1",
                name="Pathway one",
                p_value=0.001,
                intersections=["A", "B"],
            ),
            dict(
                pattern=1,
                query="pattern_1",
                source="GO:BP",
                native="GO:1",
                name="Process two",
                p_value=0.002,
                intersections=["C"],
            ),
        ]
    )


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def test_pattern_enrichment_delegates_and_keeps_provenance(sp, monkeypatch):
    captured = {}

    def fake(self, queries, **kwargs):
        captured.update(queries=queries, **kwargs)
        return pd.DataFrame(
            {"query": list(queries), "intersections": list(queries.values())}
        )

    monkeypatch.setattr(KEGGFinder, "enrich_gene_set", fake)
    result = sp.enrich_patterns(organism="drerio", labels=[1, 0], sources=["KEGG"])
    assert captured["queries"] == {"pattern_0": ["C"], "pattern_1": ["A", "B"]}
    assert captured["background"] == ["A", "B", "C"]
    assert captured["organism"] == "drerio"
    assert result["pattern"].tolist() == [1, 0]
    assert result.attrs["pattern_gene_sets"] == {1: ["C"], 0: ["A", "B"]}


@pytest.mark.parametrize("sparse", [True, False])
def test_spatial_scores_use_hits_not_the_whole_pattern_and_preserve_data(
    sp, enrichment, sparse
):
    if not sparse:
        sp.adata.X = sp.adata.X.toarray()
    before = sp.adata.X.copy()
    original = enrichment.copy(deep=True)
    scores, terms = sp.score_pathways(enrichment, min_genes=1)
    np.testing.assert_array_equal(scores.iloc[:, 0], [0.25, 1, 3])
    np.testing.assert_array_equal(scores.iloc[:, 1], [9, 6, 3])
    assert scores.index.equals(sp.adata.obs_names)
    assert terms.loc["pathway_0", "genes"] == ("A", "B")
    assert sp.adata.X.dtype == before.dtype
    if sparse:
        np.testing.assert_array_equal(sp.adata.X.data, before.data)
    else:
        np.testing.assert_array_equal(sp.adata.X, before)
    pd.testing.assert_frame_equal(original, enrichment)


def test_explicit_mapping_filters_out_other_pattern_genes(sp, enrichment):
    enrichment.at[0, "intersections"] = ["ENS_A", "B", "C", "unknown"]
    with pytest.warns(UserWarning, match="unmapped"):
        scores, terms = sp.score_pathways(
            enrichment, labels=[0], gene_map={"ENS_A": "A"}
        )
    np.testing.assert_array_equal(scores.iloc[:, 0], [0.25, 1, 3])
    assert terms.iloc[0]["unmapped_genes"] == ("C", "unknown")
    assert terms.iloc[0]["n_genes"] == 2
    assert terms.iloc[0]["n_hits"] == 4


def test_existing_results_require_explicit_pattern_mapping(sp, enrichment):
    enrichment = enrichment.drop(columns="pattern")
    with pytest.raises(ValueError, match="pattern_map"):
        sp.score_pathways(enrichment)
    scores, terms = sp.score_pathways(
        enrichment, pattern_map={"pattern_0": 0, "pattern_1": 1}, min_genes=1
    )
    assert terms["pattern"].tolist() == [0, 1]


def test_layer_is_used_without_rounding(sp, enrichment):
    sp.adata.layers["normalized"] = sp.adata.X.astype(np.float32) * 0.5
    scores, terms = sp.score_pathways(enrichment, labels=[0], layer="normalized")
    np.testing.assert_array_equal(scores.iloc[:, 0], [0.125, 0.5, 1.5])
    assert scores.iloc[:, 0].dtype == np.float32
    assert terms.iloc[0]["layer"] == "normalized"


def test_distinct_terms_and_same_term_in_other_patterns_use_their_own_hits(
    sp, enrichment
):
    enrichment.at[1, "source"] = "KEGG"
    enrichment.at[1, "native"] = "KEGG:1"
    another = enrichment.iloc[[0]].copy(deep=True)
    another.at[0, "native"] = "KEGG:2"
    another.at[0, "intersections"] = ["A"]
    enrichment = pd.concat([enrichment, another], ignore_index=True)
    scores, terms = sp.score_pathways(enrichment, min_genes=1, top_n=2)
    np.testing.assert_array_equal(scores.iloc[:, 0], [0.25, 1, 3])
    np.testing.assert_array_equal(scores.iloc[:, 1], [0.125, 0.5, 2])
    np.testing.assert_array_equal(scores.iloc[:, 2], [9, 6, 3])
    assert terms["pattern"].tolist() == [0, 0, 1]


def test_pathway_rendering_exports_and_does_not_overwrite_patterns(
    sp, enrichment, tmp_path
):
    original_patterns = sp.patterns_matrix_dict
    original_labels = sp.genes_labels.copy(deep=True)
    fig, axes, scores, terms = sp.plot.plot_pathways(
        enrichment,
        min_genes=1,
        display="expression",
        heatmap=True,
        output_path=tmp_path / "paths.svg",
        show=False,
    )
    assert axes.shape == (1, 2)
    assert "KEGG:1" in axes[0, 0].get_title(loc="left")
    assert "genes 2/2" in axes[0, 0].get_title(loc="left")
    assert axes[0, 0].images[0].get_array()[0, 0] == 0.25
    assert axes[0, 0].images[0].norm is axes[0, 1].images[0].norm
    assert sp.patterns_matrix_dict is original_patterns
    pd.testing.assert_frame_equal(original_labels, sp.genes_labels)
    assert (tmp_path / "paths.svg").stat().st_size > 0


def test_duplicate_grid_positions_average_for_display_only(sp, enrichment):
    sp.adata.obs.loc["s2", ["x", "y"]] = [0, 0]
    _, axes, scores, _ = sp.plot.plot_pathways(
        enrichment, labels=[0], display="expression", heatmap=True, show=False
    )
    assert axes[0, 0].images[0].get_array()[0, 0] == 0.625
    np.testing.assert_array_equal(scores.iloc[:, 0], [0.25, 1, 3])


def test_no_hits_and_stale_membership_fail_explicitly(sp, enrichment):
    with pytest.raises(ValueError, match="enough mapped"):
        sp.score_pathways(enrichment, min_genes=3)
    enrichment.attrs["pattern_gene_sets"] = {0: ["A"], 1: ["C"]}
    with pytest.raises(ValueError, match="membership changed"):
        sp.score_pathways(enrichment)


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"top_n": 0}, "top_n"),
        ({"p_value_cutoff": 0}, "p_value_cutoff"),
        ({"layer": "absent"}, "layer"),
        ({"term_ids": ["bad"]}, "Unknown term_ids"),
        ({"sources": ["GO:CC"]}, "No enriched"),
    ],
)
def test_invalid_selection_is_actionable(sp, enrichment, kwargs, match):
    with pytest.raises(ValueError, match=match):
        sp.score_pathways(enrichment, **kwargs)


def test_no_silent_serialized_list_or_negative_expression(sp, enrichment):
    enrichment.at[0, "intersections"] = "['A', 'B']"
    with pytest.raises(ValueError, match="serialized"):
        sp.score_pathways(enrichment)
    enrichment.at[0, "intersections"] = ["A", "B"]
    sp.adata.X.data[0] = -1
    with pytest.raises(ValueError, match="nonnegative"):
        sp.score_pathways(enrichment)


def test_continuous_density_has_unit_mass_without_altering_spot_scores(sp, enrichment):
    expected, _ = sp.score_pathways(enrichment, min_genes=1)
    _, axes, scores, terms = sp.plot.plot_pathways(
        enrichment, min_genes=1, contours=0, show=False
    )
    pd.testing.assert_frame_equal(scores, expected)
    for ax in axes.flat:
        assert len(ax.images) == 1
        assert not ax.collections
        assert ax.images[0].get_interpolation() == "bilinear"
        density = ax.images[0].get_array()
        assert density.sum() == pytest.approx(1)
        assert density[1, 0] > 0
    assert terms.attrs["display"] == "density"


def test_density_obeys_explicit_tissue_mask_and_bandwidth(sp, enrichment):
    mask = np.array([[True, False], [True, True]])
    _, axes, _, terms = sp.plot.plot_pathways(
        enrichment, labels=[0], tissue_mask=mask, bandwidth=0.5, show=False
    )
    density = axes[0, 0].images[0].get_array()
    assert np.ma.getmaskarray(density)[1, 0]
    assert density.sum() == pytest.approx(1)
    assert terms.attrs["bandwidth"] == 0.5
    with pytest.raises(ValueError, match="tissue_mask"):
        sp.plot.plot_pathways(
            enrichment, tissue_mask=np.ones((3, 3), dtype=bool), show=False
        )
    with pytest.raises(ValueError, match="bandwidth"):
        sp.plot.plot_pathways(enrichment, bandwidth=0, show=False)
