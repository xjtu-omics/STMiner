"""Regression coverage for repeated analysis and upstream state changes."""

import importlib

import matplotlib

matplotlib.use("Agg")
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from STMiner import SPFinder
from STMiner.services import spfinder_analysis, spfinder_io


spfinder_module = importlib.import_module("STMiner.SPFinder")
CLUSTER_FIELDS = ("genes_labels", "kmeans_fit_result", "mds_features")
DATA_FIELDS = (
    "patterns", "image_gmm", "global_distance", "custom_pattern",
    "genes_distance_array", "filtered_distance_array", "all_labels",
) + CLUSTER_FIELDS


def make_adata(prefix="A"):
    return ad.AnnData(
        X=np.array([[1.25, 0], [0, 2.5], [3.75, 5]], dtype=np.float64),
        obs=pd.DataFrame({"x": [0, 2, 4], "y": [0, 1, 2]}, index=["s0", "s1", "s2"]),
        var=pd.DataFrame(index=[f"{prefix}0", f"{prefix}1"]),
    )


def seed_results(sp):
    for name in DATA_FIELDS:
        setattr(sp, name, object())
    sp.csr_dict = {"old_gene": object()}
    sp.patterns_matrix_dict = {0: np.ones((5, 3)), 9: np.ones((5, 3))}
    sp.patterns_binary_matrix_dict = {0: np.ones((5, 3)), 9: np.ones((5, 3))}
    sp._highly_variable_genes = ["old_gene"]
    sp._gene_expression_edge = {"old_gene": object()}


def assert_no_arrays(sp):
    assert sp.patterns_matrix_dict == {}
    assert sp.patterns_binary_matrix_dict == {}
    assert sp.all_labels is None


def assert_no_clusters(sp):
    for name in CLUSTER_FIELDS:
        assert getattr(sp, name) is None
    assert sp.filtered_distance_array is None
    assert_no_arrays(sp)


def assert_no_data_results(sp):
    for name in DATA_FIELDS:
        assert getattr(sp, name) is None
    assert sp.csr_dict == {}
    assert sp._gene_expression_edge == {}
    assert_no_arrays(sp)


@pytest.fixture
def sp():
    obj = SPFinder(make_adata())
    seed_results(obj)
    return obj


@pytest.fixture
def hvg_stub(monkeypatch):
    # Keep real filtering/normalization, replacing only HVG fitting on tiny data.
    def select(adata, **kwargs):
        adata.var["highly_variable"] = True

    monkeypatch.setattr(spfinder_module.sc.pp, "highly_variable_genes", select)


@pytest.mark.parametrize("same_object", [False, True])
def test_set_adata_discards_results_even_for_same_object(sp, same_object):
    data = sp.adata if same_object else make_adata("B")
    before = data.X.copy()
    plot = sp.plot
    sp.set_adata(data)
    assert sp.adata is data
    assert sp.plot is plot and plot.sp is sp
    assert_no_data_results(sp)
    assert sp._highly_variable_genes == []
    assert sp._scope == (0, 4)
    assert sp.adata.X.dtype == before.dtype
    np.testing.assert_array_equal(sp.adata.X, before)


@pytest.mark.parametrize("method,reader", [
    ("read_h5ad", "read_h5ad"),
    ("read_bmk_dir", "read_bmk"),
    ("read_gem", "read_gem_file"),
])
def test_readers_invalidate_previous_dataset(sp, monkeypatch, method, reader):
    data = make_adata("B")
    monkeypatch.setattr(spfinder_io, reader, lambda *args, **kwargs: data)
    getattr(sp, method)("unused")
    assert sp.adata is data
    assert_no_data_results(sp)


def test_invalid_coordinates_cannot_leave_old_results(sp):
    invalid = make_adata("B")
    invalid.obs.drop(columns=["x", "y"], inplace=True)
    with pytest.raises(ValueError, match="must contain either"):
        sp.set_adata(invalid)
    assert_no_data_results(sp)
    assert sp._scope is None


def test_merge_bin_invalidates_models_and_preserves_expression(sp):
    before = sp.adata.X.copy()
    sp.merge_bin(2)
    assert_no_data_results(sp)
    assert sp._highly_variable_genes == []
    assert sp._scope == (0, 2)
    np.testing.assert_array_equal(sp.adata.obs["x"], [0, 1, 2])
    np.testing.assert_array_equal(sp.adata.X, before)
    assert sp.adata.X.dtype == before.dtype


def test_preprocess_invalidates_results_and_updates_scope(sp, hvg_stub):
    sp.adata.X[2] = 0
    sp.preprocess(False, False, False, min_cells=1, min_genes=1)
    assert_no_data_results(sp)
    assert sp._highly_variable_genes == ["A0", "A1"]
    assert sp._scope == (0, 2)
    assert sp.adata.shape == (2, 2)
    assert sp.adata.X.dtype == np.float64


def test_failed_preprocess_clears_results_before_partial_mutation(sp, monkeypatch):
    def fail(adata, **kwargs):
        adata.X[0, 0] = 99
        raise RuntimeError("filter failed")

    monkeypatch.setattr(spfinder_module.sc.pp, "filter_genes", fail)
    with pytest.raises(RuntimeError, match="filter failed"):
        sp.preprocess(False, False, False, 1, 1)
    assert sp.adata.X[0, 0] == 99
    assert_no_data_results(sp)
    assert sp._highly_variable_genes == []


def test_rebuilding_csr_invalidates_previous_analysis(sp, hvg_stub):
    sp.get_genes_csr_array(min_cells=1, normalize=False)
    assert set(sp.csr_dict) == {"A0", "A1"}
    assert sp.patterns is None and sp.global_distance is None
    assert sp.genes_distance_array is None
    assert_no_clusters(sp)
    assert sp.csr_dict["A0"].dtype == np.float64
    assert sp.csr_dict["A0"][0, 0] == 1.25


def test_new_dataset_ot_never_uses_previous_gene_matrices(sp, hvg_stub, monkeypatch):
    # Use enough observations for the existing automatic min_cells=50 default.
    data = ad.concat([make_adata("B")] * 30, index_unique="-")
    sp.set_adata(data)
    observed = []

    def distance(background, gene):
        observed.append((background.copy(), gene.copy()))
        return float(gene.sum())

    monkeypatch.setattr(spfinder_module, "calculate_ot_distance", distance)
    sp.spatial_high_variable_genes()
    assert set(sp.global_distance["Gene"]) == {"B0", "B1"}
    assert set(sp.csr_dict) == {"B0", "B1"}
    assert len(observed) == 2
    np.testing.assert_allclose(
        observed[0][0].toarray(), (observed[0][1] + observed[1][1]).toarray()
    )


@pytest.mark.parametrize("fail", [False, True])
def test_refit_never_keeps_old_models_or_downstream_results(sp, hvg_stub, monkeypatch, fail):
    fitted = {"A0": object()}

    def fit(*args, **kwargs):
        if fail:
            raise RuntimeError("fit failed")
        return fitted

    monkeypatch.setattr(spfinder_module, "fit_gmms", fit)
    sp.fit_pattern(min_cells=1)
    assert sp.patterns is (None if fail else fitted)
    assert sp.csr_dict == {} and sp.global_distance is None
    assert sp.genes_distance_array is None
    assert_no_clusters(sp)


@pytest.mark.parametrize("method,builder", [
    ("gmm", "build_gmm_distance_array"),
    ("mse", "build_mse_distance_array"),
    ("cs", "build_cosine_similarity_array"),
    ("ot", "build_ot_distance_array"),
])
@pytest.mark.parametrize("fail", [False, True])
def test_distance_rebuild_invalidates_only_downstream(sp, monkeypatch, method, builder, fail):
    upstream = (sp.patterns, sp.csr_dict, sp.global_distance, sp.custom_pattern)
    matrix = pd.DataFrame([[0.0, 1.0], [1.0, 0.0]], index=["A0", "A1"], columns=["A0", "A1"])

    def build(*args):
        if fail:
            raise RuntimeError("distance failed")
        return matrix

    monkeypatch.setattr(spfinder_analysis, builder, build)
    if fail:
        with pytest.raises(RuntimeError, match="distance failed"):
            sp.build_distance_array(method)
        assert sp.genes_distance_array is None
    else:
        sp.build_distance_array(method)
        assert sp.genes_distance_array is matrix
    assert_no_clusters(sp)
    assert (sp.patterns, sp.csr_dict, sp.global_distance, sp.custom_pattern) == upstream


@pytest.mark.parametrize("filtered", [False, True])
@pytest.mark.parametrize("fail", [False, True])
def test_recluster_discards_old_maps_and_labels(sp, monkeypatch, filtered, fail):
    matrix = pd.DataFrame([[0.0, 1.0], [1.0, 0.0]], index=["A0", "A1"], columns=["A0", "A1"])
    sp.genes_distance_array = matrix
    labels = pd.DataFrame({"gene_id": ["A0", "A1"], "labels": [0, 0]})
    model, features = object(), object()

    def cluster(*args, **kwargs):
        if fail:
            raise RuntimeError("cluster failed")
        return labels, model, features

    monkeypatch.setattr(spfinder_analysis, "cluster", cluster)
    if fail:
        with pytest.raises(RuntimeError, match="cluster failed"):
            sp.cluster_gene(1, use_highly_variable_gene=filtered)
        for name in CLUSTER_FIELDS:
            assert getattr(sp, name) is None
    else:
        sp.cluster_gene(1, use_highly_variable_gene=filtered)
        assert sp.genes_labels is labels
        assert sp.kmeans_fit_result is model and sp.mds_features is features
    assert sp.genes_distance_array is matrix
    if not filtered:
        assert sp.filtered_distance_array is None
    assert_no_arrays(sp)


def test_pattern_rebuild_removes_obsolete_categories(sp):
    sp.genes_labels = pd.DataFrame({"gene_id": ["A0", "A1"], "labels": [0, 0]})
    sp.get_pattern_array()
    assert set(sp.patterns_matrix_dict) == {0}
    assert set(sp.patterns_binary_matrix_dict) == {0}
    assert sp.all_labels is None


def test_failed_pattern_rebuild_does_not_publish_partial_maps(sp, monkeypatch):
    sp.genes_labels = pd.DataFrame({"gene_id": ["A0", "A1"], "labels": [0, 1]})
    calls = []

    def build(*args):
        calls.append(args)
        if len(calls) == 2:
            raise RuntimeError("pattern failed")
        return np.ones((5, 3)), np.ones((5, 3))

    monkeypatch.setattr(sp, "_genes_to_pattern", build)
    with pytest.raises(RuntimeError, match="pattern failed"):
        sp.get_pattern_array()
    assert len(calls) == 2
    assert_no_arrays(sp)
