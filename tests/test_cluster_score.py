"""Clustering evaluation must not overwrite the active analysis."""

from copy import deepcopy
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import cdist

from STMiner import SPFinder
from STMiner.services import spfinder_analysis as analysis


@pytest.fixture
def sp():
    obj = SPFinder()
    points = np.array([[0, 0], [0, 1], [5, 5], [5, 6], [10, 0], [10, 1]], dtype=np.float64)
    names = list("ABCDEF")
    obj.genes_distance_array = pd.DataFrame(cdist(points, points), index=names, columns=names)
    obj.filtered_distance_array = obj.genes_distance_array.iloc[:4, :4].copy()
    obj.genes_labels = pd.DataFrame({"gene_id": names, "labels": [0, 0, 1, 1, 2, 2]})
    obj.kmeans_fit_result = SimpleNamespace(labels_=np.array([0, 0, 1, 1, 2, 2]))
    obj.mds_features = points.copy()
    obj.patterns_matrix_dict = {0: np.array([[0.125, 1.25]])}
    obj.patterns_binary_matrix_dict = {0: np.array([[1, 1]])}
    obj.all_labels = pd.DataFrame({"label": [0], "value": [0.25]}, index=["A"])
    yield obj
    plt.close("all")


def capture(sp):
    references = vars(sp).copy()
    values = deepcopy({key: value for key, value in references.items() if key != "plot"})
    return references, values


def assert_unchanged(sp, snapshot):
    references, values = snapshot
    assert vars(sp).keys() == references.keys()
    for key, reference in references.items():
        assert getattr(sp, key) is reference
    for key in ("genes_distance_array", "filtered_distance_array", "genes_labels", "all_labels"):
        pd.testing.assert_frame_equal(getattr(sp, key), values[key])
    np.testing.assert_array_equal(sp.mds_features, values["mds_features"])
    np.testing.assert_array_equal(sp.kmeans_fit_result.labels_, values["kmeans_fit_result"].labels_)
    for key in ("patterns_matrix_dict", "patterns_binary_matrix_dict"):
        assert getattr(sp, key).keys() == values[key].keys()
        for label, arr in values[key].items():
            np.testing.assert_array_equal(getattr(sp, key)[label], arr)
            assert getattr(sp, key)[label].dtype == arr.dtype


def test_real_clustering_plot_preserves_existing_results(sp, monkeypatch):
    snapshot = capture(sp)
    shown = []
    monkeypatch.setattr(plt, "show", lambda: shown.append(True))
    sp.plot.plot_cluster_score(2, 2, 3)
    assert shown == [True]
    assert len(plt.gca().lines) >= 3
    assert_unchanged(sp, snapshot)


@pytest.mark.parametrize("stage", ["cluster", "metric", "render"])
def test_failed_evaluation_or_plot_preserves_existing_results(sp, monkeypatch, stage):
    snapshot = capture(sp)
    calls = []

    def fail(*args, **kwargs):
        raise RuntimeError("evaluation failed")

    def cluster(*args, **kwargs):
        calls.append(kwargs["n_clusters"])
        if stage == "cluster" and len(calls) == 2:
            fail()
        return None, SimpleNamespace(labels_=np.array([0, 0, 1, 1, 2, 2])), None

    monkeypatch.setattr(analysis, "cluster", cluster)
    if stage == "metric":
        monkeypatch.setattr(analysis, "silhouette_score", fail)
    elif stage == "render":
        monkeypatch.setattr(plt, "show", fail)
    with pytest.raises(RuntimeError, match="evaluation failed"):
        sp.plot.plot_cluster_score(2, 2, 3)
    assert calls == ([2] if stage == "metric" else [2, 3])
    assert_unchanged(sp, snapshot)


def test_evaluation_preserves_metric_inputs_and_raw_scores(sp, monkeypatch):
    calls = []
    labels = np.array([0, 0, 1, 1, 2, 2])

    def cluster(matrix, **kwargs):
        assert matrix is sp.genes_distance_array
        calls.append(kwargs)
        return None, SimpleNamespace(labels_=labels), None

    def metric(value):
        def score(matrix, assigned):
            assert matrix is sp.genes_distance_array
            assert assigned is labels
            return value
        return score

    monkeypatch.setattr(analysis, "cluster", cluster)
    monkeypatch.setattr(analysis, "davies_bouldin_score", metric(2.0))
    monkeypatch.setattr(analysis, "calinski_harabasz_score", metric(8.0))
    monkeypatch.setattr(analysis, "silhouette_score", metric(0.25))
    result = analysis.evaluate_cluster_scores(sp.genes_distance_array, 4, 2, 3)
    assert calls == [dict(n_clusters=2, mds_components=4), dict(n_clusters=3, mds_components=4)]
    expected = pd.DataFrame(
        {"1/Davies-Bouldin": [0.5, 0.5], "Silhouette": [0.25, 0.25], "Calinski-Harabasz": [8.0, 8.0]},
        index=[2, 3],
    )
    pd.testing.assert_frame_equal(result, expected)
