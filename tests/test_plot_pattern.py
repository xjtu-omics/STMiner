from itertools import product
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from STMiner.Plot.plot import Plot


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


@pytest.fixture
def plot():
    a = np.zeros((4, 7), dtype=np.float64)
    a[0, 1], a[3, 5] = 1.23456789012345, 6.98765432109876
    sp = SimpleNamespace(
        genes_labels=pd.DataFrame({"labels": [2, 0, 2], "gene_id": ["A", "B", "C"]}),
        patterns_matrix_dict={2: a, 0: a * 2},
        adata=SimpleNamespace(obs=pd.DataFrame({"x": [0, 3, 2], "y": [1, 5, 6]})),
    )
    return Plot(sp)


def test_sorted_panels_auto_layout_no_global_side_effects(plot, monkeypatch):
    previous = matplotlib.rcParams.copy()
    monkeypatch.setattr(plt, "show", lambda: pytest.fail("show=False must not display"))
    fig, axes = plot.plot_pattern(show=False)
    assert axes.shape == (1, 2)
    assert axes[0, 0].get_title(loc="left") == "Pattern 0 · 1 gene"
    assert axes[0, 1].get_title(loc="left") == "Pattern 2 · 2 genes"
    assert axes[0, 0].get_aspect() == 1
    for key in ["font.family", "font.size", "pdf.fonttype", "axes.grid"]:
        assert matplotlib.rcParams[key] == previous[key]
    assert len(fig.axes) == 4


@pytest.mark.parametrize(
    "rotate,reverse_x,reverse_y", list(product([False, True], repeat=3))
)
def test_spots_heatmap_and_background_share_coordinates(
    plot, rotate, reverse_x, reverse_y
):
    kwargs = dict(
        rotate=rotate,
        reverse_x=reverse_x,
        reverse_y=reverse_y,
        colorbar=False,
        show=False,
    )
    _, spot_axes = plot.plot_pattern(plot_bg=True, **kwargs)
    _, heat_axes = plot.plot_pattern(heatmap=True, **kwargs)
    for spot_ax, heat_ax in zip(spot_axes.flat, heat_axes.flat):
        # Every panel must receive its own measured-tissue background.
        assert len(spot_ax.collections) == 2
        background = set(map(tuple, spot_ax.collections[0].get_offsets()))
        spots = spot_ax.collections[1]
        coordinates = np.asarray(spots.get_offsets(), dtype=int)
        displayed = heat_ax.images[0].get_array()
        np.testing.assert_array_equal(
            displayed[coordinates[:, 1], coordinates[:, 0]], spots.get_array()
        )
        assert set(map(tuple, coordinates)).issubset(background)
        assert spot_ax.get_xlim() == heat_ax.get_xlim()
        assert spot_ax.get_ylim() == heat_ax.get_ylim()
        assert spots.norm.vmax == heat_ax.images[0].norm.vmax


def test_sparse_precision_order_and_shared_scale(plot):
    original = plot.sp.patterns_matrix_dict[2].copy()
    plot.sp.patterns_matrix_dict[2] = csr_matrix(original)
    _, axes = plot.plot_pattern(labels=[2, 0], shared_scale=True, show=False)
    norm = axes[0, 0].collections[0].norm
    assert norm is axes[0, 1].collections[0].norm
    np.testing.assert_array_equal(
        axes[0, 0].collections[0].get_array(), original[original > 0]
    )
    np.testing.assert_array_equal(plot.sp.patterns_matrix_dict[2].toarray(), original)
    assert plot.sp.patterns_matrix_dict[2].dtype == original.dtype


@pytest.mark.parametrize("heatmap", [False, True])
def test_empty_and_constant_signal(plot, heatmap):
    plot.sp.patterns_matrix_dict[0][:] = 0
    plot.sp.patterns_matrix_dict[2][:] = 7
    _, axes = plot.plot_pattern(heatmap=heatmap, show=False)
    assert any(text.get_text() == "No signal" for text in axes[0, 0].texts)
    artist = axes[0, 1].images[0] if heatmap else axes[0, 1].collections[0]
    assert artist.norm.vmin == 0
    assert artist.norm.vmax == 7
    assert np.isfinite(artist.norm(7))


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"labels": []}, "empty"),
        ({"labels": [99]}, "Unknown"),
        ({"num_cols": 0}, "num_cols"),
        ({"num_cols": 1.5}, "num_cols"),
        ({"s": float("nan")}, "s must"),
        ({"vmax": 101}, "percentile"),
        ({"vmax": 0}, "vmax"),
        ({"aspect": 0}, "aspect"),
        ({"figsize": (1, float("inf"))}, "figsize"),
        ({"cmap": "not-a-colormap"}, "colormap"),
        ({"heatmap": "yes"}, "boolean"),
    ],
)
def test_invalid_parameters_fail_before_creating_figure(plot, kwargs, match):
    with pytest.raises(ValueError, match=match):
        plot.plot_pattern(show=False, **kwargs)
    assert not plt.get_fignums()


@pytest.mark.parametrize(
    "value,match",
    [
        (np.nan, "finite"),
        (np.inf, "finite"),
        (-1, "negative"),
    ],
)
def test_invalid_intensity_has_actionable_error(plot, value, match):
    plot.sp.patterns_matrix_dict[0][0, 0] = value
    with pytest.raises(ValueError, match=match):
        plot.plot_pattern(show=False)
    assert not plt.get_fignums()


def test_missing_state_and_coordinate_bounds(plot):
    plot.sp.adata.obs.loc[0, "x"] = 100
    with pytest.raises(ValueError, match="outside"):
        plot.plot_pattern(plot_bg=True, show=False)
    # Background validation must not be imposed on a standalone pattern plot.
    plot.plot_pattern(show=False)
    plot.sp.patterns_matrix_dict = {}
    with pytest.raises(ValueError, match="get_pattern_array"):
        plot.plot_pattern(show=False)
    plot.sp.genes_labels = None
    with pytest.raises(ValueError, match="cluster_gene"):
        plot.plot_pattern(show=False)


@pytest.mark.parametrize(
    "rotate,reverse_x,reverse_y", list(product([False, True], repeat=3))
)
def test_image_and_pattern_landmark_stay_aligned(
    plot, tmp_path, monkeypatch, rotate, reverse_x, reverse_y
):
    image = np.zeros((7, 4, 3), dtype=np.float32)
    image[7 - 1 - 1, 0] = [1, 0, 0]
    path = tmp_path / "aligned.png"
    plt.imsave(path, image)
    import STMiner.Plot._pattern as renderer

    original_reader = renderer.mpimg.imread
    reads = []

    def read_once(path):
        reads.append(path)
        return original_reader(path)

    monkeypatch.setattr(renderer.mpimg, "imread", read_once)
    _, axes = plot.plot_pattern(
        image_path=path,
        rotate=rotate,
        reverse_x=reverse_x,
        reverse_y=reverse_y,
        show=False,
    )
    assert len(reads) == 1
    for ax in axes.flat:
        spots = ax.collections[0]
        # The lower-intensity spot marks the same location as the red pixel.
        point = np.asarray(spots.get_offsets()[np.argmin(spots.get_array())], dtype=int)
        np.testing.assert_array_equal(
            ax.images[0].get_array()[point[1], point[0], :3], [1, 0, 0]
        )


def test_export_format_follows_suffix_and_hides_spare_axes(plot, tmp_path):
    plot.sp.genes_labels.loc[3] = [3, "D"]
    plot.sp.patterns_matrix_dict[3] = plot.sp.patterns_matrix_dict[2].copy()
    path = tmp_path / "nested" / "patterns.png"
    fig, axes = plot.plot_pattern(num_cols=2, output_path=path, dpi=100, show=False)
    assert axes.shape == (2, 2)
    assert not axes[1, 1].get_visible()
    assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    fig.canvas.draw()
    for ax in list(axes.flat)[:3]:
        assert ax.get_position().width > 0
    plot.plot_pattern(labels=[0], output_path=tmp_path / "vector.svg", show=False)
    assert "<svg" in (tmp_path / "vector.svg").read_text(encoding="utf-8")


def test_bad_image_paths_fail_explicitly(plot, tmp_path):
    with pytest.raises(ValueError, match="not a directory"):
        plot.plot_pattern(output_path=tmp_path, show=False)
    with pytest.raises(FileNotFoundError, match="Background image"):
        plot.plot_pattern(image_path=tmp_path, show=False)
    with pytest.raises(FileNotFoundError):
        plot.plot_pattern(image_path=tmp_path / "missing.png", show=False)


def test_legacy_positional_arguments_still_work(plot):
    _, axes = plot.plot_pattern(
        "viridis",
        95,
        1,
        True,
        True,
        False,
        False,
        8,
        None,
        False,
        1,
        0.55,
        None,
        True,
        show=False,
        colorbar=False,
    )
    assert axes.shape == (2, 1)
    assert axes[0, 0].get_aspect() == 0.55
    assert axes[0, 0].collections[1].get_sizes()[0] == 8
