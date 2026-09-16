"""Validated spatial pattern rendering; coordinates follow STMiner's [x, y] convention."""

from collections.abc import Mapping
from numbers import Integral, Real
from pathlib import Path

import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import issparse

from .enrichment import _NATURE_STYLE


def _adjust_arr(arr, rotate, reverse_x, reverse_y):
    if rotate:
        arr = np.rot90(arr)
    if reverse_y:
        arr = np.flipud(arr)
    if reverse_x:
        arr = np.fliplr(arr)
    return arr


def _positive_number(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite positive number.")
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive number.")


def _pattern_data(sp, labels):
    table = getattr(sp, "genes_labels", None)
    if not isinstance(table, pd.DataFrame) or "labels" not in table or table.empty:
        raise ValueError("No gene clusters available. Run cluster_gene() first.")
    if table["labels"].isna().any():
        raise ValueError("genes_labels contains missing cluster labels.")
    available = list(pd.unique(table["labels"]))
    if labels is None:
        try:
            selected = sorted(available)
        except TypeError:
            selected = available
    else:
        if isinstance(labels, (str, bytes)):
            raise ValueError("labels must be a non-empty sequence, e.g. labels=[0].")
        try:
            selected = list(dict.fromkeys(labels))
        except TypeError as exc:
            raise ValueError("labels must be a sequence of cluster labels.") from exc
        if not selected:
            raise ValueError("labels must not be empty.")
        unknown = [label for label in selected if label not in available]
        if unknown:
            raise ValueError(f"Unknown pattern labels: {unknown}.")

    matrices = getattr(sp, "patterns_matrix_dict", None)
    if not isinstance(matrices, Mapping) or any(
        label not in matrices for label in selected
    ):
        raise ValueError("Pattern matrices are missing. Run get_pattern_array() first.")
    arrays = []
    for label in selected:
        matrix = matrices[label]
        array = matrix.toarray() if issparse(matrix) else np.asarray(matrix)
        if array.ndim != 2 or 0 in array.shape:
            raise ValueError(f"Pattern {label!r} must be a non-empty 2-D matrix.")
        if array.dtype.kind not in "biuf" or not np.isfinite(array).all():
            raise ValueError(f"Pattern {label!r} must contain finite real values.")
        if (array < 0).any():
            raise ValueError(f"Pattern {label!r} contains negative intensity values.")
        if arrays and array.shape != arrays[0].shape:
            raise ValueError(
                "All selected pattern matrices must have the same spatial shape."
            )
        arrays.append(array)
    counts = table.groupby("labels", sort=False).size()
    return selected, arrays, counts


def _background_points(sp, shape, rotate, reverse_x, reverse_y):
    adata = getattr(sp, "adata", None)
    if adata is None or not {"x", "y"}.issubset(adata.obs.columns):
        raise ValueError("plot_bg=True requires adata.obs['x'] and adata.obs['y'].")
    points = adata.obs[["x", "y"]].to_numpy(copy=True)
    if points.dtype.kind not in "iuf" or not np.isfinite(points).all():
        raise ValueError("Background coordinates must be finite numbers.")
    if (points < 0).any() or not np.equal(points, np.floor(points)).all():
        raise ValueError("Background coordinates must be non-negative grid integers.")
    if (points >= np.asarray(shape)).any():
        raise ValueError("Background coordinates fall outside the pattern matrix.")
    x, y = points.T
    nx, ny = shape
    if rotate:
        x, y = ny - 1 - y, x
        nx, ny = ny, nx
    if reverse_y:
        x = nx - 1 - x
    if reverse_x:
        y = ny - 1 - y
    return x, y


def plot_patterns(
    sp,
    *,
    cmap,
    vmax,
    num_cols,
    rotate,
    reverse_y,
    reverse_x,
    heatmap,
    s,
    image_path,
    rotate_img,
    k,
    aspect,
    output_path,
    plot_bg,
    labels,
    colorbar,
    shared_scale,
    figsize,
    show,
    dpi,
    panel_titles=None,
    colorbar_label="Pattern intensity",
    interpolation="nearest",
    contours=0,
):
    for name, value in {
        "rotate": rotate,
        "reverse_y": reverse_y,
        "reverse_x": reverse_x,
        "heatmap": heatmap,
        "rotate_img": rotate_img,
        "plot_bg": plot_bg,
        "colorbar": colorbar,
        "shared_scale": shared_scale,
        "show": show,
    }.items():
        if not isinstance(value, (bool, np.bool_)):
            raise ValueError(f"{name} must be a boolean.")
    _positive_number(vmax, "vmax")
    if vmax > 100:
        raise ValueError("vmax is a percentile and must be in (0, 100].")
    _positive_number(aspect, "aspect")
    _positive_number(dpi, "dpi")
    if s is not None:
        _positive_number(s, "s")
    if num_cols is not None and (
        isinstance(num_cols, bool)
        or not isinstance(num_cols, Integral)
        or num_cols <= 0
    ):
        raise ValueError("num_cols must be a positive integer or None (automatic).")
    if isinstance(k, bool) or not isinstance(k, Integral):
        raise ValueError("k must be an integer number of image quarter-turns.")
    if figsize is not None:
        if not isinstance(figsize, (tuple, list)) or len(figsize) != 2:
            raise ValueError("figsize must contain width and height in inches.")
        for value in figsize:
            _positive_number(value, "figsize")
    try:
        colormap = mpl.colormaps.get_cmap(
            cmap if cmap is not None else "viridis"
        ).copy()
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Unknown colormap: {cmap!r}.") from exc
    colormap.set_bad((0, 0, 0, 0))

    output = None
    if output_path is not None:
        output = Path(output_path)
        if output.is_dir():
            raise ValueError("output_path must be a file, not a directory.")
        if not output.suffix:
            output = output.with_suffix(".pdf")
        if output.suffix.lower() not in {
            ".pdf",
            ".svg",
            ".png",
            ".eps",
            ".tif",
            ".tiff",
        }:
            raise ValueError("output_path must use PDF, SVG, PNG, EPS, or TIFF.")

    selected, arrays, counts = _pattern_data(sp, labels)
    background = (
        _background_points(sp, arrays[0].shape, rotate, reverse_x, reverse_y)
        if plot_bg
        else None
    )
    arrays = [_adjust_arr(a, rotate, reverse_x, reverse_y) for a in arrays]
    image = None
    if image_path is not None:
        path = Path(image_path)
        if not path.is_file():
            raise FileNotFoundError(f"Background image is not a file: {path}")
        try:
            image = mpimg.imread(path)
        except (OSError, ValueError, SyntaxError) as exc:
            raise ValueError(f"Cannot read background image: {path}") from exc
        if rotate_img:
            image = np.rot90(image, k=k)
        # Convert the top-origin image into [x, y] before the common transform.
        image = np.flip(np.swapaxes(image, 0, 1), axis=1)
        image = np.swapaxes(_adjust_arr(image, rotate, reverse_x, reverse_y), 0, 1)

    positive_values = [a[a > 0] for a in arrays]
    if shared_scale:
        pooled = np.concatenate(positive_values)
        upper = np.percentile(pooled, vmax) if pooled.size else 1
        norms = [mpl.colors.Normalize(0, upper)] * len(arrays)
    else:
        norms = [
            mpl.colors.Normalize(0, np.percentile(v, vmax) if v.size else 1)
            for v in positive_values
        ]
    auto_cols = (
        len(arrays) if len(arrays) <= 3 else min(3, int(np.ceil(np.sqrt(len(arrays)))))
    )
    ncols = min(len(arrays), num_cols or auto_cols)
    nrows = int(np.ceil(len(arrays) / ncols))
    nx, ny = arrays[0].shape
    if figsize is None:
        width = 3.5 if ncols == 1 else 7.2
        panel_height = np.clip((width / ncols) * ny / nx * aspect, 1.6, 4.0)
        title_height = (
            0.6
            if panel_titles is None
            else 0.35
            + 0.15 * max(str(title).count("\n") + 1 for title in panel_titles.values())
        )
        figsize = (width, nrows * (panel_height + title_height))
    extent = (-0.5, nx - 0.5, -0.5, ny - 0.5)

    with mpl.rc_context(_NATURE_STYLE):
        fig, axes = plt.subplots(
            nrows, ncols, figsize=figsize, squeeze=False, layout="constrained"
        )
        used_axes = list(axes.flat)[: len(arrays)]
        artists = []
        for index, (ax, label, array, norm) in enumerate(
            zip(used_axes, selected, arrays, norms)
        ):
            # Point area follows grid density; explicit s retains its old meaning.
            size = (
                s
                if s is not None
                else np.clip(
                    (0.65 * 72 * min(figsize[0] / ncols / nx, figsize[1] / nrows / ny))
                    ** 2,
                    0.15,
                    36,
                )
            )
            if image is not None:
                ax.imshow(image, origin="lower", extent=extent, interpolation="nearest")
            if background is not None:
                ax.scatter(
                    *background, s=size, c="#DEDEDE", edgecolors="none", zorder=1
                )
            if heatmap:
                artist = ax.imshow(
                    np.ma.masked_equal(array.T, 0),
                    origin="lower",
                    extent=extent,
                    interpolation=interpolation,
                    cmap=colormap,
                    norm=norm,
                    zorder=2,
                )
            else:
                x, y = np.nonzero(array)
                artist = ax.scatter(
                    x,
                    y,
                    c=array[x, y],
                    s=size,
                    cmap=colormap,
                    norm=norm,
                    edgecolors="none",
                    zorder=2,
                )
            artists.append(artist)
            if contours and min(array.shape) > 1 and positive_values[index].size:
                surface = np.ma.masked_less_equal(array.T, 0)
                levels = np.linspace(0, surface.max(), contours + 2)[1:-1]
                levels = levels[(levels > surface.min()) & (levels < surface.max())]
                if levels.size:
                    ax.contour(
                        np.arange(nx),
                        np.arange(ny),
                        surface,
                        levels=levels,
                        colors="#365167",
                        linewidths=0.4,
                        alpha=0.45,
                        zorder=3,
                    )
            ax.set(xlim=extent[:2], ylim=extent[2:], aspect=aspect)
            ax.set_axis_off()
            noun = "gene" if counts[label] == 1 else "genes"
            title = (
                f"Pattern {label} · {counts[label]} {noun}"
                if panel_titles is None
                else panel_titles[label]
            )
            ax.set_title(title, loc="left", pad=6)
            if not positive_values[index].size:
                ax.text(
                    0.5,
                    0.5,
                    "No signal",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    color="#666666",
                )
            elif colorbar and not shared_scale:
                bar = fig.colorbar(
                    artist, ax=ax, shrink=0.75, pad=0.025, fraction=0.045
                )
                bar.set_label(colorbar_label, fontsize=7)
                bar.outline.set_visible(False)
                bar.ax.tick_params(labelsize=6, length=2, width=0.5)
        for ax in list(axes.flat)[len(arrays) :]:
            ax.set_visible(False)
        if colorbar and shared_scale and any(v.size for v in positive_values):
            bar = fig.colorbar(
                artists[0], ax=used_axes, shrink=0.8, pad=0.025, fraction=0.035
            )
            bar.set_label(f"{colorbar_label} (shared scale)", fontsize=7)
            bar.outline.set_visible(False)
            bar.ax.tick_params(labelsize=6, length=2, width=0.5)
        if output is not None:
            output.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output, dpi=dpi, facecolor="white")
        if show:
            plt.show()
    return fig, axes
