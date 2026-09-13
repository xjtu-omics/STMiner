from __future__ import annotations

import textwrap
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D


_SOURCE_ORDER = ("GO:BP", "GO:MF", "GO:CC", "KEGG")
_SOURCE_TITLES = {
    "GO:BP": "GO biological process",
    "GO:MF": "GO molecular function",
    "GO:CC": "GO cellular component",
    "KEGG": "KEGG pathway",
}
_NATURE_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.facecolor": "white",
}


def _validate_positive_integer(value: Any, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise ValueError(f"{name} must be a positive integer.")
    return value


def _normalise_sources(sources: Iterable[str] | str | None) -> list[str] | None:
    if sources is None:
        return None
    if isinstance(sources, str):
        sources = [sources]
    try:
        values = [str(source).strip() for source in sources]
    except TypeError as exc:
        raise TypeError("sources must be an iterable of source names.") from exc
    values = list(dict.fromkeys(source for source in values if source))
    if not values:
        raise ValueError("sources must contain at least one source name.")
    return values


def _prepare_enrichment_data(
    result: pd.DataFrame,
    *,
    top_n: int,
    sources: Iterable[str] | str | None,
    p_value_cutoff: float | None,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    if not isinstance(result, pd.DataFrame):
        raise TypeError("result must be a pandas DataFrame.")
    required = {"source", "native", "name", "p_value", "intersection_size"}
    missing = sorted(required.difference(result.columns))
    if missing:
        raise ValueError("result is missing required columns: " + ", ".join(missing))
    if result.empty:
        raise ValueError("result contains no enrichment terms to plot.")

    data = result.copy(deep=True)
    data["source"] = data["source"].astype(str)
    data["native"] = data["native"].astype(str)
    data["name"] = data["name"].astype(str)
    data["query"] = (
        data["query"].fillna("Query").astype(str)
        if "query" in data.columns
        else "Query"
    )

    selected_sources = _normalise_sources(sources)
    if selected_sources is not None:
        data = data[data["source"].isin(selected_sources)].copy()
        if data.empty:
            raise ValueError("no enrichment terms match the requested sources.")

    data["p_value"] = pd.to_numeric(data["p_value"], errors="coerce")
    data["intersection_size"] = pd.to_numeric(
        data["intersection_size"], errors="coerce"
    )
    if data[["p_value", "intersection_size"]].isna().any().any():
        raise ValueError("p_value and intersection_size must contain numeric values.")
    if (~np.isfinite(data["p_value"])).any() or not data["p_value"].between(0, 1).all():
        raise ValueError("p_value values must be finite and between 0 and 1.")
    if (~np.isfinite(data["intersection_size"])).any() or (
        data["intersection_size"] <= 0
    ).any():
        raise ValueError("intersection_size values must be finite and positive.")

    if "precision" in data.columns:
        data["gene_ratio"] = pd.to_numeric(data["precision"], errors="coerce")
    elif "query_size" in data.columns:
        query_size = pd.to_numeric(data["query_size"], errors="coerce")
        if (
            query_size.isna().any()
            or (~np.isfinite(query_size)).any()
            or (query_size <= 0).any()
        ):
            raise ValueError("query_size values must be finite and positive.")
        data["gene_ratio"] = data["intersection_size"] / query_size
    else:
        raise ValueError("result must contain either precision or query_size.")
    if data["gene_ratio"].isna().any() or (~np.isfinite(data["gene_ratio"])).any():
        raise ValueError("precision values must be finite.")
    if not data["gene_ratio"].between(0, 1).all():
        raise ValueError("precision values must be between 0 and 1.")

    if p_value_cutoff is not None:
        if (
            not isinstance(p_value_cutoff, (int, float))
            or isinstance(p_value_cutoff, bool)
            or not 0 < p_value_cutoff <= 1
        ):
            raise ValueError("p_value_cutoff must be greater than 0 and at most 1.")
        data = data[data["p_value"] <= p_value_cutoff].copy()
    if data.empty:
        raise ValueError("no enrichment terms remain after filtering.")

    positive_floor = np.nextafter(0.0, 1.0)
    data["minus_log10_p"] = -np.log10(data["p_value"].clip(lower=positive_floor))
    data = data.sort_values(
        ["source", "query", "p_value", "intersection_size", "native"],
        ascending=[True, True, True, False, True],
        kind="mergesort",
    )
    data = data.groupby(["source", "query"], sort=False).head(top_n).copy()

    available_sources = list(dict.fromkeys(data["source"]))
    if selected_sources is None:
        source_order = [
            source for source in _SOURCE_ORDER if source in available_sources
        ]
        source_order.extend(sorted(set(available_sources).difference(source_order)))
    else:
        source_order = [
            source for source in selected_sources if source in available_sources
        ]
    query_order = list(dict.fromkeys(data["query"]))
    return data, source_order, query_order


def _marker_areas(values: pd.Series) -> np.ndarray:
    minimum = float(values.min())
    maximum = float(values.max())
    if minimum == maximum:
        return np.full(len(values), 55.0)
    return 24.0 + 92.0 * (values.to_numpy(dtype=float) - minimum) / (maximum - minimum)


def _size_legend_handles(values: pd.Series) -> list[Line2D]:
    candidates = np.quantile(values.to_numpy(dtype=float), [0, 0.5, 1])
    labels = sorted(set(int(round(value)) for value in candidates))
    reference = pd.Series(labels, dtype=float)
    areas = _marker_areas(reference)
    return [
        Line2D(
            [],
            [],
            linestyle="none",
            marker="o",
            markersize=float(np.sqrt(area)),
            markerfacecolor="#8C8C8C",
            markeredgecolor="#262626",
            markeredgewidth=0.35,
            label=str(label),
        )
        for label, area in zip(labels, areas)
    ]


def _wrap_term(term: str, width: int) -> str:
    return "\n".join(
        textwrap.wrap(term, width=width, break_long_words=False, break_on_hyphens=False)
    )


def _save_figure(fig: Figure, save_path: str | Path, dpi: int) -> None:
    path = Path(save_path)
    supported = {".pdf", ".svg", ".png", ".tif", ".tiff", ".eps"}
    if path.suffix.lower() not in supported:
        raise ValueError("save_path must use PDF, SVG, PNG, TIFF, or EPS format.")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")


def plot_enrichment(
    result: pd.DataFrame,
    *,
    top_n: int = 5,
    sources: Iterable[str] | str | None = None,
    p_value_cutoff: float | None = None,
    figsize: tuple[float, float] | None = None,
    cmap: str = "viridis",
    title: str | None = None,
    save_path: str | Path | None = None,
    dpi: int = 600,
    show: bool = False,
) -> tuple[Figure, list[mpl.axes.Axes]]:
    """Plot g:Profiler enrichment results as a Nature-style bubble plot."""
    top_n = _validate_positive_integer(top_n, "top_n")
    dpi = _validate_positive_integer(dpi, "dpi")
    if not isinstance(show, bool):
        raise TypeError("show must be a boolean.")
    if title is not None and not isinstance(title, str):
        raise TypeError("title must be a string or None.")
    try:
        colormap = mpl.colormaps[cmap]
    except (KeyError, TypeError) as exc:
        raise ValueError(f"Unknown matplotlib colormap: {cmap!r}.") from exc

    data, source_order, query_order = _prepare_enrichment_data(
        result,
        top_n=top_n,
        sources=sources,
        p_value_cutoff=p_value_cutoff,
    )
    panel_terms = {
        source: data.loc[data["source"] == source, ["native", "name"]]
        .drop_duplicates()
        .shape[0]
        for source in source_order
    }
    if figsize is None:
        width = 7.2
        height = max(
            2.6,
            0.24 * sum(panel_terms.values()) + 0.7 * len(source_order) + 0.55,
        )
        figsize = (width, height)
    elif (
        not isinstance(figsize, tuple)
        or len(figsize) != 2
        or any(not isinstance(value, (int, float)) or value <= 0 for value in figsize)
    ):
        raise ValueError("figsize must be a tuple of two positive numbers.")

    color_min = float(data["minus_log10_p"].min())
    color_max = float(data["minus_log10_p"].max())
    if color_min == color_max:
        color_min = max(0.0, color_min - 0.5)
        color_max += 0.5
    normalization = mpl.colors.Normalize(vmin=color_min, vmax=color_max)
    multiple_queries = len(query_order) > 1

    with mpl.rc_context(_NATURE_STYLE):
        fig = plt.figure(figsize=figsize, facecolor="white")
        bottom = 0.15
        top = 0.89 if title else 0.96
        grid = fig.add_gridspec(
            len(source_order),
            2,
            width_ratios=[1, 0.035],
            height_ratios=[panel_terms[source] + 1.5 for source in source_order],
            left=0.34,
            right=0.93,
            bottom=bottom,
            top=top,
            hspace=0.38 if multiple_queries else 0.62,
            wspace=0.08,
        )
        axes: list[mpl.axes.Axes] = []
        all_sizes = _marker_areas(data["intersection_size"])
        data = data.assign(marker_area=all_sizes)

        for panel_index, source in enumerate(source_order):
            ax = fig.add_subplot(grid[panel_index, 0])
            axes.append(ax)
            panel = data[data["source"] == source].copy()
            term_order = (
                panel.groupby(["native", "name"], sort=False)["p_value"]
                .min()
                .sort_values(ascending=False, kind="mergesort")
                .index.tolist()
            )
            term_positions = {term: index for index, term in enumerate(term_order)}
            panel["term_key"] = list(zip(panel["native"], panel["name"]))
            panel["y_position"] = panel["term_key"].map(term_positions)

            if multiple_queries:
                query_positions = {
                    query: index for index, query in enumerate(query_order)
                }
                x_values = panel["query"].map(query_positions)
            else:
                x_values = panel["gene_ratio"]

            ax.scatter(
                x_values,
                panel["y_position"],
                s=panel["marker_area"],
                c=panel["minus_log10_p"],
                cmap=colormap,
                norm=normalization,
                edgecolors="#202020",
                linewidths=0.35,
                alpha=0.92,
                zorder=3,
            )
            ax.set_yticks(range(len(term_order)))
            ax.set_yticklabels([_wrap_term(term[1], 45) for term in term_order])
            ax.set_ylim(-0.65, len(term_order) - 0.35)
            if multiple_queries:
                ax.set_xticks(range(len(query_order)))
                ax.set_xlim(-0.55, len(query_order) - 0.45)
                if panel_index == len(source_order) - 1:
                    ax.set_xticklabels(query_order, rotation=30, ha="right")
                    ax.set_xlabel("Gene cluster")
                else:
                    ax.tick_params(axis="x", labelbottom=False)
            else:
                maximum_ratio = float(panel["gene_ratio"].max())
                ax.set_xlim(0, min(1.0, maximum_ratio * 1.15 + 0.01))
                ax.xaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1.0))
                ax.set_xlabel("Gene ratio")
            ax.set_title(
                _SOURCE_TITLES.get(source, source),
                loc="left",
                fontweight="bold",
                pad=5,
            )
            ax.grid(axis="x", color="#D9D9D9", linewidth=0.45, alpha=0.75)
            ax.set_axisbelow(True)
            ax.tick_params(axis="y", length=0, pad=3)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)

        color_axis = fig.add_subplot(grid[:, 1])
        scalar_mappable = mpl.cm.ScalarMappable(norm=normalization, cmap=colormap)
        colorbar = fig.colorbar(scalar_mappable, cax=color_axis)
        colorbar.set_label(r"$-\log_{10}$(adjusted $P$)", labelpad=5)
        colorbar.outline.set_linewidth(0.5)

        handles = _size_legend_handles(data["intersection_size"])
        fig.legend(
            handles=handles,
            title="Intersecting genes",
            loc="lower right",
            bbox_to_anchor=(0.93, 0.01),
            frameon=False,
            ncol=len(handles),
            handletextpad=0.3,
            columnspacing=1.0,
        )
        if title:
            fig.suptitle(
                title,
                x=0.34,
                y=0.975,
                ha="left",
                fontsize=9,
                fontweight="bold",
            )
        if save_path is not None:
            _save_figure(fig, save_path, dpi)
        if show:
            plt.show()
    return fig, axes
