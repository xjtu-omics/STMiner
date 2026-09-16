"""Render pattern-specific enrichment hits in their observed spatial context."""

import textwrap
from types import SimpleNamespace
from numbers import Integral

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.ndimage import distance_transform_edt, gaussian_filter

from STMiner.services.pathways import score_pathways
from ._pattern import plot_patterns, _positive_number


def plot_pathways(
    sp,
    enrichment,
    *,
    labels=None,
    sources=None,
    term_ids=None,
    top_n=2,
    p_value_cutoff=0.05,
    min_genes=2,
    layer=None,
    gene_map=None,
    pattern_map=None,
    display="density",
    bandwidth=1.5,
    contours=3,
    tissue_mask=None,
    support_radius=1.5,
    **plot_kwargs,
):
    """Return (fig, axes, spot_scores, terms) without changing the SPFinder state.

    Scores are mean expression of the pattern's mapped enrichment-hit genes at
    each observation. Repeated grid positions are averaged for display only.
    The default display is a Gaussian-smoothed, expression-weighted spatial
    density, renormalized to sum to one within its support. Bandwidth is in grid
    units. By default support extends at most support_radius from observed cells;
    pass a registered boolean [x, y] tissue_mask to define it explicitly. Contours
    are fractions of peak density, not confidence regions. This is a descriptive
    visualization, not an inferred pathway-activation probability.
    Returned spot_scores always retain unsmoothed mean expression. Use
    display='expression' for the original spot/heatmap rendering.
    Shared scaling is the default; use shared_scale=False for per-panel scaling.
    Rendering accepts the existing plot_pattern options, including background
    image alignment, rotations, heatmap, output_path, and show.
    """
    if not isinstance(display, str) or display not in {"density", "expression"}:
        raise ValueError("display must be 'density' or 'expression'.")
    _positive_number(bandwidth, "bandwidth")
    _positive_number(support_radius, "support_radius")
    if (
        isinstance(contours, bool)
        or not isinstance(contours, Integral)
        or not 0 <= contours <= 10
    ):
        raise ValueError("contours must be an integer from 0 to 10.")
    if display == "density" and plot_kwargs.get("heatmap") is False:
        raise ValueError(
            "Density display is continuous; use display='expression' for scatter."
        )
    options = dict(
        cmap="YlGnBu" if display == "density" else None,
        vmax=99,
        num_cols=2,
        rotate=False,
        reverse_y=False,
        reverse_x=False,
        heatmap=display == "density",
        s=None,
        image_path=None,
        rotate_img=False,
        k=1,
        aspect=1,
        output_path=None,
        plot_bg=display != "density",
        colorbar=True,
        shared_scale=True,
        figsize=None,
        show=True,
        dpi=600,
    )
    unknown = set(plot_kwargs).difference(options)
    if unknown:
        raise TypeError(f"Unknown pathway plotting options: {sorted(unknown)}.")
    options.update(plot_kwargs)
    adata = getattr(sp, "adata", None)
    if adata is None or not {"x", "y"}.issubset(adata.obs.columns) or adata.n_obs == 0:
        raise ValueError(
            "Pathway plots require non-empty adata.obs x/y grid coordinates."
        )
    coordinates = adata.obs[["x", "y"]].to_numpy()
    if coordinates.dtype.kind not in "iuf" or not np.isfinite(coordinates).all():
        raise ValueError("Spatial coordinates must be finite grid numbers.")
    if (coordinates < 0).any() or not np.equal(
        coordinates, np.floor(coordinates)
    ).all():
        raise ValueError("Spatial coordinates must be non-negative grid integers.")
    shape = tuple(int(value) + 1 for value in coordinates.max(axis=0))
    if shape[0] * shape[1] > 25_000_000:
        raise ValueError(
            "The spatial grid exceeds 25 million cells; bin coordinates before plotting."
        )
    if tissue_mask is not None:
        tissue_mask = np.asarray(tissue_mask)
        if (
            tissue_mask.shape != shape
            or tissue_mask.dtype != bool
            or not tissue_mask.any()
        ):
            raise ValueError(
                "tissue_mask must be a non-empty boolean array matching the [x, y] grid."
            )
    scores, terms = score_pathways(
        sp,
        enrichment,
        labels=labels,
        sources=sources,
        term_ids=term_ids,
        top_n=top_n,
        p_value_cutoff=p_value_cutoff,
        min_genes=min_genes,
        layer=layer,
        gene_map=gene_map,
        pattern_map=pattern_map,
    )
    x, y = coordinates.T.astype(np.intp)
    cell_counts = csr_matrix(
        (np.ones(len(x), dtype=np.intp), (x, y)), shape=shape
    ).toarray()
    support = (
        tissue_mask
        if tissue_mask is not None
        else distance_transform_edt(cell_counts == 0) <= support_radius
    )
    matrices, membership, titles = {}, [], {}
    for panel, term in terms.iterrows():
        values = scores[panel].to_numpy()
        grid = csr_matrix((values, (x, y)), shape=shape).toarray()
        np.divide(grid, cell_counts, out=grid, where=cell_counts > 0)
        if display == "density":
            grid[~support] = 0
            grid = gaussian_filter(grid, sigma=bandwidth, mode="constant", truncate=3.0)
            grid[~support] = 0
            total = grid.sum()
            if total > 0:
                grid /= total
        matrices[panel] = grid
        membership.extend((panel, gene) for gene in term["genes"])
        name = textwrap.fill(term["name"], width=34, break_long_words=False)
        titles[panel] = (
            f"Pattern {term['pattern']} · {term['native']}\n{name}\n"
            f"adj. P = {term['p_value']:.2g} · genes {term['n_genes']}/{term['n_hits']}"
        )
    state = SimpleNamespace(
        adata=adata,
        genes_labels=pd.DataFrame(membership, columns=["labels", "gene_id"]),
        patterns_matrix_dict=matrices,
    )
    fig, axes = plot_patterns(
        state,
        labels=list(terms.index),
        panel_titles=titles,
        colorbar_label=(
            "Spatial density (sum = 1)"
            if display == "density"
            else "Mean hit-gene expression"
        ),
        interpolation="bilinear" if display == "density" else "nearest",
        contours=contours if display == "density" else 0,
        **options,
    )
    terms.attrs["display"] = display
    if display == "density":
        terms.attrs["bandwidth"] = bandwidth
        terms.attrs["contours"] = contours
        terms.attrs["support"] = (
            "supplied tissue_mask"
            if tissue_mask is not None
            else f"within {support_radius} grid units of observed positions"
        )
        terms.attrs["density_definition"] = (
            "Gaussian-smoothed mean expression, masked and normalized to unit grid sum per term (unless no signal)."
        )
    return fig, axes, scores, terms
