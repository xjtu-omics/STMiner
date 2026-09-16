"""Pattern enrichment and auditable spatial expression summaries of term hits."""

from collections.abc import Mapping
from numbers import Integral, Real
import warnings

import numpy as np
import pandas as pd
from scipy.sparse import issparse


def pattern_gene_sets(sp, labels=None):
    adata = getattr(sp, "adata", None)
    if adata is None or adata.n_obs == 0 or not adata.var_names.is_unique:
        raise ValueError("A non-empty AnnData with unique var_names is required.")
    table = getattr(sp, "genes_labels", None)
    if (
        not isinstance(table, pd.DataFrame)
        or not {"labels", "gene_id"}.issubset(table.columns)
        or table.empty
    ):
        raise ValueError(
            "Run cluster_gene() first; genes_labels needs labels and gene_id."
        )
    if table[["labels", "gene_id"]].isna().any().any():
        raise ValueError("genes_labels contains missing labels or gene identifiers.")
    if not table["gene_id"].map(lambda x: isinstance(x, str)).all():
        raise ValueError("gene_id must contain identifiers matching adata.var_names.")
    if (table.groupby("gene_id")["labels"].nunique() > 1).any():
        raise ValueError("A gene is assigned to more than one pattern.")
    available = list(pd.unique(table["labels"]))
    if labels is None:
        try:
            labels = sorted(available)
        except TypeError:
            labels = available
    else:
        if isinstance(labels, (str, bytes)):
            raise ValueError("labels must be a sequence, e.g. [0, 1].")
        try:
            labels = list(dict.fromkeys(labels))
        except TypeError as exc:
            raise ValueError("labels must be a sequence of pattern labels.") from exc
        if not labels or any(label not in available for label in labels):
            raise ValueError("labels is empty or contains an unknown pattern.")
    groups = {}
    for label in labels:
        genes = list(dict.fromkeys(table.loc[table["labels"] == label, "gene_id"]))
        missing = [gene for gene in genes if gene not in adata.var_names]
        if missing:
            raise ValueError(
                f"Pattern {label!r} has genes absent from adata.var_names: {missing[:5]}."
            )
        groups[label] = genes
    return groups


def enrich_patterns(sp, *, organism, labels=None, background=None, **kwargs):
    """Delegate named pattern gene sets to the existing g:Profiler wrapper."""
    from STMiner.KEGGFinder import KEGGFinder

    groups = pattern_gene_sets(sp, labels)
    query_to_pattern = {f"pattern_{i}": label for i, label in enumerate(groups)}
    queries = {query: groups[label] for query, label in query_to_pattern.items()}
    if background is None:
        background = list(sp.adata.var_names)
    result = (
        KEGGFinder()
        .enrich_gene_set(queries, organism=organism, background=background, **kwargs)
        .copy(deep=True)
    )
    if "query" not in result or not result["query"].isin(query_to_pattern).all():
        raise ValueError(
            "The enrichment response contains an unknown query identifier."
        )
    result.insert(0, "pattern", result["query"].map(query_to_pattern))
    result.attrs["pattern_gene_sets"] = groups
    result.attrs["organism"] = organism
    result.attrs["correction_method"] = kwargs.get("correction_method", "g_SCS")
    return result


def _string_list(value, name):
    if isinstance(value, str):
        return [value]
    try:
        result = list(dict.fromkeys(value))
    except TypeError as exc:
        raise ValueError(f"{name} must be a sequence of identifiers.") from exc
    if not result or any(
        not isinstance(item, str) or not item.strip() for item in result
    ):
        raise ValueError(f"{name} must contain non-empty string identifiers.")
    return result


def score_pathways(
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
):
    """Return observation-level mean hit expression and one metadata row per panel.

    Only genes in the term intersection, the current pattern, and var_names are
    used. No normalization, imputation, or integer conversion is applied to X.
    Scores use the current X or requested layer; normalized nonnegative expression
    is recommended. The selected list can include GO terms as well as pathways.
    """
    groups = pattern_gene_sets(sp, labels)
    for name, value in [("top_n", top_n), ("min_genes", min_genes)]:
        if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
            raise ValueError(f"{name} must be a positive integer.")
    if (
        isinstance(p_value_cutoff, bool)
        or not isinstance(p_value_cutoff, Real)
        or not 0 < p_value_cutoff <= 1
    ):
        raise ValueError("p_value_cutoff must be in (0, 1].")
    required = {"source", "native", "name", "p_value", "intersections"}
    if not isinstance(enrichment, pd.DataFrame) or not required.issubset(
        enrichment.columns
    ):
        raise ValueError(
            "enrichment requires source, native, name, p_value, and intersections columns."
        )
    data = enrichment.copy(deep=True).reset_index(drop=True)
    if "pattern" not in data:
        if not isinstance(pattern_map, Mapping) or "query" not in data:
            raise ValueError(
                "Use sp.enrich_patterns() or provide pattern_map={query: pattern_label}."
            )
        data["pattern"] = data["query"].map(pattern_map)
    elif pattern_map is not None:
        raise ValueError("pattern_map is only used when the pattern column is absent.")
    if data["pattern"].isna().any():
        raise ValueError("Some enrichment queries have no pattern mapping.")
    all_groups = pattern_gene_sets(sp)
    if not data["pattern"].isin(all_groups).all():
        raise ValueError(
            "The enrichment result refers to patterns not in genes_labels."
        )
    snapshot = enrichment.attrs.get("pattern_gene_sets")
    if snapshot is not None:
        for label, genes in groups.items():
            if label in snapshot and set(genes) != set(snapshot[label]):
                raise ValueError(
                    "Pattern membership changed after enrichment; rerun enrich_patterns()."
                )
    data = data[data["pattern"].isin(groups)].copy()
    if sources is not None:
        data = data[data["source"].isin(_string_list(sources, "sources"))].copy()
    if term_ids is not None:
        selected_ids = _string_list(term_ids, "term_ids")
        unknown = set(selected_ids).difference(data["native"])
        if unknown:
            raise ValueError(
                f"Unknown term_ids in the selected patterns/sources: {sorted(unknown)}."
            )
        data = data[data["native"].isin(selected_ids)].copy()
    for column in ["source", "native", "name"]:
        if not data[column].map(lambda x: isinstance(x, str) and bool(x.strip())).all():
            raise ValueError(f"{column} must contain non-empty strings.")
    p_values = pd.to_numeric(data["p_value"], errors="coerce")
    if not np.isfinite(p_values).all() or not p_values.between(0, 1).all():
        raise ValueError("p_value must contain finite adjusted P values in [0, 1].")
    data = data[p_values <= p_value_cutoff].copy()
    data["p_value"] = p_values.loc[data.index]
    if data.empty:
        raise ValueError(
            "No enriched terms pass the requested pattern/source/P-value filters."
        )
    if data.duplicated(["pattern", "source", "native"]).any():
        raise ValueError(
            "Duplicate pattern/source/term rows; use one enrichment result per query."
        )
    if gene_map is not None and not isinstance(gene_map, Mapping):
        raise ValueError("gene_map must map enrichment identifiers to adata.var_names.")
    gene_map = {} if gene_map is None else gene_map
    if any(
        not isinstance(k, str) or not isinstance(v, str) for k, v in gene_map.items()
    ):
        raise ValueError(
            "gene_map keys and values must be strings (one target per identifier)."
        )
    adata = sp.adata
    if layer is not None and layer not in adata.layers:
        raise ValueError(f"Expression layer {layer!r} does not exist.")
    expression = adata.X if layer is None else adata.layers[layer]
    if expression is None:
        raise ValueError("The selected expression matrix is empty.")

    records, vectors, skipped = [], {}, []
    for pattern, members in groups.items():
        rows = data[data["pattern"] == pattern].sort_values(
            ["p_value", "source", "native"], kind="mergesort"
        )
        accepted = 0
        for _, row in rows.iterrows():
            hits = row["intersections"]
            if isinstance(hits, (str, bytes)) or not isinstance(
                hits, (list, tuple, np.ndarray)
            ):
                raise ValueError(
                    "intersections must contain gene lists, not serialized strings or counts."
                )
            if any(not isinstance(gene, str) or not gene for gene in hits):
                raise ValueError(
                    "Each intersection must contain string gene identifiers."
                )
            hits = list(dict.fromkeys(hits))
            used, missing = [], []
            for hit in hits:
                target = gene_map.get(hit, hit)
                if target in members and target in adata.var_names:
                    if target not in used:
                        used.append(target)
                else:
                    missing.append(hit)
            if len(used) < min_genes:
                skipped.append((pattern, row["native"], len(used)))
                continue
            subset = expression[:, adata.var_names.get_indexer(used)]
            values = subset.data if issparse(subset) else np.asarray(subset)
            if (
                values.dtype.kind not in "biuf"
                or not np.isfinite(values).all()
                or (values < 0).any()
            ):
                raise ValueError(
                    "Pathway maps require finite nonnegative expression; choose a normalized nonnegative layer."
                )
            score = np.asarray(subset.mean(axis=1)).ravel()
            if not np.isfinite(score).all():
                raise ValueError(
                    "Mean expression is not finite; check expression magnitudes."
                )
            panel = f"pathway_{len(records)}"
            vectors[panel] = score
            records.append(
                {
                    "panel": panel,
                    "pattern": pattern,
                    "source": row["source"],
                    "native": row["native"],
                    "name": row["name"],
                    "p_value": row["p_value"],
                    "n_hits": len(hits),
                    "n_genes": len(used),
                    "genes": tuple(used),
                    "unmapped_genes": tuple(missing),
                    "layer": "X" if layer is None else layer,
                }
            )
            accepted += 1
            if accepted == top_n:
                break
    if not records:
        raise ValueError(
            "No terms have enough mapped pattern genes. Check gene_map or min_genes."
        )
    terms = pd.DataFrame.from_records(records).set_index("panel")
    terms.attrs["skipped_terms"] = skipped
    terms.attrs["unrepresented_patterns"] = [
        label for label in groups if label not in set(terms["pattern"])
    ]
    terms.attrs["score_definition"] = (
        "Arithmetic mean expression of mapped intersection genes per observation."
    )
    if (
        skipped
        or terms["unmapped_genes"].map(bool).any()
        or terms.attrs["unrepresented_patterns"]
    ):
        warnings.warn(
            "Some hits were unmapped, below min_genes, or some patterns had no qualifying terms; "
            "inspect unmapped_genes and terms.attrs['skipped_terms'/'unrepresented_patterns'].",
            UserWarning,
            stacklevel=2,
        )
    scores = pd.DataFrame(vectors, index=adata.obs_names.copy())
    return scores, terms
