import warnings

import numpy as np
import pandas as pd
import anndata as ad

from scipy.sparse import issparse


class STMinerIOError(ValueError):
    """IO layer error with actionable context for users."""


def format_io_error(
    message: str,
    *,
    path: str | None = None,
    operation: str | None = None,
    expected: str | None = None,
    actual: str | None = None,
    hint: str | None = None,
) -> str:
    details = [message]
    if operation:
        details.append(f"operation={operation}")
    if path:
        details.append(f"path={path}")
    if expected:
        details.append(f"expected={expected}")
    if actual:
        details.append(f"actual={actual}")
    if hint:
        details.append(f"hint={hint}")
    return " | ".join(details)


def require_positive_int(value, name: str, *, path: str | None = None) -> int:
    if isinstance(value, bool):
        raise STMinerIOError(
            format_io_error(
                f"Invalid parameter '{name}'",
                operation="validate parameter",
                path=path,
                expected="positive integer",
                actual=f"{value} ({type(value).__name__})",
            )
        )

    if isinstance(value, (float, np.floating)):
        if not np.isfinite(value):
            raise STMinerIOError(
                format_io_error(
                    f"Invalid parameter '{name}'",
                    operation="validate parameter",
                    path=path,
                    expected="finite positive integer",
                    actual=str(value),
                )
            )
        converted_value = int(value)
        warnings.warn(
            format_io_error(
                f"Float parameter '{name}' converted to int",
                operation="validate parameter",
                path=path,
                expected="integer input",
                actual=f"{value} -> {converted_value}",
                hint="Pass int explicitly to avoid implicit conversion",
            ),
            UserWarning,
            stacklevel=2,
        )
        value = converted_value
    elif not isinstance(value, (int, np.integer)):
        raise STMinerIOError(
            format_io_error(
                f"Invalid parameter '{name}'",
                operation="validate parameter",
                path=path,
                expected="> 0",
                actual=f"{value} ({type(value).__name__})",
            )
        )

    value = int(value)
    if value <= 0:
        raise STMinerIOError(
            format_io_error(
                f"Invalid parameter '{name}'",
                operation="validate parameter",
                path=path,
                expected="> 0",
                actual=str(value),
            )
        )
    return value


def validate_spatial_array(spatial, *, path: str | None = None, source: str = "adata.obsm['spatial']") -> np.ndarray:
    arr = np.asarray(spatial)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise STMinerIOError(
            format_io_error(
                "Invalid spatial coordinates",
                operation="validate spatial",
                path=path,
                expected="2D array with at least 2 columns",
                actual=f"shape={arr.shape}",
                hint=f"Ensure {source} exists and stores x/y coordinates",
            )
        )
    if not np.issubdtype(arr.dtype, np.number):
        raise STMinerIOError(
            format_io_error(
                "Spatial coordinates must be numeric",
                operation="validate spatial",
                path=path,
                expected="numeric array",
                actual=str(arr.dtype),
                hint=f"Convert {source} to numeric values before loading",
            )
        )
    if np.isnan(arr[:, :2]).any() or np.isinf(arr[:, :2]).any():
        raise STMinerIOError(
            format_io_error(
                "Spatial coordinates contain NaN/Inf",
                operation="validate spatial",
                path=path,
                expected="finite numeric x/y values",
                actual="contains NaN or Inf",
            )
        )
    return arr


def merge_bin_coordinate(coordinate: np.ndarray, coordinate_min: int, bin_size: int):
    bin_size = require_positive_int(bin_size, "bin_size")
    return np.floor((coordinate - coordinate_min) / bin_size).astype(int)


def get_bin_center(bin_coordinate: np.ndarray, coordinate_min: int, bin_size: int):
    bin_size = require_positive_int(bin_size, "bin_size")
    return bin_coordinate * bin_size + coordinate_min + int(bin_size / 2)


def bin_spatial_adata(adata, bin_size=50, method="sum"):
    bin_size = require_positive_int(bin_size, "bin_size")
    if method not in {"sum", "mean"}:
        raise STMinerIOError(
            format_io_error(
                "Invalid aggregation method",
                operation="bin_spatial_adata",
                expected="'sum' or 'mean'",
                actual=repr(method),
            )
        )
    if "spatial" not in adata.obsm:
        raise STMinerIOError(
            format_io_error(
                "Missing spatial coordinates",
                operation="bin_spatial_adata",
                expected="adata.obsm['spatial']",
                actual=f"obsm_keys={list(adata.obsm.keys())}",
                hint="Load data with spatial coordinates before merge_bin=True",
            )
        )

    coords = validate_spatial_array(adata.obsm["spatial"], source="adata.obsm['spatial']")
    if adata.n_obs == 0:
        raise STMinerIOError(
            format_io_error(
                "Cannot bin empty AnnData",
                operation="bin_spatial_adata",
                expected="adata.n_obs > 0",
                actual=str(adata.n_obs),
            )
        )

    grid_x = (coords[:, 0] // bin_size).astype(int)
    grid_y = (coords[:, 1] // bin_size).astype(int)

    grid_df = pd.DataFrame({"grid_x": grid_x, "grid_y": grid_y, "index": np.arange(len(adata))})
    grouped = grid_df.groupby(["grid_x", "grid_y"])["index"].apply(list)

    new_X = []
    new_coords = []
    new_x = []
    new_y = []
    new_obs_names = []
    for (gx, gy), indices in grouped.items():
        sub_X = adata.X[indices]
        if issparse(sub_X):
            sub_X = sub_X.toarray()
        if method == "sum":
            aggregated_X = np.sum(sub_X, axis=0)
        else:
            aggregated_X = np.mean(sub_X, axis=0)

        new_X.append(aggregated_X.reshape(1, -1))
        new_coords.append([(gx + 0.5) * bin_size, (gy + 0.5) * bin_size])
        new_x.append(gx)
        new_y.append(gy)
        new_obs_names.append(f"bin_{gx}_{gy}")

    if len(new_X) == 0:
        raise STMinerIOError(
            format_io_error(
                "No valid bins generated from spatial coordinates",
                operation="bin_spatial_adata",
                expected="at least one non-empty spatial bin",
                actual="0 bins",
            )
        )

    new_X = np.vstack(new_X)
    new_coords = np.array(new_coords, dtype=int)

    new_adata = ad.AnnData(X=new_X, var=adata.var.copy(), uns=adata.uns.copy())
    new_adata.obsm["spatial"] = new_coords
    new_adata.obs_names = new_obs_names
    new_adata.obs['x'] = np.array(new_x, dtype=int)
    new_adata.obs['y'] = np.array(new_y, dtype=int)

    return new_adata
