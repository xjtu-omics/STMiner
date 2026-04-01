import pandas as pd
import scanpy as sc
import os

from STMiner.IO.IOUtil import (
    STMinerIOError,
    format_io_error,
    merge_bin_coordinate,
    require_positive_int,
    validate_spatial_array,
)


def read_bmk(file_path, bin_size=1):
    path = str(file_path)
    require_positive_int(bin_size, "bin_size", path=path)

    if not os.path.isdir(file_path):
        raise NotADirectoryError(
            format_io_error(
                "BMK input path is not a directory",
                operation="read_bmk",
                path=path,
                expected="directory containing 10x mtx files",
                hint="Pass the directory path, not a single file",
            )
        )

    matrix_exists = any(os.path.exists(os.path.join(path, f)) for f in ["matrix.mtx", "matrix.mtx.gz"])
    barcode_exists = any(os.path.exists(os.path.join(path, f)) for f in ["barcodes.tsv", "barcodes.tsv.gz"])
    feature_exists = any(os.path.exists(os.path.join(path, f)) for f in ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"])
    if not (matrix_exists and barcode_exists and feature_exists):
        raise FileNotFoundError(
            format_io_error(
                "Missing 10x matrix files in BMK directory",
                operation="read_bmk",
                path=path,
                expected=(
                    "matrix.mtx(.gz), barcodes.tsv(.gz), and features.tsv(.gz) or genes.tsv(.gz)"
                ),
                actual=f"found={sorted(os.listdir(path))[:20]}",
            )
        )

    barcode_pos_path = os.path.join(path, "barcodes_pos.tsv.gz")
    if not os.path.isfile(barcode_pos_path):
        raise FileNotFoundError(
            format_io_error(
                "Required file barcodes_pos.tsv.gz not found",
                operation="read_bmk",
                path=barcode_pos_path,
                expected="gzip TSV with columns: barcode, array_col, array_row",
            )
        )

    try:
        adata = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)
    except Exception as e:
        raise STMinerIOError(
            format_io_error(
                "Failed to read 10x matrix directory",
                operation="scanpy.read_10x_mtx",
                path=path,
                actual=f"{type(e).__name__}: {e}",
                hint="Verify matrix/features/barcodes files are complete and consistent",
            )
        ) from e

    try:
        barcode_pos = pd.read_csv(
            barcode_pos_path,
            compression='gzip',
            sep='\t',
            names=["barcode", "array_col", "array_row"],
            header=None,
        )
    except Exception as e:
        raise STMinerIOError(
            format_io_error(
                "Failed to parse barcodes_pos.tsv.gz",
                operation="pandas.read_csv",
                path=barcode_pos_path,
                actual=f"{type(e).__name__}: {e}",
                hint="Ensure the file is tab-separated gzip with 3 columns",
            )
        ) from e

    if barcode_pos.empty:
        raise STMinerIOError(
            format_io_error(
                "barcodes_pos.tsv.gz is empty",
                operation="validate barcode positions",
                path=barcode_pos_path,
                expected="at least one row",
            )
        )

    if barcode_pos.shape[1] != 3:
        raise STMinerIOError(
            format_io_error(
                "Invalid barcodes_pos.tsv.gz schema",
                operation="validate barcode positions",
                path=barcode_pos_path,
                expected="3 columns: barcode, array_col, array_row",
                actual=f"columns={barcode_pos.shape[1]}",
            )
        )

    if adata.n_obs != len(barcode_pos):
        raise STMinerIOError(
            format_io_error(
                "Barcode count mismatch between matrix and barcodes_pos",
                operation="validate barcode alignment",
                path=path,
                expected=f"{adata.n_obs} rows in barcodes_pos.tsv.gz",
                actual=f"{len(barcode_pos)}",
                hint="Use barcodes_pos file generated from the same matrix",
            )
        )

    if barcode_pos["barcode"].isnull().any():
        raise STMinerIOError(
            format_io_error(
                "Barcode column contains missing values",
                operation="validate barcode alignment",
                path=barcode_pos_path,
                expected="no missing barcode",
            )
        )

    barcode_pos["barcode"] = barcode_pos["barcode"].astype(str)
    duplicated = barcode_pos["barcode"].duplicated()
    if duplicated.any():
        duplicated_values = barcode_pos.loc[duplicated, "barcode"].head(5).tolist()
        raise STMinerIOError(
            format_io_error(
                "Duplicate barcodes found in barcodes_pos.tsv.gz",
                operation="validate barcode alignment",
                path=barcode_pos_path,
                expected="unique barcode values",
                actual=f"examples={duplicated_values}",
            )
        )

    missing_in_pos = adata.obs_names.difference(barcode_pos["barcode"])
    if len(missing_in_pos) > 0:
        raise STMinerIOError(
            format_io_error(
                "Some matrix barcodes are missing in barcodes_pos.tsv.gz",
                operation="validate barcode alignment",
                path=barcode_pos_path,
                expected="all matrix barcodes present in barcode position file",
                actual=f"missing_count={len(missing_in_pos)}, examples={missing_in_pos[:5].tolist()}",
            )
        )

    barcode_pos = barcode_pos.set_index("barcode", drop=False)
    barcode_pos = barcode_pos.loc[adata.obs_names]

    barcode_pos['in_tissue'] = 1
    barcode_pos = barcode_pos[['barcode', 'in_tissue', 'array_row', 'array_col']]
    adata.obs = barcode_pos[['in_tissue', 'array_row', 'array_col']].copy()

    obsm = barcode_pos[['array_col', 'array_row']].to_numpy()
    validate_spatial_array(obsm, path=barcode_pos_path, source="barcodes_pos.tsv.gz")
    adata.obsm["spatial"] = obsm
    if (('x' not in adata.obs.keys()) | ('y' not in adata.obs.keys())) and 'spatial' in adata.obsm.keys():
        position = adata.obsm['spatial']
        x_min = position[:, 0].min()
        y_min = position[:, 1].min()
        adata.obs['x'] = merge_bin_coordinate(position[:, 0], x_min, bin_size=bin_size)
        adata.obs['y'] = merge_bin_coordinate(position[:, 1], y_min, bin_size=bin_size)
        adata.misc = {'up': position[:, 1].min(),
                      'down': position[:, 1].max(),
                      'left': position[:, 0].min(),
                      'right': position[:, 0].max()}
    return adata
