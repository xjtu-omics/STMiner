import os

import scanpy as sc

from STMiner.IO.IOUtil import *


def read_h5ad(file, bin_size=1, merge_bin=False, amplification=1):
    """
    Reads an h5ad file and processes spatial data.
    Parameters:
    -----------
    file : str
        Path to the h5ad file to be read.
    bin_size : int, optional
        Size of the bin for spatial data aggregation. Default is 1.
    merge_bin : bool, optional
        If True, spatial data will be binned using the specified bin_size. Default is False.
    amplification : int, optional
        Factor by which to amplify the spatial coordinates. Default is 1.
    Returns:
    --------
    adata : AnnData
        An AnnData object containing the processed spatial data.
    """
    path = str(file)
    require_positive_int(bin_size, "bin_size", path=path)
    require_positive_int(amplification, "amplification", path=path)

    if not os.path.isfile(path):
        raise FileNotFoundError(
            format_io_error(
                "H5AD file does not exist",
                operation="read_h5ad",
                path=path,
                expected="existing .h5ad file",
                hint="Check file path and access permission",
            )
        )

    try:
        adata = sc.read_h5ad(path)
    except Exception as e:
        raise STMinerIOError(
            format_io_error(
                "Failed to read h5ad file",
                operation="scanpy.read_h5ad",
                path=path,
                expected="valid AnnData .h5ad file",
                actual=f"{type(e).__name__}: {e}",
                hint="File may be corrupted or not in AnnData format",
            )
        ) from e

    if "spatial" not in adata.obsm:
        raise STMinerIOError(
            format_io_error(
                "Missing required spatial coordinates",
                operation="read_h5ad",
                path=path,
                expected="adata.obsm['spatial']",
                actual=f"obsm_keys={list(adata.obsm.keys())}",
                hint="Provide h5ad with spatial coordinates in adata.obsm['spatial']",
            )
        )

    position = validate_spatial_array(adata.obsm["spatial"], path=path)

    if merge_bin:
        if int(amplification) != 1:
            adata = adata.copy()
            adata.obsm["spatial"] = position * np.int32(amplification)
        try:
            return bin_spatial_adata(adata, bin_size=bin_size)
        except Exception as e:
            if isinstance(e, STMinerIOError):
                raise
            raise STMinerIOError(
                format_io_error(
                    "Failed while binning spatial AnnData",
                    operation="bin_spatial_adata",
                    path=path,
                    actual=f"{type(e).__name__}: {e}",
                    hint="Check that adata.obsm['spatial'] exists and contains valid coordinates",
                )
            ) from e

    amplification = np.int32(amplification)
    x_min = position[:, 0].min() * amplification
    y_min = position[:, 1].min() * amplification
    adata.obs["x"] = merge_bin_coordinate(
        position[:, 0] * amplification, x_min, bin_size=bin_size
    )
    adata.obs["y"] = merge_bin_coordinate(
        position[:, 1] * amplification, y_min, bin_size=bin_size
    )
    adata.misc = {
        "up": position[:, 1].min(),
        "down": position[:, 1].max(),
        "left": position[:, 0].min(),
        "right": position[:, 0].max(),
    }
    return adata
