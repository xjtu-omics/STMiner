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
        
    if merge_bin:
        return bin_spatial_adata(sc.read_h5ad(file), bin_size=bin_size)
    else:
        amplification = np.int32(amplification)
        adata = sc.read_h5ad(file)

        position = adata.obsm["spatial"]
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
