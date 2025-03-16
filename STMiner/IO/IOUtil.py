import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from scipy.sparse import issparse

def merge_bin_coordinate(coordinate: np.ndarray, coordinate_min: int, bin_size: int):
    return np.floor((coordinate - coordinate_min) / bin_size).astype(int)

def get_bin_center(bin_coordinate: np.ndarray, coordinate_min: int, bin_size: int):
    return bin_coordinate * bin_size + coordinate_min + int(bin_size / 2)

def bin_spatial_adata(adata, bin_size=50, method="sum"):
    coords = adata.obsm["spatial"]
    grid_x = (coords[:, 0] // bin_size).astype(int)
    grid_y = (coords[:, 1] // bin_size).astype(int)

    grid_df = pd.DataFrame({"grid_x": grid_x, "grid_y": grid_y, "index": np.arange(len(adata))})
    grouped = grid_df.groupby(["grid_x", "grid_y"])["index"].apply(list)
    
    new_X = []
    new_coords = []
    new_obs_names = []
    for (gx, gy), indices in grouped.items():
        sub_X = adata.X[indices] 
        if issparse(sub_X):
            sub_X = sub_X.toarray()
        if method == "sum":
            aggregated_X = np.sum(sub_X, axis=0)
        elif method == "mean":
            aggregated_X = np.mean(sub_X, axis=0)
        
        new_X.append(aggregated_X.reshape(1, -1))
        new_coords.append([(gx + 0.5) * bin_size, (gy + 0.5) * bin_size])
        new_obs_names.append(f"bin_{gx}_{gy}")

    new_X = np.vstack(new_X)
    new_coords = np.array(new_coords, dtype=int)

    new_adata = ad.AnnData(X=new_X, var=adata.var.copy(), uns=adata.uns.copy())
    new_adata.obsm["spatial"] = new_coords
    new_adata.obs_names = new_obs_names 
    new_adata.obs['x'] = new_coords[:, 0]
    new_adata.obs['y'] = new_coords[:, 1]
    
    return new_adata
