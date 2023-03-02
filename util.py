import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

from tqdm import tqdm
from scipy.signal import convolve2d
from scipy.signal import convolve
from scipy.sparse import csr_matrix


def get_3D_matrix(adata):
    x_max = int(adata.obs['x'].max())
    y_max = int(adata.obs['y'].max())
    # the spatial coordinates should be in adata.obs
    threeD_array = np.zeros((int(x_max), int(y_max), int(adata.var.shape[0])))
    print('Transfer anndate to 3D matrix...')
    for spot in tqdm(adata, bar_format='{l_bar}{bar:20}{r_bar}{percentage:3.0f}%'):
        x = int(spot.obs['x'])-1
        y = int(spot.obs['y'])-1
        threeD_array[x, y] = spot.X.toarray()
    return threeD_array


def convolve(array, big_kernel=True):
    kernel_A = np.array([[-1, -1, -1, -1, -1],
                         [-1,  2,  2,  2, -1],
                         [-1,  2,  8,  2, -1],
                         [-1,  2,  2,  2, -1],
                         [-1, -1, -1, -1, -1]])

    kernel_B = np.array([[0, -1, 0],
                         [-1, 5, -1],
                         [0, -1, 0]])

    n, m, k = array.shape
    # convolve each 2D layer
    output_array = np.zeros((n, m, k))
    print('Convolve each 2D layer...')
    
    if big_kernel:
        for i in tqdm(range(k), bar_format='{l_bar}{bar:20}{r_bar}{percentage:3.0f}%'):
            output_array[:, :, i] = convolve2d(
                array[:, :, i], kernel_A, mode='same')
    else:
        for i in tqdm(range(k), bar_format='{l_bar}{bar:20}{r_bar}{percentage:3.0f}%'):
            output_array[:, :, i] = convolve2d(
                array[:, :, i], kernel_B, mode='same')
            
    output_array = np.where(output_array < 0, 0, output_array)
    return output_array


def update_anndata(array, adata):
    print('Update anndata...')
    for spot in tqdm(adata, bar_format='{l_bar}{bar:20}{r_bar}{percentage:3.0f}%'):
        x = int(spot.obs['x'])-1
        y = int(spot.obs['y'])-1
        spot.X = csr_matrix(array[x, y])


def run_convolve(adata):
    result = convolve(get_3D_matrix(adata))
    update_anndata(result, adata)


def add_spatial_position(adata, position_file):
    position_df = pd.read_csv(position_file,
                              sep=',',
                              header=None,
                              index_col=0)
    # set the column names
    position_df.columns = ['in_tissue', 'array_row',
                           'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    # get the cell positions
    cell_info_df = position_df.loc[adata.obs.index]
    # set the matrix position
    adata.obs['x'] = cell_info_df['array_row'].to_numpy()
    adata.obs['y'] = cell_info_df['array_col'].to_numpy()
    # set the figure position
    adata.obs['fig_x'] = cell_info_df['pxl_row_in_fullres'].to_numpy()
    adata.obs['fig_y'] = cell_info_df['pxl_col_in_fullres'].to_numpy()
    adata.obsm['spatial'] = np.array([cell_info_df['pxl_row_in_fullres'].to_numpy(), cell_info_df['pxl_col_in_fullres'].to_numpy()]).T
