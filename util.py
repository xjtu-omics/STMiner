import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import tifffile as tiff

from tqdm import tqdm
from PIL import Image
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
    pixel_array = np.array([cell_info_df['pxl_row_in_fullres'].to_numpy(),
                            cell_info_df['pxl_col_in_fullres'].to_numpy()]).T
    adata.uns['spatial'] = pixel_array
    adata.obsm['spatial'] = pixel_array


def add_image(adata, image, spatial_key='spatial', library_id='tissue', tissue_hires_scalef=1, spot_diameter_fullres=89):
    # spot_diameter_fullres:
    # this is the diameter of the capture area for each observation.
    # In the case of Visium, we usually call them “spots” and this value is set to ~89.
    # Ref: https://squidpy.readthedocs.io/en/stable/auto_tutorials/tutorial_read_spatial.html
    adata.uns[spatial_key] = {library_id: {}}
    adata.uns[spatial_key][library_id]["images"] = {}
    adata.uns[spatial_key][library_id]["images"] = {"hires": image}
    adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": tissue_hires_scalef,
                                                          "spot_diameter_fullres": spot_diameter_fullres}


class TissueImage:
    def __init__(self, img_path):
        self.img_path = img_path
        self.img = tiff.imread(img_path)[:, :, :3]

    def rotate(self):
        self.img = np.rot90(self.img, k=1)

    def flip_lr(self):
        self.img = np.fliplr(self.img)

    def flip_ud(self):
        self.img = np.flipud(self.img)

    def get_image(self):
        return self.img

    def preview(self):
        pil_img = Image.fromarray(self.img)
        thumbnail_size = (self.img.shape[1]//10, self.img.shape[0]//10)
        pil_img.thumbnail(thumbnail_size)
        pil_img.show()
