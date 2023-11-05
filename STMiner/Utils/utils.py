import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.signal as signal
import tifffile as tiff
from PIL import Image
from scipy.signal import convolve2d
from scipy.sparse import csr_matrix, issparse
from tqdm import tqdm


def get_3d_matrix(adata: anndata):
    """
    get the 3D matrix for adata
    """
    x_max = int(adata.obs['x'].max())
    y_max = int(adata.obs['y'].max())
    # the spatial coordinates should be in adata.obs
    three_d_array = np.zeros((int(x_max), int(y_max), int(adata.var.shape[0])), dtype=np.int32)
    print('Transfer anndata to 3D matrix...')
    for spot in tqdm(adata, bar_format='{l_bar}{bar:20}{r_bar}{percentage:3.0f}%'):
        x = int(spot.obs['x']) - 1
        y = int(spot.obs['y']) - 1
        three_d_array[x, y] = spot.X.toarray()
    return three_d_array


def get_gaussian_kernel(size=5, sigma=1):
    kernel = np.zeros((size, size))
    center = size // 2
    for i in range(size):
        for j in range(size):
            x = i - center
            y = j - center
            kernel[i, j] = np.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2))
    kernel /= kernel.sum()
    return kernel


def get_laplacian_kernel():
    return np.array([[0, -1, 0],
                     [-1, 5, -1],
                     [0, -1, 0]])


def get_mean_filter_kernel(size=3):
    return np.ones([size, size]) * (1 / size ** 2)


def convolve(array, method, kernel_size=3):
    x, y, gene_index = array.shape
    # convolve each 2D layer
    output_array = np.zeros((x, y, gene_index))
    print('Convolve each 2D layer...')
    if method == 'gaussian':
        for i in tqdm(range(gene_index)):
            output_array[:, :, i] = convolve2d(array[:, :, i],
                                               get_gaussian_kernel(size=kernel_size),
                                               mode='same')
    elif method == 'mid':
        for i in tqdm(range(gene_index)):
            output_array[:, :, i] = signal.medfilt2d(array[:, :, i], kernel_size=kernel_size)
    output_array = np.where(output_array < 0, 0, output_array)
    return output_array


def update_anndata(array: np.array, adata: anndata):
    print('Update anndata...')
    if issparse(adata[0].X):
        for spot in tqdm(adata):
            x = int(spot.obs['x']) - 1
            y = int(spot.obs['y']) - 1
            spot.X = csr_matrix(array[x, y])
    else:
        for spot in tqdm(adata):
            x = int(spot.obs['x']) - 1
            y = int(spot.obs['y']) - 1
            spot.X = array[x, y]


def is_sparse(array):
    return issparse(array)


def run_sharp(adata: anndata):
    result = convolve(get_3d_matrix(adata), get_laplacian_kernel())
    update_anndata(result, adata)


def run_gaussian(adata: anndata, shape=5, sigma=1):
    result = convolve(get_3d_matrix(adata), get_gaussian_kernel(shape, sigma))
    update_anndata(result, adata)


def add_spatial_position(adata: anndata, position_file: str):
    position_df = pd.read_csv(position_file, sep=',', header=None, index_col=0)
    # set the column names
    position_df.columns = ['in_tissue',
                           'array_row',
                           'array_col',
                           'pxl_row_in_fullres',
                           'pxl_col_in_fullres']
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
    # adata.uns['spatial'] = pixel_array
    adata.obsm['spatial'] = pixel_array


def add_image(adata, image, spatial_key='spatial', library_id='tissue',
              tissue_hires_scalef=1, spot_diameter_fullres=89):
    # spot_diameter_fullres:
    # this is the diameter of the capture area for each observation.
    # In the case of Visium, we usually call them “spots” and this value is set to ~89.
    # Ref: https://squidpy.readthedocs.io/en/stable/auto_tutorials/tutorial_read_spatial.html
    adata.uns[spatial_key] = {library_id: {}}
    adata.uns[spatial_key][library_id]["images"] = {}
    adata.uns[spatial_key][library_id]["images"] = {"hires": image}
    scalefactors = {"tissue_hires_scalef": tissue_hires_scalef,
                    "spot_diameter_fullres": spot_diameter_fullres}
    adata.uns[spatial_key][library_id]["scalefactors"] = scalefactors


def plot_scatter(adata, size=5):
    sc.pl.scatter(adata,
                  x='x',
                  y='y',
                  color='log1p_n_genes_by_counts',
                  size=size,
                  palette=['#FF0000', '#00FF00', '#0000FF'])


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
        thumbnail_size = (self.img.shape[1] // 10, self.img.shape[0] // 10)
        pil_img.thumbnail(thumbnail_size)
        pil_img.show()
