import anndata
import numpy as np
import pandas as pd
import tifffile as tiff

from Algorithm import *
from tqdm import tqdm
from PIL import Image
from scipy.signal import convolve2d
from scipy.sparse import csr_matrix


def get_3d_matrix(adata):
    adata: anndata
    """
    get the 3D matrix for adata
    """
    x_max = int(adata.obs['x'].max())
    y_max = int(adata.obs['y'].max())
    # the spatial coordinates should be in adata.obs
    three_d_array = np.zeros((int(x_max), int(y_max), int(adata.var.shape[0])))
    print('Transfer anndate to 3D matrix...')
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


def convolve(array, kernel):
    n, m, k = array.shape
    # convolve each 2D layer
    output_array = np.zeros((n, m, k))
    print('Convolve each 2D layer...')
    for i in tqdm(range(k)):
        output_array[:, :, i] = convolve2d(array[:, :, i],
                                           kernel,
                                           mode='same')
    output_array = np.where(output_array < 0, 0, output_array)
    return output_array


def update_anndata(array, adata):
    array: np.array
    adata: anndata

    print('Update anndata...')
    for spot in tqdm(adata):
        x = int(spot.obs['x']) - 1
        y = int(spot.obs['y']) - 1
        spot.X = csr_matrix(array[x, y])


def run_sharp(adata):
    result = convolve(get_3d_matrix(adata), get_laplacian_kernel())
    update_anndata(result, adata)


def run_gaussian(adata, shape=5, sigma=1):
    result = convolve(get_3d_matrix(adata), get_gaussian_kernel(shape, sigma))
    update_anndata(result, adata)


def add_spatial_position(adata, position_file):
    adata: anndata
    position_file: str

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
    adata: anndata
    image: str

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


def bhat_distance(gmm1, gmm2):
    gmm1_weights = gmm1.weights_
    gmm1_means = gmm1.means_
    gmm1_covs = gmm1.covariances_
    gmm2_weights = gmm2.weights_
    gmm2_means = gmm2.means_
    gmm2_covs = gmm2.covariances_
    n_components = gmm1_weights.size
    # calculate the Bhattacharyya distance
    bhat_dist = np.zeros((n_components, n_components))
    for i in range(n_components):
        for j in range(n_components):
            bhat_dist[i, j] = get_hellinger_distance(gmm1_covs[i], gmm1_means[i], gmm2_covs[j], gmm2_means[j])
    # TODO: consider the weight of each component
    # min_bd = bh_dist * gmm1_weights.reshape(n_components, 1)
    min_cost = linear_sum(bhat_dist)
    return min_cost


def get_bh_distance(gmm1_covs, gmm1_means, gmm2_covs, gmm2_means):
    mean_cov = (gmm1_covs + gmm2_covs) / 2
    mean_cov_det = np.linalg.det(mean_cov)
    mean_cov_inv = np.linalg.inv(mean_cov)
    means_diff = gmm1_means - gmm2_means
    first_term = (means_diff.T @ mean_cov_inv @ means_diff) / 8
    second_term = np.log(
            mean_cov_det / (np.sqrt(np.linalg.det(gmm1_covs) * np.linalg.det(gmm1_covs)))) / 2
    result = first_term + second_term
    return result


def get_hellinger_distance(gmm1_covs, gmm1_means, gmm2_covs, gmm2_means):
    mean_cov = (gmm1_covs + gmm2_covs) / 2
    mean_cov_det = np.linalg.det(mean_cov)
    mean_cov_inv = np.linalg.inv(mean_cov)
    means_diff = gmm1_means - gmm2_means
    first_term = np.exp(-(means_diff.T @ mean_cov_inv @ means_diff) / 8)
    second_term = (np.linalg.det(gmm1_covs) ** (1 / 4) * np.linalg.det(gmm1_covs) ** (1 / 4)) / np.sqrt(mean_cov_det)
    result = 1 - first_term * second_term
    return result


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
