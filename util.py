import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import tifffile as tiff
import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm
from PIL import Image
from scipy.signal import convolve2d
from scipy.sparse import csr_matrix
from scipy.stats import multivariate_normal

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


def get_gaussian_kernel(size=5, sigma=1):
    kernel = np.zeros((size, size))
    center = size // 2
    for i in range(size):
        for j in range(size):
            x = i - center
            y = j - center
            kernel[i, j] = np.exp(-(x**2 + y**2) / (2 * sigma**2))

    kernel /= kernel.sum()
    return kernel


def get_laplacian_kernel():
    return np.array([[0, -1, 0],
                     [-1, 5, -1],
                     [0, -1, 0]])


def get_mean_filter_kernel(size=3):
    return np.ones([size, size])*(1/size**2)


def convolve(array, kernel):
    n, m, k = array.shape
    # convolve each 2D layer
    output_array = np.zeros((n, m, k))
    print('Convolve each 2D layer...')
    for i in tqdm(range(k), bar_format='{l_bar}{bar:20}{r_bar}{percentage:3.0f}%'):
        output_array[:, :, i] = convolve2d(array[:, :, i],
                                           kernel,
                                           mode='same')

    output_array = np.where(output_array < 0, 0, output_array)
    return output_array


def update_anndata(array, adata):
    print('Update anndata...')
    for spot in tqdm(adata, bar_format='{l_bar}{bar:20}{r_bar}{percentage:3.0f}%'):
        x = int(spot.obs['x'])-1
        y = int(spot.obs['y'])-1
        spot.X = csr_matrix(array[x, y])


def run_sharp(adata):
    result = convolve(get_3D_matrix(adata), get_laplacian_kernel())
    update_anndata(result, adata)


def run_gaussian(adata, shape=5, sigma=1):
    result = convolve(get_3D_matrix(adata), get_gaussian_kernel(shape, sigma))
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
    # adata.uns['spatial'] = pixel_array
    adata.obsm['spatial'] = pixel_array


def add_image(adata, image, spatial_key='spatial', library_id='tissue', tissue_hires_scalef=1, spot_diameter_fullres=89):
    adata: anndata
    image: str
    # spot_diameter_fullres:
    # this is the diameter of the capture area for each observation.
    # In the case of Visium, we usually call them “spots” and this value is set to ~89.
    # Ref: https://squidpy.readthedocs.io/en/stable/auto_tutorials/tutorial_read_spatial.html
    adata.uns[spatial_key] = {library_id: {}}
    adata.uns[spatial_key][library_id]["images"] = {}
    adata.uns[spatial_key][library_id]["images"] = {"hires": image}
    adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": tissue_hires_scalef,
                                                          "spot_diameter_fullres": spot_diameter_fullres}


def calculate_bhattacharyya_distances(gmm1, gmm2):
    mu1 = gmm1.means_
    cov1 = gmm1.covariances_
    w1 = gmm1.weights_
    mu2 = gmm2.means_
    cov2 = gmm2.covariances_
    w2 = gmm2.weights_

    cov = (cov1.mean(axis=0) + cov2.mean(axis=0)) / 2
    diff_transpose_cov = (mu1.mean(axis=0) - mu2.mean(axis=0)).T.dot(np.linalg.inv(cov))
    sqrt_bc = np.sqrt((w1 * w2).sum()) * np.exp(-1/8 * diff_transpose_cov.dot(mu1.mean(axis=0) - mu2.mean(axis=0)))
    bd = -np.log(sqrt_bc)
    return bd

def bd(gmm1, gmm2):
    # 定义两个GMM
    gmm1_weights = gmm1.weights_
    gmm1_means = gmm1.means_
    gmm1_covs = gmm1.covariances_

    gmm2_weights = gmm.weights_
    gmm2_means = gmm.means_
    gmm2_covs = gmm.covariances_
    n_components = gmm2_weights.size
    # 计算所有组件之间的Bhattacharyya距离
    bhat_dist = np.zeros((n_components, n_components))
    for i in range(n_components):
        for j in range(n_components):
            bhat_dist[i, j] = 0.25 * (np.log(np.linalg.det(0.5*(gmm1_covs[i]+gmm2_covs[j]))) 
                                    - 0.5*np.log(np.linalg.det(gmm1_covs[i])) 
                                    - 0.5*np.log(np.linalg.det(gmm2_covs[j]))
                                    + 0.5*((gmm1_means[i] - gmm2_means[j]).T @ np.linalg.inv(0.5*(gmm1_covs[i]+gmm2_covs[j]))
                                            @ (gmm1_means[i] - gmm2_means[j]))
                                    + 0.125*np.trace(np.linalg.inv(0.5*(gmm1_covs[i]+gmm2_covs[j])) @ gmm1_covs[i])
                                    + 0.125*np.trace(np.linalg.inv(0.5*(gmm1_covs[i]+gmm2_covs[j])) @ gmm2_covs[j])
                                    )
    weight_bhat_dist = bhat_dist * gmm1_weights.reshape(n_components, 1)
    min_sum_dist = np.sum(np.amin(min_bd, axis=1))
    return min_sum_dist


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
        sc.pl.spatial()
