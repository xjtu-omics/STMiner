import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

from anndata import AnnData
from scipy.sparse import csr_matrix
from scipy.stats import wasserstein_distance


def read_gem_file(gem_file):
    gem_dataframe = pd.read_csv(str(gem_file), sep='\t', comment='#', header=0)
    return gem_dataframe


def cos_similarity_distance(x, y):
    """Calculate the cosin similarity distance between two vectors"""
    if x.shape != y.shape:
        raise ValueError('x and y must have same dimensions.')
    x = np.array(x)
    y = np.array(y)
    return 1 - np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y))


def sub_index(center_index, distance, edge_index=0):
    return center_index - distance if center_index - distance > edge_index else edge_index


def add_index(center_index, distance, edge_index):
    return center_index + distance + 1 if center_index + distance <= edge_index else edge_index


def get_surround_matrix(matrix, index, distance):
    '''Get the surrounding matrix by index'''
    if distance < 1:
        raise ValueError("The value of distance must be larger than 0")
    row_edge = matrix.shape[0]
    col_edge = matrix.shape[1]
    row_index = index[0]
    col_index = index[1]
    return matrix[sub_index(row_index, distance): add_index(row_index, distance, row_edge),
           sub_index(col_index, distance): add_index(col_index, distance, col_edge)]


def merge_bin_coor(coor: np.ndarray, coor_min: int, bin_size: int):
    return np.floor((coor - coor_min) / bin_size).astype(int)


def get_bin_center(bin_coor: np.ndarray, coor_min: int, bin_size: int):
    return bin_coor * bin_size + coor_min + int(bin_size / 2)


def get_coor_matrix(df, merge_count=20):
    '''Transform gene expression df to expression matrix'''
    count = df.loc[:, ['x', 'y', 'UMICount']].groupby(['x', 'y'], as_index=False).sum('UMICount')
    bins = [np.array(range(count['x'].min(), count['x'].max() + merge_count, merge_count)),
            np.array(range(count['y'].min(), count['y'].max() + merge_count, merge_count))]
    matrix, xedges, yedges = np.histogram2d(x=count['x'],
                                            y=count['y'],
                                            weights=count['UMICount'],
                                            bins=bins)
    return matrix


def find_track(df):
    H = get_coor_matrix(df)
    # Calculate by columns
    col_total = []
    for i in range(1, H.shape[1] - 1):
        dis_1 = wasserstein_distance(H[:, i], H[:, i + 1])
        dis_2 = wasserstein_distance(H[:, i], H[:, i - 1])
        wd_average = (dis_1 + dis_2) / 2
        col_total.append(wd_average)
    # Calculate by row
    row_total = []
    for i in range(1, H.shape[0] - 1):
        dis_1 = wasserstein_distance(H[i, :], H[i + 1, :])
        dis_2 = wasserstein_distance(H[i, :], H[i - 1, :])
        wd_average = (dis_1 + dis_2) / 2
        row_total.append(wd_average)
    return {'row': row_total, 'col': col_total}


def find_undersampled_pixel(H, distance, mean_percentage=0.2):
    """
    distance: The Chebyshev distance of surround_matrix to pixel location
    mean_percentage: The percentage threshold to compare against the mean value of surround_matrix
    """
    bad_matrix = np.zeros(H.shape)
    for index, value in np.ndenumerate(H):
        surround_matrix = get_surround_matrix(H, index, distance)
        is_lower_than_mean = surround_matrix.mean() * mean_percentage > value
        not_blank_region = np.median(surround_matrix) > 0
        is_square_matrix = surround_matrix.shape[0] == surround_matrix.shape[1]
        is_undersampled = is_lower_than_mean & is_square_matrix & not_blank_region & (
                value > 0)
        if is_undersampled:
            bad_matrix[index[0], index[1]] = 1
    return bad_matrix


def enhance_df_info(df, bin_size=20):
    """Calculate the position of each merged bin."""
    if 'MIDCounts' in df.columns:
        df.rename(columns={'MIDCounts': 'UMICount'}, inplace=True)
    elif 'MIDCount' in df.columns:
        df.rename(columns={'MIDCount': 'UMICount'}, inplace=True)
    x_min = df['x'].min()
    y_min = df['y'].min()
    df['bin_x'] = merge_bin_coor(df['x'].values, x_min, bin_size)
    df['bin_y'] = merge_bin_coor(df['y'].values, y_min, bin_size)
    df['x_center'] = get_bin_center(df['bin_x'], x_min, bin_size)
    df['y_center'] = get_bin_center(df['bin_y'], y_min, bin_size)
    df['cell_id'] = df['bin_x'].astype(str) + '_' + df['bin_y'].astype(str)
    return df


def get_annadata(df):
    """
    Get [bins * genes] matrix from inhanced df(enhance_df_info)\n
    The cell_id is the id of each merged bin (size = bin_size * bin_size)\n
    """
    # cells_id
    cells = df['cell_id'].unique()
    cells_dict = dict(zip(cells, range(0, len(cells))))
    rows = df['cell_id'].map(cells_dict)
    # genes_id
    genes = df['geneID'].unique()
    genes_dict = dict(zip(genes, range(0, len(genes))))
    cols = df['geneID'].map(genes_dict)
    # position
    position_df = df[['cell_id', 'bin_x', 'bin_y', 'x_center', 'y_center']]
    position_df.drop_duplicates(subset=['cell_id'])
    x_dict = dict(zip(position_df['cell_id'], position_df['bin_x']))
    y_dict = dict(zip(position_df['cell_id'], position_df['bin_y']))
    x_center_dict = dict(zip(position_df['cell_id'], position_df['x_center']))
    y_center_dict = dict(zip(position_df['cell_id'], position_df['y_center']))
    umi_arr = csr_matrix((df['UMICount'], (rows, cols)),
                         shape=(cells.shape[0], genes.shape[0]),
                         dtype=np.int32)
    adata = AnnData(X=umi_arr, dtype=np.float32)
    adata.var['gene_ids'] = genes
    adata.obs['cell_id'] = cells
    adata.obs['x'] = list(map(lambda id: x_dict[id], cells))
    adata.obs['y'] = list(map(lambda id: y_dict[id], cells))
    adata.obs['fig_x'] = list(map(lambda id: x_center_dict[id], cells))
    adata.obs['fig_y'] = list(map(lambda id: y_center_dict[id], cells))
    return adata


def view_under_sampled_matrix(matrix, dpi=200):
    plt.figure(dpi=dpi)
    sns.heatmap(matrix,
                cmap=sns.blend_palette(['white', 'black'],
                                       as_cmap=True))


def view_genes_matrix(matrix, dpi=200):
    max_value = np.mean(matrix) + np.std(matrix) * 3
    matrix[matrix > max_value] = max_value
    std = np.std(matrix)
    plt.figure(dpi=dpi)
    # disable the legend
    plt.axis('off')
    figure = sns.heatmap(matrix,
                         cbar=False,
                         cmap=sns.blend_palette(['white',
                                                 'blue'],
                                                as_cmap=True))
    return figure


def normalize_gene_count(adata, exclude_highly_expressed=True, max_fraction=0.05):
    sc.pp.normalize_total(adata,
                          exclude_highly_expressed=exclude_highly_expressed,
                          max_fraction=max_fraction,
                          inplace=True)


def clean_up_data(adata, **args):
    if len(args.items()) == 0:
        raise ValueError("None filter conditions specified.")
    for key, value in args.items():
        if key == 'min_cells_for_gene':
            sc.pp.filter_genes(adata, min_cells=value)
        elif key == 'max_cells_for_gene':
            sc.pp.filter_genes(adata, max_cells=value)
        elif key == 'min_genes_for_cell':
            sc.pp.filter_cells(adata, min_genes=value)
        elif key == 'max_genes_for_cell':
            sc.pp.filter_cells(adata, max_genes=value)
