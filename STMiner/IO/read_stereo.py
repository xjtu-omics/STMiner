import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from anndata import AnnData
from scipy.sparse import csr_matrix
from scipy.stats import wasserstein_distance

from STMiner.IO.IOUtil import *


def read_gem_file(gem_file, bin_size=40):
    gem_dataframe = pd.read_csv(str(gem_file), sep='\t', comment='#', header=0)
    df = enhance_df_info(gem_dataframe, bin_size=bin_size)
    adata = get_anndata(df)
    return adata


def cos_similarity_distance(x, y):
    """Calculate the cosine similarity distance between two vectors"""
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
    """Get the surrounding matrix by index"""
    if distance < 1:
        raise ValueError("The value of distance must be larger than 0")
    row_edge = matrix.shape[0]
    col_edge = matrix.shape[1]
    row_index = index[0]
    col_index = index[1]
    return matrix[
           sub_index(row_index, distance): add_index(row_index, distance, row_edge),
           sub_index(col_index, distance): add_index(col_index, distance, col_edge)]


def get_coordinate_matrix(df, merge_count=20):
    """Transform gene expression df to expression matrix"""
    count = df.loc[:, ['x', 'y', 'UMICount']].groupby(['x', 'y'], as_index=False).sum('UMICount')
    bins = [np.array(range(count['x'].min(), count['x'].max() + merge_count, merge_count)),
            np.array(range(count['y'].min(), count['y'].max() + merge_count, merge_count))]
    matrix, x_edges, y_edges = np.histogram2d(x=count['x'],
                                              y=count['y'],
                                              weights=count['UMICount'],
                                              bins=bins)
    return matrix


def find_track(df):
    mat = get_coordinate_matrix(df)
    # Calculate by columns
    col_total = []
    for i in range(1, mat.shape[1] - 1):
        dis_1 = wasserstein_distance(mat[:, i], mat[:, i + 1])
        dis_2 = wasserstein_distance(mat[:, i], mat[:, i - 1])
        wd_average = (dis_1 + dis_2) / 2
        col_total.append(wd_average)
    # Calculate by row
    row_total = []
    for i in range(1, mat.shape[0] - 1):
        dis_1 = wasserstein_distance(mat[i, :], mat[i + 1, :])
        dis_2 = wasserstein_distance(mat[i, :], mat[i - 1, :])
        wd_average = (dis_1 + dis_2) / 2
        row_total.append(wd_average)
    return {'row': row_total, 'col': col_total}


def find_under_sampled_pixel(arr, distance, mean_percentage=0.2):
    """
    distance: The Chebyshev distance of surround_matrix to pixel location
    mean_percentage: The percentage threshold to compare against the mean value of surround_matrix
    """
    bad_matrix = np.zeros(arr.shape)
    for index, value in np.ndenumerate(arr):
        surround_matrix = get_surround_matrix(arr, index, distance)
        is_lower_than_mean = surround_matrix.mean() * mean_percentage > value
        not_blank_region = np.median(surround_matrix) > 0
        is_square_matrix = surround_matrix.shape[0] == surround_matrix.shape[1]
        is_under_sampled = is_lower_than_mean & is_square_matrix & not_blank_region & (
                value > 0)
        if is_under_sampled:
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
    df['bin_x'] = merge_bin_coordinate(df['x'].values, x_min, bin_size)
    df['bin_y'] = merge_bin_coordinate(df['y'].values, y_min, bin_size)
    df['x_center'] = get_bin_center(df['bin_x'], x_min, bin_size)
    df['y_center'] = get_bin_center(df['bin_y'], y_min, bin_size)
    df['cell_id'] = df['bin_x'].astype(str) + '_' + df['bin_y'].astype(str)
    return df


def get_anndata(df):
    """
    Get [bins * genes] matrix from enhanced df(enhance_df_info)\n
    The cell_id is the id of each merged bin (size = bin_size * bin_size)\n
    :param df:
    :type df: pd.DataFrame
    :return:
    :rtype: AnnData
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
    adata.obs['x'] = list(map(lambda gene_id: x_dict[gene_id], cells))
    adata.obs['y'] = list(map(lambda gene_id: y_dict[gene_id], cells))
    adata.obs['fig_x'] = list(map(lambda gene_id: x_center_dict[gene_id], cells))
    adata.obs['fig_y'] = list(map(lambda gene_id: y_center_dict[gene_id], cells))
    adata.var_names = list(adata.var['gene_ids'])
    adata.obsm['spatial'] = adata.obs[['x', 'y']].values
    return adata


def view_under_sampled_matrix(matrix, dpi=200):
    plt.figure(dpi=dpi)
    sns.heatmap(matrix,
                cmap=sns.blend_palette(['white', 'black'],
                                       as_cmap=True))


def view_genes_matrix(matrix, dpi=200):
    max_value = np.mean(matrix) + np.std(matrix) * 3
    matrix[matrix > max_value] = max_value
    plt.figure(dpi=dpi)
    # disable the legend
    plt.axis('off')
    figure = sns.heatmap(matrix,
                         cbar=False,
                         cmap=sns.blend_palette(['white',
                                                 'blue'],
                                                as_cmap=True))
    return figure
