import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from scipy.sparse import csr_matrix
from Algorithm.distribution import get_exp_array


def plot_gene(adata,
              gene,
              cmap='Spectral_r',
              reverse_y=False,
              reverse_x=False,
              spot_size=None):
    arr = get_exp_array(adata, gene)
    arr = _adjust_arr(arr, reverse_x, reverse_y)
    sparse_matrix = csr_matrix(arr)
    if spot_size is not None:
        sns.scatterplot(x=sparse_matrix.nonzero()[1],
                        y=sparse_matrix.nonzero()[0],
                        c=sparse_matrix.data,
                        s=spot_size,
                        cmap=cmap)
    else:
        sns.scatterplot(x=sparse_matrix.nonzero()[1],
                        y=sparse_matrix.nonzero()[0],
                        c=sparse_matrix.data,
                        cmap=cmap)


def plot_genes(result=None,
               label=None,
               adata=None,
               gene_list=None,
               n_gene=12,
               cmap=None,
               num_cols=4,
               vmax=100,
               vmin=0,
               reverse_y=False,
               reverse_x=False,
               plot_type='heatmap'):
    if gene_list is None:
        if label is None or result is None:
            raise 'Error: Parameter [label] and [result] should not be None.'
        else:
            gene_list = list(result[result['labels'] == label]['gene_id'])[:n_gene]
    genes_count = len(gene_list)
    axes, fig = _get_figure(genes_count, num_cols)
    fig.subplots_adjust(hspace=0.5)
    for i, gene in enumerate(gene_list):
        row = i // num_cols
        col = i % num_cols
        if len(axes.shape) == 1:
            ax = axes[i]
        else:
            ax = axes[row, col]
        arr = get_exp_array(adata, gene)
        arr = _adjust_arr(arr, reverse_x, reverse_y)
        sns.set(style="white")
        if cmap is None:
            cmap = sns.color_palette("Spectral_r", as_cmap=True)
        if plot_type == 'heatmap':
            sns.heatmap(arr,
                        cbar=False,
                        ax=ax,
                        cmap=cmap,
                        vmax=np.percentile(arr, vmax),
                        vmin=np.percentile(arr, vmin))
        elif plot_type == 'scatter':
            sparse_matrix = csr_matrix(arr)
            sns.scatterplot(x=sparse_matrix.nonzero()[1],
                            y=sparse_matrix.nonzero()[0],
                            ax=ax,
                            c=sparse_matrix.data,
                            cmap=cmap)
        ax.set_axis_off()
        ax.set_title(gene)
    plt.tight_layout()
    plt.show()


def plot_pattern(result, adata, cmap=None, vmax=99, num_cols=4):
    label_list = set(result['labels'])
    plot_count = len(label_list)
    axes, fig = _get_figure(plot_count, num_cols)
    fig.subplots_adjust(hspace=0.5)
    for i, label in enumerate(label_list):
        row = i // num_cols
        col = i % num_cols
        if len(axes.shape) == 1:
            ax = axes[i]
        else:
            ax = axes[row, col]
        gene_list = list(result[result['labels'] == label]['gene_id'])
        total_count = np.zeros(get_exp_array(adata, gene_list[0]).shape)
        for gene in gene_list:
            exp_matrix = get_exp_array(adata, gene)
            total_sum = np.sum(exp_matrix)
            scale_factor = 100 / total_sum
            scaled_matrix = exp_matrix * scale_factor
            total_count += scaled_matrix

        sns.heatmap(total_count,
                    ax=ax,
                    cbar=False,
                    cmap=cmap if cmap is not None else 'viridis',
                    vmax=np.percentile(total_count, vmax))
        ax.axis('off')
        ax.set_title('Pattern ' + str(label))
    plt.tight_layout()
    plt.show()


def _get_figure(fig_count, num_cols):
    num_rows = (fig_count + num_cols - 1) // num_cols
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 3 * num_rows))
    # Disable the axis for each subplot
    for ax in axes.flat:
        ax.axis('off')
    return axes, fig


def _adjust_arr(arr, reverse_x, reverse_y):
    if reverse_y:
        arr = np.flipud(arr)
    if reverse_x:
        arr = np.fliplr(arr)
    return arr
