import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from Algorithm.distribution import get_exp_array


def plot_heatmap(result=None,
                 label=None,
                 adata=None,
                 gene_list=None,
                 cmap=None,
                 num_cols=4,
                 vmax=100,
                 vmin=0):
    """
    Plot the heatmap for the given list.
    :param gene_list:
    :type gene_list:
    :param cmap:
    :type cmap:
    :param result:
    :type result:
    :param adata:
    :type adata:
    :param label:
    :type label:
    :param num_cols:
    :type num_cols:
    :param vmax:
    :type vmax:
    :param vmin:
    :type vmin:
    """
    if gene_list is None:
        if label is None or result is None:
            raise 'Error: Parameter [label] and [result] should not be None.'
        else:
            gene_list = list(result[result['labels'] == label]['gene_id'])

    new_cmap = cmap if cmap is not None else 'viridis'
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
        sns.heatmap(arr,
                    cbar=False,
                    ax=ax,
                    cmap=new_cmap,
                    vmax=np.percentile(arr, vmax),
                    vmin=np.percentile(arr, vmin)
                    )
        ax.axis('off')
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
