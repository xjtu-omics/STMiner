import os
import matplotlib.pyplot as plt
import numpy as np
import umap
import pandas as pd
import seaborn as sns
from scipy.sparse import csr_matrix
from sklearn.metrics import davies_bouldin_score, calinski_harabasz_score
from sklearn.metrics import silhouette_score
from sklearn.manifold import TSNE
from collections import Counter
from STMiner.Algorithm.distance import get_exp_array
import matplotlib.image as mpimg


def _adjust_arr(arr, rotate, reverse_x, reverse_y):
    if rotate:
        arr = np.rot90(arr)
    if reverse_y:
        arr = np.flipud(arr)
    if reverse_x:
        arr = np.fliplr(arr)
    return arr


def _get_figure(fig_count, num_cols):
    num_rows = (fig_count + num_cols - 1) // num_cols
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 3 * num_rows))
    # Disable the axis for each subplot
    for ax in axes.flat:
        ax.axis('off')
    return axes, fig


def scale_array(exp_matrix, total_count):
    total_sum = np.sum(exp_matrix)
    scale_factor = 100 / total_sum
    scaled_matrix = exp_matrix * scale_factor
    total_count += scaled_matrix
    return total_count


class Plot:
    def __init__(self, sp):
        self.sp = sp
        self.count_array = {}

    def plot_gene(self,
                  gene,
                  cmap='Spectral_r',
                  reverse_y=False,
                  reverse_x=False,
                  rotate=False,
                  spot_size=None,
                  log1p=False,
                  save_path='',
                  dpi=400):
        arr = get_exp_array(self.sp.adata, gene)
        arr = _adjust_arr(arr, rotate, reverse_x, reverse_y)
        if log1p:
            arr = np.log1p(arr)
        sparse_matrix = csr_matrix(arr)
        if spot_size is not None:
            ax = sns.scatterplot(x=sparse_matrix.nonzero()[1],
                                 y=sparse_matrix.nonzero()[0],
                                 c=sparse_matrix.data,
                                 s=spot_size,
                                 cmap=cmap)
        else:
            ax = sns.scatterplot(x=sparse_matrix.nonzero()[1],
                                 y=sparse_matrix.nonzero()[0],
                                 c=sparse_matrix.data,
                                 cmap=cmap)
        ax.set_axis_off()
        ax.set_title(gene)
        if len(save_path) > 0:
            fig = ax.get_figure()
            if save_path[-1] != '/':
                save_path += '/'
            save_path += gene
            save_path += '.png'
            fig.savefig(fname=save_path, dpi=dpi)
        plt.show()

    def plot_genes(self,
                   label=None,
                   gene_list=None,
                   n_gene=12,
                   cmap=None,
                   num_cols=4,
                   vmax=100,
                   vmin=0,
                   rotate=False,
                   reverse_y=False,
                   reverse_x=False,
                   plot_type='heatmap',
                   s=1):
        result = self.sp.genes_labels
        adata = self.sp.adata
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
            arr = _adjust_arr(arr, rotate, reverse_x, reverse_y)
            sns.set(style="white")
            if cmap is None:
                cmap = sns.color_palette("viridis", as_cmap=True)
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
                                cmap=cmap,
                                s=s)
            ax.set_axis_off()
            ax.set_title(gene)
        plt.tight_layout()
        plt.show()

    def plot_pattern(self,
                     cmap=None,
                     vmax=99,
                     num_cols=4,
                     rotate=False,
                     reverse_y=False,
                     reverse_x=False,
                     vote_rate=0.2,
                     heatmap=True,
                     s=1,
                     image_path=None,
                     rotate_img=False,
                     k=1,
                     aspect=1,
                     output_path=None):
        # clean up count array
        self.count_array = None
        # initialize values
        result = self.sp.genes_labels
        adata = self.sp.adata
        label_list = set(result['labels'])
        plot_count = len(label_list)
        axes, fig = _get_figure(plot_count, num_cols=num_cols)
        fig.subplots_adjust(hspace=0.5)
        # initialize the count_array
        self.count_array = np.zeros(get_exp_array(adata, self.sp.adata.var.index[0]).shape)
        for i, label in enumerate(label_list):
            row = i // num_cols
            col = i % num_cols
            if len(axes.shape) == 1:
                ax = axes[i]
            else:
                ax = axes[row, col]
            gene_list = list(result[result['labels'] == label]['gene_id'])
            total_count = np.zeros(get_exp_array(adata, gene_list[0]).shape)
            total_coo_list = []
            vote_array = np.zeros(get_exp_array(adata, gene_list[0]).shape)
            for gene in gene_list:
                exp_matrix = get_exp_array(adata, gene)
                # calculate nonzero index
                non_zero_coo_list = np.vstack((np.nonzero(exp_matrix))).T.tolist()
                for coo in non_zero_coo_list:
                    total_coo_list.append(tuple(coo))
                total_count = scale_array(exp_matrix, total_count)
            count_dict = Counter(total_coo_list)
            for ele, count in count_dict.items():
                if int(count) / len(gene_list) >= vote_rate:
                    vote_array[ele] = 1
            total_count = total_count * vote_array
            total_count = _adjust_arr(total_count, rotate, reverse_x, reverse_y)
            binary_arr = np.where(total_count != 0, 1, total_count)
            self.count_array[label] = binary_arr
            if heatmap:
                sns.heatmap(total_count,
                            ax=ax,
                            cbar=False,
                            cmap=cmap if cmap is not None else 'viridis',
                            vmax=np.percentile(total_count, vmax))
            else:
                if isinstance(image_path, str) and os.path.isfile(image_path):
                    bg_img = mpimg.imread(image_path)
                    if rotate_img:
                        bg_img = np.rot90(bg_img, k=k)
                    ax.imshow(bg_img, extent=[0, total_count.shape[0], 0, total_count.shape[1]], aspect=aspect)
                sparse_matrix = csr_matrix(total_count)
                sns.scatterplot(x=sparse_matrix.nonzero()[0],
                                y=sparse_matrix.nonzero()[1],
                                c=sparse_matrix.data,
                                ax=ax,
                                cmap=cmap if cmap is not None else 'viridis',
                                s=s,
                                edgecolor='none')
                ax.set_xlim(0, total_count.shape[0])
                ax.set_ylim(0, total_count.shape[1])
            ax.set_title('Pattern ' + str(label))
        if output_path is not None and os.path.isdir(output_path):
            plt.savefig(os.path.join(output_path, "./scatterplot.eps"), dpi=1000, format='eps')
        plt.tight_layout()
        plt.show()

    def plot_count(self, pattern_list, cmap=None, s=None):

        sum_array = np.zeros(self.count_array[pattern_list[0]].shape)
        for i in pattern_list:
            sum_array += self.count_array[i]

        sparse_matrix = csr_matrix(sum_array)
        sns.scatterplot(x=sparse_matrix.nonzero()[0],
                        y=sparse_matrix.nonzero()[1],
                        c=sparse_matrix.data,
                        cmap=cmap if cmap is not None else 'viridis',
                        s=s if s is not None else 10,
                        edgecolor='none')

    def plot_cluster_score(self, mds_comp, min_cluster, max_cluster):
        db_dict = {}
        ch_dict = {}
        si_dict = {}
        for cluster_number in range(min_cluster, max_cluster + 1):
            self.sp.cluster_gene(self, n_clusters=cluster_number, mds_components=mds_comp)
            db_dict[cluster_number] = 1 / davies_bouldin_score(self.sp.genes_distance_array,
                                                               self.sp.kmeans_fit_result.labels_)
            ch_dict[cluster_number] = calinski_harabasz_score(self.sp.genes_distance_array,
                                                              self.sp.kmeans_fit_result.labels_)
            si_dict[cluster_number] = silhouette_score(self.sp.genes_distance_array, self.sp.kmeans_fit_result.labels_)
        score_df = pd.DataFrame([db_dict, si_dict, ch_dict],
                                index=['1/Davies-Bouldin', 'Silhouette', 'Calinski-Harabasz']).T
        norm_score_df = (score_df - score_df.min()) / (score_df.max() - score_df.min())
        sns.lineplot(norm_score_df, markers=True)
        plt.xticks(list(range(min_cluster, max_cluster + 1, 1)))
        plt.title("Evaluate Clustering Performance")
        plt.xlabel("Number of Clusters")
        plt.ylabel("Normalized Score")
        plt.show()

    def plot_tsne(self, method='tsne'):
        n_clusters = len(set(self.sp.kmeans_fit_result.labels_))
        if n_clusters <= 10:
            cmap = 'tab10'
        elif 10 < n_clusters & n_clusters <= 20:
            cmap = 'tab20'
        else:
            cmap = 'viridis'
        if method == 'tsne':
            tsne = TSNE(n_components=2, random_state=42)
            embedded_data = tsne.fit_transform(self.sp.mds_features)
        else:
            umap_model = umap.UMAP(n_neighbors=5, min_dist=0.3, n_components=2)
            embedded_data = umap_model.fit_transform(self.sp.mds_features)
        plt.figure(figsize=(8, 6))
        plt.scatter(embedded_data[:, 0],
                    embedded_data[:, 1],
                    c=self.sp.kmeans_fit_result.labels_,
                    cmap=cmap)
        plt.title("T-SNE")
        plt.xlabel("Dimension 1")
        plt.ylabel("Dimension 2")
        plt.grid(False)
        plt.xticks([])
        plt.yticks([])
        plt.show()
