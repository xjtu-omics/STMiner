import os

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import umap
from matplotlib.colors import ListedColormap
from scipy.sparse import csr_matrix
from sklearn.manifold import TSNE

from STMiner.Algorithm.AlgUtils import get_exp_array
from STMiner.services.spfinder_analysis import evaluate_cluster_scores
from ._pattern import _adjust_arr, plot_patterns


def _get_figure(fig_count, num_cols):
    if fig_count <= 0:
        raise ValueError("fig_count must be a positive integer.")
    if num_cols <= 0:
        raise ValueError("num_cols must be a positive integer.")
    num_rows = (fig_count + num_cols - 1) // num_cols
    # squeeze=False keeps axes as a 2-D array even for a single subplot,
    # so axes.flat and axes[row, col] always work.
    fig, axes = plt.subplots(
        num_rows, num_cols, figsize=(12, 3 * num_rows), squeeze=False
    )
    # Disable the axis for each subplot
    for ax in axes.flat:
        ax.axis("off")
    return axes, fig


def is_path(image_path):
    return isinstance(image_path, str) and (
        os.path.isfile(image_path) or os.path.isdir(image_path)
    )


class Plot:
    def __init__(self, sp):
        self.sp = sp

    def plot_gene(
        self,
        gene,
        cmap="Spectral_r",
        reverse_y=False,
        reverse_x=False,
        rotate=False,
        figsize=(8, 6),
        s=5,
        log1p=False,
        save_path="",
        format="eps",
        dpi=400,
        vmax=99,
    ):
        """
        Plots the spatial expression of a given gene on a scatter plot.
        
        Args:
            gene (str): The name of the gene to plot.
            cmap (str, optional): Colormap to use for gene expression. Default is 'Spectral_r'.
            reverse_y (bool, optional): If True, reverse the y-axis. Default is False.
            reverse_x (bool, optional): If True, reverse the x-axis. Default is False.
            rotate (bool, optional): If True, rotate the plot by 90 degrees. Default is False.
            figsize (tuple, optional): Figure size in inches (width, height). Default is (8, 6).
            s (int or None, optional): Size of the scatter plot points. If None, uses default size. Default is 5.
            log1p (bool, optional): If True, apply log1p transformation to expression values. Default is False.
            save_path (str, optional): Directory path to save the plot. If empty, the plot is not saved. Default is ''.
            format (str, optional): File format for saving the plot (e.g., 'eps', 'png'). Default is 'eps'.
            dpi (int, optional): Dots per inch for saved figure. Default is 400.
            vmax (float, optional): Percentile value to set the upper limit of the color scale. Default is 99.
        Returns:
            None
        Side Effects:
            - Displays the plot using matplotlib.
            - Saves the plot to the specified path if `save_path` is provided.
        """
        global_matrix = self.get_global_matrix(reverse_x, reverse_y, rotate)
        plt.figure(figsize=figsize)
        sns.scatterplot(
            x=global_matrix.nonzero()[0],
            y=global_matrix.nonzero()[1],
            s=s,
            color="#ced4da",
            edgecolor="none",
        )

        arr = get_exp_array(self.sp.adata, gene)
        arr = _adjust_arr(arr, rotate, reverse_x, reverse_y)
        if log1p:
            arr = np.log1p(arr)
        sparse_matrix = csr_matrix(arr)
        if s is not None:
            ax = sns.scatterplot(
                x=sparse_matrix.nonzero()[0],
                y=sparse_matrix.nonzero()[1],
                c=sparse_matrix.data,
                s=s,
                edgecolor="none",
                vmax=np.percentile(sparse_matrix.data, vmax),
                cmap=cmap,
            )
        else:
            ax = sns.scatterplot(
                x=sparse_matrix.nonzero()[0],
                y=sparse_matrix.nonzero()[1],
                c=sparse_matrix.data,
                edgecolor="none",
                vmax=np.percentile(sparse_matrix.data, vmax),
                cmap=cmap,
            )
        ax.set_axis_off()
        ax.set_title(gene)
        if is_path(save_path):
            fig = ax.get_figure()
            if save_path[-1] != "/":
                save_path += "/"
            save_path += gene
            save_path += "." + format
            fig.savefig(fname=save_path, dpi=dpi, format=format, bbox_inches="tight")
        plt.show()

    def get_global_matrix(self, reverse_x, reverse_y, rotate):

        adata = self.sp.adata
        expression_data = np.array(adata.X.sum(axis=1)).flatten()
        row_indices = np.array(adata.obs["x"].values).flatten()
        column_indices = np.array(adata.obs["y"].values).flatten()
        global_matrix = csr_matrix((expression_data, (row_indices, column_indices)))
        global_matrix = _adjust_arr(
            global_matrix.todense(), rotate, reverse_x, reverse_y
        )
        global_matrix = csr_matrix(global_matrix)
        return global_matrix

    def plot_genes(
        self,
        label=None,
        gene_list=None,
        n_gene=12,
        cmap=None,
        num_cols=4,
        vmax=99,
        vmin=0,
        rotate=False,
        reverse_y=False,
        reverse_x=False,
        plot_type="scatter",
        s=1,
    ):
        result = self.sp.genes_labels
        adata = self.sp.adata
        if gene_list is None:
            if label is None or result is None:
                raise ValueError(
                    "Parameter [label] and [result] should not be None."
                )
            else:
                gene_list = list(result[result["labels"] == label]["gene_id"])[:n_gene]
        genes_count = len(gene_list)
        if genes_count == 0:
            raise ValueError(
                "No genes to plot: [gene_list] is empty or [label] matched no genes."
            )
        axes, fig = _get_figure(genes_count, num_cols)
        fig.subplots_adjust(hspace=0.5)
        for i, gene in enumerate(gene_list):
            row = i // num_cols
            col = i % num_cols
            ax = axes[row, col]
            arr = get_exp_array(adata, gene)
            arr = _adjust_arr(arr, rotate, reverse_x, reverse_y)
            sns.set(style="white")
            if cmap is None:
                cmap = sns.color_palette("viridis", as_cmap=True)
            if plot_type == "heatmap":
                sns.heatmap(
                    arr,
                    cbar=False,
                    ax=ax,
                    cmap=cmap,
                    vmax=np.percentile(arr, vmax),
                    vmin=np.percentile(arr, vmin),
                )
            elif plot_type == "scatter":
                sparse_matrix = csr_matrix(arr)
                sns.scatterplot(
                    x=sparse_matrix.nonzero()[0],
                    y=sparse_matrix.nonzero()[1],
                    ax=ax,
                    c=sparse_matrix.data,
                    cmap=cmap,
                    s=s,
                )
            ax.set_axis_off()
            ax.set_title(gene)
        plt.tight_layout()
        plt.show()

    def plot_pattern(
        self,
        cmap=None,
        vmax=99,
        num_cols=None,
        rotate=False,
        reverse_y=False,
        reverse_x=False,
        heatmap=False,
        s=None,
        image_path=None,
        rotate_img=False,
        k=1,
        aspect=1,
        output_path=None,
        plot_bg=False,
        *,
        labels=None,
        colorbar=True,
        shared_scale=False,
        figsize=None,
        show=True,
        dpi=600,
    ):
        """
        Plot spatial patterns with automatic layout and Nature-style typography.

        Most calls need only ``labels`` and optionally ``output_path``. Existing
        positional arguments and keyword names remain supported.

        Parameters
        ----------
        labels : sequence, optional
            Selected cluster labels, in display order. Defaults to sorted labels
            (first-seen order for mixed types). Requires get_pattern_array().
        heatmap : bool
            False for spots, True for a pixel heatmap. Both use matrix[x, y],
            x increasing right and y increasing up, with equal geometry.
        vmax : float
            Upper color percentile in (0, 100], computed on positive intensities.
            This clips display colors only; the original values are unchanged.
        shared_scale : bool
            Use one pooled color scale across panels (default: per-panel scales).
            Values are composite pattern intensities, not normalized gene counts.
        colorbar : bool
            Display numeric color scales; enabled by default.
        num_cols, s, figsize : optional
            Columns, spot area in points squared, and figure size in inches.
            None selects automatically. One column is 89 mm, otherwise 183 mm.
        cmap : str or Colormap, optional
            Defaults to viridis.
        image_path : str or Path, optional
            Image already aligned to the full untransformed pattern grid. The
            image is stretched to that grid; this is not image registration.
            It follows all subsequent pattern rotations and flips.
        rotate, reverse_x, reverse_y : bool
            Legacy matrix transforms: np.rot90, np.fliplr, and np.flipud,
            respectively (rotation is applied first). With [x, y], rotate
            appears clockwise, reverse_x flips displayed y, reverse_y flips x.
        rotate_img, k : bool, int
            Optional image-only pre-rotation (k quarter-turns) for aligning a
            source image, before applying the common pattern transforms.
        plot_bg : bool
            Show all measured adata.obs x/y grid positions in light gray.
        aspect : positive float
            Displayed y/x unit ratio, applied to every layer; default 1.
        output_path : str or Path, optional
            Export path; suffix selects PDF, SVG, PNG, TIFF, or EPS. A path
            without a suffix gets .pdf. Layout is applied before saving.
        show : bool
            Display the figure (default True); use False in batch workflows.
        dpi : positive number
            Raster export resolution (default 600).

        Returns
        -------
        fig, axes
            Matplotlib Figure and a 2-D axes array, including hidden spare slots.

        Notes
        -----
        Empty patterns are marked 'No signal'. Malformed, negative, or nonfinite
        matrices raise ValueError before a figure is created. Heatmap orientation
        is now consistent with scatter; older heatmaps used row/column axes.
        """
        return plot_patterns(
            self.sp,
            cmap=cmap,
            vmax=vmax,
            num_cols=num_cols,
            rotate=rotate,
            reverse_y=reverse_y,
            reverse_x=reverse_x,
            heatmap=heatmap,
            s=s,
            image_path=image_path,
            rotate_img=rotate_img,
            k=k,
            aspect=aspect,
            output_path=output_path,
            plot_bg=plot_bg,
            labels=labels,
            colorbar=colorbar,
            shared_scale=shared_scale,
            figsize=figsize,
            show=show,
            dpi=dpi,
        )

    def plot_pathways(self, enrichment, **kwargs):
        """Map enriched pattern-specific gene intersections back into space.

        Use sp.enrich_patterns() to obtain enrichment with pattern identities.
        For an existing g:Profiler table, pass pattern_map={query: pattern_label}.
        Select sources, labels, term_ids, top_n (2), p_value_cutoff (0.05), and
        min_genes (2). layer selects expression input; gene_map explicitly maps
        unmatched hit identifiers to var_names. The default display='density'
        renders Gaussian-smoothed spatial densities (bandwidth=1.5 grid units)
        with three relative-peak contours and no spot markers. shared_scale=True
        by default; display='expression' retains the original mean-expression view.
        tissue_mask and support_radius restrict density support. Returned spot
        scores remain unsmoothed, irrespective of the chosen display mode.

        Returns (fig, axes, spot_scores, terms). Scores are mean expression of
        mapped hit genes, not spatial P values or a calibrated activity score.
        Terms record genes used, missing hits, source, pattern, and adjusted P.
        """
        from .pathways import plot_pathways

        return plot_pathways(self.sp, enrichment, **kwargs)

    def plot_intersection(
        self,
        pattern_list,
        cmap=None,
        s=None,
        rotate=False,
        reverse_x=False,
        reverse_y=False,
        figsize=(12, 8),
        image_path=None,
        rotate_img=False,
        plot_bg=True,
        k=1,
        bgs=5,
        aspect=1,
    ):
        sum_array = np.zeros(self.sp.patterns_binary_matrix_dict[pattern_list[0]].shape)
        flag = 1
        for i in pattern_list:
            sum_array += np.where(self.sp.patterns_binary_matrix_dict[i] > 0, flag, 0)
            flag += 1
        sum_array = _adjust_arr(
            sum_array, rotate=rotate, reverse_x=reverse_x, reverse_y=reverse_y
        )
        sparse_matrix = csr_matrix(sum_array)
        plt.figure(figsize=figsize)
        sns.set_style("white")
        if is_path(image_path):
            bg_img = mpimg.imread(image_path)
            if rotate_img:
                bg_img = np.rot90(bg_img, k=k)
            plt.imshow(
                bg_img,
                extent=[0, sum_array.shape[0], 0, sum_array.shape[1]],
                aspect=aspect,
            )
        if plot_bg:
            global_matrix = self.get_global_matrix(reverse_x, reverse_y, rotate)
            sns.scatterplot(
                x=global_matrix.nonzero()[0],
                y=global_matrix.nonzero()[1],
                s=bgs,
                color="#ced4da",
                edgecolor="none",
            )
        default_cmap = ListedColormap(["#06d6a0", "#fb8500", "#ff006e"])
        sns.scatterplot(
            x=sparse_matrix.nonzero()[0],
            y=sparse_matrix.nonzero()[1],
            c=sparse_matrix.data,
            cmap=cmap if cmap is not None else default_cmap,
            s=s if s is not None else 10,
            edgecolor="none",
        )
        plt.axis("off")

    def plot_cluster_score(self, mds_comp, min_cluster, max_cluster):
        """Plot candidate clustering scores without changing analysis results."""
        score_df = evaluate_cluster_scores(
            self.sp.genes_distance_array,
            mds_components=mds_comp,
            min_cluster=min_cluster,
            max_cluster=max_cluster,
        )
        norm_score_df = (score_df - score_df.min()) / (score_df.max() - score_df.min())
        sns.lineplot(norm_score_df, markers=True)
        plt.xticks(list(range(min_cluster, max_cluster + 1, 1)))
        plt.title("Evaluate Clustering Performance")
        plt.xlabel("Number of Clusters")
        plt.ylabel("Normalized Score")
        plt.show()

    def plot_tsne(self, method="tsne", s=10, show_bar=False):
        n_clusters = len(set(self.sp.kmeans_fit_result.labels_))
        if n_clusters <= 10:
            cmap = "tab10"
        elif 10 < n_clusters & n_clusters <= 20:
            cmap = "tab20"
        else:
            cmap = "viridis"
        if method == "tsne":
            tsne = TSNE(n_components=2)
            embedded_data = tsne.fit_transform(self.sp.mds_features)
        else:
            umap_model = umap.UMAP(n_neighbors=5, min_dist=0.3, n_components=2)
            embedded_data = umap_model.fit_transform(self.sp.mds_features)
        plt.figure(figsize=(8, 6))
        scatter = plt.scatter(
            embedded_data[:, 0],
            embedded_data[:, 1],
            c=self.sp.kmeans_fit_result.labels_,
            cmap=cmap,
            s=s,
        )

        plt.xlabel("Dimension 1")
        plt.ylabel("Dimension 2")
        plt.grid(False)
        plt.xticks([])
        plt.yticks([])
        if show_bar:
            plt.colorbar(scatter, label="Labels")
        plt.show()
