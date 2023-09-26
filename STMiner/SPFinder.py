from typing import Optional

import scanpy as sc
from anndata import AnnData

from STMiner.CustomApp.App import App
from STMiner.Plot.plot import Plot
from STMiner.Algorithm.algorithm import cluster
from STMiner.Algorithm.distance import build_gmm_distance_array, compare_gmm_distance
from STMiner.Algorithm.distribution import view_gmm, fit_gmms, get_gmm_from_image
from STMiner.IO.read_h5ad import read_h5ad
from STMiner.IO.IOUtil import merge_bin_coordinate
from STMiner.IO.read_stereo import read_gem_file


class SPFinder:
    def __init__(self, adata: Optional[AnnData] = None):
        self.adata = None
        self.genes_patterns = None
        self.genes_distance_array = None
        self.genes_labels = None
        self.kmeans_fit_result = None
        self.mds_features = None
        self.image_gmm = None
        self._gene_expression_edge = {}
        self._highly_variable_genes = []
        self._scope = ()
        self.plot = Plot(self)
        self.app = App()
        if adata is not None:
            self.set_adata(adata)

    def set_adata(self, adata):
        self.adata = adata
        self._scope = (0, max(adata.obs['y'].max(), adata.obs['x'].max()))

    def read_h5ad(self, file, amplification=1, bin_size=1):
        self.set_adata(read_h5ad(file, amplification=amplification, bin_size=bin_size))

    def read_gem(self, file, bin_size=40):
        self.set_adata(read_gem_file(file, bin_size=bin_size))

    def merge_bin(self, bin_width):
        self.adata.obs['x'] = merge_bin_coordinate(self.adata.obs['x'],
                                                   self.adata.obs['x'].min(),
                                                   bin_size=bin_width)
        self.adata.obs['y'] = merge_bin_coordinate(self.adata.obs['y'],
                                                   self.adata.obs['y'].min(),
                                                   bin_size=bin_width)

    def load_marked_image(self, file):
        self.image_gmm = get_gmm_from_image(file, self.adata)

    def compare_image_to_genes(self):
        """
        Compares the GMM between the marked image and the gene expression matrix.
        :return: pd.DataFrame
        """
        return compare_gmm_distance(self.image_gmm, self.genes_patterns)

    def compare_gene_to_genes(self, gene_name):
        gene_gmm = self.genes_patterns[gene_name]
        return compare_gmm_distance(gene_gmm, self.genes_patterns)

    def fit_pattern(self,
                    n_top_genes,
                    n_comp=20,
                    normalize=True,
                    exclude_highly_expressed=False,
                    log1p=False,
                    min_cells=200):
        """
        Given a distance matrix with the distances between each pair of objects in a set, and a chosen number of
    dimensions, N, an MDS algorithm places each object into N-dimensional space (a lower-dimensional representation)
    such that the between-object distances are preserved as well as possible.
    After that, run **K-Means** clustering and get the labels.

    Ref:
     - https://scikit-learn.org/stable/modules/manifold.html#multidimensional-scaling
     - Multidimensional scaling. (2023, March 28). In Wikipedia. https://en.wikipedia.org/wiki/Multidimensional_scaling
        :param n_top_genes: number of top-genes to fit pattern
        :type n_top_genes: int
        :param n_comp: number of components to fit GMM
        :type n_comp: int
        :param normalize: Run normalize or not, default: True
        :type normalize: boolean
        :param exclude_highly_expressed:
        :type exclude_highly_expressed: boolean
        :param log1p: Run log1p or not, default: False
        :type log1p: boolean
        :param min_cells: minimum number of cells for gene
        :type min_cells: int
        """
        self.preprocess(normalize, exclude_highly_expressed, log1p, min_cells, n_top_genes)
        self.genes_patterns = fit_gmms(self.adata,
                                       self._highly_variable_genes,
                                       n_comp=n_comp)

    def preprocess(self, normalize, exclude_highly_expressed, log1p, min_cells, n_top_genes):
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
        sc.pp.highly_variable_genes(self.adata,
                                    flavor='seurat_v3',
                                    n_top_genes=n_top_genes)
        self._highly_variable_genes = list(self.adata.var[self.adata.var['highly_variable']].index)
        if normalize:
            sc.pp.normalize_total(self.adata, exclude_highly_expressed=exclude_highly_expressed)
        if log1p:
            sc.pp.log1p(self.adata)

    def build_distance_array(self):
        self.genes_distance_array = build_gmm_distance_array(self.genes_patterns)

    def cluster_gene(self, n_clusters, mds_components=20):
        self.genes_labels, self.kmeans_fit_result, self.mds_features = cluster(self.genes_distance_array,
                                                                               n_clusters=n_clusters,
                                                                               mds_components=mds_components)

    def plot_gmm(self, gene_name, cmap=None):
        gmm = self.genes_patterns[gene_name]
        view_gmm(gmm, scope=self._scope, cmap=cmap)

    def flush_app(self):
        self.app = App()
