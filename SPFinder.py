from typing import Optional

from anndata import AnnData

from Algorithm.distance import build_gmm_distance_array
from Algorithm.graph import *
from IO.IOUtil import merge_bin_coordinate
from IO.read_10X import read_10x_h5ad
from IO.read_stereo import read_gem_file
from Utils.plot import *


class SPFinder:
    def __init__(self, adata: Optional[AnnData] = None):
        self.adata = None
        self.genes_patterns = None
        self.genes_distance_array = None
        self.genes_labels = None
        self._gene_expression_edge = {}
        self._highly_variable_genes = []
        self._kmeans_fit_result = None
        self._scope = ()
        self._old_adata = None
        if adata is not None:
            self.set_adata(adata)

    def set_adata(self, adata):
        self.adata = adata
        self._scope = (0, max(adata.obs['y'].max(), adata.obs['x'].max()))
        self._old_adata = self.adata.copy()

    def read_10x(self, file, amplification=1, bin_size=1):
        self.set_adata(read_10x_h5ad(file, amplification=amplification, bin_size=bin_size))

    def read_gem(self, file, bin_size=40):
        self.set_adata(read_gem_file(file, bin_size=bin_size))

    def merge_bin(self, bin_width):
        self.adata.obs['x'] = merge_bin_coordinate(self.adata.obs['x'],
                                                   self.adata.obs['x'].min(),
                                                   bin_size=bin_width)
        self.adata.obs['y'] = merge_bin_coordinate(self.adata.obs['y'],
                                                   self.adata.obs['y'].min(),
                                                   bin_size=bin_width)

    def fit_pattern(self,
                    n_top_genes,
                    n_comp,
                    normalize=True,
                    exclude_highly_expressed=False,
                    log1p=False,
                    min_cells=200):
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
        self.genes_labels, self._kmeans_fit_result = cluster(self.genes_distance_array,
                                                             n_clusters=n_clusters,
                                                             mds_components=mds_components)

    def plot_pattern(self, vmax=100):
        plot_pattern(self.genes_labels, self.adata, vmax=vmax)

    def plot_heatmap(self, label, vmax=100, num_cols=4):
        plot_heatmap(result=self.genes_labels,
                     label=label,
                     adata=self.adata,
                     num_cols=num_cols,
                     vmax=vmax)

    def plot_gmm(self, gene_name, cmap=None):
        gmm = self.genes_patterns[gene_name]
        view_gmm(gmm, scope=self._scope, cmap=cmap)

    def clean_up(self):
        self.__init__(adata=self._old_adata)
