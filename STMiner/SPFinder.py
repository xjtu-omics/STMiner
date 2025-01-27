import multiprocessing
from collections import Counter
from typing import Optional

import numpy as np
import scanpy as sc
from anndata import AnnData
from scipy.stats import zscore

from STMiner.Algorithm.algorithm import cluster
from STMiner.Algorithm.distance import *
from STMiner.Algorithm.distance import compare_gmm_distance
from STMiner.Algorithm.distribution import get_gmm, array_to_list
from STMiner.Algorithm.distribution import view_gmm, fit_gmms, get_gmm_from_image
from STMiner.IO.IOUtil import merge_bin_coordinate
from STMiner.IO.read_bmk import read_bmk
from STMiner.IO.read_h5ad import read_h5ad
from STMiner.IO.read_stereo import read_gem_file
from STMiner.Plot import Plot


def scale_array(exp_matrix, total_count):
    total_sum = np.sum(exp_matrix)
    scale_factor = 100 / total_sum
    scaled_matrix = exp_matrix * scale_factor
    total_count += scaled_matrix
    return total_count


class SPFinder:
    def __init__(self, adata: Optional[AnnData] = None):
        self.adata = None
        self.patterns = None
        self.genes_distance_array = None
        self.filtered_distance_array = None
        self.genes_labels = None
        self.kmeans_fit_result = None
        self.mds_features = None
        self.image_gmm = None
        self.global_distance = None
        self.all_labels = None
        self.custom_pattern = None
        self.csr_dict = {}
        self.patterns_matrix_dict = {}
        self.patterns_binary_matrix_dict = {}
        self.plot = Plot(self)
        self._highly_variable_genes = []
        self._gene_expression_edge = {}
        self._scope = ()
        # self.app = App()
        if adata is not None:
            self.set_adata(adata)

    def set_adata(self, adata):
        self.adata = adata
        self._scope = (0, max(adata.obs["y"].max(), adata.obs["x"].max()))

    def read_h5ad(self, file, amplification=1, bin_size=1):
        self.set_adata(read_h5ad(file, amplification=amplification, bin_size=bin_size))

    def read_bmk_dir(self, file, bin_size=1):
        self.set_adata(read_bmk(file, bin_size=bin_size))

    def read_gem(self, file, bin_size=40):
        self.set_adata(read_gem_file(file, bin_size=bin_size))

    def merge_bin(self, bin_width):
        self.adata.obs["x"] = merge_bin_coordinate(
            self.adata.obs["x"], self.adata.obs["x"].min(), bin_size=bin_width
        )
        self.adata.obs["y"] = merge_bin_coordinate(
            self.adata.obs["y"], self.adata.obs["y"].min(), bin_size=bin_width
        )

    def load_marked_image(self, file):
        self.image_gmm = get_gmm_from_image(file, self.adata)

    def compare_image_to_genes(self):
        """
        Compares the GMM between the marked image and the gene expression matrix.
        :return: pd.DataFrame
        """
        return compare_gmm_distance(self.image_gmm, self.patterns)

    def compare_gene_to_genes(self, gene_name):
        gene_gmm = self.patterns[gene_name]
        return compare_gmm_distance(gene_gmm, self.patterns)

    def get_genes_csr_array(
            self,
            min_cells: int,
            normalize: bool = True,
            exclude_highly_expressed: bool = False,
            log1p: bool = False,
            vmax: int = 99,
            gene_list: list = None,
    ):
        error_gene_list = []
        self.csr_dict = {}
        self.preprocess(normalize, exclude_highly_expressed, log1p, min_cells)
        if gene_list is not None:
            arr_list = gene_list
        else:
            arr_list = list(self.adata.var.index)
        for gene in tqdm(arr_list, desc="Parsing distance array..."):
            gene_adata = self.adata[:, gene]
            # Sometimes more than 1 genes has same name.
            if gene_adata.var.shape[0] != 1:
                if gene in error_gene_list:
                    pass
                else:
                    error_gene_list.append(gene)
                    print("Gene [" + gene + "] has more than one index, skip.")
                continue
            row_indices = np.array(gene_adata.obs["x"].values).flatten()
            column_indices = np.array(gene_adata.obs["y"].values).flatten()
            try:
                data = np.array(gene_adata.X.todense()).flatten()
            except AttributeError as e:
                data = np.array(gene_adata.X).flatten()
            try:
                if np.percentile(data, vmax) > 0:
                    data[data > np.percentile(data, vmax)] = np.percentile(data, vmax)
                gene_csr = csr_matrix((data, (row_indices, column_indices)))
                self.csr_dict[gene] = gene_csr
            except Exception as e:
                print("Error when parse gene " + gene + "\nError: ")
                print(e)

    def spatial_high_variable_genes(self, vmax: int = 99, thread: int = 1):
        """
        Compute the optimal transport (OT) distance matrix for high variable genes.

        :param vmax: The percentile threshold to clip data values.
        :type vmax: int, default 99
        :param thread: Number of threads to use for parallel processing. If <= 1, single-threaded execution is used.
        :type thread: int, default 1

        This function first ensures that `csr_dict` is populated. It then creates a global matrix using the sum of
        expression values across cells and computes OT distances between this matrix and each gene's CSR matrix from `csr_dict`.

        The resulting OT distances are stored in a DataFrame sorted by increasing distance,
        with genes and their corresponding distances as columns.

        Note: Exceptions during calculation are logged with the gene key and error message.
        """
        if len(self.csr_dict) == 0:
            self.get_genes_csr_array(min_cells=1000, vmax=vmax, normalize=True)
        # Process data and create global sparse matrix
        data = np.array(self.adata.X.sum(axis=1)).flatten()
        data[data > np.percentile(data, vmax)] = np.percentile(data, vmax)
        row_indices = np.array(self.adata.obs["x"].values).flatten()
        column_indices = np.array(self.adata.obs["y"].values).flatten()
        global_matrix = csr_matrix((data, (row_indices, column_indices)))
        # Compare ot distance
        if (not isinstance(thread, int)) or (thread <= 1):
            distance_dict = {}
            for key in tqdm(
                    list(self.csr_dict.keys()), desc="Computing ot distances..."
            ):
                try:
                    distance_dict[key] = calculate_ot_distance(
                        global_matrix, self.csr_dict[key]
                    )
                except Exception as e:
                    print(key)
                    print(e)
            self.global_distance = pd.DataFrame(
                list(distance_dict.items()), columns=["Gene", "Distance"]
            ).sort_values(by="Distance", ascending=False)
            # Calculating the z-score for the log-transformed values of the
            # 'Distance' column in the 'global_distance' DataFrame using the zscore function from the
            # scipy.stats library.
            self.global_distance["z-score"] = zscore(
                np.log1p(self.global_distance["Distance"])
            )
        else:
            result_dict = multiprocessing.Manager().dict()
            pool = multiprocessing.Pool(processes=thread)
            for key in tqdm(list(self.csr_dict.keys())):
                pool.apply_async(
                    self._mpl_worker, args=(global_matrix, key, result_dict)
                )
            pool.close()
            pool.join()
            self.global_distance = pd.DataFrame(dict(result_dict), index=[0]).T
            self.global_distance["z-score"] = zscore(
                np.log1p(self.global_distance["Distance"])
            )

    def _mpl_worker(self, global_matrix, key, result_dict):
        res = calculate_ot_distance(global_matrix, self.csr_dict[key])
        result_dict[key] = res

    def fit_pattern(
            self,
            n_top_genes: int = -1,
            n_comp: int = 20,
            normalize: bool = True,
            exclude_highly_expressed: bool = False,
            log1p: bool = False,
            min_cells: int = 20,
            gene_list: list = None,
            remove_low_exp_spots: bool = False,
    ):
        """
            Given a distance matrix with the distances between each pair of objects in a set, and a chosen number of
        dimensions, N, an MDS algorithm places each object into N-dimensional space (a lower-dimensional representation)
        such that the between-object distances are preserved as well as possible.
        After that, run **K-Means** clustering and get the labels.

        Ref:
         - https://scikit-learn.org/stable/modules/manifold.html#multidimensional-scaling
         - Multidimensional scaling. (2023, March 28). In Wikipedia. https://en.wikipedia.org/wiki/Multidimensional_scaling
            :param n_top_genes: number of top high variable genes to fit pattern, if n_top_genes <= 0, fit all the genes.
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
            :param gene_list:
            :type min_cells: list
            :param remove_low_exp_spots:
        """
        self.preprocess(
            normalize, exclude_highly_expressed, log1p, min_cells, n_top_genes
        )
        if gene_list is not None and isinstance(gene_list, list) and len(gene_list) > 0:
            gene_to_fit = gene_list
        else:
            if n_top_genes > 0:
                gene_to_fit = self._highly_variable_genes
            else:
                gene_to_fit = list(self.adata.var.index)
        try:
            self.patterns = fit_gmms(
                self.adata,
                gene_to_fit,
                n_comp=n_comp,
                remove_low_exp_spots=remove_low_exp_spots,
            )
        except Exception as e:
            print(e)

    def preprocess(
            self,
            normalize: bool,
            exclude_highly_expressed: bool,
            log1p: bool,
            min_cells: int,
            n_top_genes: int = 2000
    ):
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
        sc.pp.highly_variable_genes(
            self.adata, flavor="seurat_v3", n_top_genes=n_top_genes
        )
        self._highly_variable_genes = list(
            self.adata.var[self.adata.var["highly_variable"]].index
        )
        if normalize:
            sc.pp.normalize_total(
                self.adata, exclude_highly_expressed=exclude_highly_expressed
            )
        if log1p:
            sc.pp.log1p(self.adata)

    def build_distance_array(self, method="gmm", gene_list=None):
        """
        Build a distance array for genes based on the specified method.

        This function supports four distance calculation methods: "gmm" (Gaussian Mixture Model),
        "mse" (Mean Squared Error), "cs" (Cosine Similarity), and "ot" (Optimal Transport).
        If no gene list is provided, all genes are used.

        Parameters:
        - method (str): The distance calculation method to use. Default is "gmm".
        - gene_list (list): A list of specific genes to use. If not provided, all genes are used.

        Returns:
        No direct return value, but updates the `self.genes_distance_array` attribute with the calculated distances.
        """
        # If no gene list is provided, use all genes
        if gene_list is None:
            gene_list = list(self.adata.var.index)

        # Build the distance array based on the specified method
        if method == "gmm":
            self.genes_distance_array = build_gmm_distance_array(self.patterns)
        elif method == "mse":
            self.genes_distance_array = build_mse_distance_array(self.adata, gene_list)
        elif method == "cs":
            self.genes_distance_array = build_cosine_similarity_array(self.adata, gene_list)
        elif method == "ot":
            self.genes_distance_array = build_ot_distance_array(self.csr_dict, gene_list)
        else:
            # Raise an error if the method is unknown
            raise ValueError("Unknown method, method should be one of gmm, mse, cs, ot")

    def get_pattern_array(self, vote_rate: int = 0, mode: str = "vote"):
        self.patterns_binary_matrix_dict = {}
        if mode == "vote":
            label_list = set(self.genes_labels["labels"])
            for label in label_list:
                gene_list = list(self.genes_labels[self.genes_labels["labels"] == label]["gene_id"])
                binary_arr, total_count = self._genes_to_pattern(gene_list, vote_rate)
                self.patterns_matrix_dict[label] = total_count
                self.patterns_binary_matrix_dict[label] = binary_arr
        elif mode == "test":
            p_value_threshold = 0.05
            # TODO: rewrite test mode, improve run time.
            pass
        else:
            raise ValueError("mode should be vote or test")

    def _genes_to_pattern(self, gene_list, vote_rate):

        total_count = np.zeros(get_exp_array(self.adata, gene_list[0]).shape)
        total_coo_list = []
        vote_array = np.zeros(get_exp_array(self.adata, gene_list[0]).shape)
        for gene in gene_list:
            exp_matrix = get_exp_array(self.adata, gene)
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
        binary_arr = np.where(total_count != 0, 1, total_count)
        return binary_arr, total_count

    def get_custom_pattern(self, gene_list, n_components=20, vote_rate: int = 0, mode: str = "vote"):
        if mode == "vote":
            _, total_count = self._genes_to_pattern(gene_list, vote_rate)
            from sklearn import mixture
            _gmm = mixture.GaussianMixture(n_components=n_components)
            _gmm.fit(array_to_list(np.round(total_count).astype(np.int32)))
            self.custom_pattern = _gmm
        elif mode == "test":
            p_value_threshold = 0.05
            # TODO: rewrite test mode, improve run time.
            pass
        else:
            raise ValueError("mode should be vote or test")

    def cluster_gene(
            self,
            n_clusters: int,
            mds_components=20,
            use_highly_variable_gene=False,
            n_top_genes=500):
        # TODO: genes_labels should be int not float
        if use_highly_variable_gene:
            df = pd.DataFrame(self.genes_distance_array.mean(axis=1), columns=["mean"])
            df = df.sort_values(by="mean", ascending=False)
            hv_gene_list = list(df[:n_top_genes].index)
            self.filtered_distance_array = self.genes_distance_array.loc[
                hv_gene_list, hv_gene_list
            ]
            self.genes_labels, self.kmeans_fit_result, self.mds_features = cluster(
                self.filtered_distance_array,
                n_clusters=n_clusters,
                mds_components=mds_components,
            )
        else:
            self.genes_labels, self.kmeans_fit_result, self.mds_features = cluster(
                self.genes_distance_array,
                n_clusters=n_clusters,
                mds_components=mds_components,
            )

    def plot_gmm(self, gene_name, cmap=None):
        gmm = self.patterns[gene_name]
        view_gmm(gmm, scope=self._scope, cmap=cmap)

    def get_all_labels(self):
        df_list = []
        for i in self.patterns_binary_matrix_dict:
            gmm = get_gmm(self.patterns_binary_matrix_dict[i], n_comp=20)
            df = compare_gmm_distance(gmm, self.patterns)
            df.rename(columns={0: i}, inplace=True)
            df_list.append(df)
        total = pd.concat(df_list, axis=1)
        label = total.idxmin(axis=1)
        all_labels = pd.DataFrame(label, columns=["label"])
        all_labels["value"] = [
            total.iloc[i, col_index]
            for i, col_index in enumerate(all_labels["label"].to_list())
        ]
        self.all_labels = all_labels

    def get_pattern_of_given_genes(self, gene_list, n_comp=20):
        _genes = []
        if self.adata is None:
            raise ValueError("Please load ST data first.")
        for i in gene_list:
            if i in list(self.adata.var.index):
                _genes.append(i)

        # Get expression patterns of interested gene set
        self.fit_pattern(n_comp=n_comp, gene_list=_genes)
        self.cluster_gene(n_clusters=1, mds_components=2)
        self.patterns_matrix_dict = None
        self.get_custom_pattern(gene_list=gene_list, n_components=n_comp, vote_rate=0)

    # def flush_app(self):
    #     self.app = App()
