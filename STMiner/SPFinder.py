import multiprocessing
from collections import Counter
from typing import Optional

import scanpy as sc
from anndata import AnnData

from STMiner.Algorithm.algorithm import cluster
from STMiner.Algorithm.distance import *
from STMiner.Algorithm.distance import compare_gmm_distance
from STMiner.Algorithm.distribution import get_gmm
from STMiner.Algorithm.distribution import view_gmm, fit_gmms, get_gmm_from_image
from STMiner.IO.IOUtil import merge_bin_coordinate
from STMiner.IO.read_h5ad import read_h5ad
from STMiner.IO.read_stereo import read_gem_file
from STMiner.Plot import Plot


def scale_array(exp_matrix, total_count):
    """
    Scales the given matrix so that the sum of its elements equals 100 and adds the scaled matrix to the total count.

    Args:
    exp_matrix: numpy array, representing the matrix to be scaled.
    total_count: numpy array, representing the cumulative sum of the scaled matrices.

    Returns:
    total_count: The updated cumulative sum after adding the scaled matrix.
    """
    # Calculate the sum of all elements in the matrix
    total_sum = np.sum(exp_matrix)
    # Compute the scaling factor to make the sum of matrix elements equal to 100
    scale_factor = 100 / total_sum
    # Apply the scaling factor to the matrix
    scaled_matrix = exp_matrix * scale_factor
    # Add the scaled matrix to the total count
    total_count += scaled_matrix
    return total_count


class SPFinder:
    def __init__(self, adata: Optional[AnnData] = None):
        self.adata = None
        self.patterns = None
        self.patterns_matrix_dict = {}
        self.patterns_binary_matrix_dict = {}
        self.genes_distance_array = None
        self.filtered_distance_array = None
        self.genes_labels = None
        self.kmeans_fit_result = None
        self.mds_features = None
        self.image_gmm = None
        self.csr_dict = {}
        self.global_distance = None
        self.all_labels = None
        self._gene_expression_edge = {}
        self._highly_variable_genes = []
        self._scope = ()
        self.plot = Plot(self)
        # self.app = App()
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
        return compare_gmm_distance(self.image_gmm, self.patterns)

    def compare_gene_to_genes(self, gene_name):
        gene_gmm = self.patterns[gene_name]
        return compare_gmm_distance(gene_gmm, self.patterns)

    def get_genes_csr_array(self,
                            min_cells,
                            normalize=True,
                            exclude_highly_expressed=False,
                            log1p=False,
                            vmax=99,
                            gene_list=None):
        error_gene_list = []
        self.csr_dict = {}
        self.preprocess(normalize, exclude_highly_expressed, log1p, min_cells)
        if gene_list is not None:
            arr_list = gene_list
        else:
            arr_list = list(self.adata.var.index)
        for gene in tqdm(arr_list, desc='Parsing distance array...'):
            gene_adata = self.adata[:, gene]
            # Sometimes more than 1 genes has same name.
            if gene_adata.var.shape[0] != 1:
                if gene in error_gene_list:
                    pass
                else:
                    error_gene_list.append(gene)
                    print('Gene [' + gene + '] has more than one index, skip.')
                continue
            row_indices = np.array(gene_adata.obs['x'].values).flatten()
            column_indices = np.array(gene_adata.obs['y'].values).flatten()
            try:
                data = np.array(gene_adata.X.todense()).flatten()
            except AttributeError as e:
                data = np.array(gene_adata.X).flatten()
            try:
                data[data > np.percentile(data, vmax)] = np.percentile(data, vmax)
                gene_csr = csr_matrix((data, (row_indices, column_indices)))
                self.csr_dict[gene] = gene_csr
            except Exception as e:
                print('Error when parse gene ' + gene + '\nError: ')
                print(e)

    def spatial_high_variable_genes(self, vmax=99, thread=1):
        if len(self.csr_dict) == 0:
            self.get_genes_csr_array(min_cells=1000, vmax=vmax, normalize=True)
        data = np.array(self.adata.X.sum(axis=1)).flatten()
        data[data > np.percentile(data, vmax)] = np.percentile(data, vmax)
        row_indices = np.array(self.adata.obs['x'].values).flatten()
        column_indices = np.array(self.adata.obs['y'].values).flatten()
        global_matrix = csr_matrix((data, (row_indices, column_indices)))
        # Compare ot distance
        if (not isinstance(thread, int)) or (thread <= 1):
            distance_dict = {}
            for key in tqdm(list(self.csr_dict.keys()), desc='Computing ot distances...'):
                try:
                    distance_dict[key] = calculate_ot_distance(global_matrix, self.csr_dict[key])
                except Exception as e:
                    print(key)
                    print(e)
            self.global_distance = pd.DataFrame(list(distance_dict.items()),
                                                columns=['Gene', 'Distance']).sort_values(by='Distance',
                                                                                          ascending=False)
        else:
            result_dict = multiprocessing.Manager().dict()
            pool = multiprocessing.Pool(processes=thread)
            for key in tqdm(list(self.csr_dict.keys())):
                pool.apply_async(self._mpl_worker, args=(global_matrix, self.csr_dict[key], result_dict))
            pool.close()
            pool.join()
            self.global_distance = dict(result_dict)

    def _mpl_worker(self, global_matrix, key, result_dict):
        result_dict[key] = calculate_ot_distance(global_matrix, self.csr_dict[key])

    def fit_pattern(self,
                    n_top_genes=-1,
                    n_comp=20,
                    normalize=True,
                    exclude_highly_expressed=False,
                    log1p=False,
                    min_cells=20,
                    gene_list=None,
                    remove_low_exp_spots=False):
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
        self.preprocess(normalize, exclude_highly_expressed, log1p, min_cells, n_top_genes)
        if gene_list is not None and isinstance(gene_list, list) and len(gene_list) > 0:
            gene_to_fit = gene_list
        else:
            if n_top_genes > 0:
                gene_to_fit = self._highly_variable_genes
            else:
                gene_to_fit = list(self.adata.var.index)
        try:
            self.patterns = fit_gmms(self.adata,
                                     gene_to_fit,
                                     n_comp=n_comp,
                                     remove_low_exp_spots=remove_low_exp_spots)
        except Exception as e:
            print(e)

    def preprocess(self, normalize, exclude_highly_expressed, log1p, min_cells, n_top_genes=2000):
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
        sc.pp.highly_variable_genes(self.adata,
                                    flavor='seurat_v3',
                                    n_top_genes=n_top_genes)
        self._highly_variable_genes = list(self.adata.var[self.adata.var['highly_variable']].index)
        if normalize:
            sc.pp.normalize_total(self.adata, exclude_highly_expressed=exclude_highly_expressed)
        if log1p:
            sc.pp.log1p(self.adata)

    def build_distance_array(self, method='gmm', gene_list=None):
        """
        Constructs a distance array between genes using the specified method.

        Args:

        1. method (str, optional): The method used to calculate distances. Defaults to 'gmm'.
            - 'gmm': Gaussian Mixture Model distance
            - 'mse': Mean Squared Error distance
            - 'cs': Cosine Similarity
            - 'ot': Optimal Transport distance
        2. gene_list (list, optional): A list of gene names to consider. If None, uses all genes in `self.adata.var.index`. Defaults to None.

        Attributes:
        self.genes_distance_array (numpy array): The constructed distance array based on the chosen method.

        """
        if gene_list is None:
            gene_list = list(self.adata.var.index)

        if method == 'gmm':
            self.genes_distance_array = build_gmm_distance_array(self.patterns)
        elif method == 'mse':
            self.genes_distance_array = build_mse_distance_array(self.adata, gene_list)
        elif method == 'cs':
            self.genes_distance_array = build_cosine_similarity_array(self.adata, gene_list)
        elif method == 'ot':
            self.genes_distance_array = build_ot_distance_array(self.csr_dict, gene_list)
        else:
            raise ValueError(f"Unsupported method '{method}'. Available methods: 'gmm', 'mse', 'cs', 'ot'")

    def get_pattern_array(self, vote_rate=0.2):
        """
        Generates pattern arrays by aggregating expression data based on gene labels, applying a voting scheme,
        and converting to array.

        Args:

        - vote_rate (float, optional): The minimum proportion of genes that must express at a position for it to be considered significant. Defaults to 0.2.

        Attributes Updated:

        - self.patterns_binary_matrix_dict: A dictionary where keys are gene labels and values are binary matrices representing patterns.
        - self.patterns_matrix_dict: A dictionary similar to `patterns_binary_matrix_dict`, but with non-binary aggregated expression.

        Procedure:

        1. Initializes empty dictionaries for storing pattern matrices.
        2. Iterates over unique gene labels.
        3. For each label, collects gene IDs and their expression matrices.
        4. Aggregates expression counts across genes, also tracking positions with non-zero expression.
        5. Applies a voting mechanism to determine significant expression positions based on `vote_rate`.
        6. Converts the scaled aggregate matrix to a binary form based on presence of expression.
        7. Stores the resulting matrices in respective dictionaries.
        """
        self.patterns_binary_matrix_dict = {}
        label_list = set(self.genes_labels['labels'])

        for label in label_list:
            gene_list = list(self.genes_labels[self.genes_labels['labels'] == label]['gene_id'])
            total_count = np.zeros(get_exp_array(self.adata, gene_list[0]).shape)
            total_coo_list = []  # List to collect coordinates of non-zero elements
            vote_array = np.zeros(get_exp_array(self.adata, gene_list[0]).shape)

            for gene in gene_list:
                exp_matrix = get_exp_array(self.adata, gene)
                # Get non-zero element coordinates
                non_zero_coo_list = np.vstack(np.nonzero(exp_matrix)).T.tolist()
                # Collect all non-zero element coordinates
                total_coo_list.extend(non_zero_coo_list)
                # Scale and accumulate expression
                total_count = scale_array(exp_matrix, total_count)

            # Count occurrences of each coordinate and apply voting threshold
            count_dict = Counter(total_coo_list)
            for ele, count in count_dict.items():
                if int(count) / len(gene_list) >= vote_rate:
                    vote_array[tuple(ele)] = 1

            # Apply voting mask, convert to binary, and store results
            total_count *= vote_array
            binary_arr = np.where(total_count != 0, 1, 0)  # Convert to binary
            self.patterns_matrix_dict[label] = total_count
            self.patterns_binary_matrix_dict[label] = binary_arr

    def cluster_gene(self,
                     n_clusters,
                     mds_components=20,
                     use_highly_variable_gene=False,
                     n_top_genes=500,
                     gene_list=None):
        # TODO: genes_labels should be int not float
        if use_highly_variable_gene:
            df = pd.DataFrame(self.genes_distance_array.mean(axis=1), columns=['mean'])
            df = df.sort_values(by='mean', ascending=False)
            hv_gene_list = list(df[:n_top_genes].index)
            self.filtered_distance_array = self.genes_distance_array.loc[hv_gene_list, hv_gene_list]
            self.genes_labels, self.kmeans_fit_result, self.mds_features = cluster(self.filtered_distance_array,
                                                                                   n_clusters=n_clusters,
                                                                                   mds_components=mds_components)
        else:
            self.genes_labels, self.kmeans_fit_result, self.mds_features = cluster(self.genes_distance_array,
                                                                                   n_clusters=n_clusters,
                                                                                   mds_components=mds_components)

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
        all_labels = pd.DataFrame(label, columns=['label'])
        all_labels['value'] = [total.iloc[i, col_index] for i, col_index in enumerate(all_labels['label'].to_list())]
        self.all_labels = all_labels

    # def flush_app(self):
    #     self.app = App()
