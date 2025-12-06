import importlib
import multiprocessing
from collections import Counter
from typing import Optional

import scanpy as sc
from anndata import AnnData
from scipy.stats import zscore
from sklearn import mixture

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
    """
    SPFinder is a class for spatial pattern discovery and analysis in spatial transcriptomics data.
    This class provides methods for reading, preprocessing, and analyzing spatial transcriptomics data,
    including gene expression matrix handling, spatial binning, pattern extraction using Gaussian Mixture Models (GMMs),
    distance calculations (e.g., optimal transport, cosine similarity, mean squared error), clustering,
    and visualization of spatial gene expression patterns.

    - This class is designed for spatial transcriptomics data analysis and requires AnnData and related dependencies.
    - Some methods rely on external utility functions and classes (e.g., Plot, fit_gmms, calculate_ot_distance).
    - Multiprocessing is supported for some distance calculations.
    """

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

    def set_adata(self, adata) -> None:
        """
        Assigns the provided AnnData object to the instance.

        Args:
            adata (AnnData): The annotated data matrix to be set for the instance.
        """
        self.adata = adata        
        position = self.adata.obsm["spatial"]
        x_min = position[:, 0].min()
        y_min = position[:, 1].min()
        self.adata.obs["x"] = merge_bin_coordinate(position[:, 0], x_min, bin_size=1)
        self.adata.obs["y"] = merge_bin_coordinate(position[:, 1], y_min, bin_size=1)
        self._scope = (0, max(self.adata.obs["y"].max(), self.adata.obs["x"].max()))

    def read_h5ad(self, file, amplification=1, bin_size=1, merge_bin=False):
        """
        Reads an h5ad file and sets the object's adata attribute with the loaded data.

        Args:
            file (str): Path to the h5ad file to be read.
            amplification (int, optional): Amplification factor to apply to the data. Defaults to 1.
            bin_size (int, optional): Size of the bin for binning the data. Defaults to 1.
            merge_bin (bool, optional): Whether to merge bins during reading. Defaults to False.

        Returns:
            None
        """
        self.set_adata(
            read_h5ad(
                file,
                amplification=amplification,
                bin_size=bin_size,
                merge_bin=merge_bin,
            )
        )
        self._scope = (0, max(self.adata.obs["y"].max(), self.adata.obs["x"].max()))

    def read_bmk_dir(self, file, bin_size=1):
        self.set_adata(read_bmk(file, bin_size=bin_size))

    def read_gem(self, file, bin_size=40):
        self.set_adata(read_gem_file(file, bin_size=bin_size))

    def merge_bin(self, bin_width):
        """
        Merge spatial coordinates into bins of a specified width.

        This method updates the 'x' and 'y' columns in `self.adata.obs` by grouping their values into bins of size `bin_width`.
        The binning is performed using the `merge_bin_coordinate` function, starting from the minimum value of each coordinate.

        Args:
            bin_width (int or float): The width of each bin to merge coordinates into.

        Returns:
            None

        Notes:
            - Assumes `self.adata.obs` is a pandas DataFrame with 'x' and 'y' columns.
            - The `merge_bin_coordinate` function should accept a coordinate array, a minimum value, and a bin size.
        """
        self.adata.obs["x"] = merge_bin_coordinate(
            self.adata.obs["x"], self.adata.obs["x"].min(), bin_size=bin_width
        )
        self.adata.obs["y"] = merge_bin_coordinate(
            self.adata.obs["y"], self.adata.obs["y"].min(), bin_size=bin_width
        )

    def load_marked_image(self, file, n_components=10):
        self.image_gmm = get_gmm_from_image(file, self.adata, n_components=n_components)

    # def compare_image_to_genes(self):
    #     """
    #     Compares the GMM between the marked image and the gene expression matrix.

    #     Returns:
    #         pd.DataFrame
    #     """
    #     return compare_gmm_distance(self.image_gmm, self.patterns)

    def compare_gene_to_genes(self, gene_name):
        """
        Compares the Gaussian Mixture Model (GMM) of a specified gene to the GMMs of all genes in the current patterns.

        Args:
            gene_name (str): The name of the gene whose GMM will be compared to others.

        Returns:
            dict: A dictionary containing the distances between the specified gene's GMM and the GMMs of all genes in the patterns.
        """
        gene_gmm = self.patterns[gene_name]
        return compare_gmm_distance(gene_gmm, self.patterns)

    def get_genes_csr_array(
        self,
        min_cells: int,
        min_genes: int = 1,
        normalize: bool = True,
        exclude_highly_expressed: bool = False,
        log1p: bool = False,
        vmax: int = 100,
        gene_list: list = None,
    ):
        """
        Generates a dictionary of compressed sparse row (CSR) matrices for gene expression data.
        This method processes the AnnData object (`self.adata`) to extract gene expression arrays
        for each gene, optionally normalizing, excluding highly expressed genes, and applying log1p
        transformation. The resulting matrices are stored in `self.csr_dict` with gene names as keys.

        Args:
            min_cells : int
                Minimum number of cells a gene must be expressed in to be included.
            min_genes : int, optional
                Minimum number of genes a cell must express to be included (default: 1).
            normalize : bool, optional
                Whether to normalize the data before processing (default: True).
            exclude_highly_expressed : bool, optional
                Whether to exclude highly expressed genes (default: False).
            log1p : bool, optional
                Whether to apply log1p transformation to the data (default: False).
            vmax : int, optional
                Percentile value to cap gene expression values (default: 100).
            gene_list : list, optional
                List of gene names to process. If None, all genes in `self.adata` are used (default: None).

        Returns:
            None
        """

        error_gene_list = []
        self.csr_dict = {}
        self.preprocess(
            normalize, exclude_highly_expressed, log1p, min_cells, min_genes
        )
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

    def spatial_high_variable_genes(self, vmax: int = 100, thread: int = 1):
        """
        Identifies spatially high variable genes by comparing each gene's spatial expression pattern
        to a global expression matrix using optimal transport (OT) distance.
        This method processes the spatial transcriptomics data to create a global sparse matrix,
        then computes the OT distance between the global matrix and each gene-specific matrix.
        The results are stored in a DataFrame with gene names, distances, and z-scores of the
        log-transformed distances. Optionally, multiprocessing can be used to speed up computation.

        Args:
            vmax: int, optional
                The upper percentile threshold for capping expression values in the global matrix (default is 100).
            thread: int, optional
                The number of threads to use for multiprocessing. If thread <= 1, computation is done serially (default is 1).

        Returns:
            None

        Notes:
            - Requires `self.csr_dict` to be populated; otherwise, it will be generated.
            - Uses tqdm for progress display and multiprocessing for parallel computation if `thread > 1`.
            - Handles exceptions during OT distance calculation and prints the gene key and error.
            - The results are stored in the `self.global_distance` attribute as a pandas DataFrame with columns: "Gene": Gene names. | "Distance": OT distance to the global matrix. | "z-score": Z-score of the log-transformed distances.
        """

        if len(self.csr_dict) == 0:
            self.get_genes_csr_array(min_cells=50, vmax=vmax, normalize=True)
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
            print("Using multiprocessing, thread is {thread}.".format(thread=thread))
            with multiprocessing.Pool(processes=thread) as pool:
                results = pool.starmap(
                    self._mpl_worker,
                    [
                        (global_matrix, key, self.csr_dict)
                        for key in self.csr_dict.keys()
                    ],
                )

            normal_dict = dict(results)
            self.global_distance = pd.DataFrame(
                list(normal_dict.items()), columns=["Gene", "Distance"]
            ).sort_values(by="Distance", ascending=False)
            self.global_distance["z-score"] = zscore(
                np.log1p(self.global_distance["Distance"])
            )

    def _mpl_worker(self, global_matrix, key, csr_dict):
        res = calculate_ot_distance(global_matrix, csr_dict[key])
        return key, res

    def fit_pattern(
        self,
        n_top_genes: int = -1,
        n_comp: int = 20,
        normalize: bool = False,
        exclude_highly_expressed: bool = False,
        log1p: bool = False,
        min_cells: int = 20,
        gene_list: list = None,
        remove_low_exp_spots: bool = False,
    ):
        """
        Fits gene expression patterns using Gaussian Mixture Models (GMMs) on selected genes.
        This method preprocesses the data and fits GMMs to the expression profiles of specified genes,
        allowing for the identification of spatial patterns in gene expression.

        Args:
            n_top_genes : int, optional (default: -1)
                Number of top highly variable genes to use. If -1, use all genes.
            n_comp : int, optional (default: 20)
                Number of GMM components to fit for each gene.
            normalize : bool, optional (default: False)
                Whether to normalize the data before fitting.
            exclude_highly_expressed : bool, optional (default: False)
                Whether to exclude highly expressed genes during preprocessing.
            log1p : bool, optional (default: False)
                Whether to apply log1p transformation to the data.
            min_cells : int, optional (default: 20)
                Minimum number of cells a gene must be expressed in to be included.
            gene_list : list, optional (default: None)
                List of gene names to fit patterns for. If None, uses top genes or all genes.
            remove_low_exp_spots : bool, optional (default: False)
                Whether to remove spots with low expression before fitting.

        Notes:
            The fitted patterns are stored in the `self.patterns` attribute.
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
        min_genes: int,
        n_top_genes: int = 2000,
    ):
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
        sc.pp.filter_cells(self.adata, min_genes=min_genes)
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

        Args:
            method (str): The distance calculation method to use. Default is "gmm".
            gene_list (list): A list of specific genes to use. If not provided, all genes are used.

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
            self.genes_distance_array = build_cosine_similarity_array(
                self.adata, gene_list
            )
        elif method == "ot":
            self.genes_distance_array = build_ot_distance_array(
                self.csr_dict, gene_list
            )
        else:
            # Raise an error if the method is unknown
            raise ValueError("Unknown method, method should be one of gmm, mse, cs, ot")

    def get_pattern_array(self, vote_rate: int = 0, mode: str = "vote"):
        self.patterns_binary_matrix_dict = {}
        if mode == "vote":
            label_list = set(self.genes_labels["labels"])
            for label in label_list:
                gene_list = list(
                    self.genes_labels[self.genes_labels["labels"] == label]["gene_id"]
                )
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
        # Initialize matrices
        exp_shape = get_exp_array(self.adata, gene_list[0]).shape
        total_coo_list = []
        total_count = np.zeros(exp_shape)
        vote_array = np.zeros(exp_shape)

        # Efficiently accumulate coordinates and counts
        for gene in tqdm(gene_list, desc="Accumulating gene expression..."):
            exp_matrix = get_exp_array(self.adata, gene)
            # Directly collect non-zero coordinates in a more efficient manner
            non_zero_coords = np.nonzero(exp_matrix)
            total_coo_list.extend(zip(*non_zero_coords))

            # Scale and accumulate the exp_matrix
            total_count += scale_array(exp_matrix, total_count)

        # Count occurrences of each coordinate
        count_dict = Counter(total_coo_list)

        # Vote-based update
        for ele, count in count_dict.items():
            if count / len(gene_list) >= vote_rate:
                vote_array[ele] = 1

        # Apply vote mask to the total count
        total_count *= vote_array

        # Return binary array and scaled total count
        binary_arr = np.where(total_count != 0, 1, total_count)
        return binary_arr, total_count

    def get_custom_pattern(
        self, gene_list, n_components=20, vote_rate: int = 0, mode: str = "vote"
    ):
        """
        Generates a custom pattern model based on a list of genes using either a voting mechanism or a test mode.

        Args:
            gene_list (list): List of gene identifiers to be used for pattern extraction.
            n_components (int, optional): Number of components for the Gaussian Mixture Model (GMM). Defaults to 20.
            vote_rate (int, optional): Threshold for voting mechanism in pattern extraction. Defaults to 0.
            mode (str, optional): Mode of operation, either "vote" for GMM-based pattern extraction or "test" for statistical testing. Defaults to "vote".

        Raises:
            ValueError: If the mode is not "vote" or "test".

        Notes:
            - In "vote" mode, fits a Gaussian Mixture Model to the gene pattern data.
            - In "test" mode, statistical testing is intended but not yet implemented.
        """
        if mode == "vote":
            _, total_count = self._genes_to_pattern(gene_list, vote_rate)
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
        n_top_genes=500,
    ):
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
        """
        Plots the Gaussian Mixture Model (GMM) for a specified gene.

        Args:
            gene_name (str): The name of the gene whose GMM is to be plotted.
            cmap (str or matplotlib.colors.Colormap, optional): Colormap to use for plotting. Defaults to None.

        Returns:
            None
        """

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
        for gene in gene_list:
            if gene in list(self.adata.var.index):
                _genes.append(gene)

        # Get expression patterns of interested gene set
        self.get_custom_pattern(gene_list=_genes, n_components=n_comp, vote_rate=0)

    
    def gui(self):
        """Launch a Gradio interface for the main SPFinder workflow."""

        if importlib.util.find_spec("gradio") is None:
            raise ImportError(
                "Gradio is required for the GUI. Please install it with `pip install gradio`."
            )

        import gradio as gr

        def _check_adata_loaded():
            if self.adata is None:
                raise ValueError("Please load spatial transcriptomics data first.")

        def load_data(file, file_path, amplification, bin_size, merge_bin):
            try:
                path = None
                if file is not None:
                    path = file.name
                elif file_path:
                    path = file_path
                if path is None:
                    return "Please provide an h5ad file path or upload a file."
                amplification = int(amplification)
                bin_size = int(bin_size)
                self.read_h5ad(
                    file=path,
                    amplification=amplification,
                    bin_size=bin_size,
                    merge_bin=merge_bin,
                )
                obs_num, var_num = self.adata.n_obs, self.adata.n_vars
                return f"Loaded data from {path}. Spots: {obs_num}, Genes: {var_num}."
            except Exception as exc:
                return f"Failed to load data: {exc}"

        def find_spatial_hvg(min_cells, min_genes, vmax, thread):
            try:
                _check_adata_loaded()
                min_cells = int(min_cells)
                min_genes = int(min_genes)
                vmax = int(vmax)
                thread = int(thread)
                self.get_genes_csr_array(
                    min_cells=min_cells,
                    min_genes=min_genes,
                    vmax=vmax,
                )
                self.spatial_high_variable_genes(vmax=vmax, thread=thread)
                preview = self.global_distance.head(20)
                return "Computed spatially variable genes.", preview
            except Exception as exc:
                return f"Failed to compute: {exc}", None

        def fit_patterns(n_top_genes, n_comp, min_cells, gene_list_text):
            try:
                _check_adata_loaded()
                n_top_genes = int(n_top_genes)
                n_comp = int(n_comp)
                min_cells = int(min_cells)
                genes = None
                if gene_list_text:
                    genes = [g.strip() for g in gene_list_text.split("\n") if g.strip()]
                self.fit_pattern(
                    n_top_genes=n_top_genes,
                    n_comp=n_comp,
                    min_cells=min_cells,
                    gene_list=genes,
                )
                return "Finished fitting GMM patterns."
            except Exception as exc:
                return f"Failed to fit patterns: {exc}"

        def cluster_genes(method, n_clusters, mds_components, use_highly_variable_gene, n_top_genes, gene_list_text):
            try:
                _check_adata_loaded()
                n_clusters = int(n_clusters)
                mds_components = int(mds_components)
                n_top_genes = int(n_top_genes)
                genes = None
                if gene_list_text:
                    genes = [g.strip() for g in gene_list_text.split("\n") if g.strip()]
                self.build_distance_array(method=method, gene_list=genes)
                self.cluster_gene(
                    n_clusters=n_clusters,
                    mds_components=mds_components,
                    use_highly_variable_gene=use_highly_variable_gene,
                    n_top_genes=n_top_genes,
                )
                preview = None
                if self.genes_labels is not None:
                    preview = self.genes_labels.head(20)
                return "Clustering completed.", preview
            except Exception as exc:
                return f"Failed to cluster genes: {exc}", None

        with gr.Blocks(title="STMiner SPFinder GUI") as demo:
            gr.Markdown(
                """
                # STMiner.SPFinder GUI

                load spatial data, find spatially variable genes,
                fit Gaussian Mixture Models, build distance arrays, and cluster genes.
                """
            )

            with gr.Tab("Step.1-Load data"):
                gr.Markdown(
                    "Upload or provide an h5ad path, then configure binning options to load the dataset."
                )
                with gr.Row():
                    file_input = gr.File(label="Upload h5ad", file_types=[".h5ad"])
                    path_input = gr.Textbox(label="Or enter h5ad path", placeholder="/path/to/file.h5ad")
                with gr.Row():
                    amplification = gr.Number(label="Amplification", value=1, precision=0)
                    bin_size = gr.Number(label="Bin size", value=1, precision=0)
                    merge_bin = gr.Checkbox(label="Merge bin", value=False)
                load_btn = gr.Button("Load dataset")
                load_status = gr.Textbox(label="Status")
                load_btn.click(
                    load_data,
                    inputs=[file_input, path_input, amplification, bin_size, merge_bin],
                    outputs=load_status,
                )

            with gr.Tab("Step.2-Spatial HVG"):
                gr.Markdown("Compute spatially variable genes following the README workflow.")
                with gr.Row():
                    min_cells = gr.Number(label="Min cells", value=100, precision=0)
                    min_genes = gr.Number(label="Min genes", value=1, precision=0)
                    vmax = gr.Number(label="Vmax percentile", value=100, precision=0)
                    thread = gr.Number(label="Threads", value=1, precision=0)
                hvg_btn = gr.Button("Find spatial HVGs")
                hvg_status = gr.Textbox(label="Status")
                hvg_preview = gr.DataFrame(label="Top spatially variable genes", interactive=False)
                hvg_btn.click(
                    find_spatial_hvg,
                    inputs=[min_cells, min_genes, vmax, thread],
                    outputs=[hvg_status, hvg_preview],
                )

            with gr.Tab("Step.3-Fit patterns"):
                gr.Markdown("Fit GMM patterns for genes.")
                with gr.Row():
                    n_top_genes = gr.Number(label="Top highly variable genes (-1 for all)", value=2000, precision=0)
                    n_comp = gr.Number(label="GMM components", value=20, precision=0)
                    min_cells_fit = gr.Number(label="Min cells", value=20, precision=0)
                gene_list_fit = gr.Textbox(
                    label="Optional gene list (one per line)",
                    placeholder="GeneA\nGeneB",
                )
                fit_btn = gr.Button("Fit patterns")
                fit_status = gr.Textbox(label="Status")
                fit_btn.click(
                    fit_patterns,
                    inputs=[n_top_genes, n_comp, min_cells_fit, gene_list_fit],
                    outputs=fit_status,
                )

            with gr.Tab("Step.4-Distance & clustering"):
                gr.Markdown("Build distance arrays and cluster genes into spatial patterns.")
                with gr.Row():
                    method = gr.Radio(
                        label="Distance method",
                        choices=["gmm", "mse", "cs", "ot"],
                        value="gmm",
                    )
                    n_clusters = gr.Number(label="Number of clusters", value=6, precision=0)
                    mds_components = gr.Number(label="MDS components", value=20, precision=0)
                use_hvg = gr.Checkbox(label="Use highly variable genes for clustering", value=False)
                n_top_cluster = gr.Number(label="Top genes when using HVG", value=500, precision=0)
                gene_list_cluster = gr.Textbox(
                    label="Optional gene list for distance array (one per line)",
                    placeholder="GeneA\nGeneB",
                )
                cluster_btn = gr.Button("Build & cluster")
                cluster_status = gr.Textbox(label="Status")
                cluster_preview = gr.DataFrame(label="Gene labels", interactive=False)
                cluster_btn.click(
                    cluster_genes,
                    inputs=[
                        method,
                        n_clusters,
                        mds_components,
                        use_hvg,
                        n_top_cluster,
                        gene_list_cluster,
                    ],
                    outputs=[cluster_status, cluster_preview],
                )

        demo.launch()
