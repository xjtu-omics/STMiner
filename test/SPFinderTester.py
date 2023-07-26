from SPFinder import SPFinder
from test.testUtils import *
from Algorithm.graph import *
from scipy import sparse


class SPFinderTester(SPFinder):
    def __init__(self):
        super().__init__()
        self.noise_dict = {}
        self.noise_type = 'uniform'
        self.noise_argument = 1

    def set_noise_type(self, noise_type, noise_argument):
        self.noise_type = noise_type
        self.noise_argument = noise_argument

    def fit_pattern(self, n_top_genes, n_comp):
        sc.pp.highly_variable_genes(self.adata,
                                    flavor='seurat_v3',
                                    n_top_genes=n_top_genes)
        self._highly_variable_genes = list(self.adata.var[self.adata.var['highly_variable']].index)
        self.genes_patterns = self.fit_gmms(
            self._highly_variable_genes,
            n_comp=n_comp)

    def fit_gmms(self,

                 gene_name_list,
                 n_comp=5,
                 binary=False,
                 threshold=95,
                 max_iter=100,
                 reg_covar=1e-3,
                 cut=False):
        gmm_dict = {}
        dropped_genes_count = 0
        for gene_id in tqdm(gene_name_list, desc='Fitting GMM...'):
            try:
                fit_result = self.fit_gmm(
                    gene_id,
                    n_comp=n_comp,
                    binary=binary,
                    threshold=threshold,
                    max_iter=max_iter,
                    reg_covar=reg_covar,
                    cut=cut)
                if fit_result is not None:
                    gmm_dict[gene_id] = fit_result
                else:
                    dropped_genes_count += 1
            except Exception as e:
                error_msg = str(e)
                raise ValueError("Gene id: " + gene_id + "\nError: " + error_msg)
        if dropped_genes_count > 0:
            print('Number of dropped genes: ' + str(dropped_genes_count))
        return gmm_dict

    def fit_gmm(self,
                gene_name: str,
                n_comp: int = 2,
                cut: bool = False,
                binary: bool = False,
                threshold: int = 95,
                max_iter: int = 200,
                reg_covar=1e-3):
        result = self.preprocess_array(binary, cut, gene_name, threshold)
        # Add noise
        if self.noise_type == 'gauss':
            noise = get_gauss_noise(result, self.noise_argument)
            noised = result + noise
        elif self.noise_type == 'periodicity':
            noise = get_periodicity_noise(result, self.noise_argument)
            noised = result * noise
        elif self.noise_type == 'sp':
            noise = get_salt_pepper_noise(result, self.noise_argument)
            noised = result * noise
        else:
            noise = get_uniform_noise(result, self.noise_argument)
            noised = result + noise
        self.noise_dict[gene_name] = noise
        # Number of unique center must be larger than the number of components.
        if len(set(map(tuple, noised))) >= n_comp:
            gmm = mixture.GaussianMixture(n_components=n_comp, max_iter=max_iter, reg_covar=reg_covar)
            gmm.fit(noised)
            return gmm
        else:
            return None

    def preprocess_array(self, binary, cut, gene_name, threshold):
        dense_array = self.get_exp_array(gene_name)
        if binary:
            if threshold > 100 | threshold < 0:
                print('Warning: the threshold is illegal, the value in [0, 100] is accepted.')
                threshold = 100
            binary_arr = np.where(dense_array > np.percentile(dense_array, threshold), 1, 0)
            result = array_to_list(binary_arr)
        else:
            if cut:
                dense_array[dense_array < np.percentile(dense_array, threshold)] = 0
            result = array_to_list(dense_array)
        return result

    def get_exp_array(self, gene_name):
        exp_array = self.adata[:, self.adata.var_names == gene_name].X
        if sparse.issparse(exp_array):
            data = np.array(exp_array.todense())
        else:
            data = np.array(exp_array)
        sparse_matrix = sparse.coo_matrix((data[:, 0], (np.array(self.adata.obs['x']), np.array(self.adata.obs['y']))))
        dense_array = np.array(np.round(sparse_matrix.todense()), dtype=np.int32)
        return dense_array
