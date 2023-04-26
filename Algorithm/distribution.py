import anndata
import multiprocessing
import seaborn as sns
import matplotlib.pyplot as plt

from Algorithm.Algorithm import *
from sklearn import mixture


def distribution_distance(gmm1, gmm2):
    gmm1_weights = gmm1.weights_
    gmm1_means = gmm1.means_
    gmm1_covs = gmm1.covariances_
    gmm2_weights = gmm2.weights_
    gmm2_means = gmm2.means_
    gmm2_covs = gmm2.covariances_
    n_components = gmm1_weights.size
    # calculate the distance
    distance_array = np.zeros((n_components, n_components))
    for i in range(n_components):
        for j in range(n_components):
            distance_array[i, j] = get_hellinger_distance(gmm1_covs[i], gmm1_means[i], gmm2_covs[j], gmm2_means[j])
    # TODO: consider the weight of each component
    min_cost = linear_sum(distance_array)
    return min_cost


def get_bh_distance(gmm1_covs, gmm1_means, gmm2_covs, gmm2_means):
    mean_cov = (gmm1_covs + gmm2_covs) / 2
    mean_cov_det = np.linalg.det(mean_cov)
    mean_cov_inv = np.linalg.inv(mean_cov)
    means_diff = gmm1_means - gmm2_means
    first_term = (means_diff.T @ mean_cov_inv @ means_diff) / 8
    second_term = np.log(
        mean_cov_det / (np.sqrt(np.linalg.det(gmm1_covs) * np.linalg.det(gmm1_covs)))) / 2
    result = first_term + second_term
    return result


def get_hellinger_distance(gmm1_covs, gmm1_means, gmm2_covs, gmm2_means):
    mean_cov = (gmm1_covs + gmm2_covs) / 2
    mean_cov_det = np.linalg.det(mean_cov)
    mean_cov_inv = np.linalg.inv(mean_cov)
    means_diff = gmm1_means - gmm2_means
    first_term = np.exp(-(means_diff.T @ mean_cov_inv @ means_diff) / 8)
    second_term = (np.linalg.det(gmm1_covs) ** (1 / 4) * np.linalg.det(gmm1_covs) ** (1 / 4)) / np.sqrt(mean_cov_det)
    result = 1 - first_term * second_term
    return result


def fit_gmm(adata, gene_name, n_comp=5, max_iter=1000):
    adata: anndata
    gene_name: str
    n_comp: int
    max_iter: int

    sample = []
    exp_array = adata[:, adata.var_names == gene_name]
    for index, value in enumerate(exp_array.X.todense()):
        exp_count = int(value * 10)
        if exp_count > 0:
            for i in range(exp_count):
                sample.append([int(exp_array[index].obs.fig_x), int(exp_array[index].obs.fig_y)])

    gmm = mixture.GaussianMixture(n_components=n_comp, max_iter=max_iter)
    gmm.fit(sample)
    return gmm


def fit_gmms(adata: anndata,
             gene_name_list: list,
             n_comp: int = 5,
             max_iter: int = 1000,
             thread: int = 4):
    manager = multiprocessing.Manager()
    shared_dict = manager.dict()
    pool = multiprocessing.Pool(processes=thread)
    for i in gene_name_list:
        pool.apply_async(_fit_worker, args=(shared_dict, adata, i, n_comp, max_iter))
    pool.close()
    pool.join()
    return shared_dict


def _fit_worker(shared_dict, adata, gene_name, n_comp, max_iter):
    shared_dict[gene_name] = fit_gmm(adata, gene_name, n_comp, max_iter)


def view_gmm(gmm, plot_type: str = '3d'):
    x = np.linspace(0, 30000, 100)
    y = np.linspace(0, 30000, 100)
    x_range, y_range = np.meshgrid(x, y)
    x_y = np.column_stack([x_range.flat, y_range.flat])
    # calculate the density
    density = gmm.score_samples(x_y)
    density = density.reshape(x_range.shape)
    if plot_type == '3d':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x_range, y_range, density, cmap='viridis')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('pdf')
        ax.set_box_aspect([2, 2, 1])
        ax.set_title('Probability Density Function Surface')
        # ax.grid(False)
        ax.view_init(elev=30, azim=235)
        plt.show()
    if plot_type == '2d':
        sns.heatmap(np.exp(density))
        plt.show()
