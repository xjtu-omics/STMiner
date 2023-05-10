import numpy as np

from sklearn.manifold import MDS
from sklearn.cluster import KMeans
from scipy.optimize import linear_sum_assignment


def linear_sum(cost: np.array):
    """
    Solve the linear sum assignment problem. The goal is to find a complete assignment of workers to jobs of minimal cost.

    \min_{x}\sum_{i=1}^{n}\sum_{j=1}^{m} c_{i, j} x_{i, j}

    Ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html
    :param cost:
    :type cost:
    :return:
    :rtype:
    """
    row_ind, col_ind = linear_sum_assignment(cost)
    min_cost = cost[row_ind, col_ind].sum()
    return min_cost


def cluster(distance_array: np.array,
            mds_components: int = 20,
            n_clusters: int = 10) -> np.array:
    """
    After that, run K-Means clustering and get the labels.
    :param distance_array: Distance array
    :type distance_array: numpy.Array
    :param mds_components: Number of dimensions
    :type mds_components: int
    :param n_clusters: Number of clusters
    :type n_clusters: int
    :return: A list of labels for each element
    :rtype: numpy.Array
    """
    mds = MDS(n_components=mds_components, dissimilarity='precomputed')
    embedding = mds.fit_transform(distance_array)
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(embedding)
    return kmeans.labels_
