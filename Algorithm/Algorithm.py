import numpy as np
import pandas as pd

from sklearn.manifold import MDS
from sklearn.cluster import KMeans
from scipy.optimize import linear_sum_assignment


def linear_sum(cost: np.array) -> np.float:
    """
    Solve the linear sum assignment problem. The goal is to find a complete assignment of workers to jobs of minimal
    cost.

    LATEX: \min_{x}\sum_{i=1}^{n}\sum_{j=1}^{m} c_{i, j} x_{i, j}

    Ref:
     - https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html
    :param cost:
    :type cost:
    :return:
    :rtype:
    """
    row_ind, col_ind = linear_sum_assignment(cost)
    min_cost = cost[row_ind, col_ind].sum()
    return min_cost


def cluster(distance_array: pd.DataFrame,
            mds_components: int = 20,
            n_clusters: int = 10) -> np.array:
    """
    Given a distance matrix with the distances between each pair of objects in a set, and a chosen number of
    dimensions, N, an MDS algorithm places each object into N-dimensional space (a lower-dimensional representation)
    such that the between-object distances are preserved as well as possible.
    After that, run **K-Means** clustering and get the labels.

    Ref:
     - https://scikit-learn.org/stable/modules/manifold.html#multidimensional-scaling
     - Multidimensional scaling. (2023, March 28). In Wikipedia. https://en.wikipedia.org/wiki/Multidimensional_scaling
    :param distance_array: Distance array dataframe
    :type distance_array: pd.DataFrame
    :param mds_components: Number of dimensions
    :type mds_components: int
    :param n_clusters: Number of clusters
    :type n_clusters: int
    :return: A list of labels for each element.
    :rtype: numpy.Array
    """
    index = distance_array.index
    mds = MDS(n_components=mds_components, dissimilarity='precomputed')
    print('Embedding...')
    embedding = mds.fit_transform(distance_array)
    print('Clustering...')
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(embedding)
    result = pd.DataFrame({'gene_id': list(index), 'labels': list(kmeans.labels_)})
    return result
