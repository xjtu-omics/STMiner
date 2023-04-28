import numpy as np

from scipy.optimize import linear_sum_assignment


def linear_sum(cost: np.array):
    row_ind, col_ind = linear_sum_assignment(cost)
    min_cost = cost[row_ind, col_ind].sum()
    return min_cost
