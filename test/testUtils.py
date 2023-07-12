import numpy as np
from skimage.util import random_noise


def add_salt_noise(matrix, percentage: int) -> np.array:
    amount = percentage / 100
    added = random_noise(matrix,
                         mode='salt',
                         amount=amount)
    return added


def add_pepper_noise(matrix, percentage: int) -> np.array:
    amount = percentage / 100
    added = random_noise(matrix,
                         mode='pepper',
                         amount=amount)
    return added


def add_sp_noise(matrix, percentage: int) -> np.array:
    amount = percentage / 100
    added = random_noise(matrix,
                         mode='s&p',
                         amount=amount)
    return added


def add_gauss_noise(matrix, mean: float) -> np.array:
    added = random_noise(matrix,
                         mode='gaussian',
                         mean=mean)
    return added


def add_periodicity_noise(matrix,
                          interval: int,
                          multiplier: int = 10,
                          axis=1) -> np.array:
    n_row = matrix.shape[0]
    n_col = matrix.shape[1]
    for i in range(0, n_row, interval):
        matrix[i, :] = matrix[i, :] * multiplier
    if axis != 1:
        for i in range(0, n_col, interval):
            matrix[:, i] = matrix[:, i] * multiplier
    return matrix
