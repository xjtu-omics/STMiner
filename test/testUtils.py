from typing import Union

import numpy as np


def get_salt_pepper_noise(matrix, noise_percentage: float) -> np.array:
    probability_array = np.random.random(matrix.shape)
    probability_array[probability_array <= noise_percentage] = 0
    probability_array[probability_array > noise_percentage] = 1
    return probability_array


def get_gauss_noise(matrix, mean: Union[float, int]) -> np.array:
    noise = np.random.normal(loc=mean, scale=.1, size=matrix.shape)
    noise[noise < 0] = 0
    matrix += noise
    return noise


def get_periodicity_noise(matrix,
                          interval: int,
                          multiplier: int = 2,
                          axis=1) -> np.array:
    n_row = matrix.shape[0]
    n_col = matrix.shape[1]
    for i in range(0, n_row, interval):
        matrix[i, :] = matrix[i, :] * multiplier
    if axis != 1:
        for i in range(0, n_col, interval):
            matrix[:, i] = matrix[:, i] * multiplier
    return matrix


def get_uniform_noise(matrix, mean):
    base = np.random.random(matrix.shape)
    noise = base * mean
    matrix += noise
    return noise
