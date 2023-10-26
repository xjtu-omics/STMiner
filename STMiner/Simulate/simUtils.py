from typing import Union

import numpy as np


def under_sampling(matrix, percentage: float) -> np.array:
    probability_array = np.random.random(matrix.shape)
    probability_array[probability_array <= percentage] = 0
    probability_array[probability_array > percentage] = 1
    return probability_array


def get_gauss_noise(matrix, mean: Union[float, int], scale) -> np.array:
    noise = np.random.normal(loc=mean * scale, size=matrix.shape)
    noise[noise < 0] = 0
    return noise


def get_periodicity_noise(matrix,
                          interval: int,
                          scale_factor: int = 0.1) -> np.array:
    noise_arr = np.ones(matrix.shape)
    n_row = matrix.shape[0]
    n_col = matrix.shape[1]
    for i in range(0, n_row, interval):
        noise_arr[i, :] = scale_factor
    for i in range(0, n_col, interval):
        noise_arr[:, i] = scale_factor
    return noise_arr


def get_uniform_noise(matrix, mean):
    base = np.random.random(matrix.shape)
    noise = base * mean
    return noise
