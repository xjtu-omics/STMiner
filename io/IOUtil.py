import numpy as np


def merge_bin_coordinate(coordinate: np.ndarray, coordinate_min: int, bin_size: int):
    return np.floor((coordinate - coordinate_min) / bin_size).astype(int)


def get_bin_center(bin_coordinate: np.ndarray, coordinate_min: int, bin_size: int):
    return bin_coordinate * bin_size + coordinate_min + int(bin_size / 2)
