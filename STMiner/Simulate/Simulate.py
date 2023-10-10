from STMiner.Simulate.simUtils import *
import numpy as np
from random import choice
from anndata import AnnData


class Simulator:
    def __init__(self, gene_exp_array):
        self.gene_exp_array = gene_exp_array
        self.noise_type = None
        self.noise_argument = None

    def set_noise_type(self, noise_type, noise_argument):
        if noise_type == 'gauss':
            if noise_argument <= 0:
                raise ValueError('noise_argument should larger than 0')
        elif noise_type == 'periodicity':
            if noise_argument <= 1:
                raise ValueError('noise_argument should larger than 1')
        elif noise_type == 'sp':
            if noise_argument >= 1:
                raise ValueError('noise_argument should smaller than 1')
        self.noise_type = noise_type
        self.noise_argument = noise_argument

    def generate(self, offset_radius, count, add_noise=True, offset_probability=0.5):
        indices = np.column_stack(np.where(self.gene_exp_array))
        for i in range(count):
            scale_factor = np.random.uniform(0.5, 2)
            sim_array = np.zeros(self.gene_exp_array.shape)
            row_count = sim_array.shape[0]
            col_count = sim_array.shape[1]

            scaled_arr = self.gene_exp_array * scale_factor
            for index, value in np.ndenumerate(scaled_arr):
                if value != 0:
                    if np.random.rand() > (1 - offset_probability):
                        row_index = choice(
                            list(range(max(index[0] - offset_radius, 0), min(index[0] + offset_radius + 1, row_count))))
                        col_index = choice(
                            list(range(max(index[1] - offset_radius, 0), min(index[1] + offset_radius + 1, col_count))))
                        sim_array[row_index, col_index] = value
                    else:
                        sim_array[index] = value

            if add_noise:
                if self.noise_type == 'gauss':
                    noise = get_gauss_noise(sim_array, self.noise_argument)
                    sim_array = sim_array + noise
                elif self.noise_type == 'periodicity':
                    noise = get_periodicity_noise(sim_array, self.noise_argument)
                    sim_array = np.array(sim_array) * noise
                elif self.noise_type == 'sp':
                    noise = get_salt_pepper_noise(sim_array, self.noise_argument)
                    sim_array = np.array(sim_array) * noise
                else:
                    noise = get_uniform_noise(sim_array, self.noise_argument)
                    sim_array = sim_array + noise


            values = sim_array[indices[:, 0], indices[:, 1]]
            adata = AnnData()
        return adata
