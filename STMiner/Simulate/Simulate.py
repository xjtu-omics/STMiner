from random import choice
from typing import List

import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

from STMiner.Simulate.simUtils import *


class Simulator:
    def __init__(self, gene_exp_array: Union[np.array, List[np.array]]):
        self.gene_exp_array_list = gene_exp_array
        self.noise_type = None
        self.noise_argument = None

    def set_noise_type(self, noise_type, noise_argument):
        if noise_type == 'gauss':
            if noise_argument <= 0 or noise_argument > 1:
                raise ValueError('noise_argument should be between 0 and 1')
        elif noise_type == 'periodicity':
            if noise_argument < 0:
                raise ValueError('noise_argument should larger than 0')
        elif noise_type == 'undersampling':
            if noise_argument >= 1:
                raise ValueError('noise_argument should smaller than 1')
        self.noise_type = noise_type
        self.noise_argument = noise_argument

    def generate(self, offset_radius, count, add_noise=True, offset_probability=0.5):
        df = None
        if not isinstance(self.gene_exp_array_list, list):
            self.gene_exp_array_list = [self.gene_exp_array_list]
        for gene_index, gene_arr in enumerate(self.gene_exp_array_list):
            for i in range(count):
                sim_array = np.zeros(gene_arr.shape)
                row_count = sim_array.shape[0]
                col_count = sim_array.shape[1]

                # Random adjust the expression level, from 1/2 to 2 times.
                scale_factor = np.random.uniform(0.5, 2)
                scaled_arr = gene_arr * scale_factor

                for index, value in np.ndenumerate(scaled_arr):
                    if value != 0:
                        if np.random.rand() > (1 - offset_probability):
                            row_index = choice(
                                list(range(max(index[0] - offset_radius, 0),
                                           min(index[0] + offset_radius + 1, row_count))))
                            col_index = choice(
                                list(range(max(index[1] - offset_radius, 0),
                                           min(index[1] + offset_radius + 1, col_count))))
                            sim_array[row_index, col_index] = value
                        else:
                            sim_array[index] = value

                if add_noise:
                    if self.noise_type == 'gauss':
                        noise = get_gauss_noise(sim_array, sim_array.mean(), self.noise_argument)
                        sim_array = sim_array + noise
                    elif self.noise_type == 'periodicity':
                        noise = get_periodicity_noise(sim_array, self.noise_argument)
                        sim_array = np.array(sim_array) * noise
                    elif self.noise_type == 'undersampling':
                        noise = under_sampling(sim_array, self.noise_argument)
                        sim_array = np.array(sim_array) * noise
                    else:
                        noise = get_uniform_noise(sim_array, self.noise_argument)
                        sim_array = sim_array + noise
                sparse_array = csr_matrix(sim_array)
                row_indices, col_indices = sparse_array.nonzero()
                index = list(zip(row_indices, col_indices))
                gene_name = 'gene_P' + str(gene_index) + '_N' + str(i)
                if df is not None:
                    tmp = pd.DataFrame(sparse_array.data, index=index, columns=[gene_name])
                    df = pd.concat([df, tmp], axis=1)
                else:
                    df = pd.DataFrame(sparse_array.data, index=index, columns=[gene_name])
        arr_no_nan = np.nan_to_num(df.values)
        adata = AnnData(arr_no_nan)
        adata.obs['cell_id'] = list(df.index)
        adata.obs['x'] = [x[0] for x in list(df.index)]
        adata.obs['y'] = [x[1] for x in list(df.index)]
        adata.var['gene_ids'] = list(df.columns)
        adata.var_names = list(df.columns)
        return adata
