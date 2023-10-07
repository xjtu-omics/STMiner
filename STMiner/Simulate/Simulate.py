from simUtils import *
import numpy as np
from random import random

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

    def generate(self, count):
        sim_array = np.zeros(self.gene_exp_array.shape)
        for index, value in np.ndenumerate(self.gene_exp_array):
            if value != 0:
                if random() > 0.1:
                    sim_array[index] = value

        if self.noise_type == 'gauss':
            noise = get_gauss_noise(sim_array, self.noise_argument)
            noised = sim_array + noise
        elif self.noise_type == 'periodicity':
            noise = get_periodicity_noise(sim_array, self.noise_argument)
            noised = np.array(sim_array) * noise
        elif self.noise_type == 'sp':
            noise = get_salt_pepper_noise(sim_array, self.noise_argument)
            noised = np.array(sim_array) * noise
        else:
            noise = get_uniform_noise(sim_array, self.noise_argument)
            noised = sim_array + noise
