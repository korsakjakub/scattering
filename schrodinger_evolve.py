from config import Config
import numpy as np


class SchrodingerEvolve:
    def __init__(self, params=None):
        if params is None:
            params = {}
        self.reduced_mass = params["reduced_mass"] if "reduced_mass" in params else Config.get["reduced_mass"]
        self.energy = params["energy"] if "energy" in params else Config.get["energy"]
        self.De = params["De"] if "De" in params else Config.get["De"]
        self.Re = params["Re"] if "Re" in params else Config.get["Re"]
        self.R0 = params["R0"] if "R0" in params else Config.get["R0"]
        self.grid_step = 2 * np.pi / self.k(self.Re) / 50

    def k(self, r):
        return np.sqrt(2 * self.reduced_mass * (self.energy - self.effective_potential(r)))

    def eval_func_series(self, n, grid_step, func):
        return func(self.R0 + n * grid_step)

    def evolve_wave_function_ratio(self, wave_function_ratio, n):
        """
        Evolve the wave function using the Numerov's method
        :param wave_function_ratio: The ratio phi[m+1] / phi[m] at m = n-1
        :param n: the index of series
        :return: wavefunction ratio at step m + 1
        """

        def k_squared(r):
            return self.k(r)**2

        kn_sq = self.eval_func_series(n, self.grid_step, k_squared) ** 2
        knn_sq = self.eval_func_series(n - 1, self.grid_step, k_squared) ** 2

        return (2 * (1 - 5 * self.grid_step * kn_sq / 12) -
                (1 + self.grid_step ** 2 * knn_sq / 12) / wave_function_ratio) / \
               (1 + self.grid_step ** 2 * knn_sq / 12)

    def effective_potential(self, x):
        return self.potential(x)

    def potential(self, x):
        return self.De * ((self.Re / x) ** 12 - 2 * (self.Re / x) ** 6)

    def wave_function_ratio_to_wave_function(self, initial_value):
        return []