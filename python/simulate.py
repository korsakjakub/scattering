import numpy as np

import matplotlib.pyplot as plt
import config

De = 0.0011141
Re = 10.98
u = 1822.888
m39k = 38.963703 * u
m40k = 39.963999 * u
m41k = 40.961825 * u


def k_sq(position, mu, col_energy, angular_momentum):
    def potential(pos):
        return De * (pow(Re / pos, 12) - 2 * pow(Re / pos, 6))

    def total_potential(pos, ang_momentum):
        return potential(pos) + ang_momentum * (ang_momentum + 1) / (2 * mu * pow(pos, 2))

    return 2 * mu * (col_energy - total_potential(position, angular_momentum))


def wave_function_ratio_get_n(prev_ratio, index, grid_step, k_sq_array):
    gridStepSq = pow(grid_step, 2)
    return (2 * (1 - 5 * gridStepSq * k_sq_array[index] / 12) - (
            1 + gridStepSq * k_sq_array[index - 1] / 12) / prev_ratio) / (1 + gridStepSq * k_sq_array[index + 1] / 12)


def reduced_mass(m1, m2):
    return m1 * m2 / (m1 + m2)


def simulate(cnf, collision_energy, angular_momentum) -> list:
    ratio = float(cnf['ratio'])
    initial_position = float(cnf['initialposition'])
    initial_phi = float(cnf['initialphi'])
    index_range = int(float(cnf['indexrange']))

    mu = reduced_mass(m39k, m40k)
    grid_step = 2 * np.pi / (50 * np.sqrt(k_sq(Re, mu, collision_energy, angular_momentum)))
    k = [k_sq(initial_position + i * grid_step, mu, collision_energy, angular_momentum) for i in range(index_range + 1)]
    ratios = []
    for i in range(index_range):
        ratio = wave_function_ratio_get_n(ratio, i, grid_step, k)
        ratios.append(ratio)

    phi = [initial_phi]
    for i in range(index_range):
        phi.append(phi[i] * ratios[i])
    return phi


if __name__ == '__main__':
    col_en = 1e-8
    ang_mom = 0.0
    conf = config.load()
    phi = simulate(conf, col_en, ang_mom)

    plt.figure(figsize=(10, 10))
    plt.plot(phi)
    plt.savefig("./python.png")
