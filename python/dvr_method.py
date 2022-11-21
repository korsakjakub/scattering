import numpy as np
from simulate import reduced_mass
from scipy import linalg
import matplotlib.pyplot as plt

De = 0.0011141
Re = 10.98
u = 1822.888
m39k = 38.963703 * u
m40k = 39.963999 * u
m41k = 40.961825 * u


def collocation_grid(min, max, size) -> np.ndarray:
    return np.array([min + i * (max - min)/size for i in range(size-1)])

def kinetic_energy(min, max, size, mu) -> np.ndarray:
    diagonal = np.diag([np.pi**2 * ((2*size**2+1)/3.0 - pow(np.sin(np.pi * i/size), -2))/(4.0*mu*pow(max-min, 2)) for i in range(1, size)])
    off_diagonal = np.array([[0.0 if i == j else np.pi**2*pow(-1, i-j) * (pow(np.sin(np.pi*(i-j)/(2*size)), -2)-pow(np.sin(np.pi*(i+j)/(2*size)), -2))/(4*mu*pow(max-min, 2)) for i in range(1, size)] for j in range(1, size)])
    return diagonal + off_diagonal

def potential_energy(min, max, size, ang_momentum, mu, coll_grid) -> np.ndarray:
    def potential(pos):
        return De * (pow(Re / pos, 12) - 2 * pow(Re / pos, 6))
    return np.diag([ang_momentum*(ang_momentum+1)/(2*mu*pow(coll_grid[i-1], 2)) + potential(coll_grid[i-1]) for i in range(1, size)])

def plot_spectrum(grid, vecs):
    [plt.scatter(grid, v, s=2) for i, v in enumerate(vecs)]
    plt.savefig("./python/whole_spectrum.png")
    """
    plt.scatter(grid, vecs[0], s=2)
    plt.scatter(grid, vecs[1], s=2)
    plt.scatter(grid, vecs[2], s=2)
    plt.scatter(grid, vecs[3], s=2)
    plt.savefig("./python/first_four.png")
    """

def plot_means(vals, means):
    plt.scatter(vals, means)
    plt.savefig("./python/mean_positions.png")

def plot_anharmonicity(diff):
    m = [n for n, _ in enumerate(diff)]
    plt.scatter(m, diff)
    plt.savefig("./python/anharmonicity.png")

if __name__ == '__main__':
    R0 = 9.5
    RN = 14.5
    grid_size = 1000
    mu = reduced_mass(m39k, m39k)
    angular_momentum = 0.0

    grid = collocation_grid(R0, RN, grid_size)
    kinetic_en = kinetic_energy(R0, RN, grid_size, mu)
    potential_en = potential_energy(R0, RN, grid_size, angular_momentum, mu, grid)
    hamiltonian = kinetic_en + potential_en
    vals, vecs = linalg.eigh(hamiltonian)
    vecs = vecs.T
    solutions = [[val, vecs[i]] for i, val in enumerate(vals)]
    v_solutions = []
    for i, el in enumerate(solutions):
        if el[0] < 0.0:
            v_solutions.append(el)

    vals = -np.array(v_solutions, dtype=object)[:, 0]
    vecs = np.array(v_solutions, dtype=object)[:, 1]
    # plot_spectrum(grid, vecs)
    means = [np.sum([grid[i] * abs(vec[i])**2 for i, _ in enumerate(vec)]) for _, vec in enumerate(vecs)]
    # plot_means(vals, means)
    omega0 = 0.5 * vals[0]
    harmonic_energies = [omega0 * (0.5 + v) for v, _ in enumerate(vals)]
    diffs = harmonic_energies - vals
    plot_anharmonicity(diffs)