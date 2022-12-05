import numpy as np
import matplotlib.pyplot as plt

hertz_to_hartree = 1.519828500716e-16
au_to_gauss = 2.350518086e9
g_e = 2.002319
m_p = 1836.15
bohr_magneton = 1 / 2
nuclear_magneton = bohr_magneton / m_p


def i_ps_p(state, spin, nuclear):
    result = np.sqrt((spin * (spin + 1) - state[0] * (state[0] - 1)) *
                     (nuclear * (nuclear + 1) - state[1] * (state[1] + 1)))
    return result, [state[0] - 1, state[1] + 1]


def i_ms_p(state, spin, nuclear):
    result = np.sqrt((spin * (spin + 1) - state[0] * (state[0] + 1)) *
                     (nuclear * (nuclear + 1) - state[1] * (state[1] - 1)))
    return result, [state[0] + 1, state[1] - 1]


def hamiltonian(params):
    mag_field = params['B']
    states = params['states']
    g_n = params['g_n']
    hyperfine_coupling = params['hyperfine_coupling']
    spin = params['spin']
    nuclear = params['nuclear']

    dim = len(states)

    diagonal, off_diagonal = np.zeros((dim, dim)), np.zeros((dim, dim))

    for i in range(dim):
        diagonal[i][i] = g_e * bohr_magneton * mag_field * states[i][0] + \
                     g_n * nuclear_magneton * mag_field * states[i][1] + \
                     hyperfine_coupling * states[i][0] * states[i][1]

    for state in states:
        result_pm, state_pm = i_ps_p(state, spin, nuclear)[0], i_ps_p(state, spin, nuclear)[1]
        result_mp, state_mp = i_ms_p(state, spin, nuclear)[0], i_ms_p(state, spin, nuclear)[1]

        if result_pm != 0:
            off_diagonal[states.index(state)][states.index(state_pm)] = \
                1 / 2 * hyperfine_coupling * result_pm
        elif result_mp != 0:
            off_diagonal[states.index(state)][states.index(state_mp)] = \
                1 / 2 * hyperfine_coupling * result_mp

    ham = diagonal + off_diagonal
    energy, vec = np.linalg.eig(ham)
    return np.sort(energy/hertz_to_hartree)


def splitting(element):
    states = []
    g_n, hyperfine_coupling, spin, nuclear = np.zeros(4)
    if element == 'Li':
        spin, nuclear = 1 / 2, 1
        electron_spins = np.array([-spin, spin])
        nuclear_spins = np.arange(-nuclear, nuclear + 1, 1)
        for i in range(electron_spins.size):
            for j in range(nuclear_spins.size):
                states.append([electron_spins[i], nuclear_spins[j]])
        g_n = -0.000447
        hyperfine_coupling = 152.1e6 * hertz_to_hartree

    elif element == 'Rb':
        spin, nuclear = 1 / 2, 5 / 2
        electron_spins = np.array([-spin, spin])
        nuclear_spins = np.arange(-nuclear, nuclear + 1, 1)
        for i in range(electron_spins.size):
            for j in range(nuclear_spins.size):
                states.append([electron_spins[i], nuclear_spins[j]])
        g_n = -0.000294
        hyperfine_coupling = 1011.9e6 * hertz_to_hartree

    mag_field = np.linspace(0, 2000, 500) / au_to_gauss

    for i, _ in enumerate(states):
        energies = []
        for B in mag_field:
            energies.append(hamiltonian({
                'B': B,
                'states': states,
                'g_n': g_n,
                'hyperfine_coupling': hyperfine_coupling,
                'spin': spin,
                'nuclear': nuclear
            })[i])
        plt.plot(mag_field * au_to_gauss, energies)

    plt.xlabel("B/Gauss")
    plt.ylabel("E/GHz")
    plt.title('Atomic hyperfine spectra of %s' % element)
    plt.savefig(f'{element}.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    splitting('Li')
    splitting('Rb')
