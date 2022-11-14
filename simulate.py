from schrodinger_evolve import SchrodingerEvolve


def simulate():
    se = SchrodingerEvolve()
    wave_function_ratio = 1.0
    for n in range(2, 10):
        wave_function_ratio = se.evolve_wave_function_ratio(wave_function_ratio, n)
        print(wave_function_ratio)


if __name__ == '__main__':
    simulate()
