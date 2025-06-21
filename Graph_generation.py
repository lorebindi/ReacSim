import numpy as np

from Constants import *

import matplotlib.pyplot as plt

def plot_gillepsie(evolution, t_max, save_dir=None):
    if evolution is None: raise Exception("Evolution is None")
    if TIME not in evolution: raise Exception("Evolution does not contain time")
    if TIME in evolution.keys() and len(evolution) == 1: raise Exception("Evolution does not contain data.")

    time_values = evolution[TIME]
    specie_names = [k for k in evolution if k != TIME]

    plt.figure(figsize=(10, 6))

    for specie in specie_names:
        valori = evolution[specie]
        plt.plot(time_values, valori, label=specie)

    plt.xlabel("Time")
    plt.ylabel("Molecule count")
    plt.title("Evolution of species over time")
    plt.legend(loc="upper right")
    plt.grid(True)
    plt.xlim(0, t_max)
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.show()

def ode_plot(result):

    for i, name in enumerate(result.odesys.names):  # nomi delle specie
        value = result.yout[:, i]
        plt.plot(result.xout, value, label=name)

    plt.xlabel("Tempo")
    plt.ylabel("Concentrazione")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend()
    plt.grid(True)
    plt.title("Evoluzione delle specie nel tempo")
    plt.show()

'''
def kinetic_constant_plot (x_vals, y_vals, k):

    plt.figure(figsize=(8, 6))
    plt.scatter(x_vals, y_vals, color='blue', label='Data points (v, x)')

    # Create line for fitted model y = k*x
    x_line = np.linspace(min(x_vals), max(x_vals), 100)
    y_line = k_opt[0] * x_line
    plt.plot(x_line, y_line, color='red', label=f'Fitted line: v = {k_opt[0]:.4f} * x')

    plt.xlabel('Rate Expression (x)')
    plt.ylabel('Reaction Rate (v)')
    plt.title(f'Kinetic constant inference for {constant}')
    plt.legend()
    plt.grid(True)
    plt.show()
'''