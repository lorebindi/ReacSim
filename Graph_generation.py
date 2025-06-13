from Constants import *

import matplotlib.pyplot as plt

def plot(evolution, save_dir=None):
    if evolution is None: return None
    if TIME not in evolution: return None

    time_values = evolution[TIME]
    specie_names = [k for k in evolution if k != TIME]

    plt.figure(figsize=(10, 6))

    for specie in specie_names:
        valori = evolution[specie]
        plt.plot(time_values, valori, label=specie, marker='o')

    plt.xlabel("Time")
    plt.ylabel("Molecule count")
    plt.title("Evolution of species over time")
    plt.legend(loc="upper right")
    plt.grid(True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.show()