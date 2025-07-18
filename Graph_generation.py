import os
import numpy as np
from Constants import *
import matplotlib.pyplot as plt

def plot_gillepsie(evolution, t_max, filename=None, dfs_csv= None):
    if evolution is None: raise Exception("Evolution is None")
    if TIME not in evolution: raise Exception("Evolution does not contain time")
    if TIME in evolution.keys() and len(evolution) == 1: raise Exception("Evolution does not contain data.")

    time_values = evolution[TIME]
    specie_names = [k for k in evolution if k != TIME]

    plt.figure(figsize=(8, 6))

    custom_colors = {
        'Prey': 'blue',
        'Predator': 'green',
        'Susceptible': 'darkcyan',
        'Infected': 'red',
        'Recovered': 'purple'
    }

    for specie in specie_names:
        valori = evolution[specie]
        plt.plot(time_values, valori, label=specie, color=custom_colors[specie])

    # Scatter plot dai punti del DataFrame, se disponibile
    for df in dfs_csv.values():
        for specie in specie_names:
            if TIME in df.columns and specie in df.columns:
                plt.scatter(df[TIME], df[specie], s=40, marker='o', label=f"{specie} (data)", edgecolors='black', alpha=0.6, color=custom_colors[specie])

    plt.xlabel("Time")
    plt.ylabel("Molecule count")
    plt.title(f"Gillespie: {filename}")
    plt.legend(loc="upper left")
    plt.grid(True)
    plt.xlim(0, t_max)
    plt.ylim(bottom=0)
    plt.tight_layout()

    '''if GRAPH_FOLDER:
        os.makedirs(GRAPH_FOLDER, exist_ok=True)
        file_path = os.path.join(GRAPH_FOLDER, "gillepsie_plot.png")
        if os.path.exists(file_path):
            os.remove(file_path)  # Elimina il file precedente
        plt.savefig(file_path)
        print(f"Grafico salvato in: {file_path}")'''
    #plt.savefig("/home/lorenzo/Desktop/ReacSim/Graphs/gillespie.png", dpi=900)
    plt.show()

def ode_plot(rr, t_max, filename, show=True):

    import pylab as p

    result = rr.getSimulationData()

    if result is None:
        raise Exception("no simulation result")

    # assume result is a standard numpy array

    selections = rr.timeCourseSelections

    if len(result.shape) != 2 or result.shape[1] != len(selections):
        raise Exception("simulation result columns not equal to number of selections,"
                        "likely a simulation has not been run")

    times = result[:,0]

    for i in range(1, len(selections)):
        series = result[:,i]
        name = selections[i]
        p.plot(times, series, label=str(name))

        p.legend()

    p.grid(True)
    p.xlim(0, t_max)
    p.title(f"ODE: {filename}")
    p.tight_layout()

    #p.savefig("/home/lorenzo/Desktop/ReacSim/Graphs/ode.png", dpi=900)

    if show:
        p.show()


def stochastic_rate_constant_plot (x_vals, y_vals, k, constant):

    plt.figure(figsize=(8, 6))
    plt.scatter(x_vals, y_vals, color='blue', label='Data points $(v(\\Delta t),\ x(\\Delta t))$')

    # Create line for fitted model y = k*x
    x_line = np.linspace(min(x_vals), max(x_vals), 100)
    y_line = k * x_line
    plt.plot(x_line, y_line, color='red', label=f'Fitted line: v = {k:.4f} * x')

    plt.xlabel('Mass Action Term (x)')
    plt.ylabel('Reaction Rate (v)')
    plt.title(f'Stochastic rate constant inference for {constant}')
    plt.legend()
    plt.grid(True)
    #plt.savefig("/home/lorenzo/Desktop/ReacSim/Graphs/constant.png" , dpi=900)
    plt.show()

    '''if GRAPH_FOLDER:
        os.makedirs(GRAPH_FOLDER, exist_ok=True)
        file_path = os.path.join(GRAPH_FOLDER, "kinetic_costant_plot.png")
        if os.path.exists(file_path):
            os.remove(file_path)  # Elimina il file precedente
        plt.savefig(file_path)
        print(f"Grafico salvato in: {file_path}")'''



