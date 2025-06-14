# This is a sample Python script.
import os
import time

from Parser import *
from Simulation_methods import *
from Graph_generation import *
T_MAX=100

if __name__ == '__main__':
    t=time.time()
    current_dir = os.path.dirname(__file__)
    file_path = os.path.join(current_dir, "Example", "Generated/lotka-volterra.xml")
    model = read_sbml_file(file_path)

    evolution = gillespie_ssa(model, T_MAX)

    plot(evolution, T_MAX) #TODO add memorization plot
    print(f"Execution time: {time.time() - t:.3f}")
