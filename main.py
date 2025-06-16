# This is a sample Python script.
import os
import time

from Parser import *
from Simulation_methods import *
from Graph_generation import *
T_MAX=1

if __name__ == '__main__':
    t=time.time()
    current_dir = os.path.dirname(__file__)
    file_path = os.path.join(current_dir, "Example", "Generated/lotka-volterra.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v3/GastricSlowWaveActivity.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v4/Ciliberto2003.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v4/You2010_General_Yeast_mRNA_Translation.xml")
    model = read_sbml_file(file_path)

    evolution = gillespie_ssa(model, T_MAX)

    plot(evolution, T_MAX) #TODO add memorization plot
    print(f"Execution time: {time.time() - t:.3f}")
