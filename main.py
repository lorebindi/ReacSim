# This is a sample Python script.
import os
import threading
import time

from Parser import *
from Gillespie_events import *
from Graph_generation import *
from ODE_simulation import *
from Kinetic_constant_inference import  *
T_MAX = 1
T_MAX_EVAL = 10


if __name__ == '__main__':

    current_dir = os.path.dirname(__file__)
    file_path = os.path.join(current_dir, "Example", "Generated/lotka-volterra-event.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v3/GastricSlowWaveActivity.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v4/Ciliberto2003.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v4/You2010_General_Yeast_mRNA_Translation.xml")
    model = read_sbml_file(file_path)

    file_path_csv = os.path.join(current_dir, "Inference_of_kinetic_laws", "1.csv")


    t = time.time()
    evolution = gillespie_ssa(model, T_MAX)
    plot_gillepsie(evolution, T_MAX) #TODO add memorization plot
    print(f"Execution time with Gillespie: {time.time() - t:.3f}")

    t = time.time()
    simulate_odes_road_runner(file_path, T_MAX)
    print(f"Execution time ODEs: {time.time() - t:.3f}")
