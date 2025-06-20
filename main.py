# This is a sample Python script.
import os
import threading
import time

from Parser import *
from Gillespie_events import *
from Graph_generation import *
from ODE_simulation import *
T_MAX = 1
T_MAX_EVAL = 10

def delayed_simulation():
    simulate_odes_road_runner(file_path, T_MAX)

if __name__ == '__main__':
    t=time.time()
    current_dir = os.path.dirname(__file__)
    file_path = os.path.join(current_dir, "Example", "Generated/lotka-volterra-event.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v3/GastricSlowWaveActivity.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v4/Ciliberto2003.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v4/You2010_General_Yeast_mRNA_Translation.xml")
    model = read_sbml_file(file_path)

    '''evolution = gillespie_ssa(model, T_MAX)
    plot_gillepsie(evolution, T_MAX) #TODO add memorization plot

    #simulate_odes_road_runner(file_path, T_MAX)

    print(f"Execution time: {time.time() - t:.3f}")'''

    evolution_result = {}
    def run_gillespie():
        evolution_result['data'] = gillespie_ssa(model, T_MAX)

    thread = threading.Thread(target=run_gillespie)
    thread.start()
    thread.join(timeout=T_MAX_EVAL)

    if thread.is_alive():
        print(f"⚠️ Gillespie SSA did not finish in {T_MAX_EVAL} seconds.")
        t = time.time()
    else:
        plot_gillepsie(evolution_result['data'], T_MAX)

    simulate_odes_road_runner(file_path, T_MAX)

    print(f"Execution time: {time.time() - t:.3f}")
