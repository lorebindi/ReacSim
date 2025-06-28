# This is a sample Python script.
import sys
import time

from Gillespie_events import *
from Graph_generation import *
from ODE_simulation import *

T_MAX = 10

'''if __name__ == '__main__':

    current_dir = os.path.dirname(__file__)
    file_path = os.path.join(current_dir, "Example", "Generated/lotka-volterra-event.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v3/GastricSlowWaveActivity.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v4/Ciliberto2003.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v4/You2010_General_Yeast_mRNA_Translation.xml")
    file_path_csv = os.path.join(current_dir, "Example", "Inference_of_kinetic_laws")

    parser=Parser(file_path, file_path_csv  )




    t = time.time()
    Gillespie = Gillespie(parser, T_MAX)
    Gillespie.gillespie_ssa()
    evolution=Gillespie.get_evolution()
    plot_gillepsie(evolution, T_MAX) #TODO add memorization plot
    print(f"Execution time with Gillespie: {time.time() - t:.3f}")'''

'''t = time.time()
    simulate_odes_road_runner(file_path, T_MAX)
    print(f"Execution time ODEs: {time.time() - t:.3f}")'''

if __name__ == '__main__':
    # , "lotka_volterra_p_2" "SIR_p"
    file_to_valuate = ["testBeta"]

    current_dir = os.path.dirname(__file__)
    test_folder = os.path.join(current_dir, "Test")
    file_path_csv = os.path.join(current_dir, "Example", "Inference_of_kinetic_laws")

    for filename in os.listdir(test_folder):
        if filename.endswith(".xml") and filename[:-4] in file_to_valuate:
            file_path = os.path.join(test_folder, filename)
            print(f"\nProcessing file: {filename}")
            parser = Parser(file_path, file_path_csv)
            try:

                t = time.time()
                Gillespie_sim = Gillespie(parser, T_MAX)
                Gillespie_sim.gillespie_ssa()
                evolution = Gillespie_sim.get_evolution()
                print(plot_gillepsie(evolution, T_MAX, filename[:-4], parser.dfs))
                print(f"Execution time with Gillespie for {filename}: {time.time() - t:.3f}")
            except Exception as e:
                print(f"{RED}[Gillespie ERROR] {filename}: {e}{RESET}")

            '''try:
                # Optional ODE Simulation:
                t = time.time()
                simulate_odes_road_runner(file_path, T_MAX, filename[:-4])
                print(f"Execution time ODEs for {filename}: {time.time() - t:.3f}")
            except Exception as e:
                print(f"{RED}[ODEs ERROR] {filename}: {e}{RESET}")'''

