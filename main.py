# This is a sample Python script.
import argparse
import time
from multiprocessing import Process
from pathlib import Path

from Gillespie_events import *
from Graph_generation import *
from ODE_simulation import *
import Parser

T_MAX = 10

'''if __name__ == '__main__':

    current_dir = os.path.dirname(__file__)
    file_path = os.path.join(current_dir, "Example", "Generated/lotka-volterra-event.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v3/GastricSlowWaveActivity.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v4/Ciliberto2003.xml")
    #file_path = os.path.join(current_dir, "Example", "BioModels/l2v4/You2010_Ge neral_Yeast_mRNA_Translation.xml")
    file_path_csv = os.path.join(current_dir, "Example", "Inference_of_stochastic_rate_constant")

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


def run_gillespie(parser, filename):
    try:
        t = time.time()
        Gillespie_sim = Gillespie(parser, T_MAX)
        Gillespie_sim.gillespie_ssa()
        evolution = Gillespie_sim.get_evolution()
        plot_gillepsie(evolution, T_MAX, filename[:-4], parser.dfs)
        print(f"Execution time with Gillespie for {filename}: {time.time() - t:.3f}")
    except Exception as e:
        print(f"{RED}[Gillespie ERROR] {filename}: {e}{RESET}")

if __name__ == '__main__':
    # Single parser with all command-line options
    parser = argparse.ArgumentParser(
        description="Run Gillespie with optional timeout and ODE fallback. Process SBML files."
    )

    # Optional timeout for Gillespie simulation
    parser.add_argument(
        "--max-time-gillespie",
        type=float,
        default=0.0,
        help="Max time (in seconds) for Gillespie simulation. If 0, run without timeout."
    )

    # One or more SBML files to process
    parser.add_argument(
        '--filesbml',
        type=Path,
        nargs='+',
        required=True,
        help="One or more SBML file paths to process"
    )

    # Parse all arguments
    args = parser.parse_args()

    # Read values
    MAX_TIME_GILLESPIE = args.max_time_gillespie
    files_to_process = [p.resolve() for p in args.filesbml]  # Absolute paths

    # Path to save CSVs
    current_dir = os.path.dirname(__file__)
    file_path_csv = os.path.join(current_dir, "Example", "Inference_of_stochastic_rate_constant")

    for file_path in files_to_process:
        if file_path.exists() and file_path.suffix == ".xml":
            filename = file_path.name
            print(f"\nProcessing file: {filename}")
            parser_instance = Parser.Parser(str(file_path), file_path_csv)

            if MAX_TIME_GILLESPIE > 0.0:
                # Run Gillespie in a separate process with timeout
                p = Process(target=run_gillespie, args=(parser_instance, filename))
                p.start()
                p.join(timeout=MAX_TIME_GILLESPIE)

                if p.is_alive():
                    print(f"{RED}[TIMEOUT] Gillespie took too long for {filename}, terminating...{RESET}")
                    p.terminate()
                    p.join()
                else:
                    print(f"Gillespie simulation for {filename} completed within allowed time.")
            else:
                # Run Gillespie directly (no timeout)
                run_gillespie(parser_instance, filename)

            # Always run ODE simulation
            try:
                t1 = time.time()
                simulate_odes_road_runner(str(file_path), T_MAX, filename[:-4])
                print(f"Execution time ODEs for {filename}: {time.time() - t1:.3f}")
            except Exception as e:
                print(f"{RED}[ODEs ERROR] {filename}: {e}{RESET}")
        else:
            print(f"{RED}[SKIP] File does not exist or is not .xml: {file_path}{RESET}")