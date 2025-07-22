import roadrunner
from Graph_generation import *

def simulate_odes_road_runner(file_path, t_max, filename):
    rr = roadrunner.RoadRunner(file_path)
    #rr.setIntegrator('gillespie')
    result = rr.simulate(0, t_max, 100)

    result[result < 0] = 0
    #data = rr.getSimulationData()
    ode_plot(rr,t_max, filename)







































