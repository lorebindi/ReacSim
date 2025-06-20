import roadrunner

def plot(rr, t_max, show=True):

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
    p.tight_layout()

    if show:
        p.show()

def simulate_odes_road_runner(file_path, t_max):
    rr = roadrunner.RoadRunner(file_path)
    #rr.setIntegrator('gillespie')
    result = rr.simulate(0, t_max, 100)
    result[result < 0] = 0
    #data = rr.getSimulationData()
    plot(rr,t_max)







































