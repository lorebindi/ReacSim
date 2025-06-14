from Constants import *

import Parser
import math
import random

def evaluate_rate(formula, species, parameters):
    local_dict = {}
    local_dict.update(species)
    local_dict.update(parameters)

    safe_globals = {"__builtins__": None, "pow": pow, "math": math}
    try:
        val = eval(formula, safe_globals, local_dict)
        return val
    except:
        return 0.0

def gillespie_ssa (model, t_max):
    t = 0.0

    state = Parser.extract_species(model)
    if state is None:
        raise Exception("Model has no species")

    reactions = Parser.extract_reactions(model)
    if reactions is None:
        raise Exception("Model has no reactions")

    parameters = Parser.extract_parameters(model)
    if parameters is None:
        raise Exception("Model has no parameters")

    evolution = {TIME: [t]}
    evolution.update({id: [amount] for id,amount in state.items()})

    while t < t_max:
        propensities = []
        for r in reactions:
            propensities.append(evaluate_rate(r[RATE], state, parameters))

        a0 = sum(propensities)
        if a0 == 0:
            break

        # New time
        y = random.random()
        tau = -math.log(y) / a0
        t += tau

        #Reaction choose
        n = random.random() * a0
        cumulative = 0.0
        for i, a in enumerate(propensities):
            cumulative += a
            if n < cumulative:
                chosen = i
                break

        # Update of state
        r = reactions[chosen]
        for s, stoich in r[REACTANTS].items():
            state[s] -= stoich
        for p, stoich in r[PRODUCTS].items():
            state[p] += stoich

        evolution[TIME].append(t)
        for s_id in state:
            evolution[s_id].append(state[s_id])

    return evolution


#TODO 1. trovare libreria che traduca un sistema di reazioni in un sistema di ODE
#    2. trovare  la libreria per eseguire il sistema di ODE (es. )