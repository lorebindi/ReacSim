from Constants import *

import Parser
import math
import random

def evaluate_rate(formula, species, parameters):
    local_dict = {}
    local_dict.update(species)
    local_dict.update(parameters)

    try:
        return eval(formula, {"__builtins__": None, "math": math}, local_dict)
    except:
        return 0.0

def gillespie_ssa (model, t_max):
    t = 0.0

    state = Parser.extract_species(model)
    if state is None:
        return None

    reactions = Parser.extract_reactions(model)
    if reactions is None:
        return None

    evolution = {TIME: [t]}
    evolution.update({id: [amount] for id,amount in state.items()})

    while t < t_max:
        propensities = []
        for r in reactions:
            propensities.append(evaluate_rate(r[RATE], state, r[PARAMETERS]))

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
