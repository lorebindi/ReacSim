from Constants import *

import Parser
import math
import random

def gillespie_ssa (model, t_max):
    def evaluate_expr(expr, error_message, safe_globals = SAFE_GLOBALS_BASE):
        # Merge state and parameters in a single dictionary for expression evaluation
        local_scope = {**state, **parameters, "time": t}

        try:
            return eval(expr, safe_globals, local_scope)
        except:
            raise Exception(
                f"{error_message} — Evaluation failed.\n"
                f"Expression: {expr}\n"
            )

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

    events = Parser.extract_events(model)

    evolution = {TIME: [t]}
    evolution.update({id: [amount] for id,amount in state.items()})

    while t < t_max:
        propensities = []
        for r in reactions:
            propensities.append(evaluate_expr(r[RATE], ERROR_KINETIC_LAW, SAFE_GLOBALS_RATE))

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
            if state[s] < 0:
                state[s] = 0
        for p, stoich in r[PRODUCTS].items():
            state[p] += stoich
            if state[p] < 0:
                state[p] = 0

        # (L2V3) An important question is whether an event can be triggered at or
        # before the initial simulation time (t ≤ 0). The answer is no: events
        # can only be triggered after the simulation has started, that is, for t > 0.
        for event in events:
            trigger_expr = event[TRIGGER]
            val_expr = evaluate_expr(trigger_expr, ERROR_TRIGGER)
            # A trigger is active when its expr evaluates to True and its PREVIOUS value is False.
            if not event[PREVIOUS] and val_expr:
                for ea in event[LIST_OF_EVENT_ASSIGMENT]:
                    var_id = ea.getVariable()
                    value = evaluate_expr(libsbml.formulaToString(ea.getMath()),ERROR_EVENT_ASSIGNMENTS)

                    # Apply to the right variable
                    if var_id in state:
                        state[var_id] = value
                    elif var_id in parameters:
                        parameters[var_id] = value
                    else:
                        raise Exception(f"Cannot apply assignment to unknown variable: {var_id}")
            event[PREVIOUS] = val_expr

        evolution[TIME].append(t)
        for s_id in state:
            evolution[s_id].append(state[s_id])

    return evolution


#TODO 1. trovare libreria che traduca un sistema di reazioni in un sistema di ODE
#    2. trovare  la libreria per eseguire il sistema di ODE (es. )