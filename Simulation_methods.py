from Constants import *

import Parser
import math
import random


'''
1. Calculate τ (tau).
2. Iterate over the scheduled events and add them to to_eval_events if:
                t < delay_time < t + τ (where min_delay_time = min{delay_time}) and
                evaluate_expr(trigger_expr, ERROR_TRIGGER) at t = delay_time returns True
   Then
    a. If min_delay_time is not None:
        I. Iterate again over the scheduled events and remove those for which:
                evaluate_expr(trigger_expr, ERROR_TRIGGER) at t = delay_time returns False
        II. Set t = delay_time
    b. Otherwise, set t = t + τ
3. Iterate over all events, evaluate them at time t (which may be t + τ or min_delay_time),
    and if their evaluation changes from False → True, add them to to_eval_events.
4. Evaluate all events in to_eval_events, in no particular order.
    Before evaluating each event, always verify that its trigger condition remains True.
    (e.g., 1 - (check secondary trigger) - 2 - ...)
5. Only if to_eval_events is empty, execute the reaction.
'''
def gillespie_ssa (model, t_max):

    def evaluate_expr(expr, error_message, time_value, safe_globals = SAFE_GLOBALS_BASE):
        # Merge state and parameters in a single dictionary for expression evaluation
        local_scope = {**state, **parameters, "time": time_value}

        try:
            return eval(expr, safe_globals, local_scope)
        except:
            raise Exception(
                f"{error_message} — Evaluation failed.\n"
                f"Expression: {expr}\n"
            )

    def apply_events_assigment(event_assigment):
        for ea in event_assigment:
            var_id = ea.getVariable()
            value = evaluate_expr(libsbml.formulaToString(ea.getMath()), ERROR_EVENT_ASSIGNMENTS, t)

            # Apply to the right variable
            if var_id in state:
                state[var_id] = value
            elif var_id in parameters:
                parameters[var_id] = value
            else:
                raise Exception(f"Cannot apply assignment to unknown variable: {var_id}")

    t = 0.0

    pending_event_delay = {} # The set of all events with active delay

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
    evolution.update({id_specie: [amount] for id_specie,amount in state.items()})

    while t < t_max:
        propensities = []
        for r in reactions:
            propensities.append(evaluate_expr(r[RATE_FORMULA], ERROR_KINETIC_LAW, t, SAFE_GLOBALS_RATE))
        a0 = sum(propensities)
        if a0 == 0:
            break

        # New time
        y = random.random()
        tau = -math.log(y) / a0

        # Select delay events
        min_delay_time = t + tau
        to_eval_events = []

        for delay_time, event in pending_event_delay.values():
            # evaluation of the trigger at delay_time
            if evaluate_expr(event[TRIGGER_FORMULA], ERROR_TRIGGER, delay_time):
                # take only the events with the least delay_time
                if t < delay_time <= t + tau:
                    if delay_time < min_delay_time:
                        to_eval_events=[event]
                        min_delay_time = delay_time
                    elif delay_time == min_delay_time:
                        to_eval_events.append(event)

        t = min_delay_time

        # There are some events to evaluate
        if to_eval_events != []:
            for id_event, (_, event) in pending_event_delay.items():
                # Delete the events whose trigger is False
                if not evaluate_expr(event[TRIGGER_FORMULA], ERROR_TRIGGER, t):
                    del pending_event_delay[id_event]

        # Adding all the True-valuated events without delay or delay=0 at time t
        for event in events:
            trigger_expr = event[TRIGGER_FORMULA]
            val_trigger = evaluate_expr(trigger_expr, ERROR_TRIGGER,t) # at time t
            # The event isn't triggered
            if not(not event[PREVIOUS] and val_trigger):
                event[PREVIOUS] = val_trigger
                continue

            # Triggered event: its expr evaluates to True and its PREVIOUS value is False.

            delay = event[DELAY_FORMULA]
            if delay is None:
                # Triggered event without delay
                to_eval_events.append(event)
            else:
                delay_time = evaluate_expr(delay, ERROR_DELAY, t)

                if delay_time == 0:
                    # Triggered event's delay tag contains zero
                    to_eval_events.append(event)
                elif delay_time > 0:
                    # Triggered event's delay tag contains a value > 0
                    pending_event_delay[event[ID]] = (t + delay_time, event) # Adding new event with delay
                else:
                    raise Exception("The delay time is lesser than 0.")

            event[PREVIOUS] = True

        # There is one or more event(s) to evaluate and apply its EventAssigment
        if to_eval_events != []:
            for event in to_eval_events:
                # Check that the trigger is still True
                if evaluate_expr(event[TRIGGER_FORMULA], ERROR_TRIGGER, t):
                    apply_events_assigment(event[LIST_OF_EVENT_ASSIGMENT_FORMULA])

            # Updating the evolution of species
            evolution[TIME].append(t)
            for s_id in state:
                evolution[s_id].append(state[s_id])
            continue

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
        # Updating the evolution of species
        evolution[TIME].append(t)
        for s_id in state:
            evolution[s_id].append(state[s_id])

    return evolution


#TODO 1. trovare libreria che traduca un sistema di reazioni in un sistema di ODE
#    2. trovare  la libreria per eseguire il sistema di ODE (es. )


# (L2V3) An important question is whether an event can be triggered at or
# before the initial simulation time (t ≤ 0). The answer is no: events
# can only be triggered after the simulation has started, that is, for t > 0.
# If the trigger condition becomes false before the delay expires, the event
# is cancelled (As if the persistent attribute of the trigger were set to True).