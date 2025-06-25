# Librerie standard
import math
import random

# Moduli locali
from Constants import *
from Parser import *

'''
1. Calculate τ (tau).
2. Iterate over the scheduled events and add them to to_eval_events if:
                t < delay_time < t + τ (where min_delay_time = min{delay_time}) and
                evaluate_expr(trigger_expr, ERROR_TRIGGER) at t = delay_time returns True
   Then
    a. If min_delay_time is not None:
        I. Iterate again over the scheduled events and remove those for which:
                                  raise Exception(f"Variable '{string_of_variable}' not found in model.")
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

# (L2V3) An important question is whether an event can be triggered at or
# before the initial simulation time (t ≤ 0). The answer is no: events
# can only be triggered after the simulation has started, that is, for t > 0.
# If the trigger condition becomes false before the delay expires, the event
# is cancelled (As if the persistent attribute of the trigger were set to True).

class Gillespie:
    def __init__(self, parser, t_max):
        self.t_max = t_max
        self.t = 0.0
        self.pending_event_delay = {} # The set of all events with active delay
        self.parser = parser
        
        if self.parser.species is None:
            raise Exception("Model has no species")

        if self.parser.reactions is None:
            raise Exception("Model has no reactions")

        if self.parser.parameters is None:
            raise Exception("Model has no parameters")

        self.evolution = {TIME: [self.t]}
        self.evolution.update({id_specie: [amount] for id_specie,amount in self.parser.species.items()})

    def evaluate_expr(self, expr, error_message, time_value, safe_globals = SAFE_GLOBALS_BASE):
        # Merge state and parameters in a single dictionary for expression evaluation
        local_scope = {**self.parser.species, **self.parser.parameters, "time": time_value}

        try:
            return eval(expr, safe_globals, local_scope)
        except:
            raise Exception(
                f"{error_message} — Evaluation failed.\n"
                f"Expression: {expr}\n")
        
    def apply_events_assigment(self, event_id, events_assigments, values_from_trigger_time):
        # Update the current state and parameters with the values captured at the trigger time
        for var_ea, value_ea in values_from_trigger_time.items():
            if var_ea in self.parser.species:
                self.parser.species[var_ea] = value_ea
            elif var_ea in self.parser.parameters:
                self.parser.parameters[var_ea] = value_ea

        for ea in events_assigments:
            var_id = ea.getVariable()
            value = self.evaluate_expr(libsbml.formulaToString(ea.getMath()), ERROR_EVENT_ASSIGNMENTS, self.t)

            if value < 0:
                raise Exception(f"Impossible to do apply the event {event_id}")

            # Apply to the right variable
            if var_id in self.parser.species:
                self.parser.species[var_id] = value
            elif var_id in self.parser.parameters:
                self.parser.parameters[var_id] = value
            else:
                raise Exception(f"Cannot apply assignment to unknown variable: {var_id}")

    def gillespie_ssa (self):
        while self.t < self.t_max:
            propensities = []
            for r in self.parser.reactions:
                propensity = self.evaluate_expr(r[RATE_FORMULA], ERROR_KINETIC_LAW, self.t, SAFE_GLOBALS_RATE)
                if propensity < 0:
                    raise Exception("Negative propensity is denied")
                propensities.append(propensity)
            a0 = sum(propensities)
            if a0 == 0:
                break

            # New time
            y = random.random()
            tau = -math.log(y) / a0

            # Select delay events
            min_delay_time = self.t + tau

            to_eval_events = []
            to_remove = []

            for delay_time, event in self.pending_event_delay.values():
                # evaluation of the trigger at delay_time
                if self.evaluate_expr(event[TRIGGER_FORMULA], ERROR_TRIGGER, delay_time):
                    # take only the events with the least delay_time
                    if self.t < delay_time <= self.t + tau:
                        if delay_time < min_delay_time:
                            to_eval_events=[event]
                            min_delay_time = delay_time
                        elif delay_time == min_delay_time:
                            to_eval_events.append(event)
                else:
                    to_remove.append(event[ID])

            for event_id in to_remove:
                self.pending_event_delay.pop(event_id, None)

            self.t = min_delay_time

            # There are some events to evaluate
            if to_eval_events != []:
                to_remove = []
                for id_event, (_, event) in self.pending_event_delay.items():
                    # Delete the events whose trigger is False
                    if not self.evaluate_expr(event[TRIGGER_FORMULA], ERROR_TRIGGER, self.t):
                        event[VALUES_FROM_TRIGGER_TIME] = {}
                        to_remove.append(id_event)

                for event_id in to_remove:
                    self.pending_event_delay.pop(event_id, None)

            # Adding all the True-valuated events without delay or delay=0 at time t
            for event in self.parser.events:
                trigger_expr = event[TRIGGER_FORMULA]
                val_trigger = self.evaluate_expr(trigger_expr, ERROR_TRIGGER, self.t) # at time t
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
                    delay_time = self.evaluate_expr(delay, ERROR_DELAY, self.t)

                    if delay_time == 0:
                        # Triggered event's delay tag contains zero
                        to_eval_events.append(event)
                    elif delay_time > 0: # Triggered event's delay tag contains a value > 0
                        # Storing the values at trigger time for the input variables used in event assignments
                        if event[USE_VALUES_FROM_TRIGGER_TIME]:
                            for key_ea in event[VALUES_FROM_TRIGGER_TIME].keys():
                                if key_ea in self.parser.species:
                                    event[VALUES_FROM_TRIGGER_TIME][key_ea] = self.parser.species[key_ea]
                                elif key_ea in self.parser.parameters:
                                    event[VALUES_FROM_TRIGGER_TIME][key_ea] = self.parser.parameters[key_ea]

                        # Adding new event with delay
                        self.pending_event_delay[event[ID]] = (self.t + delay_time, event)

                    else:
                        raise Exception("The delay time is lesser than 0.")

                event[PREVIOUS] = True

            # There is one or more event(s) to evaluate and apply its EventAssigment
            if to_eval_events != []:
                for event in to_eval_events:
                    # Check that the trigger is still True
                    if self.evaluate_expr(event[TRIGGER_FORMULA], ERROR_TRIGGER, self.t):
                        # event[VALUES_FROM_TRIGGER_TIME] may be {} if eventAssigment use only constant or if
                        # useValuesFromTriggerTime="false"
                        self.apply_events_assigment(event[ID], event[LIST_OF_EVENT_ASSIGMENT], event[VALUES_FROM_TRIGGER_TIME])

                # Updating the evolution of species
                self.evolution[TIME].append(self.t)
                for s_id in self.parser.species:
                    self.evolution[s_id].append(self.parser.species[s_id])
                continue

            #Reaction choose between 0 and a0
            n = random.random() * a0
            cumulative = 0.0
            for i, a in enumerate(propensities):
                cumulative += a
                if n < cumulative:
                    chosen = i
                    break

            # Update of state
            r = self.parser.reactions[chosen]
            for s, stoich in r[REACTANTS].items():
                self.parser.species[s] -= stoich
                if self.parser.species[s] < 0:
                    self.parser.species[s] = 0
            for p, stoich in r[PRODUCTS].items():
                self.parser.species[p] += stoich
                if self.parser.species[p] < 0:
                    self.parser.species[p] = 0
            # Updating the evolution of species
            self.evolution[TIME].append(self.t)
            for s_id in self.parser.species:
                self.evolution[s_id].append(self.parser.species[s_id])

    def get_evolution(self):
        return self.evolution