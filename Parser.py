import math
import random
import libsbml
import re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

TIME = '_time_'
ID = "id"
REACTANTS = "reactants"
PRODUCTS = "products"
RATE = "rate"
PARAMETERS = "parameters"

def readSBMLfile(filename):
    reader = libsbml.SBMLReader()
    SBMLdocument = reader.readSBML(filename)

    # Check the validity of the document
    if SBMLdocument.getNumErrors() > 0:
        print("Reading error")
        SBMLdocument.printErrors()
        return None

    # Make sure the file uses level 3
    if SBMLdocument.getLevel() != 2:
        print("The file isn't level 2.")
        return None

    # Error checking
    # TODO eventualmente aggiungere nei parametri della funzione il campo "inConversion=True", per abilitare la
    #  conversione se non compatibile
    if SBMLdocument.checkL2v3Compatibility() > 0 and SBMLdocument.checkL2v4Compatibility() > 0:
        print("The file is not compatible with level 2 version 3/4")
        return None

    model = SBMLdocument.getModel()
    if model is None:
        print("Unable to create Model object.")
        return None

    return model

def extract_species(model):
    species = {}
    for s in model.getListOfSpecies():
        if s.isSetInitialAmount(): # We have the absolute initial amount.
            #species[s.getId()] = (s.getName(), s.getInitialAmount())
            species[s.getId()] = s.getInitialAmount()
        elif s.isSetInitialConcentration(): # We have to do: amount = volume * InitialConcetrantion.
            compartment = model.getCompartment(s.getCompartment())
            volume = compartment.getSize()
            try:
                # species[s.getId()] = (s.getName(), volume * s.getInitialConcentration())
                species[s.getId()] = volume * s.getInitialConcentration()
            except TypeError as e: # Example: value * None
                # species[s.getId()] = (s.getName(), None)
                # species[s.getId()] = None
                return None
        else:
            return None
    return species

def extract_reactions(model):
    reactions = []

    for r in model.getListOfReactions():
        reaction = {
            ID: r.getId(),
            REACTANTS: {},
            PRODUCTS: {},
            RATE: None,
            PARAMETERS: {}
        }

        # Reactans' coefficients
        reactants_ids = []
        for sr in r.getListOfReactants():
            reactants_ids.append(sr.getId()) #usefull to check if the kinectik law is mass action law
            reaction[REACTANTS][sr.getSpecies()] = sr.getStoichiometry()

        # Products' coefficients
        for sp in r.getListOfProducts():
            reaction[PRODUCTS][sp.getSpecies()] = sp.getStoichiometry()

        # Kinectic law
        kinetic_law = r.getKineticLaw()
        if kinetic_law is None:
            return None

        formula = kinetic_law.getFormula()
        '''if not is_mass_action(formula, reactants_ids):
            return None'''

        reaction[RATE] = formula

        #TODO: supporto ad altre o qualsiasi legge cineticha

        #TODO: inferimento della legge cinetica da dati forniti tramite file csv.
        # 1. Come diciamo quale reazione deve passare da l'inferimento?
        # 2. Generazione punti dati da cui inferire?

        # Parametri
        params = {p.getId(): p.getValue() for p in model.getListOfParameters()}
        reaction[PARAMETERS] = params

        reactions.append(reaction)

    return reactions

def is_mass_action(formula, reactant_ids):
    if formula is None:
        return False

    # Rimuovi spazi
    formula = formula.replace(" ", "")

    # La formula deve contenere solo una costante e i reagenti (es. "k*A*B")
    pattern = r'^[\w\.]+(\*[\w\.]+)*$'
    if not re.fullmatch(pattern, formula):
        return False

    # Tutti i reagenti devono essere presenti
    for rid in reactant_ids:
        if rid not in formula:
            return False

    return True

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

    state = extract_species(model)
    if state is None:
        return None

    reactions = extract_reactions(model)
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

def plot(evolution, save_dir=None):
    if evolution is None: return None
    if TIME not in evolution: return None

    time_values = evolution[TIME]
    specie_names = [k for k in evolution if k != TIME]

    plt.figure(figsize=(10, 6))

    for specie in specie_names:
        valori = evolution[specie]
        plt.plot(time_values, valori, label=specie, marker='o')

    plt.xlabel("Time")
    plt.ylabel("Molecule count")
    plt.title("Evolution of species over time")
    plt.legend(loc="upper right")
    plt.grid(True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.show()

    # Horizontal line graph with letters on the Y-axis



























