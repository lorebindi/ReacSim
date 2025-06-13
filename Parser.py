from Constants import *

import libsbml
import matplotlib.pyplot as plt

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

        reaction[RATE] = formula

        #TODO: inferimento della legge cinetica da dati forniti tramite file csv.
        # 1. Come diciamo quale reazione deve passare da l'inferimento?
        # 2. Generazione punti dati da cui inferire?

        # Parametri
        params = {p.getId(): p.getValue() for p in model.getListOfParameters()}
        reaction[PARAMETERS] = params

        reactions.append(reaction)

        #TODO: gestione degli eventi.

    return reactions


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



























