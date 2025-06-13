from Constants import *

import libsbml
import re

def read_SBML_file(filename):
    reader = libsbml.SBMLReader()
    SBMLdocument = reader.readSBML(filename)

    # Check the validity of the document
    if SBMLdocument.getNumErrors() > 0:
        raise Exception("Reading error")

    # Make sure the file uses level 3
    if SBMLdocument.getLevel() != 2:
        raise Exception("The file isn't level 2.")

    # Error checking
    # TODO eventualmente aggiungere nei parametri della funzione il campo "inConversion=True", per abilitare la
    #  conversione se non compatibile
    if SBMLdocument.checkL2v3Compatibility() > 0 and SBMLdocument.checkL2v4Compatibility() > 0:
        raise Exception("The file is not compatible with level 2 version 3/4")

    model = SBMLdocument.getModel()
    if model is None:
        raise Exception("Unable to create Model object.")

    return model

def extract_species(model):
    species = {}
    for s in model.getListOfSpecies():
        if s.isSetInitialAmount(): # We have the absolute initial amount.
            species[s.getId()] = s.getInitialAmount()
        elif s.isSetInitialConcentration(): # We have to do: amount = volume * InitialConcetrantion.
            compartment = model.getCompartment(s.getCompartment())
            volume = compartment.getSize()
            try:
                species[s.getId()] = volume * s.getInitialConcentration()
            except TypeError as e: # Example: value * None
                raise Exception("Unable to extract InitialAmount from InitialConcentration.")
        else:
            raise Exception("Neither InitialAmount nor InitialConcentration are set.")
    return species

def extract_parameters(model):
    parameters = {p.getId(): p.getValue() for p in model.getListOfParameters()}
    parameters.update({c.getId(): 1.0 for c in model.getListOfCompartments()})
    return parameters

def contains_identifier(math_ast, target_id):
    if math_ast is None:
        return False

    # Recursively search for a node with name target_id
    def traverse(node):
        if node is None:
            return False
        if node.getType() == libsbml.AST_NAME:
            return node.getName() == target_id
        for i in range(node.getNumChildren()):
            if traverse(node.getChild(i)):
                return True
        return False

    return traverse(math_ast)

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

def extract_reactions(model):
    reactions = []

    for r in model.getListOfReactions():
        reaction = {
            ID: r.getId(),
            REACTANTS: {},
            PRODUCTS: {},
            RATE: None
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
            raise Exception("Kinetic law not found.")

        '''
        We support ONLY mass-action laws as kinetic laws. If the kinetic law
        is not a mass-action law, the program will abort.
        
        If the kinetic law contains a compartment, then all species involved
        in the reaction MUST belong to that compartment.
        
        If a species does not specify an InitialAmount, we ALWAYS compute it
        as: volume * InitialConcentration.
        
        If the kinetic law contains a compartment and the species are initially
        defined only with InitialConcentration, we simply set the compartment 
        value to 1 and use the previously computed InitialAmount values
        '''
        if not is_mass_action(kinetic_law.getFormula(), reactants_ids):
            raise Exception("Reaction is not mass action.")

        for compartment in model.getListOfCompartments():
            if contains_identifier(kinetic_law.getMath(), compartment.getId()):
                    raise Exception("Kinetic law invalid for us.")

        formula = kinetic_law.getFormula()

        reaction[RATE] = formula

        #TODO: inferimento della legge cinetica da dati forniti tramite file csv.
        # 1. Come diciamo quale reazione deve passare da l'inferimento?
        # 2. Generazione punti dati da cui inferire?

        reactions.append(reaction)

        #TODO: gestione degli eventi.

    return reactions
