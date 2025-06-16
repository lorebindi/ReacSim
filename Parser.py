from Constants import *
import libsbml

def read_sbml_file(filename):
    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML(filename)

    # Check the validity of the document
    if sbml_document.getNumErrors() > 0:
        sbml_document.printErrors()
        raise Exception("Reading error")

    # Make sure the file uses level 3
    if sbml_document.getLevel() != 2:
        raise Exception("The file isn't level 2.")

    # Error checking
    # TODO eventualmente aggiungere nei parametri della funzione il campo "inConversion=True", per abilitare la
    #  conversione se non compatibile
    if sbml_document.checkL2v3Compatibility() > 0 and sbml_document.checkL2v4Compatibility() > 0:
        raise Exception("The file is not compatible with level 2 version 3/4")

    model = sbml_document.getModel()
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
            except TypeError: # Example: value * None
                raise Exception("Unable to extract InitialAmount from InitialConcentration.")
        else:
            raise Exception("Neither InitialAmount nor InitialConcentration are set.")
    return species

def extract_parameters(model):
    parameters = {p.getId(): p.getValue() for p in model.getListOfParameters()}
    parameters.update({c.getId(): 1.0 for c in model.getListOfCompartments()})
    return parameters

def print_ast(node, indent=0):
    if node is None:
        return

    print('  ' * indent + f"type: {node.getType()}, name: {node.getName()}, value: {node.getValue()}")

    for i in range(node.getNumChildren()):
        print_ast(node.getChild(i), indent + 1)

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

def is_mass_action_kinetic_law(reactants):
    if kinetic_law is None:
        raise Exception("Unable to determine mass action kinetic law.")

    ast_root = kinetic_law.getMath() # get the Abstract Syntax Tree (AST) of the kinetic law

    if ast_root is None:
        raise Exception("AST is None.")

    # If the reaction is a synthesis (  -> something) then the kinetic law is equal
    # to the kinetic constant
    if ast_root.getType() == libsbml.AST_NAME:
        return True

    def validate_mass_action_structure(ast, kinetic_constant_found=0):
        if kinetic_law is None:
            raise Exception("Kinetic law not set.")

        if ast is None:
            raise Exception("AST is None.")

        # Always * at the root
        if ast.getType() != libsbml.AST_TIMES:
            raise Exception("AST's root is not TIMES(*).")

        #  all child nodes (multiplicands)
        for child in [ast.getChild(i) for i in range(ast.getNumChildren())]:
            if child.getType() == libsbml.AST_TIMES:
                kinetic_constant_found = validate_mass_action_structure(child, kinetic_constant_found)
                # There must be at least one parameter (rate constant)
                if kinetic_constant_found > 1:
                    raise Exception("Too many kinetic constant.")
            elif child.getType() == libsbml.AST_NAME:  # kinetic constant
                kinetic_constant_found += 1
            elif child.getType() == libsbml.AST_FUNCTION_POWER:
                base = child.getChild(0)
                exponent = child.getChild(1)
                if (
                        base.getType() != libsbml.AST_NAME or
                        base.getName() not in reactants or(
                        exponent.getType() != libsbml.AST_INTEGER and
                        exponent.getType() != libsbml.AST_REAL and
                        exponent.getType() != libsbml.AST_REAL_E and
                        exponent.getType() != libsbml.AST_RATIONAL
                        ) or
                        exponent.getValue() != reactants[base.getName()]
                ):
                    raise Exception("AST_FUNCTION_POWER node is written in the wrong way.")
            else:
                raise Exception(f"AST node not expected (type: {child.getType()}, name: {child.getName()}, "
                                f"value: {child.getValue()}) in: {kinetic_law.getFormula()}")

        return kinetic_constant_found

    if validate_mass_action_structure(ast_root) == 0:
        raise Exception("Kinetic constant absent.")
    return True

'''    
    1 - reactant
    AST_TIMES
    ├── AST_NAME("beta")
    └── AST_POWER
        ├── AST_NAME("Prey")
        └── AST_CN(1)

    2 - reactants
    AST_TIMES
    ├── AST_TIMES
    │   ├── AST_NAME("beta")
    │   └── AST_POWER
    │       ├── AST_NAME("Prey")
    │       └── AST_CN(1)
    └── AST_POWER
        ├── AST_NAME("Predator")
        └── AST_CN(1)

    3 - reactants
    (times)
    ├── (times)
    │    ├── (times)
    │    │     ├── beta
    │    │     └── pow(Prey, 1)
    │    └── pow(Predator, 1)
    └── pow(SomeOther, 1)

    ...
    '''

def extract_reactions(model):
    global kinetic_law
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

        #TODO: proviamo a non scorrere l'albero più volte ma una volta sola? Unendo:
        # - is_mass_action_kinetic_law
        # - contains_identifier
        is_mass_action_kinetic_law(reaction[REACTANTS])

        for compartment in model.getListOfCompartments():
            if contains_identifier(kinetic_law.getMath(), compartment.getId()):
                raise Exception("Kinetic law invalid for us.")

        reaction[RATE] = kinetic_law.getFormula()

        #TODO: inferenza della legge cinetica da dati forniti tramite file csv.
        # 1. Come diciamo quale reazione deve passare da l'inferenza?
        # 2. Generazione punti dati da cui inferire?

        reactions.append(reaction)

        #TODO: gestione degli eventi.

    if len(reactions) == 0: raise Exception("No reactions found.")

    return reactions
