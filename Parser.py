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
    def validate_mass_action_kinetic_law(reactants):
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
                            base.getName() not in reactants or
                            exponent.getType() not in TYPE_NUMBER or
                            exponent.getValue() != reactants[base.getName()]
                    ):
                        raise Exception("AST_FUNCTION_POWER node is written in the wrong way.")
                else:
                    raise Exception(f"AST node not expected (type: {child.getType()}, name: {child.getName()}, "
                                    f"value: {child.getValue()}) in: {kinetic_law.getFormula()}")

            return kinetic_constant_found

        if kinetic_law is None:
            raise Exception("Unable to determine mass action kinetic law.")

        ast_root = kinetic_law.getMath()  # get the Abstract Syntax Tree (AST) of the kinetic law

        if ast_root is None:
            raise Exception("AST is None.")

        # If the reaction is a synthesis (  -> something) then the kinetic law is equal
        # to the kinetic constant
        if ast_root.getType() == libsbml.AST_NAME:
            return

        if validate_mass_action_structure(ast_root) == 0:
            raise Exception("Kinetic constant absent.")

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
        # - validate_mass_action_kinetic_law
        # - contains_identifier
        validate_mass_action_kinetic_law(reaction[REACTANTS])

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

def extract_events(model):
    def is_valid_identifier(name):
        return (
                model.getSpecies(name) is not None or
                model.getParameter(name) is not None or
                model.getCompartment(name) is not None or
                model.getSpeciesReference(name) is not None
        )

    def validate_trigger_boolean_expr(ast):
        if ast is None:
            raise Exception("AST is None.")

        if ast.getType() not in TYPE_TRIGGER:
            raise Exception("AST must be logical or relational operator.")

        for child in [ast.getChild(i) for i in range(ast.getNumChildren())]:
            if child.getType() in TYPE_TRIGGER:
                validate_trigger_boolean_expr(ast)
            elif child.getType() == libsbml.AST_NAME:
                if not is_valid_identifier(child.getName()):
                    raise Exception("Invalid identifier.")
                return
            elif child.getType() in TYPE_NUMBER:
                return
            else:
                raise Exception(f"AST node not expected (type: {child.getType()}, name: {child.getName()}, "
                                f"value: {child.getValue()}) in: {ast.getFormula()}")

    '''
    Supported operations:
    1. PLUS and MINUS between variables and costants.
    2. PLUS and MINUS between variables. 
    '''
    # TODO: estendere alle altre operazioni
    def validate_event_assigment(ast):
        if ast is None:
            raise Exception("AST is None.")
        if ast.getType() in TYPE_NUMBER: # variable = constant
            return
        elif ast.getType() == libsbml.AST_NAME:
            string_of_element = ast.getName()
            element = model.getElementBySId(string_of_element)

            if element is None:
                raise Exception(f"Symbol '{string_of_element}' not found in model.")

            if element.getTypeCode() not in TYPE_CODE: # element != (Species, Parameter, Compartment))
                raise Exception(f"Unsupported element type for symbol '{string_of_element}'")

            if not element.isSetUnits():
                raise Exception(f"Unsupported element type for symbol '{string_of_element}'")

            # Check if the units are different
            if not libsbml.UnitDefinition.areEquivalent(model.getUnitDefinition(element.getUnits()),
                                                        unit_definition_variable):
                raise Exception(f"Unit mismatch: '{string_of_element}' has units different from assigned variable.")
        else:
            if ast.getType() not in TYPE_OP:
                raise Exception("Only AST_PLUS and AST_MINUS are supported in the EventAssigment.")

            for i in range(ast.getNumChildren()):
                validate_event_assigment(ast.getChild(i))

    events = []
    for e in model.getListOfEvents():
        event = {
            TRIGGER: {},
            PREVIOUS: False,
            LIST_OF_EVENT_ASSIGMENT: [],
            DELAY: None,
            PRIORITY: None
        }

        trigger = e.getTrigger()
        if trigger is None:
            raise Exception("Trigger not found.")
        # Check if trigger condition is valid
        ast_trigger = trigger.getMath()
        validate_trigger_boolean_expr(ast_trigger)
        event[TRIGGER] = libsbml.formulaToString(ast_trigger)

        for event_assignment in e.getListOfEventAssignments():
            string_of_variable = event_assignment.getVariable()
            variable = model.getElementBySId(string_of_variable)
            if variable is None:
                raise Exception(f"Variable '{string_of_variable}' not found in model.")

            # Ottieni l'attributo 'units' dalla variabile (solo se esiste)
            unit_id = variable.getUnits() if variable.isSetUnits() else None
            if unit_id is None:
                raise Exception("Undefined unit definition for variable.")
            unit_definition_variable = model.getUnitDefinition(unit_id)

            validate_event_assigment(event_assignment.getMath())

            event[LIST_OF_EVENT_ASSIGMENT].append(event_assignment)

        # TODO: aggiungere delay e priority

        events.append(event)

        '''
        1. le unità di misura.
        2. variable contiene riferimenti a un Compartment, Species, SpeciesReference, Parameter
            Permettiamo: 
            a. costanti (<cn>) singole, variabili (species, ..., compartment) (<ci>) singole.
            b. somme, differenze, moltiplicazioni
        3. il blocco math restituisca un valore con la solità unità di misura di variable
        '''
    return events
