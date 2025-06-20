from Constants import *
import libsbml

class Parser:
    model = None
    def __init__(self, file_name):
        self.file_name = file_name
        Parser.model = self.read_sbml_file()
        self.species = self.extract_species()
        self.parameters = self.extract_parameters()
        self.reactions = self.extract_reactions()
        self.events = self.extract_events()

    def read_sbml_file(self):
        reader = libsbml.SBMLReader()
        sbml_document = reader.readSBML(self.file_name)

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

    def extract_species(self):
        species = {}
        for s in self.model.getListOfSpecies():
            if s.isSetInitialAmount():  # We have the absolute initial amount.
                species[s.getId()] = s.getInitialAmount()
            elif s.isSetInitialConcentration():  # We have to do: amount = volume * InitialConcetrantion.
                compartment = self.model.getCompartment(s.getCompartment())
                volume = compartment.getSize()
                try:
                    species[s.getId()] = volume * s.getInitialConcentration()
                except TypeError:  # Example: value * None
                    raise Exception("Unable to extract InitialAmount from InitialConcentration.")
            else:
                raise Exception("Neither InitialAmount nor InitialConcentration are set.")
        return species

    def extract_parameters(self):
        parameters = {p.getId(): p.getValue() for p in self.model.getListOfParameters()}
        parameters.update({c.getId(): 1.0 for c in self.model.getListOfCompartments()})
        return parameters

    def extract_reactions(self):
        reactions = []
        for r in self.model.getListOfReactions():
            reactions.append(Parser.Reaction(r).get_reaction_as_dict())
        if len(reactions) == 0: raise Exception("No reactions found.")
        return reactions

    def extract_events(self):
        events = []
        for e in self.model.getListOfEvents():
            events.append(Parser.Event(e).get_event_as_dict())
        return events

    class Reaction:
        def __init__(self, reaction):
            self.reaction = reaction
            self.id = self.reaction.getId()
            self.reactants = {}
            self.products = {}
            self.rate_formula = None
            self.kinetic_law = None
            self.extract_reaction()

        def extract_reaction(self):
            # Reactans' coefficients
            reactants_ids = []
            for sr in self.reaction.getListOfReactants():
                reactants_ids.append(sr.getId())  # usefull to check if the kinetic law is mass action law
                self.reactants[sr.getSpecies()] = sr.getStoichiometry()

            # Products' coefficients
            for sp in self.reaction.getListOfProducts():
                self.products[sp.getSpecies()] = sp.getStoichiometry()

            # Kinectic law
            self.kinetic_law = self.reaction.getKineticLaw()
            if self.kinetic_law is None:
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

            '''
            Quando si valuta la costante cinetica:
            1. Si controlla se è presente nei parametri:
                a. se è presente nei parametri si continua.
                b. se non è presente nei parametri si cerca il file. 
                    I. Se il file non c'è schianti
                    II. Altrimenti il valore inferito lo aggiungi a model
            '''

            # TODO: proviamo a non scorrere l'albero più volte ma una volta sola? Unendo:
            # - validate_mass_action_kinetic_law
            # - contains_identifier
            self.validate_mass_action_kinetic_law()

            for compartment in Parser.model.getListOfCompartments():
                if self.contains_identifier(compartment.getId()):
                    raise Exception("Kinetic law invalid for us.")

            self.rate_formula = self.kinetic_law.getFormula()

            # TODO: inferenza della legge cinetica da dati forniti tramite file csv.
            # 1. Come diciamo quale reazione deve passare da l'inferenza?
            # 2. Generazione punti dati da cui inferire?

        def validate_mass_action_structure(self, ast, kinetic_constant_found=0):
            if self.kinetic_law is None:
                raise Exception("Kinetic law not set.")

            if ast is None:
                raise Exception("AST is None.")

            # Always * at the root
            if ast.getType() != libsbml.AST_TIMES:
                raise Exception("AST's root is not TIMES(*).")

            #  all child nodes (multiplicands)
            for child in [ast.getChild(i) for i in range(ast.getNumChildren())]:
                if child.getType() == libsbml.AST_TIMES:
                    kinetic_constant_found = self.validate_mass_action_structure(child, kinetic_constant_found)
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
                            base.getName() not in self.reactants or
                            exponent.getType() not in TYPE_NUMBER or
                            exponent.getValue() != self.reactants[base.getName()]
                    ):
                        raise Exception("AST_FUNCTION_POWER node is written in the wrong way.")
                else:
                    raise Exception(f"AST node not expected (type: {child.getType()}, name: {child.getName()}, "
                                    f"value: {child.getValue()}) in: {self.kinetic_law.getFormula()}")

            return kinetic_constant_found

        def validate_mass_action_kinetic_law(self):

            if self.kinetic_law is None:
                raise Exception("Unable to determine mass action kinetic law.")

            ast_root = self.kinetic_law.getMath()  # get the Abstract Syntax Tree (AST) of the kinetic law

            if ast_root is None:
                raise Exception("AST is None.")

            # If the reaction is a synthesis (  -> something) then the kinetic law is equal
            # to the kinetic constant
            if ast_root.getType() == libsbml.AST_NAME:
                return

            if self.validate_mass_action_structure(ast_root) == 0:
                raise Exception("Kinetic constant absent.")

        def contains_identifier(self, target_id):
            math_ast = self.kinetic_law.getMath()
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

        def get_reaction_as_dict(self):
            return {
                ID: self.id,
                REACTANTS: self.reactants,
                PRODUCTS: self.products,
                RATE_FORMULA: self.rate_formula
            }
    class Event:
        def __init__(self,event):
            self.event = event
            self.id=event.getId()
            self.trigger_formula = {}
            self.previous = False
            self.list_of_event_assigment = []
            self.delay_formula = None
            self.use_values_from_trigger_time = None
            self.value_from_trigger_time = {}

            self.use_trigger_values = False

            self.extract_event()

        def extract_event(self):
            # In SBML Level 2 Version 4, the useValuesFromTriggerTime attribute is only valid if a delay is also defined.
            # This block checks that constraint and sets a flag (use_trigger_values) to True only if the event has a delay
            # and set UseValuesFromTriggerTime to True or is not setted.
            if Parser.model.getVersion() == 4:
                '''
                # Events with useValuesFromTriggerTime and without delay are not allowed
                a = e.isSetDelay()
                b = e.isSetUseValuesFromTriggerTime()
                if not e.isSetDelay() and e.isSetUseValuesFromTriggerTime():
                    raise Exception("Is not possible having an event with useValuesFromTriggerTime attribute without the delay")
                '''

                if self.event.isSetDelay() and self.event.getUseValuesFromTriggerTime():
                    self.use_trigger_values = True

                self.use_values_from_trigger_time = self.use_trigger_values

            trigger = self.event.getTrigger()
            if trigger is None:
                raise Exception("Trigger not found.")
            # Check if trigger condition is valid
            ast_trigger = trigger.getMath()
            self.validate_trigger_boolean_expr(ast_trigger)
            self.trigger_formula = (libsbml.formulaToL3String(ast_trigger).replace("&&", " and ")
                                    .replace("||", " or ").replace("!", " not "))

            # TODO: controllare che le variable siano tutte diverse
            for event_assignment in self.event.getListOfEventAssignments():
                string_of_variable = event_assignment.getVariable()
                variable = Parser.model.getElementBySId(string_of_variable)
                if variable is None:
                    raise Exception(f"Variable '{string_of_variable}' not found in model.")

                # Ottieni l'attributo 'units' dalla variabile (solo se esiste)
                unit_id = variable.getUnits() if variable.isSetUnits() else None
                if unit_id is None:
                    raise Exception("Undefined unit definition for variable.")
                unit_definition_variable = Parser.model.getUnitDefinition(unit_id)

                # Check the eventAssignment
                trigger_value_dict = self.validate_event_assigment(event_assignment.getMath(),unit_definition_variable)
                if self.use_trigger_values:
                    self.value_from_trigger_time.update(trigger_value_dict)
                self.list_of_event_assigment.append(event_assignment)

            # TODO: aggiungere priority (supportata solo al livello 3)
            # if version is 3 or there isn't the Delay tag, then we append the event and terminate by evaluating the single event
            if not (Parser.model.getVersion() == 3 or not self.event.isSetDelay()):
                delay = self.event.getDelay()
                ast_delay = delay.getMath()

                bool_constant_found, bool_parameter_found = self.validate_delay(ast_delay)
                if bool_constant_found and not bool_parameter_found:
                    raise Exception("Incorrect use of a constant for the delay.")

                self.delay_formula = libsbml.formulaToString(ast_delay)



        def is_valid_identifier(self, name):
            return (
                    Parser.model.getSpecies(name) is not None or
                    Parser.model.getParameter(name) is not None or
                    Parser.model.getCompartment(name) is not None or
                    Parser.model.getSpeciesReference(name) is not None
            )

        def validate_trigger_boolean_expr(self, ast):
            if ast is None:
                raise Exception("AST is None.")

            if ast.getType() not in TYPE_TRIGGER:
                raise Exception("AST must be logical or relational operator.")

            for child in [ast.getChild(i) for i in range(ast.getNumChildren())]:
                if child.getType() in TYPE_TRIGGER:
                    self.validate_trigger_boolean_expr(child)
                elif child.getType() == libsbml.AST_NAME:
                    if not self.is_valid_identifier(child.getName()):
                        raise Exception("Invalid identifier.")
                    return
                elif child.getType() in TYPE_NUMBER or child.getType() == libsbml.AST_NAME_TIME:
                    return
                else:
                    raise Exception(f"AST node not expected (type: {child.getType()}, name: {child.getName()}, "
                                    f"value: {child.getValue()}) in: {ast.getFormula()}")
        '''
       Supported operations:
       1. PLUS and MINUS between variables and constans. 
       2. Logical comparisons and conditions using variables, constants, and operators (LT, GT, EQ, AND, OR, NOT).
       3. Use of <csymbol> time to refer to the current simulation time.
       '''
        # TODO: estendere alle altre operazioni (TIMES, DIV) e math functions
        def validate_event_assigment(self, ast, unit_definition_variable):
            event_assignment_input_vars = {} # store the set of variables used in the eventAssignment's Math expression
            if ast is None:
                raise Exception("AST is None.")
            if ast.getType() in TYPE_NUMBER: # variable = constant
                pass
            elif ast.getType() == libsbml.AST_NAME:
                string_of_element = ast.getName()
                if self.use_trigger_values:
                    event_assignment_input_vars[string_of_element] = None
                element = Parser.model.getElementBySId(string_of_element)
                if element is None:
                    raise Exception(f"Symbol '{string_of_element}' not found in model.")

                if element.getTypeCode() not in TYPE_CODE: # element != (Species, Parameter, Compartment))
                    raise Exception(f"Unsupported element type for symbol '{string_of_element}'")

                if not element.isSetUnits():
                    raise Exception(f"Unsupported element type for symbol '{string_of_element}'")

                # Check if the units are different
                if not libsbml.UnitDefinition.areEquivalent(Parser.model.getUnitDefinition(element.getUnits()),
                                                            unit_definition_variable):
                    raise Exception(f"Unit mismatch: '{string_of_element}' has units different from assigned variable.")
            else:
                if ast.getType() not in TYPE_OP:
                    raise Exception("Only AST_PLUS and AST_MINUS are supported in the EventAssigment.")

                for i in range(ast.getNumChildren()):
                    return_value = self.validate_event_assigment(ast.getChild(i),unit_definition_variable)
                    if self.use_trigger_values:
                        event_assignment_input_vars.update(return_value)
            return event_assignment_input_vars

        """
        Validates <delay> expressions to ensure they:
        1. The expression inside the <math> block evaluate to a duration in seconds, 
            which must match the model's time units without providing any unit conversion.
        2. The expression must consist only of a single parameter or mathematical 
            operations (only SUM and MINUS) involving parameters and/or numerical 
            constants—no complex expressions or unsupported constructs are allowed.
        """
        def validate_delay(self, ast, constant_found = False, parameter_found = False):
            if ast is None:
                raise Exception("AST is None.")
            if ast.getType() in TYPE_NUMBER: # variable = constant
                return True, parameter_found
            elif ast.getType() == libsbml.AST_NAME:
                string_of_element = ast.getName()
                element = Parser.model.getElementBySId(string_of_element)
                if element is None:
                    raise Exception(f"Symbol '{string_of_element}' not found in model.")

                if element.getTypeCode() != libsbml.SBML_PARAMETER:  # element can only be a parameter
                    raise Exception(f"Unsupported element type for symbol '{string_of_element}'")

                if not element.isSetUnits():
                    raise Exception(f"Unsupported element type for symbol '{string_of_element}'")

                # Check if the units are seconds
                delay_units = element.getUnits()
                delay_ud = Parser.model.getUnitDefinition(delay_units)

                if delay_units == 'second' and delay_ud is None:
                    return constant_found, True
                elif delay_units != 'second' and delay_ud is not None:
                    list = delay_ud.getListOfUnits()
                    if len(list) != 1:
                        raise Exception(f"Too many units in {string_of_element}.")
                    if list[0].getKind() != libsbml.UNIT_KIND_SECOND:
                        raise Exception(f"Unsupported unit type for symbol '{string_of_element}'. Only second supported.")
                    return constant_found, True
                #TODO: estendere anche a unitDefinition custom per millisecond e minutes
                raise Exception(f"Unsupported unit type for symbol '{string_of_element}'. Only second supported.")

            else:
                if ast.getType() not in TYPE_OP:
                    raise Exception("Only AST_PLUS and AST_MINUS are supported in the EventAssigment.")

                for i in range(ast.getNumChildren()):
                    return_constant_found, return_parameter_found = self.validate_delay(ast.getChild(i))
                    constant_found = return_constant_found or constant_found
                    parameter_found = return_parameter_found or parameter_found

                return constant_found, parameter_found


        def get_event_as_dict(self):
            return {
                ID: self.id,
                TRIGGER_FORMULA: self.trigger_formula,
                PREVIOUS: self.previous,    # Value of trigger at t-tau (previous value)
                LIST_OF_EVENT_ASSIGMENT: self.list_of_event_assigment,
                DELAY_FORMULA: self.delay_formula,
                USE_VALUES_FROM_TRIGGER_TIME: self.use_values_from_trigger_time,
                VALUES_FROM_TRIGGER_TIME: self.value_from_trigger_time
            }
