# Standard libraries
import csv
import os

# Third part libraries
import libsbml
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


# Local modul
import Constants as constants
import Graph_generation as graphgen
import Gillespie_events as Gillespie
from main import T_MAX


class Parser:
    def __init__(self, file_name, path_inference_directory = None):
        self.file_name = file_name
        self.path_inference_directory = path_inference_directory
        self.model = self.read_sbml_file()
        self.species = self.extract_species()
        self.parameters = self.extract_parameters()
        self.dfs = {}
        self.reactions = self.extract_reactions()
        self.events = self.extract_events()

    ''' This method reads a SBML file, checks its validity and 
    the compatibility with SBML L2V3 and L2V4, and store the model.'''
    def read_sbml_file(self):
        reader = libsbml.SBMLReader()
        sbml_document = reader.readSBML(self.file_name)

        # Check the validity of the document
        if sbml_document.getNumErrors() > 0:
            sbml_document.printErrors()
            raise Exception("Reading error")

        # Make sure that the file uses level 2
        if sbml_document.getLevel() != 2:
            raise Exception("The file isn't level 2.")

        # Error checking
        if sbml_document.checkL2v3Compatibility() > 0 and sbml_document.checkL2v4Compatibility() > 0:
            raise Exception("The file is not compatible with level 2 version 3/4")

        model = sbml_document.getModel()
        if model is None:
            raise Exception("Unable to create Model object.")

        return model

    ''' This method extracts the species from the model. It returns 
    a dictionary with the species IDs as keys and the respective initial 
    amounts as values. '''
    def extract_species(self):
        species = {}
        for s in self.model.getListOfSpecies():
            if s.isSetInitialAmount() and not s.isSetInitialConcentration() and s.getHasOnlySubstanceUnits():  # We have the absolute initial amount.
                species[s.getId()] = s.getInitialAmount()
            else:
                if s.isSetInitialConcentration():
                    raise Exception("Initial Concentration cannot be set.")
                if not s.hasOnlySubstanceUnits():
                    raise Exception("Substance units cannot be unset.")
                raise Exception("InitialAmount are not set.")
        return species

    ''' This method extracts the parameters from the model. It returns 
    a dictionary with the parameters IDs as keys and the respective 
    values as values. '''
    def extract_parameters(self):
        parameters = {p.getId(): p.getValue() for p in self.model.getListOfParameters()}
        #parameters.update({c.getId(): 1.0 for c in self.model.getListOfCompartments()})
        return parameters

    '''This method extracts the reactions from the model. It returns 
    a dictionary with the reactions IDs as keys and the respective 
    details as values.'''
    def extract_reactions(self):
        reactions = []
        for r in self.model.getListOfReactions():
            reaction=Reaction(r, self)
            reactions.extend(reaction.get_reaction_as_dict())
            if reaction.constant_inferred_name is not None:
                self.dfs[reaction.constant_inferred_name] = reaction.df_csv

        if len(reactions) == 0: raise Exception("No reactions found.")
        return reactions

    '''This method extracts the events from the model. It returns 
    a dictionary with the events IDs as keys and the respective 
    details as values.'''
    def extract_events(self):
        events = []
        for e in self.model.getListOfEvents():
            events.append(Event(e, self).get_event_as_dict())
        return events

    '''This method extracts and return the stochastic rate constant name from
    kinetic law's ast.'''
    def extract_stochastic_rate_constant_name(self, ast):
        if ast is None:
            raise Exception("AST is None.")

        # Always * at the root
        if ast.getType() != libsbml.AST_TIMES:
            raise Exception("AST's root is not TIMES(*).")

        #  all child nodes (multiplicands)
        for child in [ast.getChild(i) for i in range(ast.getNumChildren())]:
            if child.getType() == libsbml.AST_NAME:  # stochastic rate constant
                return child.getName()
            return self.extract_stochastic_rate_constant_name(child)

        raise Exception("Kinetic constants not found.")

    '''This method performs multiple stochastic simulations (using the Gillespie 
    algorithm) only on the specified reaction (by building a temporal model that 
    contains its associated components: species, parameters, compartments and 
    units). It interpolates the species amounts over time, averages the results 
    across the specified number of executions, and saves the mean values to a CSV file.'''
    def export_mean_species_counts_csv(self, reaction_id, executions):
        reaction = self.model.getReaction(reaction_id)
        if reaction is None:
            raise Exception(f"Reaction '{reaction_id}' not found.")

        # Create a new SBML document and model
        new_document = libsbml.SBMLDocument(self.model.getLevel(), self.model.getVersion())
        new_model = new_document.createModel()
        new_model.setId(f"single_reaction_{reaction_id}_model")

        kinetic_law = reaction.getKineticLaw()
        if kinetic_law is None:
            raise Exception("Kinetic law not found for the reaction.")

        kinetic_constant_name = self.extract_stochastic_rate_constant_name(kinetic_law.getMath())

        # Copy the units, compartments, species, parameters involved
        involved_species = set()
        for reactant in reaction.getListOfReactants():
            involved_species.add(reactant.getSpecies())
        for product in reaction.getListOfProducts():
            involved_species.add(product.getSpecies())

        # Copy compartments needed
        compartments = set()
        for species_id in involved_species:
            species = self.model.getSpecies(species_id)
            if species is None:
                continue
            comp_id = species.getCompartment()
            compartments.add(comp_id)

        for comp_id in compartments:
            comp = self.model.getCompartment(comp_id)
            if comp:
                new_model.addCompartment(comp.clone())

        # Copy of species
        unit_ids = set()
        for species_id in involved_species:
            species = self.model.getSpecies(species_id)
            if species:
                new_model.addSpecies(species.clone())
                if species.isSetUnits():
                    unit_ids.add(species.getUnits())

        for param in self.model.getListOfParameters():
            if param.getId() == kinetic_constant_name:
                new_model.addParameter(param.clone())
                if param.isSetUnits():
                    unit_ids.add(param.getUnits())

        # Copy of the reaction
        new_model.addReaction(reaction.clone())

        # Add all required unit definitions to the new model
        for unit_id in unit_ids:
            unit_def = self.model.getUnitDefinition(unit_id)
            if unit_def:
                new_model.addUnitDefinition(unit_def.clone())

        # Ml saving
        path = f"./Example/Generated/{reaction_id}.xml"
        writer = libsbml.SBMLWriter()
        writer.writeSBMLToFile(new_document, path)

        # Interpolation of the points
        sum_per_specie = {}
        species_to_time_dict = {} #final result
        #TODO if the execution ends much earlier than T_MAX, I have to deal with the fact that values are absent at those points.
        time_query = list(range(0, T_MAX + 1))

        for i in range(executions):
            print(f"{i}\n")
            # Gillespie simulation
            new_parser = Parser(path)
            gillespie_sim = Gillespie.Gillespie(new_parser, T_MAX)
            gillespie_sim.gillespie_ssa()
            evolution = gillespie_sim.get_evolution()

            # Extract time values and species names
            time_values = evolution["time"]
            specie_names = [name for name in evolution if name != "time"]

            for specie in specie_names:
                values = evolution[specie]

                # Step-wise interpolation (values stay constant until next time point)
                interp = interp1d(
                    time_values,
                    values,
                    kind='previous',
                    bounds_error=False,
                    fill_value=(values[0], values[-1])
                )

                # Calculate interpolated values at each integer second
                interpolated_values = []
                for t in time_query:
                    v = float(interp(t))
                    interpolated_values.append(v)

                # Initialize with zeros if this species is encountered for the first time
                if specie not in sum_per_specie:
                    sum_per_specie[specie] = [0.0] * (T_MAX+ 1) #list

                # Accumulate values
                for idx in range(len(time_query)):
                    sum_per_specie[specie][idx] += interpolated_values[idx]

        # Calculate the average for each species and time point
        for specie in sum_per_specie:
            mean_dict = {}
            for idx, total in enumerate(sum_per_specie[specie]):
                media = round(total / executions, 4)
                tempo = time_query[idx]
                mean_dict[tempo] = media

            species_to_time_dict[specie] = mean_dict

        # Write to CSV file
        csv_path = f"./Example/Inference_of_kinetic_laws/{kinetic_constant_name}.csv"

        with open(csv_path, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow([constants.TIME] + [f"{species}" for species in species_to_time_dict])

            for time_point in range(0, T_MAX + 1):
                row = [time_point] + [species_to_time_dict[species][time_point] for species in species_to_time_dict]
                writer.writerow(row)

class Reaction:
    def __init__(self, reaction, parser):
        self.parser = parser
        self.reaction = reaction
        self.id = self.reaction.getId()
        self.reactants = {}
        self.products = {}
        self.rate_formula = None
        self.rate_formula_rev = None
        self.kinetic_law = None
        self.isReversible = self.reaction.getReversible()
        self.df_csv = None
        self.constant_inferred_name = None

        self.extract_reaction()

    '''This method retrieves the stoichiometric coefficients of reactants and 
    products, extracts the kinetic law, and validates its mass-action structure 
    based on whether the reaction is reversible or not. '''
    def extract_reaction(self):
        # Reactans' coefficients
        for sr in self.reaction.getListOfReactants():
            self.reactants[sr.getSpecies()] = sr.getStoichiometry()

        # Products' coefficients
        for sp in self.reaction.getListOfProducts():
            self.products[sp.getSpecies()] = sp.getStoichiometry()

        # Kinectic law
        self.kinetic_law = self.reaction.getKineticLaw()
        if self.kinetic_law is None:
            raise Exception("Kinetic law not found.")

        for compartment in self.parser.model.getListOfCompartments():
            if self.contains_identifier(compartment.getId()):
                raise Exception("Kinetic law invalid for us.")

        if not self.isReversible:
            self.validate_mass_action_kinetic_law(self.kinetic_law.getMath())

            self.rate_formula = self.kinetic_law.getFormula()
        else:
            ast = self.kinetic_law.getMath()
            if ast is None:
                raise Exception("Kinetic law not found.")
            # kinetic law = forward - reverse
            if ast.getType() != libsbml.AST_MINUS:
                raise Exception("Kinetic law with reversible attribute not correct.")

            ast_forward = ast.getChild(0)
            ast_reverse = ast.getChild(1)
            self.validate_mass_action_kinetic_law(ast_forward)
            self.validate_mass_action_kinetic_law(ast_reverse, isReverse = True)

            self.rate_formula = libsbml.formulaToString(ast_forward)
            self.rate_formula_rev = libsbml.formulaToString(ast_reverse)

    '''The method searches for a `.csv` file whose name (case-insensitive, 
    without extension) matches the provided constant. If such a file is found 
    and is a valid file, its full path is returned. If no match is found or 
    the directory is not set, it returns None.'''
    def check_file_path(self, constant):
        if self.parser.path_inference_directory is None:
            raise Exception("Path inference directory not set.")

        for filename in os.listdir(self.parser.path_inference_directory):
            if filename.endswith('.csv')  and filename[:-4].lower() == constant.lower():
                file_path = os.path.join(self.parser.path_inference_directory, filename)
                if os.path.isfile(file_path):
                    return file_path
        return None

    '''This method infers the stochastic rate constant for a given constant
    from experimental or simulated data in a CSV file. At each time interval
    extracted from the csv file, it calculates:
      - rates from each species changes and averages them to reduce the
       measurement noise.
      - product of the discrete amounts of reactants raised to their 
      respective stoichiometric coefficients.
    It then fits these values to estimate the stochastic rate constant and
    updates the model parameters to use value inferred. '''
    def stochastic_rate_constant_inference(self, constant, plot = True):
        file_path = self.check_file_path(constant)
        if file_path is None:
            raise Exception("File csv not found.")

        self.df_csv = pd.read_csv(file_path)
        self.constant_inferred_name = constant

        stoich_reactants = {reactant.getSpecies(): reactant.getStoichiometry() for reactant in
                            self.reaction.getListOfReactants()}

        stoich_products = {product.getSpecies(): product.getStoichiometry() for product in
                           self.reaction.getListOfProducts()}

        # There must be at least one specie of the reaction in the file and the time.
        species_in_csv = self.df_csv.columns
        num_of_species_in_csv = 0
        if constants.TIME not in species_in_csv:
            raise Exception("Missing time column in csv.")
        for species in stoich_reactants.keys():
            if species in species_in_csv:
                num_of_species_in_csv += 1
        for species in stoich_products.keys():
            if species in species_in_csv and species not in stoich_reactants:
                num_of_species_in_csv += 1
        if num_of_species_in_csv == 0:
            raise Exception("No species found in csv.")

        # Compute time differences between rows
        delta_t = self.df_csv[constants.TIME].diff().iloc[1:]

        v_t_list = []

        for reac_name, coeff_reac in stoich_reactants.items():
            if reac_name not in self.df_csv.columns:
                continue
            delta_p = self.df_csv[reac_name].diff().iloc[1:]
            if reac_name in stoich_products.keys():
                # If a species is both reactant and product, use the difference of
                # its stoichiometric coefficients.
                coeff_prod = stoich_products[reac_name]
                v = delta_p / ((coeff_prod - coeff_reac) * delta_t)
            else:
                #Otherwise, use the negative of the stoichiometric coefficient.
                v = delta_p / ((-1) * coeff_reac * delta_t)
            v_t_list.append(v)
        for prod_name, coeff_prod in stoich_products.items():
            if prod_name not in self.df_csv.columns or prod_name in stoich_reactants.keys():
                continue
            delta_p = self.df_csv[prod_name].diff().iloc[1:]
            v = delta_p / (coeff_prod * delta_t)
            v_t_list.append(v)

        if not v_t_list:
            raise Exception("No species found in CSV.")

        # This series containing the reaction rate values (v) over time.
        v_avg = sum(v_t_list) / len(v_t_list)

        # Calcolo del termine di legge d'azione di massa: [A]^1 * [B]^2
        rate_expr = np.ones(len(self.df_csv))
        for specie, coeff_reac in stoich_reactants.items():
            rate_expr *= self.df_csv[specie] ** coeff_reac
        rate_expr = rate_expr.iloc[1:]  # delete the first

        # Fitting: v ≈ k * rate_expr
        k_opt, *_ = np.linalg.lstsq(rate_expr.values.reshape(-1, 1), v_avg.values, rcond=None)

        if k_opt[0] < 0:
            raise ValueError(
                f"Estimated stochastic rate constant is negative: {k_opt[0]}. Check experimental data or model assumptions.")

        self.parser.parameters.update({constant: k_opt[0]})

        if plot:
            graphgen.stochastic_rate_constant_plot(rate_expr.values, v_avg.values, k_opt[0], constant)

    '''This method validates that the kinetic law AST (Abstract Syntax Tree)
     follows the expected mass-action form for a reaction. It checks that 
     the kinetic law either consists only of a stochastic rate constant 
     (for synthesis reactions) or has the proper structure.'''
    def validate_mass_action_kinetic_law(self, ast_root, isReverse = False): #verso
        if ast_root is None:
            raise Exception("AST is None.")
        # If the reaction is a synthesis (  -> something) then the kinetic law is equal
        # to the stochastic rate constant
        if ast_root.getType() == libsbml.AST_NAME:
            return

        if self.validate_mass_action_structure(ast_root, isReverse = isReverse) == 0:
            raise Exception("Stochastic rate constant absent.")

    '''This method recursively traverses the AST nodes of a kinetic law 
    to verify it matches the expected mass-action structure. Ensures that
    exactly one stochastic rate constant is present, and attempts to infer
    missing rate constants if needed.'''
    def validate_mass_action_structure(self, ast, stochastic_rate_constant_found=0, isReverse = False):
        if ast is None:
            raise Exception("AST is None.")

        # Always * at the root
        if ast.getType() != libsbml.AST_TIMES:
            raise Exception("AST's root is not TIMES(*).")

        #  all child nodes (multiplicands)
        for child in [ast.getChild(i) for i in range(ast.getNumChildren())]:
            if child.getType() == libsbml.AST_TIMES:
                stochastic_rate_constant_found = self.validate_mass_action_structure(child, stochastic_rate_constant_found)
                # There must be at least one parameter (rate constant)
                if stochastic_rate_constant_found > 1:
                    raise Exception("Too many stochastic rate constant.")
            elif child.getType() == libsbml.AST_NAME:  # stochastic rate constant
                stochastic_rate_constant_found += 1
                # Check if the stochastic rate constant should be inferred
                if child.getName() not in self.parser.parameters:
                    self.stochastic_rate_constant_inference(child.getName())
            elif child.getType() == libsbml.AST_FUNCTION_POWER:
                base = child.getChild(0)
                exponent = child.getChild(1)
                if (
                        base.getType() != libsbml.AST_NAME or
                        exponent.getType() not in constants.TYPE_NUMBER
                ):
                    if (not isReverse and exponent.getValue() != self.reactants[base.getName()] and
                            base.getName() not in self.reactants):
                        raise Exception("AST_FUNCTION_POWER node is written in the wrong way.")
                    elif (isReverse and exponent.getValue() != self.products[base.getName()] and
                        base.getName() not in self.products):
                        raise Exception("AST_FUNCTION_POWER node is written in the wrong way.")
            else:
                raise Exception(f"AST node not expected (type: {child.getType()}, name: {child.getName()}, "
                                f"value: {child.getValue()}) in: {self.kinetic_law.getFormula()}")

        return stochastic_rate_constant_found

    '''This method return True if the kinetic law contains the target_id
     parameter, otherwise return False.'''
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

    '''This method returns a dictionary representation of the reaction (usefull 
    for simulations).'''
    def get_reaction_as_dict(self):
        if not self.isReversible:
            return [{
                constants.ID: self.id,
                constants.REACTANTS: self.reactants,
                constants.PRODUCTS: self.products,
                constants.RATE_FORMULA: self.rate_formula
            }]
        else:
            return [{
                constants.ID: self.id,
                constants.REACTANTS: self.reactants,
                constants.PRODUCTS: self.products,
                constants.RATE_FORMULA: self.rate_formula
            },
            {
                constants.ID: self.id+"Rev",
                constants.REACTANTS: self.products,
                constants.PRODUCTS: self.reactants,
                constants.RATE_FORMULA: self.rate_formula_rev
            } ]


class Event:
    def __init__(self, event, parser):
        self.parser = parser
        self.event = event
        self.id = event.getId()
        self.trigger_formula = {}
        self.previous = False
        self.list_of_event_assigment = []
        self.delay_formula = None
        self.use_values_from_trigger_time = None
        self.value_from_trigger_time = {}
        self.use_trigger_values = True

        self.extract_event()

    '''This method parses and validates an SBML event, including its trigger 
    condition, event assignments, and optional delay. Ensures model compliance
    with SBML Level 2 Version 3 and 4 rules, verifies units, checks for duplicate 
    or missing variables, and prepares the event for evaluation and scheduling. '''
    def extract_event(self):
        # SBML Level 2 Version 4: the useValuesFromTriggerTime attribute is only valid if a delay is also defined.
        # This block checks that constraint and sets a flag (use_trigger_values) to True only if the event has a delay
        # and set UseValuesFromTriggerTime to True or is not setted.
        # Level 2 Version 3: The identifiers occurring in the MathML ci attributes of the
        # 2 EventAssignment object represent the value of the identifier at the point when the Event is fired.
        if self.parser.model.getVersion() == 4:
            if not (self.event.isSetDelay() and self.event.getUseValuesFromTriggerTime()):
                self.use_trigger_values = False

        self.use_values_from_trigger_time = self.use_trigger_values

        trigger = self.event.getTrigger()
        if trigger is None:
            raise Exception("Trigger not found.")
        # Check if trigger condition is valid
        ast_trigger = trigger.getMath()
        self.validate_trigger_boolean_expr(ast_trigger)
        self.trigger_formula = (libsbml.formulaToL3String(ast_trigger).replace("&&", " and ")
                                .replace("||", " or ").replace("!", " not "))
        self.previous = self.evaluate_expr(self.trigger_formula, constants.ERROR_TRIGGER, 0)

        variables = []
        for event_assignment in self.event.getListOfEventAssignments():
            string_of_variable = event_assignment.getVariable()
            if string_of_variable in variables:  # In the same event each variable MUST be different
                raise Exception(f"Variable '{string_of_variable}' is duplicated in the event '{self.id}'.")
            variables.append(string_of_variable)

            variable = self.parser.model.getElementBySId(string_of_variable)
            if variable is None:
                variable = self.parser.parameters[string_of_variable]
                if variable is None:
                    raise Exception(f"Variable '{string_of_variable}' not found in model.")
                else:
                    raise Exception(f"Variable '{string_of_variable}' not found in model but it was inferred.")

            # Ottieni l'attributo 'units' dalla variabile (solo se esiste)
            unit_id = variable.getUnits() if variable.isSetUnits() else None
            if unit_id is None:
                raise Exception("Undefined unit definition for variable.")
            unit_definition_variable = self.parser.model.getUnitDefinition(unit_id)

            # Check the eventAssignment
            trigger_value_dict = self.validate_event_assigment(event_assignment.getMath(), unit_definition_variable)
            if self.use_trigger_values:
                self.value_from_trigger_time.update(trigger_value_dict)
            self.list_of_event_assigment.append(event_assignment)

        # TODO: add priority (supported only at level 3)
        # if version is 3 or there isn't the Delay tag, then we append the event and terminate by evaluating the single event
        if not (self.parser.model.getVersion() == 3 or not self.event.isSetDelay()):
            delay = self.event.getDelay()
            ast_delay = delay.getMath()

            bool_constant_found, bool_parameter_found = self.validate_delay(ast_delay)
            if bool_constant_found and not bool_parameter_found:
                raise Exception("Incorrect use of a constant for the delay.")

            self.delay_formula = libsbml.formulaToString(ast_delay)

    '''This method return True if the identifier name is valid in the model,
    otherwise return False.'''
    def is_valid_identifier(self, name):
        return (
                self.parser.model.getSpecies(name) is not None or
                self.parser.model.getParameter(name) is not None or
                self.parser.model.getCompartment(name) is not None or
                self.parser.model.getSpeciesReference(name) is not None
        )

    '''This method recursively validates the AST of a trigger expression to 
    ensure it consists of valid logical or relational operations, identifiers, 
    constants, or time references.'''
    def validate_trigger_boolean_expr(self, ast):
        if ast is None:
            raise Exception("AST is None.")

        if ast.getType() not in constants.TYPE_TRIGGER:
            raise Exception("AST must be logical or relational operator.")

        for child in [ast.getChild(i) for i in range(ast.getNumChildren())]:
            if child.getType() in constants.TYPE_TRIGGER:
                self.validate_trigger_boolean_expr(child)
            elif child.getType() == libsbml.AST_NAME:
                if not self.is_valid_identifier(child.getName()):
                    raise Exception("Invalid identifier.")
                return
            elif child.getType() in constants.TYPE_NUMBER or child.getType() == libsbml.AST_NAME_TIME:
                return
            else:
                raise Exception(f"AST node not expected (type: {child.getType()}, name: {child.getName()}, "
                                f"value: {child.getValue()}) in: {ast.getFormula()}")

    '''This method validates the AST of an EventAssigment expression to ensure
    it contains only supported operation (currently PLUS and MINUS) and only 
    valid identifier. It also checks the compatibility between the variable's 
    and EventAssigment expression's unit measure.'''
    def validate_event_assigment(self, ast, unit_definition_variable):
        event_assignment_input_vars = {} # store the set of variables used in the eventAssignment's Math expression
        if ast is None:
            raise Exception("AST is None.")
        if ast.getType() == libsbml.AST_NAME:
            string_of_element = ast.getName()
            if self.use_trigger_values:
                event_assignment_input_vars[string_of_element] = None
            element = self.parser.model.getElementBySId(string_of_element)
            if element is None:
                raise Exception(f"Symbol '{string_of_element}' not found in model.")

            if element.getTypeCode() not in constants.TYPE_CODE: # element != (Species, Parameter, Compartment))
                raise Exception(f"Unsupported element type for symbol '{string_of_element}'")

            if not element.isSetUnits():
                raise Exception(f"Unsupported element type for symbol '{string_of_element}'")

            # Check if the units are different
            if not libsbml.UnitDefinition.areEquivalent(self.parser.model.getUnitDefinition(element.getUnits()),
                                                        unit_definition_variable):
                raise Exception(f"Unit mismatch: '{string_of_element}' has units different from assigned variable.")
        elif ast.getType() not in constants.TYPE_NUMBER:
            if ast.getType() not in constants.TYPE_OP: # TODO: Extend to other operations (TIMES, DIV) and math functions
                raise Exception("Only AST_PLUS and AST_MINUS are supported in the EventAssigment.")

            for i in range(ast.getNumChildren()):
                return_value = self.validate_event_assigment(ast.getChild(i),unit_definition_variable)
                if self.use_trigger_values:
                    event_assignment_input_vars.update(return_value)
        return event_assignment_input_vars

    '''
    This method validates <delay> expressions to ensure they:
    1. The expression inside the <math> block evaluate to a duration in seconds, 
        which must match the model's time units.
    2. The expression must consist only of a single parameter or mathematical 
        operations (only SUM and MINUS) involving parameters and/or numerical 
        constants—no complex expressions or unsupported constructs are allowed.
    '''
    def validate_delay(self, ast, constant_found = False, parameter_found = False):
        if ast is None:
            raise Exception("AST is None.")
        if ast.getType() in constants.TYPE_NUMBER: # variable = constant
            return True, parameter_found
        elif ast.getType() == libsbml.AST_NAME:
            string_of_element = ast.getName()
            element = self.parser.model.getElementBySId(string_of_element)
            if element is None:
                raise Exception(f"Symbol '{string_of_element}' not found in model.")

            if element.getTypeCode() != libsbml.SBML_PARAMETER:  # element can only be a parameter
                raise Exception(f"Unsupported element type for symbol '{string_of_element}'")

            if not element.isSetUnits():
                raise Exception(f"Unsupported element type for symbol '{string_of_element}'")

            # Check if the units are seconds
            delay_units = element.getUnits()
            delay_ud = self.parser.model.getUnitDefinition(delay_units)

            if delay_units == 'second' and delay_ud is None:
                return constant_found, True
            elif delay_units != 'second' and delay_ud is not None:
                list_of_units = delay_ud.getListOfUnits()
                if len(list_of_units) != 1:
                    raise Exception(f"Too many units in {string_of_element}.")
                if list_of_units[0].getKind() != libsbml.UNIT_KIND_SECOND:
                    raise Exception(f"Unsupported unit type for symbol '{string_of_element}'. Only second supported.")
                return constant_found, True
            #TODO: extend to custom unitDefinitions for milliseconds and minutes
            raise Exception(f"Unsupported unit type for symbol '{string_of_element}'. Only second supported.")

        else:
            if ast.getType() not in constants.TYPE_OP:
                raise Exception("Only AST_PLUS and AST_MINUS are supported in the EventAssigment.")

            for i in range(ast.getNumChildren()):
                return_constant_found, return_parameter_found = self.validate_delay(ast.getChild(i))
                constant_found = return_constant_found or constant_found
                parameter_found = return_parameter_found or parameter_found

            return constant_found, parameter_found

    ''' This method evaluates a mathematical expression using the current 
    species values, parameters, and simulation time. The evaluation is done 
    in a restricted and safe environment. '''
    def evaluate_expr(self, expr, error_message, time_value, safe_globals = constants.SAFE_GLOBALS_BASE):
        # Merge state and parameters in a single dictionary for expression evaluation
        local_scope = {**self.parser.species, **self.parser.parameters, "time": time_value}

        try:
            return eval(expr, safe_globals, local_scope)
        except:
            raise Exception(
                f"{error_message} — Evaluation failed.\n"
                f"Expression: {expr}\n")

    '''This method returns a dictionary representation of the event (usefull 
    for simulations).'''
    def get_event_as_dict(self):
        return {
            constants.ID: self.id,
            constants.TRIGGER_FORMULA: self.trigger_formula,
            constants.PREVIOUS: self.previous,    # Value of trigger at t-tau (previous value)
            constants.LIST_OF_EVENT_ASSIGMENT: self.list_of_event_assigment,
            constants.DELAY_FORMULA: self.delay_formula,
            constants.USE_VALUES_FROM_TRIGGER_TIME: self.use_values_from_trigger_time,
            constants.VALUES_FROM_TRIGGER_TIME: self.value_from_trigger_time
        }
