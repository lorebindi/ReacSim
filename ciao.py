from Constants import *
import libsbml

class Parser:
    def __init__(self, file_name):
        self.file_name = file_name
        self.model = self.read_sbml_file()
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
            reactions.append(Reaction(r).get_reaction_as_dict())
        if len(reactions) == 0: raise Exception("No reactions found.")
        return reactions

    def extract_events(self):
        events = []
        for e in self.model.getListOfEvents():
            events.append(Event(self.model, e).get_event_as_dict())
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

        def validate_mass_action_kinetic_law(self, reactants):

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
    def __init__(self, model):


