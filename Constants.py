import libsbml
import math
import operator

TIME = '_time_'
# reaction
ID = "id"
REACTANTS = "reactants"
PRODUCTS = "products"
RATE_FORMULA = "rate_formula"

TYPE_NUMBER = (libsbml.AST_INTEGER, libsbml.AST_REAL,
               libsbml.AST_REAL_E, libsbml.AST_RATIONAL)

# event
TRIGGER_FORMULA = "trigger_formula"
PREVIOUS = "previous"
IS_TIME = "is_time_trigger"
LIST_OF_EVENT_ASSIGMENT_FORMULA = "list_of_event_assigment_formula"
DELAY_FORMULA = "delay_formula"
PRIORITY = "priority"

TYPE_TRIGGER_LOGICAL = (libsbml.AST_LOGICAL_AND, libsbml.AST_LOGICAL_OR,
              libsbml.AST_LOGICAL_NOT, libsbml.AST_LOGICAL_XOR)

TYPE_TRIGGER_RELATIONAL = (libsbml.AST_RELATIONAL_EQ, libsbml.AST_RELATIONAL_GEQ,
              libsbml.AST_RELATIONAL_GT, libsbml.AST_RELATIONAL_LEQ,
              libsbml.AST_RELATIONAL_LT,libsbml.AST_RELATIONAL_NEQ)

TYPE_TRIGGER = TYPE_TRIGGER_LOGICAL + TYPE_TRIGGER_RELATIONAL

TYPE_CODE = (libsbml.SBML_PARAMETER, libsbml.SBML_SPECIES, libsbml.SBML_COMPARTMENT)

TYPE_OP = (libsbml.AST_PLUS, libsbml.AST_MINUS)

ERROR_KINETIC_LAW = "Kinetic Law"
ERROR_TRIGGER = "Trigger"
ERROR_EVENT_ASSIGNMENTS = "Event assignment"
ERROR_DELAY = "Delay"
SAFE_GLOBALS_BASE = {"__builtins__": None, "math": math}
SAFE_GLOBALS_RATE = {**SAFE_GLOBALS_BASE,
                     "pow": pow}
'''
SAFE_GLOBALS_TRIGGER = {**SAFE_GLOBALS_BASE,
                        # Relational
                        "eq": operator.eq,
                        "geq": operator.ge,
                        "gt": operator.gt,
                        "leq": operator.le,
                        "lt": operator.lt,
                        "neq": operator.ne,
                        # Logical
                        "and": lambda a, b: a and b,
                        "or": lambda a, b: a or b,
                        "not": lambda a: not a,
                        "xor": operator.xor}
'''