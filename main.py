# This is a sample Python script.
from Parser import *
from Simulation_methods import *
from Graph_generation import *
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    current_dir = os.path.dirname(__file__)
    file_path = os.path.join(current_dir, "Example", "2.xml")
    model = Parser.readSBMLfile(file_path)
    evolution = Parser.gillespie_ssa(model, 1)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
