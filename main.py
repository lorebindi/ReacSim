# This is a sample Python script.
import Parser
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    model=Parser.readSBMLfile('/home/lorenzo/Desktop/ReacSim/Example/1.xml')
    Parser.gillespie_ssa(model,10)


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
