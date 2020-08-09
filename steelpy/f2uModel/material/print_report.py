# 
# Copyright (c) 2018-2020 steelpy


# Python stdlib imports
#import datetime
#from typing import List
#

# package imports

def print_plastic_material(self):
    """
    """
    output =[]
    output.append("_______________________________________________________________________________________\n")
    output.append(" "+"\n")
    output.append("                                  MATERIAL PROPERTIES\n")
    output.append(" "+"\n")
    output.append("Name           Fy [N/mm2]  Fu [N/mm2]  E  [N/mm2]  Poisson     Rho[kg/m3]    G  [N/mm2]\n")
    output.append(".......................................................................................\n")
    output.append("\n")
    output.append("{:14s} {:1.4E}  {:1.4E}  {:1.4E}  {:1.4E}  {:1.4E}   {:1.4E}\n"
                  .format("", self.Fy.convert('megapascal').value, 
                          self.Fu.convert('megapascal').value, 
                          self.E.convert('megapascal').value, 
                          self.poisson, 
                          self.rho.convert('kilogram/metre^3').value,
                          self.G.convert('megapascal').value))
    return output