# 
# Copyright (c) 2009-2021 fem2ufo
#
# Python stdlib imports
from __future__ import annotations
#from typing import NamedTuple
#import pickle

# package imports
#from steelpy.trave3D.preprocessor.assemble import Rmatrix
#from steelpy.trave3D.processor.operations import zeros, trns_3Dv, zeros_vector
#from steelpy.process.math.vector import Vector
#
# ---------------------
#
def process_load(basic_load:dict, jbc:list, 
                 nodes, elements):
    """
    """
    print("** Processing Load")
    #
    node_load = nodal_load(nodes, basic_load)
    #
    member_load = beam_load(elements, nodes, 
                            basic_load, node_load)
    #
    loads = [node_load[i][j] 
             for i in range(len(jbc)) for j in range(6) 
             if jbc[i][j] != 0]    
    return loads, member_load

