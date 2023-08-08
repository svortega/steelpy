#
# Copyright (c) 2009-2023 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
from dataclasses import dataclass
#from typing import NamedTuple
#import pickle
#from typing import NamedTuple

#
# package imports
# steelpy.trave3D.postprocessor
#from .operations import get_reactions
from .output import  Node, Beam
from steelpy.utils.dataframe.main import DBframework

#
#
#
# --------------------
# Results
# --------------------
#
#
#
@dataclass
class Results:
    __slots__ = ['node', 'beam', 'shell']
    
    def __init__(self, noderes, nodereac, beamres, plane) -> None:
        """
        """
        self.node:tuple = Node(results=noderes, 
                               reac=nodereac,
                               plane=plane)
        #
        self.beam:tuple = Beam(results=beamres, plane=plane)
    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += self.node.__str__()
        output += self.beam.__str__()
        return output
    #
    #
    def reactions(self):
        """ """
        pass
#