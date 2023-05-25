# 
# Copyright (c) 2009-2023 fem2ufo
#

# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
#from typing import NamedTuple, Tuple, List, Iterator, Dict, Iterable


# package imports
#from steelpy.process.units.main import Units
from .codecheck.codecheck import CodeCheck
from .hydrodynamic.hydrodynamic import Hydrodynamic


#
#
class Properties:
    
    __slots__ = ['_hydrodynamic', '_code_check']
    
    def __init__(self):
        """
        """
        #self._units = Units()
        self._hydrodynamic =  Hydrodynamic()
        self._code_check = CodeCheck()
    #
    #@property
    #def units(self):
    #    """
    #    """
    #    return self._units    
    #
    #@property
    def hydrodynamic(self):
        """ """
        return self._hydrodynamic
    #
    #@property
    def design_parameters(self):
        """ """
        return self._code_check
#
    