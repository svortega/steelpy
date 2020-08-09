# 
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#from dataclasses import dataclass
#from typing import NamedTuple, Tuple, List, Iterator, Dict, Iterable


# package imports
from steelpy.process.units.main import Units
from steelpy.f2uModel.properties.codecheck.codecheck import CodeCheck
from steelpy.f2uModel.properties.hydrodynamic.hydrodynamic import Hydrodynamic


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
    @property
    def hydrodynamic(self):
        """ """
        return self._hydrodynamic
    #
    @property
    def design_parameters(self):
        """ """
        return self._code_check
#
    