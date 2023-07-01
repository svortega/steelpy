# 
# Copyright (c) 2009-2023 steelpy
#
from __future__ import annotations
# Python stdlib imports
from array import array
from dataclasses import dataclass
import math
from typing import NamedTuple

# package imports
from steelpy.process.units.main import Units

class MeanCurrent:
     __slots__ = ['c_type', '_current']
    
     def __init__(self):
          """ """
          self.c_type = 1
          self._current = 0
     #
     @property
     def Euler(self):
          """ Eularian mean """
          self.c_type = 1
          
     @Euler.setter
     def Euler(self, value:Units):
          """ Eularian mean """
          self._current = value.convert('metre/second').value
          self.c_type = 1     
     #
     @property
     def Stokes(self):
          """ Stokes depth integrated mean """
          self.c_type = 2
     
     @Stokes.setter
     def Stokes(self, value:Units):
          """ Stokes depth integrated mean """
          self._current = value.convert('metre/second').value
          self.c_type = 2