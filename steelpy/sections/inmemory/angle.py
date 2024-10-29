# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
from array import array
#from collections import namedtuple
#from dataclasses import dataclass
import math
#

# package imports
from steelpy.sections.utils.shape.angle import AngleBasic
from steelpy.sections.utils.shape.main import ShapeBasic

#
#
#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------
#
class Angle(ShapeBasic):
    __slots__ = ['_labels', '_number', '_title', 
                 '_d', '_tw', '_b', '_r']
    
    def __init__(self):
        """ """
        super().__init__()
        self._d: array = array('f', [])
        self._tw: array = array('f', [])        
        self._b: array = array('f', [])
        self._r: array = array('f', [])  
    #
    #
    def __setitem__(self, shape_name: int|str, parameters: list) -> None:
        """
        parameters = [node1, node2, material, section, roll_angle]
        """
        try:
            self._labels.index(shape_name)
            raise Exception('element {:} already exist'.format(shape_name))
        except ValueError:
            self._labels.append(shape_name)
            self._title.append('NULL')
            mnumber = next(self.get_number())
            self._number.append(mnumber)
            #
            self._d.append(parameters.d)
            self._tw.append(parameters.tw)
            self._b.append(parameters.b)
            try:
                self._r.append(parameters.r*0.50)
            except TypeError:
                self._r.append(0.50 * (self._tw[-1] / math.sqrt(2.0)))
    #
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            index = self._labels.index(shape_name)
        except ValueError:
            raise Exception(f" section name {shape_name} not found")
        #
        return AngleBasic(name=self._labels[index],
                          d=self._d[index], tw=self._tw[index],
                          b=self._b[index], r=self._r[index])
#
