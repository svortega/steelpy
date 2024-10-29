# 
# Copyright (c) 2009 steelpy
#
#
# Python stdlib imports
from __future__ import annotations
from array import array
#from dataclasses import dataclass
from collections import namedtuple
#import math
#

# package imports
from steelpy.sections.utils.shape.main import ShapeBasic
from steelpy.sections.utils.shape.channel import ChannelBasic

#
#
points = namedtuple('Points', ['y', 'z'])
#
#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------
#
class Channel(ShapeBasic):
    __slots__ = ['_labels', '_number', '_title', 
                 '_d', '_tw', '_b', '_tb']
    
    def __init__(self):
        """ """
        super().__init__()
        self._d: array = array('f', [])
        self._tw: array = array('f', [])        
        self._b: array = array('f', [])
        self._tb: array = array('f', [])  
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
            self._tb.append(parameters.tb)
    #
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            index = self._labels.index(shape_name)
        except ValueError:
            raise Exception(f" section name {shape_name} not found")
        #
        return ChannelBasic(name=self._labels[index],
                            d=self._d[index], tw=self._tw[index],
                            b=self._b[index], tb=self._tb[index])
# 
#
