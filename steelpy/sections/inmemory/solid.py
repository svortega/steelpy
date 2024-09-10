# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
from array import array
from collections import namedtuple
from dataclasses import dataclass
#import math
import re
#from typing import NamedTuple

# package imports
#from steelpy.sections.process.operations import ShapeProperty
#from steelpy.sections.process.stress import BeamStress
from steelpy.sections.utils.shape.main import ShapeBasic
from steelpy.sections.utils.shape.solid import (CircleSolid,
                                                RectangleSolid,
                                                TrapeziodSolid)

# ----------------------------------------
#      Basic Solid Shapes
# ----------------------------------------
#
points = namedtuple('Points', ['y', 'z'])
#
#
@dataclass
class SolidSection(ShapeBasic):
    """ """
    #
    def __init__(self):
        """ """
        super().__init__()
        #
        self._d: array = array('f', [])
        self._wb: array = array('f', [])
        self._wt: array = array('f', [])
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
            shape_type = parameters.pop(0)
            self._labels.append(shape_name)
            self._title.append('NULL')
            mnumber = next(self.get_number())
            self._number.append(mnumber)
            #
            self._d.append(parameters[0])
            if re.match(r"\b((solid|bar(\_)?)?circular|round)\b", shape_type, re.IGNORECASE):
                self._type.append('round')
                self._wb.append(0)
                self._wt.append(0)

            elif re.match(r"\b((solid|bar(\_)?)?square|rectangle)\b", shape_type, re.IGNORECASE):
                self._type.append('rectangle')
                self._wb.append(parameters[1])
                self._wt.append(parameters[1])

            elif re.match(r"\b((solid|bar(\_)?)?trapeziod)\b", shape_type, re.IGNORECASE):
                self._type.append('trapeziod')
                self._wb.append(parameters[1])
                try:
                    self._wt.append(parameters[2])
                except IndexError:
                    raise IOError('Trapeziod section data missing')

            else:
                raise Exception(f" section type {shape_type} not recognized")
    #
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            index = self._labels.index(shape_name)
        except ValueError:
            raise Exception(f" section name {shape_name} not found")


        shape_type = self._type[index]

        if re.match(r"\b((solid|bar(\_)?)?circular|round)\b", shape_type, re.IGNORECASE):
            return CircleSolid(d=self._d[index])

        elif re.match(r"\b((solid|bar(\_)?)?square|rectangle)\b", shape_type, re.IGNORECASE):
            return RectangleSolid(depth=self._d[index], width=self._wb[index])

        elif re.match(r"\b((solid|bar(\_)?)?trapeziod)\b", shape_type, re.IGNORECASE):
            c = abs(self._wt[index] - self._wb[index]) / 2.0
            return TrapeziodSolid(depth=self._d[index],
                                  width=self._wb[index],
                                  a=self._wt[index], c=c)

        else:
            raise Exception(f" section type {shape_type} not recognized")
    #
    #
    #def set_default(self):
    #    """ """
    #    self.cls._default = self.name
    #
    #@property
    #def height(self):
    #    """
    #    d : height of rectangular bar
    #    """
    #    return self.depth
    ##
    #@property
    #def d(self):
    #    """
    #    d : Section height of rectangular bar
    #    """
    #    return self.depth
    #@d.setter
    #def d(self, value):
    #    """
    #    d : Section height of rectangular bar
    #    """
    #    self.depth = value
    ##
    #@property
    #def w(self):
    #    """
    #    w : width of rectangular bar
    #    """
    #    return self.width
    #@w.setter
    #def w(self, value):
    #    """
    #    w : width of rectangular bar
    #    """
    #    self.width = value
    ##
    #@property
    #def wb(self):
    #    """
    #    wb : width bottom of rectangular bar
    #    """
    #    return self.width
    #
    #@wb.setter
    #def wb(self, value):
    #    """
    #    wb : width bottom of rectangular bar
    #    """
    #    self.width = value
    ##
    #@property
    #def wt(self):
    #    """
    #    wt : width top of rectangular bar
    #    """
    #    return self.a
    #
    #@wt.setter
    #def wt(self, value):
    #    """
    #    wt : width top of rectangular bar
    #    """
    #    self.a = value
    #
#
#