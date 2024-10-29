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
#from steelpy.sections.utils.operations import ShapeProperty
#from steelpy.sections.utils.stress import BeamStress
from steelpy.sections.utils.shape.main import ShapeBasic
from steelpy.sections.utils.shape.solid import (CircleSolid,
                                                RectangleSolid,
                                                TrapezoidSolid)

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
            shape_type = parameters.shape
            self._labels.append(shape_name)
            self._title.append('NULL')
            mnumber = next(self.get_number())
            self._number.append(mnumber)
            #
            if re.match(r"\b((solid(_|-|\s*)?)?circular|round((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
                self._type.append("Circular Bar")
                self._d.append(parameters.d)
                self._wb.append(0)
                self._wt.append(0)

            elif re.match(r"\b((solid(_|-|\s*)?)?square|rectangle((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
                self._type.append("Rectangle Bar")
                self._d.append(parameters.h)
                self._wb.append(parameters.b)
                self._wt.append(parameters.b)

            elif re.match(r"\b((solid(_|-|\s*)?)?trapezoid((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
                self._type.append("Trapezoid Bar")
                self._d.append(parameters.h)
                self._wb.append(parameters.b)
                #try:
                self._wt.append(parameters.a)
                #except IndexError:
                #    raise IOError('Trapeziod section data missing')

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

        if re.match(r"\b((solid(_|-|\s*)?)?circular|round((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
            return CircleSolid(name=shape_name, d=self._d[index])

        elif re.match(r"\b((solid(_|-|\s*)?)?square|rectangle((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
            return RectangleSolid(name=shape_name,
                                  depth=self._d[index],
                                  width=self._wb[index])

        elif re.match(r"\b((solid(_|-|\s*)?)?trapezoid((_|-|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
            c = abs(self._wt[index] - self._wb[index]) / 2.0
            return TrapezoidSolid(name=shape_name,
                                  depth=self._d[index],
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