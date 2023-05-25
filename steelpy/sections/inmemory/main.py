#
# Copyright (c) 2019-2023 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from collections.abc import Mapping
import re
#

# package imports
from .angle import Angle
from .tubular import TubularIM
from .tee import Tee
from .channel import Channel
from .box import Box
from .ibeam import Ibeam
from .solid import SolidSection
from ..process.operations import get_sect_properties
from .operations import SectionBasic

# ---------------------------------
#
#
#
class SectionIM(SectionBasic):
    __slots__ = ['_labels', '_number', '_title', '_type',
                 '_tubular', '_solid', '_ibeam', '_box',
                 '_channel', '_tee', '_angle', '_default']

    def __init__(self):
        """
        """
        super().__init__()
        #
        self._tubular = TubularIM()
        self._solid = SolidSection()
        self._ibeam = Ibeam()
        self._box = Box()
        self._channel = Channel()
        self._tee = Tee()
        self._angle = Angle()
    #
    def __setitem__(self, shape_name: str | int,
                    properties: list[float] | dict[str, float] | str) -> None:
        """
        """
        try:
            self._labels.index(shape_name)
            raise Exception(f'Section {shape_name} already exist')
        except ValueError:
            #
            shape_type = properties[0]
            properties = get_sect_properties(properties[1:])
            #
            mnumber = next(self.get_number())
            self._labels.append(shape_name)
            self._title.append('NULL')
            self._number.append(mnumber)
            self._type.append(shape_type)
            #
            if re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe)\b",
                        shape_type, re.IGNORECASE):
                self._ibeam[shape_name] = properties

            elif re.match(r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
                self._tee[shape_name] = properties

            elif re.match(r"\b(tub(ular)?|pipe|chs)\b", shape_type, re.IGNORECASE):
                self._tubular[shape_name] = properties

            elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod|circular|round)\b",
                          shape_type, re.IGNORECASE):
                self._solid[shape_name] = [shape_type, *properties]

            elif re.match(r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE):
                self._box[shape_name] = properties
            
            elif re.match(r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE):
                self._channel[shape_name] = properties

            elif re.match(r"\b(l|angle)\b", shape_type, re.IGNORECASE):
                self._angle[shape_name] = properties

            else:
                raise Exception(" section item {:} not recognized".format(shape_type))
            #
            self._default = shape_name
    #
    def __getitem__(self, shape_name: int):
        """
        node_name : node number
        """
        try:
            index = self._labels.index(shape_name)
            shape_type = self._type[index]
        except ValueError:
            raise KeyError(f'   *** Section {shape_name} does not exist')
        #
        if re.match(r"\b(tub(ular)?|pipe)\b", shape_type, re.IGNORECASE):
            return self._tubular[shape_name]

        elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod|circular|round)\b", shape_type, re.IGNORECASE):
            return self._solid[shape_name]
        
        elif re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe)\b", shape_type, re.IGNORECASE):
            return self._ibeam[shape_name]
        
        elif re.match(r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE):
            return self._box[shape_name]
        
        elif re.match(r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE):
            return self._channel[shape_name]
        
        elif re.match(r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
            return self._tee[shape_name]
        
        elif re.match(r"\b(l|angle)\b", shape_type, re.IGNORECASE):
            return self._angle[shape_name]
        
        else:
            raise IOError(f' Section type {shape_type} not recognised')
    #
    #
    # def get_properties(self):
    #    """
    #    """
    #    summary = {}
    #    for key, item in self._sections.items():
    #        summary[key] = item._get_properties()
    #    return summary
    #
    #
    #@property
    #def tubular(self, shape_name:int|str):
    #    """ """
    #    return self._sections[shape_name]
    #
    #@tubular.setter
    #def tubualar(self, shape_name:int|str, properties):
    #    """ """
    #    self._sections[shape_name] = TubularBasic(name=shape_name,
    #                                              d=properties[0], t=properties[1])
    #
    #
    #
    #def df(self):
    #    """ """
    #    pass
#
