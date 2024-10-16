#
# Copyright (c) 2019 steelpy
#

# Python stdlib imports
from __future__ import annotations
from array import array
#from dataclasses import dataclass
from collections.abc import Mapping
#import math
import re
#

# package imports
#
from steelpy.sections.utils.operations import get_sect_properties
from steelpy.sections.utils.shape.tubular import TubularBasic
from steelpy.sections.utils.shape.solid import RectangleSolid, CircleSolid, TrapezoidSolid
from steelpy.sections.utils.shape.ibeam import IbeamBasic, get_Isection
from steelpy.sections.utils.shape.box import BoxBasic
from steelpy.sections.utils.shape.channel import ChannelBasic
from steelpy.sections.utils.shape.tee import TeeBasic
from steelpy.sections.utils.shape.angle import AngleBasic
#
from steelpy.utils.dataframe.main import DBframework
#
#
#
#-------------------------------------------------
#
class ShapeBasic(Mapping):
    """ """
    __slots__ = ['_labels', '_title', '_number', '_type']
    
    def __init__(self):
        """
        """
        self._labels:list = []
        self._title:list = []
        self._type:list = []
        self._number:array= array('i', [])        
    # -----------------------------------------------
    #
    def __len__(self):
        return len(self._labels)

    def __iter__(self):
        return iter(self._labels)

    def __contains__(self, value):
        return value in self._labels
    #
    # -----------------------------------------------   
    #
    def get_number(self, start: int = 1):
        """
        """
        try:
            n = max(self._number) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
    #
    
#
#
#
#-------------------------------------------------
#
#
class SectionMain(ShapeBasic):
    #__slots__ = ['_labels', '_number', '_title', '_type',
    #             '_tubular', '_solid', '_ibeam', '_box',
    #             '_channel', '_tee', '_angle', '_default']

    def __init__(self):
        """
        """
        super().__init__()
    #
    def __setitem__(self, shape_name: str | int,
                    properties: list[float] | dict[str, float] | str) -> None:
        """
        """
        try:
            self._labels.index(shape_name)
            raise Exception(f'Section {shape_name} already exist')
        except ValueError:
            shape_type = properties.shape
            #
            if re.match(r"\b(i((_|-|\s*)?beam)?|w|m|s|hp|ub|uc|he|ipe|pg((_|-|\s*)?section)?)\b",
                        shape_type, re.IGNORECASE):
                # [d, tw, bf, tf, bfb, tfb, r, title]
                self._ibeam[shape_name] = properties

            elif re.match(r"\b(t(ee)?((_|\s*)?bar)?)\b", shape_type, re.IGNORECASE):
                self._tee[shape_name] = properties                  

            elif re.match(r"\b(tub(ular)?|pipe|chs((_|-|\s*)?section)?)\b", shape_type, re.IGNORECASE):
                self._tubular[shape_name] = properties

            elif re.match(r"\b((solid(_|\s*)?)?square|rectangle|trapezoid|circular|round((_|\s*)?bar)?)\b",
                          shape_type, re.IGNORECASE):
                self._solid[shape_name] = properties

            elif re.match(r"\b(b(ox)?|rhs|shs((_|-|\s*)?section)?)\b", shape_type, re.IGNORECASE):
                self._box[shape_name] = properties 

            elif re.match(r"\b(c(hannel)?((_|-|\s*)?section)?)\b", shape_type, re.IGNORECASE):
                self._channel[shape_name] = properties                 

            elif re.match(r"\b(l|angle((_|-|\s*)?section)?)\b", shape_type, re.IGNORECASE):
                self._angle[shape_name] = properties

            elif re.match(r"\b(general((_|-|\s*)?section)?)\b", shape_type, re.IGNORECASE):
                #self._angle[shape_name] = properties
                raise NotImplementedError('general section error')

            else:
                raise Exception(" section item {:} not recognized".format(shape_type))    
    #
    def __getitem__(self, shape_name: int|str):
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

        elif re.match(r"\b((solid|bar(_|-|\s*)?)?square|rectangle|trapezoid|circular|round)\b", shape_type, re.IGNORECASE):
            return self._solid[shape_name]
        
        elif re.match(r"\b(i((_|-|\s*)?beam|section)?|w|m|s|hp|ub|uc|he|ipe|pg)\b", shape_type, re.IGNORECASE):
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

#
#
#-------------------------------------------------
#
def ShapeGeometry(shape_type: str, geometry: list):
        """ """
        if re.match(r"\b(tub(ular)?|pipe)\b", shape_type, re.IGNORECASE):
            return TubularBasic(name=geometry[0],
                                diameter=geometry[3], thickness=geometry[4]) 
    
        #elif re.match(r"\b((solid|bar(\_)?)?rectangle|trapeziod|circular|round)\b", shape_type, re.IGNORECASE):
        #    return self._solid[shape_name]
        elif re.match(r"\b((solid|bar(\_)?)?circular|round)\b", shape_type, re.IGNORECASE):
            d = geometry[5]
            return CircleSolid(name=geometry[0], d=d)
    
        elif re.match(r"\b((solid|bar(\_)?)?rectangle)\b", shape_type, re.IGNORECASE):
            d = geometry[5]
            wb = geometry[7]
            return RectangleSolid(name=geometry[0], depth=d, width=wb)
    
        elif re.match(r"\b((solid|bar(\_)?)?trapezoid)\b", shape_type, re.IGNORECASE):
            d = geometry[5]
            wb = geometry[7]
            wt = geometry[9]            
            c = abs(wt - wb) / 2.0
            return TrapezoidSolid(name=geometry[0], depth=d, width=wb, a=wt, c=c)
        
        elif re.match(r"\b(i((\_)?beam|section)?|w|m|s|hp|ub|uc|he|ipe|pg)\b", shape_type, re.IGNORECASE):
            return IbeamBasic(name=geometry[0],
                              d=geometry[5], tw=geometry[6],
                              bft=geometry[7], tft=geometry[8],
                              bfb=geometry[9], tfb=geometry[10])
        
        elif re.match(r"\b(b(ox)?|rhs|shs)\b", shape_type, re.IGNORECASE):
            return BoxBasic(name=geometry[0],
                            d=geometry[5], tw=geometry[6],
                            b=geometry[7], tb=geometry[8])
        
        elif re.match(r"\b(c(hannel)?)\b", shape_type, re.IGNORECASE):
            return ChannelBasic(name=geometry[0],
                                d=geometry[5], tw=geometry[6],
                                b=geometry[7], tb=geometry[8])
        
        elif re.match(r"\b(t(ee)?)\b", shape_type, re.IGNORECASE):
            return TeeBasic(name=geometry[0],
                            d=geometry[5], tw=geometry[6],
                            b=geometry[7], tb=geometry[8])
        
        elif re.match(r"\b(l|angle)\b", shape_type, re.IGNORECASE):
            return AngleBasic(name=geometry[0],
                              d=geometry[5], tw=geometry[6],
                              b=geometry[7], r=0)
        
        else:
            raise IOError(f' Section type {shape_type} not recognised')
#
