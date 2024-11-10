#
# Copyright (c) 2019 steelpy
#

# Python stdlib imports
from __future__ import annotations
from array import array
#from dataclasses import dataclass
from collections.abc import Mapping
#import math
#import re
#

# package imports
from steelpy.sections.utils.operations import find_sect_type
#
from steelpy.sections.utils.shape.tubular import TubularBasic
from steelpy.sections.utils.shape.solid import RectangleSolid, CircleSolid, TrapezoidSolid
from steelpy.sections.utils.shape.ibeam import IbeamBasic
from steelpy.sections.utils.shape.box import BoxBasic
from steelpy.sections.utils.shape.channel import ChannelBasic
from steelpy.sections.utils.shape.tee import TeeBasic
from steelpy.sections.utils.shape.angle import AngleBasic
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
class SectionMain(Mapping):
    __slots__ = [#'_labels', '_number', '_title', '_type', 
                 '_tubular', '_solid', '_ibeam', '_box',
                 '_channel', '_tee', '_angle', '_default']

    def __init__(self):
        """
        """
        self._default: str | None = None
    #
    def __setitem__(self, shape_name: str | int,
                    properties: list[float] | dict[str, float] | str) -> None:
        """
        """
        try:
            self._labels.index(shape_name)
            raise Exception(f'Section {shape_name} already exist')
        except ValueError:
            section_type = properties.shape
            match section_type:
                case 'I section':
                    # [d, tw, bf, tf, bfb, tfb, r, title]
                    self._ibeam[shape_name] = properties
                case 'Tee':
                    self._tee[shape_name] = properties
                case 'Tubular':
                    self._tubular[shape_name] = properties
                case 'Box':
                    self._box[shape_name] = properties
                case 'Channel':
                    self._channel[shape_name] = properties
                case 'Angle':
                    self._angle[shape_name] = properties
                case 'Rectangle bar'|'Trapezoid bar'|'Circular bar':
                    self._solid[shape_name] = properties
                case 'General section':
                    raise NotImplementedError('general section error')
                case _:
                    raise Exception(f" section item {section_type} not recognized")
    #
    def __getitem__(self, shape_name: int|str):
        """
        section_name :
        """
        try:
            index = self._labels.index(shape_name)
            section_type = find_sect_type(self._type[index])
        except ValueError:
            raise KeyError(f'   *** Section {shape_name} does not exist')
        #
        match section_type:
            case 'I section':
                return self._ibeam[shape_name]
            case 'Tee':
                return self._tee[shape_name]
            case 'Tubular':
                return self._tubular[shape_name]
            case 'Box':
                return self._box[shape_name]
            case 'Channel':
                return self._channel[shape_name]
            case 'Angle':
                return self._angle[shape_name]
            case 'Rectangle bar' | 'Trapezoid bar' | 'Circular bar':
                return self._solid[shape_name]
            case 'General section':
                raise NotImplementedError('general section error')
            case _:
                raise IOError(f' Section type {section_type} not recognised')
    #
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
#-------------------------------------------------
#
#
def is_section(item)->bool:
    """ """
    if isinstance(item, TubularBasic):
        return True
    elif isinstance(item, IbeamBasic):
        return True
    elif isinstance(item, BoxBasic):
        return True
    elif isinstance(item, ChannelBasic):
        return True
    elif isinstance(item, TeeBasic):
        return True
    elif isinstance(item, AngleBasic):
        return True
    elif isinstance(item, RectangleSolid):
        return True
    elif isinstance(item, CircleSolid):
        return True
    elif isinstance(item, TrapezoidSolid):
        return True
    else:
        return False
#
#-------------------------------------------------
#
