# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
from array import array
#from dataclasses import dataclass
#from collections import namedtuple
#import math
#import re
#

# package imports
from steelpy.sections.utils.shape.ibeam import IbeamBasic
from steelpy.sections.utils.shape.main import ShapeBasic

#
#
#
#-------------------------------------------------
#
#-------------------------------------------------
#
class Ibeam(ShapeBasic):
    __slots__ = ['_labels', '_number', '_title', '_type', '_default', 
                 '_d', '_tw', '_bft', '_tft', '_bfb', '_tfb',
                 '_r']
    def __init__(self):
        """ """
        super().__init__()
        self._d: array = array('f', [])
        self._tw: array = array('f', [])
        self._bft: array = array('f', [])
        self._tft: array = array('f', [])
        self._bfb: array = array('f', [])
        self._tfb: array = array('f', [])
        self._r: array = array('f', [])   # fillet radius
    #
    #
    def __setitem__(self, shape_name: int|str, parameters: list) -> None:
        """
        Parameters
        ----------
        d   : Section Height
        tw  : Web thickness
        bf  : flange base (top)
        tf  : flange thickness (top)
        bfb : Bottom flange base
        tfb : Bottom flange thickness
        r   : root radius
        title:
        """
        try:
            self._labels.index(shape_name)
            raise Exception('Section {:} already exist'.format(shape_name))
        except ValueError:
            self._labels.append(shape_name)
            mnumber = next(self.get_number())
            self._number.append(mnumber)
            #
            self._d.append(parameters.d)
            self._tw.append(parameters.tw)
            self._bft.append(parameters.bf)
            self._tft.append(parameters.tf)
            self._bfb.append(parameters.bfb)
            self._tfb.append(parameters.tfb)
            self._r.append(parameters.r)
            self._title.append(parameters.title)
            #            
    #
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            index = self._labels.index(shape_name)
        except ValueError:
            raise Exception(f" section name {shape_name} not found")
        #
        return IbeamBasic(name=self._labels[index],
                          d=self._d[index], tw=self._tw[index],
                          bft=self._bft[index], tft=self._tft[index],
                          bfb=self._bfb[index], tfb=self._tfb[index], 
                          root_radius=self._r[index])
#
#

