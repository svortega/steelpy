#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
#from array import array
from typing import NamedTuple
from collections.abc import Mapping
#from collections import defaultdict #, Counter
#import re

# package imports
#from steelpy.ufo.load.process.utils import (get_value_point,
#                                            check_list_number)
#
#import steelpy.utils.io_module.text as common
#from steelpy.utils.dataframe.main import DBframework
#
# ---------------------------------
#
class PointNode(NamedTuple):
    """
    [fx, fy, fz, mx, my, mz,
    name, comment, load_name, system, load_complex, load_type]
    """
    fx: float
    fy: float
    fz: float
    mx: float
    my: float
    mz: float
    name: int|str
    comment: str
    load_name: int|str
    system: int
    load_complex: int
    load_type:str = "force"

    #
    def __str__(self, units: str = "si") -> str:
        """ """
        output = (f"{str(self.name):12s} {10 * ' '} "
                  f"{self.fx: 1.3e} {self.fy: 1.3e} {self.fz: 1.3e} "
                  f"{self.coordinate_system.upper():6s} "
                  f"{self.load_complex}\n")
                  
        #
        if (load_title := self.comment) == "NULL":
            load_title = ""        
        step = self.load_type
        output += (f"{step:12s} {10 * ' '} "
                   f"{self.mx: 1.3e} {self.my: 1.3e} {self.mz: 1.3e} "
                   f"{str(load_title)}\n")
        return output

    #
    @property
    def coordinate_system(self) -> str:
        if self.system != 0:
            return "local"
        return "global"
    #


#
# ---------------------------------
#
class DispNode(NamedTuple):
    """
    """
    x: float
    y: float
    z: float
    rx: float
    ry: float
    rz: float
    name: int|str
    comment: str
    load_name: int|str
    system: int
    load_complex: int
    load_type:str = "Nodal Displacement"
    #
    def __str__(self, units:str="si") -> str:
        """ """
        output = (f"{str(self.name):12s} {10 * ' '} "
                  f"{self.x: 1.3e} {self.y: 1.3e} {self.z: 1.3e} "
                  f"{self.coordinate_system.upper():6s} "
                  f"{self.load_complex}\n")
                  
        #
        if (load_title := self.comment) == "NULL":
            load_title = ""        
        step = self.load_type
        output += (f"{step:12s} {10 * ' '} "
                   f"{self.rx: 1.3e} {self.ry: 1.3e} {self.rz: 1.3e} "
                   f"{str(load_title)}\n")
        return output
    #
    @property
    def coordinate_system(self):
        if self.system != 0:
            return "local"
        return "global" 
#
# ---------------------------------
#
class NodeLoadBasic(Mapping):
    __slots__ = ['_system_flag']

    def __init__(self) -> None:
        """
        """
        # 0-global/ 1-local
        self._system_flag: int = 0 # Global system
    #
    #
    def __iter__(self):
        """
        """
        items = list(set(dict.fromkeys(self._labels)))
        return iter(items)

    def __contains__(self, value) -> bool:
        return value in self._labels

    def __len__(self) -> int:
        items = list(set(dict.fromkeys(self._labels)))
        return iter(items)

    #
    def __str__(self) -> str:
        """ """
        output = ""
        nodes = list(dict.fromkeys(self._labels))
        #nodes = set(self._labels)
        for node in nodes:
            items = self.__getitem__(node)
            for item in items:
                output += item.__str__()
                # print('---')
        return output
    #
    def get_number(self, start:int=1):
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