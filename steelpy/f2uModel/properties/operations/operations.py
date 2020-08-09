# 
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#from array import array
#import logging
from typing import  List, Union, Tuple # Iterator, Dict, Union, NamedTuple, 
#from collections.abc import Mapping

# package imports
#from steelpy.process.units.units import Units


class BasicProperty:
    __slots__ = ['name', '_groups', '_elements']
    
    def __init__(self):
        """
        """
        #self._units = Units()
        self._groups:List = []
        self._elements:List = []
        #self._default:bool = False
    
    #
    @property
    def group(self):
        """
        """
        return self._groups
    
    @group.setter
    def group(self, value:Union[List, Tuple, str, int]):
        """
        """
        if isinstance(value, (list, tuple)):
            self._groups.extend(list(value))
        else:
            self._groups.append(value)
    #
    #
    @property
    def elements(self):
        """
        """
        return self._elements
    
    @elements.setter
    def elements(self, value:Union[List, Tuple, str, int]):
        """
        """
        if isinstance(value, (list, tuple)):
            self._elements.extend(list(value))
        else:
            self._elements.append(value)
    #
    #def set_default(self):
    #    """
    #    """
    #    self._default = True