# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from collections.abc import Mapping
from dataclasses import dataclass

# package imports
from steelpy.utils.sqlite.main import ClassBasicSQL



#
@dataclass
class WindBasic(ClassBasicSQL):
    """ """
    __slots__ = ['db_file', '_component']
    #
    #
    def __init__(self, component: str|int, db_file: str):
        """
        """
        super().__init__(db_file)
        self._component = component
    
#