# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
#from datetime import datetime as dt
#from typing import NamedTuple
#import time
#import os
#

# package imports
from steelpy.ufo.mesh.sqlite.nodes import NodeSQL
from steelpy.ufo.mesh.sqlite.elements import ElementsSQL
from steelpy.ufo.mesh.sqlite.boundary import BoundarySQL
from steelpy.ufo.utils.main import ufoBasicModel
#
#
#
#
class MeshSQL(ufoBasicModel):
    """
    mesh[beam_name] = [number, element1, element2, element]
    """
    __slots__ = ['db_file', '_build', '_name', 'data_type',
                 '_nodes', '_boundaries', '_elements']
    
    def __init__(self, component:str|int,
                 name:str|None = None,
                 sql_file:str|None = None):
        """
        """
        super().__init__(component)
        # TODO: sql type
        self._name = name
        self.db_file = sql_file
        self.data_type = "sqlite"
        #
        # --------------------------------------------------
        #
        self._nodes = NodeSQL(db_system=self.data_type,
                              component=self._component, 
                              db_file=self.db_file)
        #
        self._boundaries = BoundarySQL(db_system=self.data_type,
                                       component=self._component,
                                       db_file=self.db_file)
        #
        self._elements = ElementsSQL(db_system=self.data_type,
                                     component=self._component,
                                     db_file=self.db_file)
        #         
#
#