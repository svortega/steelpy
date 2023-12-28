# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping

# package imports
#from steelpy.process.units.units import Units
from steelpy.utils.sqlite.utils import create_connection, create_table



#
class WaveBasic(Mapping):
    
    #__slots__ = ['db_file']
    #
    #
    #
    #
    def __init__(self, db_file: str):
        """
        """
        self.db_file = db_file
        #
        # create table
        conn = create_connection(self.db_file)
        with conn:
            self._create_table(conn)         
    
    def __len__(self) -> int:
        return len(self._label)

    def __iter__(self):
        """
        """
        return iter(self._label)
    
    def __contains__(self, value) -> bool:
        return value in self._label
    
#