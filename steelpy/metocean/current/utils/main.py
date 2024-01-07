# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass

# package imports
#from steelpy.utils.sqlite.main import ClassBasicSQL
from steelpy.metocean.hydrodynamic.utils.main import HydroBasic
from steelpy.utils.sqlite.utils import create_connection #, create_table



#
@dataclass
class CurrentBasic(HydroBasic):
    
    #__slots__ = ['db_file']
    #
    #
    def __init__(self, db_file: str):
        """
        """
        super().__init__(db_file)
    #
    # 
    def __setitem__(self, name:str|int,
                    value: list|tuple|dict|str) -> None:
        """
        """
        # [name, type, title, criteria]
        data = (name, 'profile', None, self._criteria)
        conn = create_connection(self.db_file)
        with conn:
            self._push_data(conn, data)
        #
        item = self.__getitem__(name)
        #
        #
        if isinstance(value, list|tuple):
            if isinstance(value[0], list|tuple):
                item.profile = value
            else:
                profile = self.get_profile(value)
                item.profile = profile
        
        elif isinstance(value, dict):
            1 / 0
            
        elif isinstance(value, str):
            1 / 0
            data = self.get_data([value])
        else:
            raise IOError('Input data not valid')
    #
#