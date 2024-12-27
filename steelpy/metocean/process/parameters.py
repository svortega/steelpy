# 
# Copyright (c) 2019 steelpy
# 
# Python stdlib imports
from __future__ import annotations

# package imports
from steelpy.utils.sqlite.main import ClassBasicSQL
from steelpy.utils.sqlite.utils import create_table

#
#
class HydroDesign(ClassBasicSQL):
    """ Hydrodynamic Design """
    __slots__ = ['_name', '_db_file', '_component']
    
    def __init__(self, component: str|int, db_file:str|None = None):
        """
        """
        super().__init__(db_file)
        self._component =  component
    #
    #
    def __setitem__(self, name: int|str,
                    parameters: list) -> None:
        """
        """
        pass
    
    def __getitem__(self, name: int|str):
        """ """
        pass    
    # ----------------------------------------
    # SQL ops
    #
    def _new_table(self, conn) -> None:
        """ """
        #
        table = "CREATE TABLE IF NOT EXISTS Parameter (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    condition_id INTEGER NOT NULL REFERENCES Condition(number),\
                    design_load TEXT, \
                    buoyancy TEXT,\
                    criterion TEXT, \
                    scaling_factor DECIMAL, \
                    title TEXT);"
        create_table(conn, table)
    #    
    