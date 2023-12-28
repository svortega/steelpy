# 
# Copyright (c) 2019 steelpy
# 
# Python stdlib imports
from __future__ import annotations

# package imports
from steelpy.utils.sqlite.utils import create_connection, create_table

#
#
class HydroDesign:
    """ Hydrodynamic Design """
    __slots__ = ['_name', '_db_file']
    
    def __init__(self, db_file:str|None = None):
        """
        """
        #self.name = name
        self._db_file = db_file
        #
        # create table
        conn = create_connection(self._db_file)
        with conn:
            self._create_table(conn)
    #
    # ----------------------------------------
    # SQL ops
    #
    def _create_table(self, conn) -> None:
        """ """
        #
        table = "CREATE TABLE IF NOT EXISTS tb_Parameters (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    condition_number INTEGER NOT NULL REFERENCES tb_Condition(number),\
                    design_load TEXT, \
                    buoyancy TEXT,\
                    criterion TEXT, \
                    scaling_factor DECIMAL, \
                    title TEXT);"
        create_table(conn, table)
    #    
    