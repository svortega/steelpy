# 
# Copyright (c) 2009 fem2ufo
#

# Python stdlib imports
from __future__ import annotations

# package imports
from steelpy.metocean.hydrodynamic.utils import HydroBasic
from steelpy.utils.sqlite.utils import create_connection, create_table


class Wind(HydroBasic):
    """
    """
    __slots__ = ['_wind', 'db_file']
    
    def __init__(self, db_file: str):
        """
        """
        super().__init__(db_file)
        #
        self._wind:dict = {}
    #
    #
    def __getitem__(self, name):
        """
        """
        1 / 0
        conn = create_connection(self.db_file)
        with conn:
            data = self._get_data(conn, current_name=name)
        #
    #
    def __setitem__(self, name:str|int, value) -> None:
        """
        """
        1 / 0
        conn = create_connection(self.db_file)
        with conn:
            self._push_data(conn, data)
    #
    # ------------------
    # SQL ops
    # ------------------
    #
    def _create_table(self, conn) -> None:
        """ """
        # Main
        table = "CREATE TABLE IF NOT EXISTS tb_Wind (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    type TEXT NOT NULL,\
                    velocity DECIMAL NOT NULL,\
                    H0 DECIMAL, \
                    formula TEXT, \
                    period_ratio DECIMAL, \
                    height_exponent DECIMAL, \
                    gus_factor DECIMAL, \
                    duration DECIMAL);"
        create_table(conn, table)
        # Profile
        #table = "CREATE TABLE IF NOT EXISTS tb_CurrentProfile (\
        #                number INTEGER PRIMARY KEY NOT NULL,\
        #                current_number NOT NULL REFERENCES tb_Current(number),\
        #                elevation DECIMAL NOT NULL,\
        #                factor DECIMAL NOT NULL);"
        #create_table(conn, table)
    #    
