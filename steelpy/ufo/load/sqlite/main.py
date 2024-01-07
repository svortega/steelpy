#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from typing import NamedTuple, Dict, List, Iterable, Union

# package imports
from .load_case import BasicLoadSQL
from .combination import LoadCombSQL
from steelpy.utils.sqlite.utils import create_connection, create_table
#
#
#
class LoadingSQL:
    __slots__ = ['_case', '_combination', 'db_file']

    def __init__(self, db_file: str,
                 db_system:str="sqlite"):
        """
        """
        1 / 0
        self.db_file = db_file
        self._case = BasicLoadSQL(db_file)
        self._combination = LoadCombSQL(db_file)
        # create node table
        conn = create_connection(self.db_file)
        with conn: 
            self._create_table()

    @property
    def case(self):
        """
        """
        return self._case

    #
    @property
    def combination(self):
        """
        """
        return self._combination
    #
    def _create_table(self):
        """ """
        table_load = "CREATE TABLE IF NOT EXISTS Load(\
                      number INTEGER PRIMARY KEY NOT NULL,\
                      name NOT NULL,\
                      title TEXT NOT NULL,\
                      type TEXT NOT NULL);"

        table_comb_load = "CREATE TABLE IF NOT EXISTS LoadCombination(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            load_id INTEGER NOT NULL REFERENCES Load(number),\
                            bl_number INTEGER REFERENCES Load(number),\
                            lc_number INTEGER REFERENCES Load(number),\
                            factor DECIMAL NOT NULL);"

        #conn = create_connection(self.db_file)
        create_table(conn, table_load)
        create_table(conn, table_comb_load)
    #
    #
    def __str__(self) -> str:
        """ """
        self._case.__str__()
        self._combination.__str__()
#
#