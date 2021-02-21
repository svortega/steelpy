# 
# Copyright (c) 2009-2021 fem2ufo
#

# Python stdlib imports
#from typing import NamedTuple, Dict, List, Iterable, Union

# package imports
from steelpy.f2uModel.load.sqlite.basic_load import BasicLoadSQL
from steelpy.f2uModel.load.sqlite.combination import LoadCombSQL
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table
#
#
#
class LoadingSQL:
    __slots__ = ['_basic', '_combination', 'db_file']

    def __init__(self, db_file: str,
                 db_system:str="sqlite"):
        """
        """
        self.db_file = db_file
        self._basic = BasicLoadSQL(db_file)
        self._combination = LoadCombSQL(db_file)
        # create node table
        self._create_table()

    @property
    def basic(self):
        """
        """
        return self._basic

    #
    @property
    def combination(self):
        """
        """
        return self._combination
    #
    def _create_table(self):
        """ """
        table_load = "CREATE TABLE IF NOT EXISTS tb_Load(\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name INTEGER NOT NULL,\
                    title TEXT NOT NULL,\
                    type TEXT NOT NULL);"

        table_comb_load = "CREATE TABLE IF NOT EXISTS tb_LoadCombIndex(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            load_number INTEGER NOT NULL REFERENCES tb_Load(number),\
                            bl_number INTEGER REFERENCES tb_Load(number),\
                            lc_number INTEGER REFERENCES tb_Load(number),\
                            factor DECIMAL NOT NULL);"

        conn = create_connection(self.db_file)
        create_table(conn, table_load)
        create_table(conn, table_comb_load)
#
#