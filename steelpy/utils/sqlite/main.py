#
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
#from dataclasses import dataclass
#from operator import itemgetter
#from typing import NamedTuple
#import sqlite3 as sqlite3
#from sqlite3 import Error
#from datetime import datetime as dt
import os
#import re
#

# package imports
#
from steelpy.utils.sqlite.utils import create_connection #, create_table
from steelpy.utils.io_module.inout import check_input_file


class ClassMainSQL:
    __slots__ = ['db_file', '_build', '_name']
    
    def __init__(self, name:str|None,
                 sql_file:str|None):
        """
        """
        #
        self._build = True 
        if sql_file:
            sql_file = check_input_file(file=sql_file,
                                        file_extension="db")
            self.db_file = sql_file
            self._build = False
        elif name:
            self.db_file = self._get_file(name=name)
            self._name = name
            #conn = create_connection(self.db_file)
            #with conn:
            #    self._create_table(conn)
        else:
            IOError('Project name or SQL file missing')
    #
    #
    def _get_file(self, name: str):
        """ """
        filename = name + ".db"
        path = os.path.abspath(filename)
        try: # remove file if exist
            os.remove(path)
        except PermissionError:
            raise PermissionError(f'close sqlite {name}')
        except FileNotFoundError:
            pass
        #
        return path
#
#
#
#
class ClassBasicSQL(Mapping):
    """ """
    __slots__ = ['db_file']
    #
    def __init__(self, db_file: str):
        """
        """
        self.db_file = db_file
        # create table
        conn = create_connection(self.db_file)
        with conn:
            self._new_table(conn)
    
    def __len__(self) -> int:
        return len(self._labels)

    def __iter__(self):
        return iter(self._labels)
    
    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    # TODO : get number from database
    def get_number(self, start:int=0):
        """
        """
        try:
            n = max(self._labels)
        except ValueError:
            n = start
        #
        while True:
            n += 1
            yield n
#