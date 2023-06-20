#
# Copyright (c) 2009-2023 fem2ufo
# 

# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
#from dataclasses import dataclass
#from typing import NamedTuple
#import re

# package imports
# steelpy.f2uModel
from steelpy.f2uModel.load.process.actions import SelfWeight
from steelpy.f2uModel.load.process.basic_load import BasicLoadBasic, LoadTypeBasic
from steelpy.f2uModel.load.sqlite.beam import BeamLoadItemSQL
from steelpy.f2uModel.load.sqlite.node import  NodeLoadItemSQL
from steelpy.f2uModel.load.sqlite.wave_load import WaveLoadItemSQL
#
from steelpy.f2uModel.mesh.sqlite.beam import BeamItemSQL
# steelpy.f2uModel
from steelpy.f2uModel.mesh.process.process_sql import create_connection, create_table
#
#import pandas as pd
#
# ---------------------------------
#
#
class BasicLoadSQL(BasicLoadBasic):
    """
    FE Load Cases

    LoadType
        |_ name
        |_ number
        |_ basic
        |_ combination_level
        |_ time_series
        |_
        |_ temperature

    **Parameters**:
      :number:  integer internal number
      :name:  string node external name
    """
    __slots__ = ['_labels', '_title','_number', 'gravity',
                 '_basic', 'db_file']
    #
    def __init__(self, db_file:str):
        """
        """
        super().__init__()
        self.db_file = db_file
        self._basic: dict = {}
        conn = create_connection(self.db_file)
        with conn: 
            self._create_table(conn)
    #
    #
    def __setitem__(self, load_name:int|str, load_title:str) -> None:
        """
        """
        load_name = str(load_name)
        try:
            self._labels.index(load_name)
            raise Warning('    *** warning load name {:} already exist'
                            .format(load_name))
        except ValueError:
            self._labels.append(load_name)
            self._title.append(load_title)
            #
            conn = create_connection(self.db_file)
            with conn:
                load_number = self._push_basic_load(conn, load_name, load_title)
            self._number.append(load_number)
            #
            self._basic[load_name] = LoadTypeSQL(name=load_name,
                                                 number=load_number,
                                                 title=load_title,
                                                 bd_file=self.db_file)
    #           
    #
    def __getitem__(self, load_name: str|int):
        """
        """
        try:
            load_name = str(load_name)
            return self._basic[load_name]
        except KeyError:
            raise IOError("load case not defined")

    #
    def __delitem__(self, load_name: str|int):
        """
        """
        load_name = str(load_name)
        del self._basic[load_name]
    #
    #
    def _push_basic_load(self, conn, load_name:int|str, load_title:str):
        """ """
        load_name = str(load_name)
        project = (load_name, load_title, "basic")
        sql = 'INSERT INTO tb_Load(name, title, type) VALUES(?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
        return cur.lastrowid
    #
    def _create_table(self, conn):
        """ """
        table_load = "CREATE TABLE IF NOT EXISTS tb_Load(\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name TEXT NOT NULL,\
                    title TEXT NOT NULL,\
                    type TEXT NOT NULL);"

        table_comb_load = "CREATE TABLE IF NOT EXISTS tb_LoadCombIndex(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            load_number INTEGER NOT NULL REFERENCES tb_Load(number),\
                            bl_number INTEGER REFERENCES tb_Load(number),\
                            lc_number INTEGER REFERENCES tb_Load(number),\
                            factor DECIMAL NOT NULL);"

        #conn = create_connection(self.db_file)
        create_table(conn, table_load)
        create_table(conn, table_comb_load)
    #
    #
    #
#
#
class LoadTypeSQL(LoadTypeBasic):
    """
    """
    __slots__ = ['_node', '_beam', '_selfweight', '_wave', 
                 'name', 'number', 'title', '_bd_file']

    def __init__(self, name: str, number: int, title: str, bd_file:str):
        """
        """
        super().__init__(name, number, title)
        #self.name = name
        #self.number = number
        #self.title = title        
        self._bd_file = bd_file
        #
        self._node = NodeLoadItemSQL(load_name=self.name,
                                     bd_file=self._bd_file)
        self._beam = BeamLoadItemSQL(load_name=self.name,
                                     bd_file=self._bd_file)
        self._selfweight = SelfWeight()
        self._wave = WaveLoadItemSQL(load_name=self.name,
                                     bd_file=self._bd_file)
    #
    #
#
