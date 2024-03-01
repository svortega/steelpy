#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
from typing import NamedTuple
#import re

# package imports
# steelpy.f2uModel
from steelpy.ufo.load.process.actions import SelfWeight
from steelpy.ufo.load.process.basic_load import LoadCaseBasic, BasicLoadMain, LoadTypeBasic
from steelpy.ufo.load.sqlite.beam import BeamLoadItemSQL, BeamLoadGloabalSQL
from steelpy.ufo.load.sqlite.node import  NodeLoadItemSQL, NodeLoadGlobalSQL
#
# steelpy.f2uModel
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework
#
#
#
class BasicLoadSQL(LoadCaseBasic):
    
    __slots__ = ['db_file', '_plane', '_component']
                # '_nodes', '_beams']
    #
    def __init__(self, db_file:str, plane: NamedTuple,
                 component: int) -> None:
        """
        """
        super().__init__()
        #
        self.db_file = db_file
        self._plane = plane
        self._component = component
        #
        self._nodes = NodeLoadGlobalSQL(component=self._component,
                                        db_file=self.db_file)
        #
        self._beams = BeamLoadGloabalSQL(component=self._component,
                                         plane=self._plane, 
                                         db_file=self.db_file)
        #
        conn = create_connection(self.db_file)
        with conn: 
            self._create_table(conn)        
    #
    @property
    def _labels(self):
        """ """
        query = ('basic', self._component, )
        table = "SELECT Load.name FROM Load \
                  WHERE level = ? AND component_id = ?;"
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]
    #
    #
    # -----------------------------------------------
    #
    def __setitem__(self, load_name:int|str, load_title:str) -> None:
        """
        """
        try:
            self._labels.index(load_name)
            raise Warning(f'    *** warning load name {load_name} already exist')
        except ValueError:
            conn = create_connection(self.db_file)
            with conn:
                self._push_load(conn, load_name, load_title)
    #           
    def __getitem__(self, load_name: str|int):
        """
        """
        try:
            index = self._labels.index(load_name)
            return LoadTypeSQL(load_name=load_name,
                               plane=self._plane,
                               component=self._component, 
                               bd_file=self.db_file)
        except ValueError:
            raise IOError("load case not defined")

    #
    # -----------------------------------------------
    #
    def _push_load(self, conn, load_name:int|str, load_title:str):
        """ """
        #load_name = str(load_name)
        project = (load_name,  self._component, "basic", load_title)
        table = 'INSERT INTO Load(name, component_id, level, title) \
                 VALUES(?,?,?,?)'
        cur = conn.cursor()
        cur.execute(table, project)
        #return cur.lastrowid
    #
    #
    def _create_table(self, conn):
        """ """
        table_load = "CREATE TABLE IF NOT EXISTS Load(\
                      number INTEGER PRIMARY KEY NOT NULL,\
                      name NOT NULL,\
                      component_id INTEGER NOT NULL REFERENCES Component(number), \
                      level TEXT NOT NULL,\
                      title TEXT );"

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
#
class LoadTypeSQL(LoadTypeBasic):
    """
    """
    __slots__ = ['_node', '_beam', '_selfweight', '_component', 
                 'name', 'number', 'title', '_db_file', '_plane']

    def __init__(self, load_name: str|int,
                 plane: NamedTuple,
                 component: int, 
                 bd_file:str):
        """
        """
        #super().__init__(name, number, title)
        self.name = load_name
        self._db_file = bd_file
        self._plane = plane
        self._component = component
        #
        self.number, self.title = self._load_spec(self.name)
        #
        self._node = NodeLoadItemSQL(load_name=load_name,
                                     component=component, 
                                     db_file=self._db_file)
        #
        self._beam = BeamLoadItemSQL(load_name=load_name,
                                     plane=self._plane,
                                     component=component, 
                                     db_file=self._db_file)
        #
        #                             
        #
        self._selfweight = SelfWeight()          
    #
    #
    #def __setitem__(self, load_name: str | int,
    #                properties: list[float]) -> None:
    #    """
    #    """
    #    super().__setitem__(load_name, properties)
    #  
    #    #1 / 0
    #
    @property
    def _labels(self):
        """ """
        query = (self._component, )
        table = 'SELECT Load.name \
                 FROM Load WHERE component_id = ?;'
        conn = create_connection(self._db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]    
    #
    def _load_spec(self, load_name: str|int):
        """ """
        query = (load_name, self._component)
        table = 'SELECT Load.number, Load.title\
                 FROM Load \
                 WHERE Load.name = ? \
                 AND component_id = ?;'
        #
        conn = create_connection(self._db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchone()
        #
        return items
    #
    #  
#
#
# ---------------------------------
#
