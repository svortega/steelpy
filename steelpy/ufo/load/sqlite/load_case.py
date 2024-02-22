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
        query = (self._component)
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
#
class BasicLoadSQLX(BasicLoadMain):
    """
    FE Load Cases

    LoadType
        |_ node
        |_ beam
        |_ shell

    **Parameters**:
      :number:  integer internal number
      :name:  string node external name
    """
    __slots__ = ['db_file', '_plane', '_basic', '_hydro']
    #
    def __init__(self, db_file:str, plane: NamedTuple) -> None:
        """
        """
        super().__init__()
        #
        self.db_file = db_file
        self._plane = plane
        #
        #self._nodes = NodeLoadGlobalSQL(db_file=self.db_file)
        #self._beams = BeamLoadGloabalSQL(db_file=self.db_file)
        self._basic = BasicLoadCase(db_file=self.db_file, plane=self._plane)
        #self._hydro = MetoceanLoadSQL(db_file=self.db_file, plane=self._plane)
        #
        conn = create_connection(self.db_file)
        with conn: 
            self._create_table(conn)
    #
    @property
    def _labels(self):
        """ """
        table = "SELECT Load.name FROM Load WHERE level = 'basic'"
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table)
            items = cur.fetchall()
        return [item[0] for item in items]
    #
    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += self._basic.__str__()
        #output += self._hydro.__str__()
        return output    
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
                               bd_file=self.db_file)
        except ValueError:
            raise IOError(f"Basic load case {load_name} not defined")

    #
    # -----------------------------------------------
    #
    def _push_load(self, conn, load_name:int|str, load_title:str):
        """ """
        #load_name = str(load_name)
        project = (load_name, load_title, "basic", None)
        sql = 'INSERT INTO Load(name, title, level, type) VALUES(?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
        #return cur.lastrowid
    #
    def _create_table(self, conn):
        """ """
        table_load = "CREATE TABLE IF NOT EXISTS Load(\
                      number INTEGER PRIMARY KEY NOT NULL,\
                      name NOT NULL,\
                      title TEXT NOT NULL,\
                      level TEXT NOT NULL,\
                      type TEXT );"

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
    #@property
    #def df(self):
    #    """basic load df"""
    #    print('basic load in')
    #    #db = DBframework()
    #    #
    #    conn = create_connection(self.db_file)
    #    with conn:        
    #        nodedf = node_load(conn)
    #        memgdf = member_load(conn)
    #    #
    #    #for name in self._labels:
    #    #    bload = self.__getitem__(name)
    #    #    nodedf = bload._node.df
    #    #    nodedf
    #    1 / 0
    #
    #@df.setter
    #def df(self, value):
    #    """basic load df"""
    #    print('basic load out')
    #    1 / 0    
    #
    #       
    #
    #def nodal(self, values:tuple|list|None=None,
    #          df=None):
    #    """Basic nodal load """
    #    if isinstance(values, (list, tuple)):
    #        if isinstance(values[0], (list,tuple)):
    #            for item in values:
    #                name = item[0]
    #                try:
    #                    self._basic[name]
    #                except IOError:
    #                    self._basic[name] = name
    #                #
    #                self._basic[name].node(item[1:])
    #        else:
    #            name = values[0]
    #            try:
    #                self._basic[name]
    #            except IOError:
    #                self._basic[name] = name
    #            #
    #            self._basic[name].node(values[1:])
    #    #
    #    #
    #    # dataframe input
    #    try:
    #        columns = list(df.columns)
    #        1 / 0
    #        #self._sections.df(df)
    #    except AttributeError:
    #        pass         
    #    #
    #    return self._basic._nodes        
    #
    #def beam(self, values:tuple|list|None=None,
    #         df=None):
    #    """Basic beam load"""
    #    if isinstance(values, (list, tuple)):
    #        if isinstance(values[0], (list,tuple)):
    #            for item in values:
    #                name = item[0]
    #                try:
    #                    self._basic[name]
    #                except IOError:
    #                    self._basic[name] = name
    #                #
    #                self._basic[name].beam(item[1:])
    #        else:
    #            name = values[0]
    #            try:
    #                self._basic[name]
    #            except IOError:
    #                self._basic[name] = name
    #            #
    #            self._basic[name].beam(values[1:])
    #    #
    #    #
    #    # dataframe input
    #    try:
    #        columns = list(df.columns)
    #        1 / 0
    #        #self._sections.df(df)
    #    except AttributeError:
    #        pass         
    #    #
    #    return self._basic._beams         
    #    
    #
#
#
def node_loadX(conn, component: int):
    """ """
    query = (component, )
    table = "SELECT Load.name, Load.title, Load.level,\
            Node.name, LoadNode.* \
            FROM Load, Node, LoadNode\
            WHERE Load.number = LoadNode.load_id \
            AND Node.number = LoadNode.node_id \
            AND LoadNode.type = 'load' \
            AND Load.component_id = ? ;"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    rows = cur.fetchall()
    #
    db = DBframework()
    cols = ['load_name', 'load_title', 'load_type', 'node_name', 
            'number', 'load_id', 'element_name', 
            'node_id','load_comment', 'load_system', 'type',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
            'x', 'y', 'z', 'rx', 'ry', 'rz']
    nodedf = db.DataFrame(data=rows, columns=cols)        
    #
    nodedf = nodedf[['load_name', 'load_type', 'load_name',
                     'load_system', 'load_comment', 'element_name',
                     'node_name', 'type',
                     'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                     'x', 'y', 'z', 'rx', 'ry', 'rz']]
    #
    #
    return nodedf
#
#
def member_loadX(conn, component: int):
    """ """
    query = (component, )
    table = "SELECT Load.name, Load.title, Load.level,\
                 Element.name, LoadBeamLine.* \
                 FROM Load, Element, LoadBeamLine \
                 WHERE Load.number = LoadBeamLine.load_id \
                 AND Element.number = LoadBeamLine.element_id \
                 AND LoadBeamLine.type = 'load' \
                 AND Load.component_id = ? ;"
    #
    cur = conn.cursor()
    # line
    cur.execute(table, query)
    rows = cur.fetchall()
    #
    #
    db = DBframework()
    cols = ['load_name', 'load_title', 'load_type', 'element_name', 
            'number', 'load_id', 'element_id', 
            'load_comment', 'load_system', 'type',
            'L0', 'qx0', 'qy0', 'qz0',
            'L1', 'qx1', 'qy1', 'qz1',
            'BS', 'OTM', 'x', 'y', 'z']
    #
    membdf = db.DataFrame(data=rows, columns=cols)      
    #
    membdf = membdf[['load_name', 'load_type', 'load_name',
                     'load_system', 'load_comment',
                     'element_name', 'type', 
                     'L0', 'qx0', 'qy0', 'qz0',
                     'L1', 'qx1', 'qy1', 'qz1',
                     'BS', 'OTM', 'x', 'y', 'z']].values    
    #
    return membdf
#