# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
from datetime import datetime as dt
#from typing import NamedTuple
#import time
import os
#

# package imports
from steelpy.ufo.mesh.sqlite.nodes import NodeSQL
from steelpy.ufo.mesh.sqlite.elements import ElementsSQL
from steelpy.ufo.mesh.sqlite.boundary import BoundarySQL

from steelpy.utils.io_module.inout import check_input_file
from steelpy.ufo.process.main import ufoBasicModel #, ModelClassBasic
from steelpy.utils.sqlite.utils import create_connection, create_table
#
from steelpy.utils.sqlite.main import ClassMainSQL
#
#
#
class MeshSQL(ClassMainSQL, ufoBasicModel):
    """
    mesh[beam_name] = [number, element1, element2, elementn]
    """
    def __init__(self, component:str|int,
                 name:str|None = None,
                 sql_file:str|None = None):
        """
        """
        super().__init__(name=name, sql_file=sql_file)
        #
        self.data_type = "sqlite"
        self._component = component
        # --------------------------------------------------
        #
        self._nodes = NodeSQL(db_system=self.data_type,
                              component=self._component, 
                              db_file=self.db_file)
        #
        self._boundaries = BoundarySQL(db_system=self.data_type,
                                       component=self._component,
                                       db_file=self.db_file)
        #
        self._elements = ElementsSQL(db_system=self.data_type,
                                     component=self._component,
                                     db_file=self.db_file)
        #         
#
#
class UFOmain(ufoBasicModel):
    """
    mesh[beam_name] = [number, element1, element2, elementn]
    """
    __slots__ = ['_name', 'db_file', '_component', #'_plane', '_df', 
                 'data_type', '_build']

    def __init__(self, name:str|int|None = None,
                 component:str|int|None = None, 
                 sql_file:str|None = None):
        """
        """
        self.data_type = "sqlite"
        #
        self._build = True
        if sql_file:
            #print('--')
            sql_file = check_input_file(file=sql_file,
                                        file_extension="db")
            self.db_file = sql_file
            self._build = False
            # fixme: name
            #self._name = sql_file.split('.')[0]
        else:
            self.db_file = self._get_file(name=name)
            #self.data_type = mesh_type
            self._name = name
            # FIXME : how to handle component? 
            if not component:
                component = name
            #
            conn = create_connection(self.db_file)
            with conn:
                self._new_table(conn)
                self._component = self._push_data(conn, component)
        #
        # --------------------------------------------------
        #
        #self._nodes = NodeSQL(db_system=self.data_type,
        #                      component=self._component, 
        #                      db_file=self.db_file)
        #
        #self._boundaries = BoundarySQL(db_system=self.data_type,
        #                               component=self._component,
        #                               db_file=self.db_file)
        #
        #self._elements = ElementsSQL(db_system=self.data_type,
        #                             component=self._component,
        #                             db_file=self.db_file)
        #        
    #
    #
    #
    def _get_file(self, name: str):
        """ """
        #BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        filename = name + ".db"
        path = os.path.abspath(filename)
        #self.db_file = path
        #directory = os.path.dirname(path)
        #
        #self.db_file:str = component + "_f2u.db"
        #if mesh_type != "ttt": #"inmemory":
        try: # remove file if exist
            os.remove(path)
        except FileNotFoundError:
            pass
        #
        return path
    #    
    #
    # --------------------------------------------
    # SQL ops
    #
    def _new_table(self, conn) -> None:
        """ """
        table = "CREATE TABLE IF NOT EXISTS Component (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    plane TEXT NOT NULL, \
                    date TEXT NOT NULL,\
                    title TEXT);"
        create_table(conn, table)
    #
    #
    def _push_data(self, conn,
                   component: str|int|None,
                   plane: str = '3D', 
                   title: str|None = None):
        """ """
        table = 'INSERT INTO Component(name, plane,\
                                       date, title)\
                            VALUES(?,?,?,?)'
        #
        date = dt.now().strftime('%Y-%m-%d')
        data = (component,  plane, date, title)
        # push
        cur = conn.cursor()
        out = cur.execute(table, data)
        return out.lastrowid
    #
    #
    def _set_type(self, component: str|int,
                  comp_type: str, title: str|None):
        """ """
        
        time=dt.now().strftime('%Y-%m-%d')
        #item = 'concept'
        #
        query = (time, comp_type, title, component)
        table = f"UPDATE Component \
                 SET date = ?, \
                     type = ?, \
                     title = ? \
                 WHERE name = ?;"
        #
        conn = create_connection(self.db_file)
        with conn:          
            cur = conn.cursor()
            comp = cur.execute(table, query)
        #
        if not comp:
            raise IOError(f' component {component} not valid')
    #    
