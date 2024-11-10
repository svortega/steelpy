# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
#import os
from datetime import datetime as dt
#
# package imports
from steelpy.ufo.mesh.main import Mesh
from steelpy.ufo.concept.main import Concept
from steelpy.ufo.properties.main import Properties
from steelpy.ufo.utils.main import ufoBasicModel
from steelpy.utils.sqlite.main import get_db_file
from steelpy.utils.sqlite.utils import create_connection, create_table

# 
#

#
class UFOmain(ufoBasicModel):
    """
    mesh[beam_name] = [number, element1, element2, element]
    """
    __slots__ = ['_name', 'db_file', '_component',
                 'data_type', '_build']

    def __init__(self, name:str|int|None = None,
                 #component:str|int|None = None, 
                 sql_file:str|None = None):
        """
        """
        #
        self.db_file, self._name, self._build = get_db_file(name, sql_file)
        #
        #if not component:
        component = self._name
        #
        super().__init__(component)
        self.data_type = "sqlite"         
        #
        if self._build:
            conn = create_connection(self.db_file)
            with conn:
                self._new_table(conn)
                self._component = self._push_data(conn, component)
    #
    # --------------------------------------------------
    #       
    #
    #
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
#
#
#
#
class UFOmodel(UFOmain):
    """
    UFO Model class
    
    Parameters:
      :number: integer
      :name: string
      :nodes: list of node's class
      :elements: list of element's class
      :materials: list of material's class
      :sets: list of groups (elements, nodes, etc)
      :sections: list of section's class
      :vectors: list of guide points
      :eccentricities: list of eccentricities
      :joints: list of joint's class
      :hinges: list of hinges definitios
      :loads: list of load's class
      :data: FE model data
      :units: FE model units
      :soil: list of soil's class
      :hydrodynamics: hydrodynamic's data
      :boundaries: list of FE model boundary
    
    Parameters:  
      :number:  integer internal number 
      :name:  string node external name
    """
    __slots__ = ['_name', '_properties', 
                 '_mesh', '_concept',
                 'db_file', '_component',
                 'data_type', '_build']

    def __init__(self, name:str|int|None = None,
                 sql_file:str|None = None) -> None:
        """
        mesh_type : sqlite/inmemory
        """
        print("-- module : ufo Version 6.50dev")
        print('{:}'.format(52*'-'))
        #
        super().__init__(name=name,
                         sql_file=sql_file)        
        #
        # mesh basic
        self._mesh = Mesh(component=self._component, 
                          sql_file=self.db_file,
                          build=self._build)
        #self._name = name
        self._properties = Properties()
        #
        self._concept = Concept(component=self._component,
                                mesh=self._mesh, 
                                properties=self._properties)

    #
    #
    # -------------------
    #
    #
    def properties(self):
        """
        """
        return self._properties
    #
    # -------------------
    #
    #@property
    def concept(self):
        """
        """
        return self._concept
    #
    # -------------------
    #
    def mesh(self): #name:str|None = None,
             #sql_file:str|None = None):
        """ """
        return self._mesh

    #
    # -------------------
    #
    def build(self, name:str|None = None) -> None:
        """
        """
        #
        if not name:
            name = self._name        
        #
        #self._sections.get_properties()
        #
        if self._concept:
            self._mesh[name] = self._concept.mesh()
        #    meshing = Meshing(concept=self._concept,
        #                      component=name, 
        #                      mesh_type=self.mesh_type)
        #    self._mesh[name] = meshing.get_mesh()
        #    self._mesh[name].renumbering()
        #    #mesh._load._basic.FER()
        #    #return mesh
        #    #_sql.write_concept(self._concept)
        #
        # check wave load case
        #
        for key, item in self._mesh.items():
            item.build()
            #item._load._basic.wave_process()
            # TODO : remove second _load for simplification
            #item._load._basic.FER(elements= item._elements)        
        #
        #
        print('end meshing')
        return self._mesh[name]
    #
    #@property
    #def plot(self):
    #    """ """
    #    return self._plot
#
#
