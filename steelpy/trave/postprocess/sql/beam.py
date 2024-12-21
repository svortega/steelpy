# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
#from collections.abc import Mapping
#from typing import NamedTuple
#import re
#
# package imports
from steelpy.trave.postprocess.utils.beam import (BeamResBasic, BeamForce,
                                                  BeamDeflection, BeamStress)
from steelpy.utils.sqlite.utils import create_connection #, create_table
from steelpy.ufo.mesh.main import MeshPlane
from steelpy.utils.dataframe.main import DBframework

#
class BeamResSQL(BeamResBasic):
    __slots__ = ['_mesh', 'db_file', '_result_name',
                 '_beams']
    
    def __init__(self, mesh, result_name:int|str,
                 db_file: str)-> None:
        """ """
        super().__init__(mesh)
        self.db_file = db_file
        self._result_name = result_name
    #
    @property
    def _labels(self):
        """ """
        query = (self._mesh._name, )
        table = "SELECT Element.name \
                 FROM Element, Result, Mesh \
                 WHERE Mesh.name = ? \
                 AND Element.mesh_id = Mesh.number \
                 AND Result.mesh_id = Mesh.number ;"
        #
        # Extract data from sqlite
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            rows = cur.fetchall()
        #
        labels = set([item[0] for item in rows])
        return list(labels) 
    #
    def __getitem__(self, beam_name: int|str)-> BeamResItem|IndexError:
        """ """
        try:
            self._labels.index(beam_name)
            return BeamResItem(beam=self._beams[beam_name],
                               mesh_name=self._mesh._name,
                               result_name=self._result_name, 
                               plane=self._mesh._plane.plane2D,
                               db_file=self.db_file)
        except ValueError:
            raise IndexError(' ** element {:} does not exist'.format(beam_name))
    #
    # -----------------------------
    #
    def _plane(self)-> MeshPlane:
        """ """
        conn = create_connection(self.db_file)
        with conn:
            sql = 'SELECT * FROM Component'
            cur = conn.cursor()
            cur.execute(sql)
            data = cur.fetchone()
        #
        if data[5] == '2D':
            return MeshPlane(plane2D=True)
        return MeshPlane(plane2D=False)
    #
    # -----------------------------
    #
    def force(self, units:str='si')-> BeamForce:
        """beam integration forces"""
        #
        conn = create_connection(self.db_file)
        with conn:        
            df = get_force(conn,
                           result_name=self._result_name,
                           mesh_name=self._mesh._name)
        #
        if self._mesh._plane.plane2D:
            df.drop(['Fz', 'Mx', 'My'], axis=1, inplace=True)
        #
        return BeamForce(df, units=units)
    #
    def displacement(self, units:str='si'):
        """beam integration forces"""
        
        return self.deflection(units)
        
    #
    def deflection(self, units:str='si')-> BeamDeflection:
        """beam integration forces"""
        #
        conn = create_connection(self.db_file)
        with conn:        
            df = get_displacement(conn,
                                  result_name=self._result_name,
                                  mesh_name=self._mesh._name)
        #
        if self._mesh._plane.plane2D:
            df.drop(['z', 'rx', 'ry'], axis=1, inplace=True)        
        #
        return BeamDeflection(df, units=units)    
    #
    def stress(self, units:str='si')-> BeamStress:
        """ beam stress """
        #
        conn = create_connection(self.db_file)
        with conn:        
            df = get_stress(conn,
                            result_name=self._result_name,
                            mesh_name=self._mesh._name)
        #
        return BeamStress(df, units=units)         
    #
    #def code_check(self):
    #    """ beam code check"""
    #    1 / 0
    #
    #def pull_data(self, query: str, cols: list):
    #    """ """
    #    conn = create_connection(self.db_file)
    #    with conn:
    #        cur = conn.cursor()
    #        cur.execute(query)
    #        data = cur.fetchall()
    #    #       
    #    db = DBframework()
    #    return db.DataFrame(data=data, columns=cols)    
    #
#
#
#
@dataclass
class BeamResItem:
    __slots__ = ['_beam', '_db_file', '_plane',
                 '_mesh_name', '_result_name']
    
    def __init__(self, beam,
                 mesh_name: int|str,
                 result_name: int|str, 
                 plane: bool, db_file: str)-> None:
        """ """
        self._beam = beam
        self._plane = plane
        self._db_file = db_file
        self._mesh_name = mesh_name
        self._result_name = result_name
    #
    def force(self, units:str='si')-> BeamForce:
        """node force"""
        beam_name = self._beam.name
        conn = create_connection(self._db_file)
        with conn:        
            df = get_force(conn,
                           element_name=beam_name,
                           result_name=self._result_name,
                           mesh_name=self._mesh_name)
        #
        if self._plane:
            df.drop(['Fz', 'Mx', 'My'], axis=1, inplace=True)
        #
        return BeamForce(df, units=units)
    #
    def deflection(self, units:str='si')-> BeamDeflection:
        """ """
        beam_name = self._beam.name
        conn = create_connection(self._db_file)
        with conn:        
            df = get_displacement(conn,
                                  element_name=beam_name,
                                  result_name=self._result_name,
                                  mesh_name=self._mesh_name)
        #
        if self._plane:
            df.drop(['z', 'rx', 'ry'], axis=1, inplace=True)
        #        
        return BeamDeflection(df, units=units)
    #
    def stress(self, units:str='si')-> BeamStress:
        """node force"""
        beam_name = self._beam.name
        conn = create_connection(self._db_file)
        with conn:        
            df = get_stress(conn,
                            element_name=beam_name,
                            result_name=self._result_name,
                            mesh_name=self._mesh_name)
        #
        return BeamStress(df, units=units)
#
# --------------------
# sql operations
#
def get_force(conn, result_name: int|str,
              mesh_name: int|str,
              element_name: int|str|None=None,
              item:str='*'):
    """ """
    query = (result_name, mesh_name, )
    table = (f'SELECT Load.name, Load.level, Component.name, Element.name,\
              ResultBeamForce.{item} \
              FROM Load,  Mesh, Component, Element, Result, ResultBeamForce \
              WHERE Result.name = ? \
              AND Mesh.name = ? \
              AND Result.mesh_id = Mesh.number \
              AND Result.number = ResultBeamForce.result_id \
              AND Mesh.component_id = Component.number \
              AND Element.number = ResultBeamForce.element_id \
              AND Load.number = ResultBeamForce.load_id \
              ')
              #WHERE Load.number = ResultBeamForce.load_id\
              #   AND Result.number = ResultBeamForce.result_id\
              #   AND Element.number = ResultBeamForce.element_id\
              #   AND Result.mesh_id = Mesh.number\
              #   AND Mesh.component_id = Component.number')

    if element_name:
        table += f' AND Element.name = {element_name}'
    table += ';'
    #
    cols = ['load_name', 'load_level', 'component_name', 'element_name',
            'number', 'result_id', 'load_id',
            'element_id', 'length', 'system',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
            'theta1', 'theta2', 'theta3']
    df = pull_data(conn, table, query, cols)
    #
    cols = ['number', 'component_name', 'load_name',
            'load_level', 'element_name', 'length', 'system',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    return df[cols]
#
#
def get_displacement(conn, result_name: int|str,
                     mesh_name: int|str,
                     element_name: int|str|None=None,
                     item:str='*'):
    """ """
    query = (result_name, mesh_name, )
    table = (f'SELECT Load.name, Load.level, Component.name, Element.name,\
              ResultBeamDeflection.{item} \
              FROM Load, Mesh, Component, Element, Result, ResultBeamDeflection \
              WHERE Result.name = ? \
              AND Mesh.name = ? \
              AND Result.mesh_id = Mesh.number \
              AND Result.number = ResultBeamDeflection.result_id \
              AND Mesh.component_id = Component.number \
              AND Element.number = ResultBeamDeflection.element_id \
              AND Load.number = ResultBeamDeflection.load_id \
              ')              
              #WHERE Load.number = ResultBeamDeflection.load_id\
              #   AND Result.number = ResultBeamDeflection.result_id\
              #   AND Element.number = ResultBeamDeflection.element_id\
              #   AND Result.mesh_id = Mesh.number\
              #   AND Mesh.component_id = Component.number')

    if element_name:
        table += f' AND Element.name = {element_name}'
    table += ';'
    #
    cols = ['load_name', 'load_level', 'component_name', 'element_name',
            'number', 'result_id', 'load_id',
            'element_id', 'length', 'system',
            'x', 'y', 'z', 'rx', 'ry', 'rz']
    df = pull_data(conn, table, query, cols)
    #
    cols = ['number', 'component_name', 'load_name',
            'load_level', 'element_name', 'length', 'system',
            'x', 'y', 'z', 'rx', 'ry', 'rz']
    return df[cols]
#
#
def get_stress(conn, result_name: int|str,
               mesh_name: int|str,
               element_name: int|str|None=None,
               item:str='*'):
    """ """
    query = (result_name, mesh_name, )
    table = (f'SELECT Load.name, Load.level, Component.name, Element.name,\
               ResultBeamStress.{item} \
               FROM Load, Node, Mesh, Component, Element, Result, ResultBeamStress \
               WHERE Result.name = ? \
               AND Mesh.name = ? \
               AND Result.mesh_id = Mesh.number \
               AND Result.number = ResultBeamStress.result_id \
               AND Mesh.component_id = Component.number \
               AND Element.number = ResultBeamStress.element_id \
               AND Load.number = ResultBeamStress.load_id \
               ')                     
               #WHERE Load.number = ResultBeamStress.load_id\
               #   AND Result.number = ResultBeamStress.result_id\
               #   AND Element.number = ResultBeamStress.element_id\
               #   AND Result.mesh_id = Mesh.number\
               #   AND Mesh.component_id = Component.number')

    if element_name:
        table += f' AND Element.name = {element_name}'
    table += ';'
    #
    cols = ['load_name', 'load_level', 'component_name', 'element_name',
            'number', 'result_id', 'load_id',
            'element_id', 'length', 'system',
            'stress_point', 'y', 'z',
            'tau_x', 'tau_y', 'tau_z',
            'sigma_x', 'sigma_y', 'sigma_z']
    df = pull_data(conn, table, query, cols)
    #
    cols = ['number', 'component_name', 'load_name',
            'load_level', 'element_name', 'length', 'system',
            'stress_point', 'y', 'z',
            'tau_x', 'tau_y', 'tau_z',
            'sigma_x', 'sigma_y', 'sigma_z']
    return df[cols]
#
#
def pull_data(conn, table: str, query: tuple, cols: list):
    """ """
    cur = conn.cursor()
    cur.execute(table, query)
    data = cur.fetchall()
    #
    db = DBframework()
    return db.DataFrame(data=data, columns=cols)
#