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
from steelpy.ufo.mesh.main import MeshPlane
from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.sqlite.utils import create_connection #, create_table

from steelpy.trave.postprocess.utils.node import (NodeResultBasic, NodeForce,
                                                  NodeDeflection, NodeReaction)
#
from steelpy.ufo.mesh.sqlite.utils import (pull_node_mesh,
                                           pull_results_mesh,
                                           pull_load_mesh)

#
class NodeResultSQL(NodeResultBasic):
    __slots__ = ['_labels', '_mesh', 
                 'plane', 'db_file', '_result_name']
    
    def __init__(self, mesh, result_name:int|str,
                 db_file: str):
        """ """
        super().__init__(mesh)
        self.db_file = db_file
        self._result_name = result_name
    #
    # -----------------------------
    #
    #@property
    #def _plane(self)-> MeshPlane:
    #    """ """
    #    conn = create_connection(self.db_file)
    #    with conn:
    #        sql = 'SELECT * FROM Component'
    #        cur = conn.cursor()
    #        cur.execute(sql)
    #        data = cur.fetchone()
    #    #
    #    if data[5] == '2D':
    #        return MeshPlane(plane2D=True)
    #    return MeshPlane(plane2D=False)
    #    
    #@property
    #def _labels(self):
    #    """ """
    #    #nodes = list(self._mesh._nodes.keys())
    #    #table = 'SELECT Nodes.name FROM Nodes ORDER BY number ASC'
    #    #conn = create_connection(self.db_file)
    #    #with conn:        
    #    #    cur = conn.cursor()
    #    #    cur.execute(table)
    #    #    items = cur.fetchall()
    #    return list(self._mesh._nodes.keys())
    #
    # ---------------------------------
    #
    def __getitem__(self, node_name: int|str) -> NodeResultItem | IndexError:
        """
        node_name : node number
        """
        try:
            self._labels.index(node_name)
            return NodeResultItem(node=self._mesh._nodes[node_name],
                                  mesh_name=self._mesh._name,
                                  result_name=self._result_name,                                  
                                  plane=self._mesh._plane.plane2D,
                                  db_file=self.db_file)
        except ValueError:
            raise IndexError('   *** node {:} does not exist'.format(node_name))
    #
    # ---------------------------------
    #
    #@property
    def force(self, units:str='si')-> NodeForce:
        """ node force"""
        conn = create_connection(self.db_file)
        with conn:        
            df = get_force(conn,
                           result_name=self._result_name,
                           mesh_name=self._mesh._name)                           
                           #plane2D=self.plane.plane2D)
        #if self._mesh._plane.plane2D:
        #    df.drop(['Fz', 'Mx', 'My'], axis=1, inplace=True)
        return NodeForce(df, units=units)
    #
    #@property
    def displacement(self, units:str='si')->NodeDeflection:
        """ node displacement"""
        conn = create_connection(self.db_file)
        with conn:         
            df = get_displacement(conn,
                                  result_name=self._result_name,
                                  mesh_name=self._mesh._name)
                                  #plane2D=self.plane.plane2D)
        #
        #if self._mesh._plane.plane2D:
        #    df.drop(['z', 'rx', 'ry'], axis=1, inplace=True)
        return NodeDeflection(df, units=units)
    #
    #@property
    def reaction(self, units:str='si')->NodeReaction:
        """ Node reactions"""
        conn = create_connection(self.db_file)
        with conn:        
            df = get_reactions(conn,
                               result_name=self._result_name,
                               mesh_name=self._mesh._name)                               
                               #plane2D=self.plane.plane2D)
        #
        #if self._mesh._plane.plane2D:
        #    df.drop(['Fz', 'Mx', 'My'], axis=1, inplace=True)
        return NodeReaction(df, units=units)
    #
#        
#
#
@dataclass
class NodeResultItem:
    __slots__ = ['_node', '_db_file', '_plane',
                 '_mesh_name', '_result_name']
    
    def __init__(self, node,
                 mesh_name: int|str,
                 result_name: int|str,
                 plane, db_file: str)-> None:
        """ """
        self._node = node
        self._plane = plane
        self._db_file = db_file
        self._mesh_name = mesh_name
        self._result_name = result_name        
    #
    def force(self, units:str='si')-> NodeForce:
        """node force"""
        node_name = self._node.name
        conn = create_connection(self._db_file)
        with conn:        
            df = get_force(conn,
                           #plane2D=self._plane.plane2D,
                           node_name=node_name,
                           result_name=self._result_name,
                           mesh_name=self._mesh_name)
        #
        if self._plane:
            df.drop(['Fz', 'Mx', 'My'], axis=1, inplace=True)
        return NodeForce(df, units=units)
    #
    def displacement(self, units:str='si')-> NodeDeflection:
        """ node displacement"""
        node_name = self._node.name
        conn = create_connection(self._db_file)
        with conn:         
            df = get_displacement(conn,
                                  #plane2D=self._plane.plane2D,
                                  node_name=node_name,
                                  result_name=self._result_name,
                                  mesh_name=self._mesh_name)
        #
        if self._plane:
            df.drop(['z', 'rx', 'ry'], axis=1, inplace=True)
        return NodeDeflection(df, units=units)
#
#
# --------------------
# sql operations
#
def get_force(conn, result_name: int|str,
              mesh_name: int|str,
              node_name: int|str|None=None):
    """ """
    query = [result_name, mesh_name]
    table = f'SELECT Load.name, Load.level, Component.name, Node.name,\
              ResultNodeForce.* \
              FROM Load, Node, Mesh, Component, Result, ResultNodeForce \
              WHERE Result.name = ? \
                 AND Mesh.name = ? \
                 AND Load.number = ResultNodeForce.load_id \
                 AND Node.number = ResultNodeForce.node_id \
                 AND Result.number = ResultNodeForce.result_id \
                 AND Result.mesh_id = Mesh.number \
                 AND Mesh.component_id = Component.number'
    
    if node_name:
        table += f'AND Node.name = ?'
        query.extend([node_name])
    table += ';'
    
    cols = ['load_name', 'load_level',
            'component_name', 'node_name',
            'number', 'result_id',
            'load_id', 'node_id',  'system',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    #
    df = pull_data_df(conn, table, tuple(query), cols)
    df['mesh_name'] = mesh_name
    df['result_name'] = result_name    
    #
    cols = ['number', 'component_name',
            'mesh_name', 'result_name', 
            'load_name', 'load_level',
            'node_name', 'system',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    #if plane2D:
    #    cols.extend(['Fx', 'Fy', 'Mz'])
    #else:
    #    cols.extend(['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'])
    df = df.astype({'Fx': 'float64',
                    'Fy': 'float64',
                    'Fz': 'float64',
                    'Mx': 'float64',
                    'My': 'float64',
                    'Mz': 'float64'}).fillna(value=float(0.0))
    
    return df[cols]
#
def get_displacement(conn, result_name: int|str,
                     mesh_name: int|str,
                     node_name: int|str|None=None):
    """ """
    query = [result_name, mesh_name]
    table = f'SELECT Load.name, Load.level, Component.name, Node.name, \
              ResultNodeDisplacement.* \
              FROM Load, Node, Mesh, Component, Result, ResultNodeDisplacement \
              WHERE Result.name = ? \
                 AND Mesh.name = ? \
                 AND Load.number = ResultNodeDisplacement.load_id \
                 AND Node.number = ResultNodeDisplacement.node_id \
                 AND Result.number = ResultNodeDisplacement.result_id \
                 AND Result.mesh_id = Mesh.number \
                 AND Mesh.component_id = Component.number '
    #
    if node_name:
        table += f' AND Node.name = ?'
        query.extend([node_name])
    table += ';'
        
    cols = ['load_name', 'load_level',
            'component_name', 'node_name', 
            'number', 'result_id', 'node_id',
            'load_id', 'system',
            'x', 'y', 'z', 'rx', 'ry', 'rz']
    #
    df = pull_data_df(conn, table, tuple(query), cols)
    df['mesh_name'] = mesh_name
    df['result_name'] = result_name
    #
    cols = ['number', 'component_name',
            'mesh_name', 'result_name', 
            'load_name', 'load_level',
            'node_name', 'system',
            'x', 'y', 'z', 'rx', 'ry', 'rz']
    #
    #if plane2D:
    #    cols.extend(['x', 'y', 'rz'])
    #else:
    #    cols.extend(['x', 'y', 'z', 'rx', 'ry', 'rz'])
    df = df.astype({'x': 'float64',
                    'y': 'float64',
                    'z': 'float64',
                    'rx': 'float64',
                    'ry': 'float64',
                    'rz': 'float64'}).fillna(value=float(0.0))
    
    return df[cols]
#
def get_reactions(conn, result_name: int|str,
                  mesh_name: int|str,
                  #plane2D: bool,
                  node_name: int|str|None=None):
    """ """
    query = [result_name, mesh_name]
    table = f'SELECT Load.name, Load.level, Component.name, Node.name,\
              ResultNodeReaction.* \
              FROM Load, Node, Mesh, Component, Result, ResultNodeReaction \
              WHERE Result.name = ? \
                 AND Mesh.name = ? \
                 AND Load.number = ResultNodeReaction.load_id \
                 AND Node.number = ResultNodeReaction.node_id \
                 AND Result.number = ResultNodeReaction.result_id \
                 AND Result.mesh_id = Mesh.number \
                 AND Mesh.component_id = Component.number'
    
    if node_name:
        table += f'AND Node.name = ?'
        query.extend([node_name])
    table += ';'
    #
    cols = ['load_name', 'load_level',
            'component_name', 'node_name',
            'number', 'result_id', 'load_id', 
            'node_id', 'system',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    #
    df = pull_data_df(conn, table, tuple(query), cols)
    df['mesh_name'] = mesh_name
    df['result_name'] = result_name    
    #
    cols = ['number', 'component_name',
            'mesh_name', 'result_name', 
            'load_name', 'load_level',
            'node_name', 'system',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    #if plane2D:
    #    cols.extend(['Fx', 'Fy', 'Mz'])
    #else:
    #    cols.extend(['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'])
    df = df.astype({'Fx': 'float64',
                    'Fy': 'float64',
                    'Fz': 'float64',
                    'Mx': 'float64',
                    'My': 'float64',
                    'Mz': 'float64'}).fillna(value=float(0.0))    
    
    return df[cols]
#
def pull_data_df(conn, table: str, query: tuple, cols: list):
    """ """
    cur = conn.cursor()
    cur.execute(table, query)
    data = cur.fetchall()
    #
    db = DBframework()
    return db.DataFrame(data=data, columns=cols)
#