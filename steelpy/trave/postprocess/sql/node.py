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

from steelpy.trave.postprocess.utils.node import (NodeResBasic, NodeForce,
                                                  NodeDeflection, NodeReaction)

#
class NodeResSQL(NodeResBasic):
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
    def __getitem__(self, node_name: int|str) -> NodeResItem|IndexError:
        """
        node_name : node number
        """
        try:
            self._labels.index(node_name)
            return NodeResItem(node=self._mesh._nodes[node_name],
                               plane=self.plane,
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
            df = get_force(conn, plane2D=self.plane.plane2D)
        return NodeForce(df, units=units)
    #
    #@property
    def displacement(self, units:str='si')->NodeDeflection:
        """ node displacement"""
        conn = create_connection(self.db_file)
        with conn:         
            dispdf = get_displacement(conn,
                                      plane2D=self.plane.plane2D)
        #
        #if self.plane.plane2D:
        #    header = ['z',  'rx',  'ry']
        #    dispdf.drop(columns=header, axis=1, inplace=True)
        return NodeDeflection(dispdf, units=units)
    #
    #@property
    def reaction(self, units:str='si')->NodeReaction:
        """ Node reactions"""
        conn = create_connection(self.db_file)
        with conn:        
            df = get_reactions(conn, plane2D=self.plane.plane2D)
        #
        #if self.plane.plane2D:
        #    header = ['Fz',  'Mx',  'My']
        #    df.drop(header, axis=1, inplace=True)
        return NodeReaction(df, units=units)
#        
#
#
@dataclass
class NodeResItem:
    __slots__ = ['_node', '_db_file', '_plane']
    
    def __init__(self, node, plane, db_file: str)-> None:
        """ """
        self._node = node
        self._plane = plane
        self._db_file = db_file
    #
    def force(self, units:str='si')-> NodeForce:
        """node force"""
        node_name = self._node.name
        conn = create_connection(self._db_file)
        with conn:        
            df = get_force(conn, plane2D=self._plane.plane2D,
                           node_name=node_name)
        #
        #if self._plane.plane2D:
        #    header = ['Fz',  'Mx',  'My']
        #    df.drop(header, axis=1, inplace=True)
        return NodeForce(df, units=units)
    #
    def displacement(self, units:str='si')-> NodeDeflection:
        """ node displacement"""
        node_name = self._node.name
        conn = create_connection(self._db_file)
        with conn:         
            df = get_displacement(conn, plane2D=self._plane.plane2D,
                                  node_name=node_name)
        #
        #if self._plane.plane2D:
        #    header = ['z',  'rx',  'ry']
        #    df.drop(header, axis=1, inplace=True)
        return NodeDeflection(df, units=units)    
#
#
# --------------------
# sql operations
#
def get_force(conn, plane2D: bool,
              node_name: int|str|None=None):
    """ """
    query = f'SELECT Load.name, Load.level, Component.name, Node.name,\
              ResultNodeForce.* \
              FROM Load, Node, Mesh, Component, Result, ResultNodeForce \
              WHERE Load.number = ResultNodeForce.load_id \
                 AND Node.number = ResultNodeForce.node_id \
                 AND Result.number = ResultNodeForce.result_id \
                 AND Result.mesh_id = Mesh.number \
                 AND Mesh.component_id = Component.number'
    
    if node_name:
        query += f'AND node_name = {node_name}'
    query += ';'
    
    cols = ['load_name', 'load_level',
            'component_name', 'node_name',
            'number', 'result_id',
            'load_id', 'node_id',  'system',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    #
    fpdf = pull_data_df(conn, query, cols)
    #
    cols = ['number', 'component_name',
            'load_name', 'load_level',
            'node_name', 'system']
    if plane2D:
        cols.extend(['Fx', 'Fy', 'Mz'])
    else:
        cols.extend(['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'])
    
    return fpdf[cols]
#
def get_displacement(conn, plane2D: bool,
                     node_name: int|str|None=None):
    """ """
    query = f'SELECT Load.name, Load.level, Component.name, Node.name, \
              ResultNodeDisplacement.* \
              FROM Load, Node, Mesh, Component, Result, ResultNodeDisplacement \
              WHERE Load.number = ResultNodeDisplacement.load_id \
                 AND Node.number = ResultNodeDisplacement.node_id \
                 AND Result.number = ResultNodeDisplacement.result_id \
                 AND Result.mesh_id = Mesh.number \
                 AND Mesh.component_id = Component.number'
    #
    if node_name:
        query += f' AND Node.name = {node_name}'
    query += ';'
        
    cols = ['load_name', 'load_level', 'component_name', 'node_name', 
            'number', 'result_id', 'node_id',
            'load_id', 'system',
            'x', 'y', 'z', 'rx', 'ry', 'rz']
    #
    dispdf = pull_data_df(conn, query, cols)
    #
    cols = ['number', 'component_name',
            'load_name', 'load_level',
            'node_name', 'system']
    
    if plane2D:
        cols.extend(['x', 'y', 'rz'])
    else:
        cols.extend(['x', 'y', 'z', 'rx', 'ry', 'rz'])
    
    return dispdf[cols]
#
def get_reactions(conn, plane2D: bool,
                  node_name: int|str|None=None):
    """ """
    query = f'SELECT Load.name, Load.level, Component.name, Node.name,\
              ResultNodeReaction.* \
              FROM Load, Node, Mesh, Component, Result, ResultNodeReaction \
              WHERE Load.number = ResultNodeReaction.load_id \
                 AND Node.number = ResultNodeReaction.node_id \
                 AND Result.number = ResultNodeReaction.result_id \
                 AND Result.mesh_id = Mesh.number \
                 AND Mesh.component_id = Component.number'
    
    if node_name:
        query += f'AND Node.name = {node_name}'
    query += ';'
    #
    cols = ['load_name', 'load_level', 'component_name', 'node_name',
            'number', 'result_id', 'load_id', 
            'node_id', 'system',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    #
    recdf = pull_data_df(conn, query, cols)
    #
    cols = ['number', 'component_name',
            'load_name', 'load_level',
            'node_name', 'system']
    if plane2D:
        cols.extend(['Fx', 'Fy', 'Mz'])
    else:
        cols.extend(['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'])
    
    return recdf[cols]
#
def pull_data_df(conn, query: str, cols: list):
    """ """
    cur = conn.cursor()
    cur.execute(query)
    data = cur.fetchall()
    #
    db = DBframework()
    return db.DataFrame(data=data, columns=cols)
#