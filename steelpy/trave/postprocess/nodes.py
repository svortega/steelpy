# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Mapping
from typing import NamedTuple
import re
#
# package imports
from steelpy.ufo.mesh.main import MeshPlane
from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.sqlite.utils import create_connection #, create_table


#
#
class NodeResBasic(Mapping):
    __slots__ = ['_mesh', '_labels', 'plane']
    
    def __init__(self, mesh) -> None:
        """
        mesh : 
        """
        self._mesh = mesh
        self._labels = list(self._mesh._nodes.keys())
        self.plane = mesh._plane # self._plane()
    #    
    #
    def __len__(self) -> int:
        return len(self._labels)

    def __iter__(self):
        """
        """
        nodes = sorted(self._labels)
        return iter(nodes)

    def __contains__(self, value) -> bool:
        return value in self._labels
    
    def __str__(self) -> str:
        """ """
        #lenght = ' m'
        #space = " "
        output = "\n"
        output += self.displacement().__str__()
        output += self.force().__str__()
        output += self.reaction().__str__()
        return output   
    #
#
#
class Nodes(NodeResBasic):
    __slots__ = ['_labels', '_mesh', 'plane', 'db_file']
    
    def __init__(self, mesh, db_file: str):
        """ """
        self.db_file = db_file
        super().__init__(mesh)
    #
    #
    # -----------------------------
    #
    #@property
    def _plane(self):
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
    #
    def __getitem__(self, node_name: int|str) -> tuple:
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
    def force(self, units:str='si'):
        """ node force"""
        conn = create_connection(self.db_file)
        with conn:        
            df = get_force(conn)
        #
        if self.plane.plane2D:
            header = ['Fz',  'Mx',  'My']
            df.drop(header, axis=1, inplace=True)
        return NodeForce(df, units=units)
    #
    #@property
    def displacement(self, units:str='si'):
        """ node displacement"""
        conn = create_connection(self.db_file)
        with conn:         
            df = get_displacement(conn)
        #
        if self.plane.plane2D:
            header = ['z',  'rx',  'ry']
            df.drop(header, axis=1, inplace=True)
            #df.rename(columns={'node_name': 'node'}, inplace=True)
        return NodeDeflection(df, units=units)
    #
    #@property
    def reaction(self, units:str='si'):
        """ Node reactions"""
        conn = create_connection(self.db_file)
        with conn:        
            df = get_reactions(conn)
        #
        if self.plane.plane2D:
            header = ['Fz',  'Mx',  'My']
            df.drop(header, axis=1, inplace=True)
        return NodeReaction(df, units=units)
    #        
    #
#
#
# --------------------
# Node ops
#
class NodeDeflection(NamedTuple):
    """ """
    df: DBframework
    units: str
    #
    def __str__(self):
        """ """
        unitsout = "Units : SI [m/radians]"
        if re.match(r"\b(us|imperial)\b", self.units, re.IGNORECASE):
            unitsout = "units : US [ft/radians]"
        #
        output = "\n"
        output += "{:}\n".format(52 * '-')
        output += "** Node Displacements | "
        output += unitsout
        output += "\n"
        output += "{:}\n".format(52 * '-')
        output += print_nitems(self.df)
        return output
#
class NodeForce(NamedTuple):
    """ """
    df: DBframework
    units: str
    #
    def __str__(self):
        """ """
        unitsout = "Units : SI [N/N-m]"
        if re.match(r"\b(us|imperial)\b", self.units, re.IGNORECASE):
            unitsout = "units : US [lb/lb-ft]"
        #
        output = "\n"
        output += "{:}\n".format(52 * '-')
        output += "** Node Force | "
        output += unitsout
        output += "\n"        
        output += "{:}\n".format(52 * '-')
        output += print_nitems(self.df)
        return output
#
class NodeReaction(NamedTuple):
    """ """
    df: DBframework
    units: str
    #
    def __str__(self):
        """ """
        unitsout = "Units : SI [N/N-m]"
        if re.match(r"\b(us|imperial)\b", self.units, re.IGNORECASE):
            unitsout = "units : US [lb/lb-ft]"
        #        
        output = "\n"
        output += "{:}\n".format(52 * '-')
        output += "** Node Reactions | "
        output += unitsout
        output += "\n"         
        output += "{:}\n".format(52 * '-')
        output += print_nitems(self.df)
        return output
    #
#
#
@dataclass
class NodeResItem:
    __slots__ = ['_node', '_db_file', '_plane']
    
    def __init__(self, node, plane, db_file: str):
        """ """
        self._node = node
        self._plane = plane
        self._db_file = db_file
    #
    def force(self, units:str='si'):
        """node force"""
        node_name = self._node.name
        conn = create_connection(self._db_file)
        with conn:        
            df = get_force(conn, node_name)
        #
        if self._plane.plane2D:
            header = ['Fz',  'Mx',  'My']
            df.drop(header, axis=1, inplace=True)        
        return NodeForce(df, units=units)
    #
    def displacement(self, units:str='si'):
        """ node displacement"""
        node_name = self._node.name
        conn = create_connection(self._db_file)
        with conn:         
            df = get_displacement(conn, node_name)
        #
        if self._plane.plane2D:
            header = ['z',  'rx',  'ry']
            df.drop(header, axis=1, inplace=True)
        return NodeDeflection(df, units=units)    
#
#
# --------------------
# sql operations
#
def get_force(conn, node_name: int|str|None=None,
              item:str='*'):
    """ """
    if node_name:
        #project = (node_name,)
        query = f'SELECT {item} FROM ResultNodeForce \
                  WHERE node_name = {node_name}'
    else:
        query = f'SELECT {item} FROM ResultNodeForce'
    
    cols = ['number', 'load_name', 'component_name',
            'load_level', 'load_system', 
            'node_name', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    
    return pull_data(conn, query, cols)
#
#
def get_displacement(conn, node_name: int|str|None=None,
                     item:str='*'):
    """ """
    if node_name:
        #project = (node_name,)
        query = f'SELECT {item} FROM ResultNodeU \
                  WHERE node_name = {node_name}'
    else:
        query = f'SELECT {item} FROM ResultNodeU'
        
    cols = ['number', 'load_name', 'component_name',
            'load_level', 'load_system', 
            'node_name', 'x', 'y', 'z', 'rx', 'ry', 'rz']
    
    return pull_data(conn, query, cols)
#
#
def get_reactions(conn, node_name: int|str|None=None,
                  item:str='*'):
    """ """
    if node_name:
        #project = (node_name,)
        query = f'SELECT {item} FROM ResultNodeReaction \
                  WHERE node_name = {node_name}'
    else:
        query = f'SELECT {item} FROM ResultNodeReaction'
    
    cols = ['number', 'load_name', 'component_name',
            'load_level', 'load_system', 
            'node_name', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    
    return pull_data(conn, query, cols)
#
#
def pull_data(conn, query: str, cols: list):
    """ """
    cur = conn.cursor()
    cur.execute(query)
    data = cur.fetchall()
    #
    db = DBframework()
    return db.DataFrame(data=data, columns=cols)
#
#
# --------------------
# printing
#
#
def print_nitems(items,
                 cols: list = ['number', 'load_name', 'component_name',
                               'load_level', 'load_system']):
    """
    """
    #
    items.rename(columns={'node_name': 'node'}, inplace=True)
    #
    ndgrp = items.groupby('load_level')
    ndtype = ndgrp.get_group('basic')
    nditems = ndtype.groupby(['load_name', 'component_name', 'load_system'])
    #
    output = ""
    # Basic
    for key, wk in nditems:
        #header1 =  wk[['load_name','component_name', 'load_system']].values
        #
        output += "-- Basic Load  Name: {:}  Component: {:} System: {:}\n".format(*key)
        #
        vals = wk.drop(cols, axis=1)
        header2 = list(vals)
        header2 = get_gap(header2)
        #
        output += header2
        output += "\n"
        output += printout(vals)
        output += "\n"
    #
    # Combination
    try:
        ndtype = ndgrp.get_group('combination')
        nditems = ndtype.groupby(['load_name', 'component_name', 'load_system'])
        #
        output += '{:}\n'.format(52 * '-')
        for key, wk in nditems:
            #header1 =  wk[['load_name','component_name', 'load_system']].values
            output += "-- Load Combination  Name: {:} Component: {:} System: {:}\n".format(*key)
            #output += header2
            #
            vals = wk.drop(cols, axis=1)
            header2 = list(vals)
            #header2 = '       '.join(header2)
            header2 = get_gap(header2)
            #
            #ndsum =  wk.groupby('node_name')[plane.hdisp].sum()
            #output += printout(ndsum.index, ndsum.values)
            output += header2
            output += "\n"
            output += printout(vals)
            output += "\n"
    except KeyError:
        pass
    #
    return output
#
#
def printout(disp:list):
    """
    Report the relative displacement of the structure
    """
    output = ""
    #node_id = disp.groupby('node')
    disp.set_index('node', inplace=True)
    for item in disp.itertuples():
        output += '{:9d} {: 1.3e} {: 1.3e} {: 1.3e}\n'.format(*item)
    return output
#
#
def get_gap(header, step:int=9):
    """ """
    new = []
    for item in header:
        gap = step - len(item) 
        new.extend([item, " " * gap])
    #
    new = " ".join(new)
    return new
#
#
# --------------------
# operations
#
def get_max_displacement(dispp):
    """ """
    columns = list(zip(*dispp))
    #
    maxval = []
    nodeitem = []
    #
    output = "\n"
    output += "Maximum displacements\n"
    for column in columns:
        nmax = [max(column), min(column)]
        maxval.append(nmax[0] if nmax[0] > abs(nmax[1]) else nmax[1])
        nodeitem.append(column.index(maxval[-1]))
    #
    output += "node {:>10}/n".format("")
    for _node in nodeitem:
        output += "{:}{:>10}/n".format(_node, "")
    #
    output = "\n"
    output += "value \n"
    for _value in maxval:
        output += "{: 1.3e} \n".format(_value)
    #
    return output
#
#