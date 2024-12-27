# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
from collections.abc import Mapping
from typing import NamedTuple
import re
#
# package imports
#from steelpy.ufo.mesh.main import MeshPlane
from steelpy.utils.dataframe.main import DBframework
#from steelpy.utils.sqlite.utils import create_connection #, create_table


#
#
class NodeResultBasic(Mapping):
    __slots__ = ['_mesh', '_labels', 'plane']
    
    def __init__(self, mesh) -> None:
        """
        mesh : 
        """
        self._mesh = mesh
        self._labels = list(self._mesh._nodes.keys())
        self.plane = mesh._plane
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
        # TODO: is node force needed? 
        #output += self.force().__str__()
        output += self.reaction().__str__()
        return output
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
        output += print_node_items(self.df)
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
        output += print_node_items(self.df)
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
        output += print_node_items(self.df)
        return output
    #
#
#
#
#
# --------------------
# printing
#
#
def print_node_items(items: DBframework.DataFrame,
                     cols: list|None = None):
    """
    """
    # default heading
    if not cols:
        cols = ['number', 'load_name', 'component_name',
                'load_level', 'system', 'mesh_name', 'result_name']
    #
    output = ""
    items.rename(columns={'node_name': 'node'}, inplace=True)
    node_grp = items.groupby(['load_level'])
    #
    # Basic
    try:
        node_type = node_grp.get_group(('basic', ))
        node_items = node_type.groupby(['load_name', 'component_name', 'system'])
        for key, item in node_items:
            output += "-- Basic Load  Name: {:}  Component: {:} System: {:}\n".format(*key)
            #
            vals = item.drop(cols, axis=1)
            header2 = list(vals)
            header2 = get_gap(header2)
            #
            output += header2
            output += "\n"
            output += printout(vals)
            output += "\n"
    except KeyError:
        pass    
    #
    # Combination
    try:
        node_type = node_grp.get_group(('combination', ))
        node_items = node_type.groupby(['load_name', 'component_name', 'system'])
        output += '{:}\n'.format(52 * '-')
        for key, item in node_items:
            output += "-- Load Combination  Name: {:} Component: {:} System: {:}\n".format(*key)
            #
            vals = item.drop(cols, axis=1)
            header2 = list(vals)
            header2 = get_gap(header2)
            #
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
        output += "{:9d} ".format(item.Index)
        for val in item[1:]:
            output += "{: 1.3e} ".format(val)
        output += "\n"
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