# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from array import array
from collections.abc import Mapping
import re
from typing import NamedTuple

#from steelpy.ufo.process.elements.nodes import (check_point_list,
#                                                check_point_dic)
#
# -----------------------
# TODO: merge with slite
class BoundaryItem(NamedTuple):
    """
    """
    x: float
    y: float
    z: float
    rx: float
    ry: float
    rz: float
    number: int
    name: int|str
    node:int
    
    def __str__(self) -> str:
        #if (name := self.name) == None:
        #    name = ""
        return "{:>12s} {:12d} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f}\n"\
            .format(str(self.name), self.node, self.x, self.y, self.z, self.rx, self.ry, self.rz)
#
#
class BoundaryNode(Mapping):
    __slots__ = ['_component', '_labels']
    
    def __init__(self, component: int):
        """
        """
        #self._labels: list[int|str] = []
        self._component = component
    #
    #
    def __len__(self) -> float:
        return len(self._labels)

    def __iter__(self):
        """
        """
        return iter(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    # ----------------------------
    # Operations
    # ----------------------------
    # TODO : why two ?
    #def get_boundary(self, name:str):
    #    """ """
        #if re.match(r"\b(fix(ed)?|encastre)\b", name, re.IGNORECASE):
        #    #self._title.append('fixed')
        #    value = [1,1,1,1,1,1]
        #elif re.match(r"\b(pinn(ed)?|roll(er)?)\b", name, re.IGNORECASE):
        #    #self._title.append('pinned')
        #    value = [1,1,1,1,0,0]
        #
        #elif re.match(r"\b(guide(d)?|roll(ed)?)\b", name, re.IGNORECASE):
        #    return [0,1,1,1,0,0]
        #
        #elif re.match(r"\b(free)\b", name, re.IGNORECASE):
        #    value = [0,0,0,0,0,0]
        #    
        #else:
        #    raise IOError("boundary type {:} not implemented".format(name))
        #
        #value = get_boundary_by_name(bname=name)
        #
        #return value
    #
    def _get_fixity(self, fixity):
        """ """
        #if isinstance(fixity, str):
        #    return get_boundary_by_name(bname=fixity)
            #if re.match(r"\b(fix(ed)?)\b", fixity, re.IGNORECASE):
            #    return [1,1,1,1,1,1]
            #
            #elif re.match(r"\b(pinn(ed)?)\b", fixity, re.IGNORECASE):
            #    return [1,1,1,1,0,0]
            #
            #elif re.match(r"\b(roll(er)?)\b", fixity, re.IGNORECASE):
            #    return [0,1,1,1,0,0]
            #
            #elif re.match(r"\b(guide(d)?)\b", fixity, re.IGNORECASE):
            #    return [1,0,0,1,0,0]
            #
            #elif re.match(r"\b(free)\b", fixity, re.IGNORECASE):
            #    return [0,0,0,0,0,0]
            #
            #else:
            #    raise IOError("boundary type {:} not implemented".format(fixity))
        
        #elif isinstance(fixity, (list, tuple)):
        #    if isinstance(fixity[0], (list, tuple)):
        #        fixity = fixity[0]
        #    return fixity
        #
        #elif isinstance(fixity, dict):
        #    return [fixity['x'], fixity['y'], fixity['z'], 
        #            fixity['rx'], fixity['ry'], fixity['rz']]
        #
        #else:
        #    raise Exception('   *** Boundary input format not recognized')
        return get_node_boundary(fixity)
    #
    #def _get_coordinates(self, coordinates):
    #    """ """
    #    if isinstance(coordinates, (list, tuple)):
    #        coordinates = check_point_list(coordinates, steps=3)
    #    elif isinstance(coordinates, dict):
    #        coordinates = check_point_dic(coordinates)
    #    else:
    #        raise Exception('Node input format not valid')
    #    return coordinates
#
#
#
def get_node_boundary(fixity:list|tuple|dict|str):
    """ """
    if isinstance(fixity, str):
        return get_boundary_by_name(bname=fixity)
    
    elif isinstance(fixity, (list, tuple)):
        if isinstance(fixity[0], (list, tuple)):
            fixity = fixity[0]
        return fixity
    
    elif isinstance(fixity, dict):
        return [fixity['x'], fixity['y'], fixity['z'], 
                fixity['rx'], fixity['ry'], fixity['rz']]
    
    else:
        raise Exception('   *** Boundary input format not recognized')    
#
def get_nboundary_dict(values: dict):
    """
    return [node, x, y, z, rx, ry, rz]
    """
    nodes = []
    fixity = [0, 0, 0, 0, 0, 0]
    for key, item in values.items():
        if re.match(r"\b(node)\b", key, re.IGNORECASE):
            if isinstance(item, list):
                1 / 0
                #nodes = {x:[] for x in item}
                nodes.extend(item)
            else:
                #nodes[item] = []
                nodes.append(item)
        elif re.match(r"\b(support)\b", key, re.IGNORECASE):
            fixity =  get_node_boundary(item)
    #
    fixity = [[item, *fixity] for item in nodes]
    return fixity
#
def get_boundary_by_name(bname: str):
    """ """
    if re.match(r"\b(fix(ed)?)\b", bname, re.IGNORECASE):
        return [1,1,1,1,1,1]
    
    elif re.match(r"\b(pinn(ed)?)\b", bname, re.IGNORECASE):
        return [1,1,1,1,0,0]
    
    elif re.match(r"\b(roll(er)?)\b", bname, re.IGNORECASE):
        return [0,1,1,1,0,0]

    elif re.match(r"\b(guide(d)?)\b", bname, re.IGNORECASE):
        return [1,0,0,1,0,0]
    
    elif re.match(r"\b(free)\b", bname, re.IGNORECASE):
        return [0,0,0,0,0,0]
        
    else:
        raise IOError(f"boundary type {bname} not implemented")
    #
    #return value
#
