# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from array import array
from collections.abc import Mapping
import re
from typing import NamedTuple

# package imports
import steelpy.utils.io_module.text as common
from steelpy.utils.dataframe.main import DBframework
#from steelpy.ufo.utils.node import find_coord

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
    __slots__ = ['_name']
    
    def __init__(self, name: int|str):
        """
        """
        self._name = name
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
    #
    def _get_fixity(self, fixity):
        """ """
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
    # ----------------------------
    # Operations
    # ----------------------------
    #
    #
    def get_number(self, start:int=1):
        """
        """
        try:
            n = max(self._number) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
#
#
#
def find_boundary_type(word_in:str) -> str:
    """
    """
    key = {"support": r"\b(node(s)?|support(s)?|constrain(s)?)\b",
           "beam": r"\b(beam(s)?)\b",}
    match = common.find_keyword(word_in, key)
    return match
#
def find_support_node(word_in:str) -> str:
    """
    """
    key = {'name': r'\b((boundar(y|ies)(_|-|\s*)?)?(id|name))\b',
           "node": r"\b(node(s)?|support(s)?)\b",
           'type': r'\b(type)\b',
           'restrain': r'\b(restrain|constrain(s)?|fixit(y|ies))\b',
           'title': r'\b(title)\b'}
    match = common.find_keyword(word_in, key)
    return match
#
def find_support_coord(word_in:str) -> str:
    """
    """
    key = {'name': r'\b((boundar(y|ies)(_|-|\s*)?)?(id|name))\b',
           'type': r'\b(type)\b',
           'restrain': r'\b(restrain|constrain(s)?|fixit(y|ies))\b',
           'title': r'\b(title)\b',
           }
    try:
        match = common.find_keyword(word_in, key)
        return match
    except IOError:
        return find_coord(word_in)
#
def find_support_boundary(word_in:str) -> str:
    """
    """
    key = {'x': r'\b((i)?(_|-|\s*)?x)\b',
           'y': r'\b((i)?(_|-|\s*)?y)\b',
           'z': r'\b((i)?(_|-|\s*)?z)\b',
           'rx': r'\b(r(_|-|\s*)?x)\b',
           'ry': r'\b(r(_|-|\s*)?y)\b',
           'rz': r'\b(r(_|-|\s*)?z)\b'}
    match = common.find_keyword(word_in, key)
    return match
#
#
def find_coord(word_in: str) -> str:
    """ """
    key: dict = {
                  "x": r"\b((coord(inate)?(s)?|elev(ation)?)?(_|-|\s*)?x)\b",
                  "y": r"\b((coord(inate)?(s)?|elev(ation)?)?(_|-|\s*)?y)\b",
                  "z": r"\b((coord(inate)?(s)?|elev(ation)?)?(_|-|\s*)?z)\b",
                  #
                  "r": r"\b((coord(inate)?(s)?|elev(ation)?)?(_|-|\s*)?r)\b",
                  "theta": r"\b((coord(inate)?(s)?|elev(ation)?)?(_|-|\s*)?theta)\b",
                  "phi": r"\b((coord(inate)?(s)?|elev(ation)?)?(_|-|\s*)?phi)\b",                
                  }

    match = common.find_keyword(word_in, key)
    return match
#
#
def get_support_df(df: DBframework.DataFrame):
    """ """
    columns = list(df.columns)
    try:
        header = {item:find_support_node(item) for item in columns}
        df.rename(columns=header, inplace=True)
        #
        #columns = list(df.columns)
        #if 'restrain' in columns:
        #    df['restrain'] = df['restrain'].apply(get_node_boundary).tolist()
        #else:
        #    raise IOError (f'rastrain missing')
    except IOError:
        header = {item:find_support_coord(item) for item in columns}
        df.rename(columns=header, inplace=True)
        #1 / 0
    #
    columns = list(df.columns)
    if not 'restrain' in columns:
        raise IOError (f'rastrain missing')
    return df
#
def get_node_boundary(fixity:list|tuple|dict|str):
    """ [x,y,z,rx,ry,rz] """
    if isinstance(fixity, str):
        return get_boundary_by_name(bname=fixity)
    
    elif isinstance(fixity, (list, tuple)):
        if isinstance(fixity[0], (list, tuple)):
            fixity = fixity[0]
        return get_fixity(fixity)
    
    elif isinstance(fixity, dict):
        fixity = {find_support_boundary(key): item
                  for key, item in fixity.items()}
        return get_fixity([fixity['x'], fixity['y'], fixity['z'], 
                           fixity['rx'], fixity['ry'], fixity['rz']])
    
    else:
        raise Exception('   *** Boundary input format not recognized')    
#
def get_fixity(data: list) -> list:
    """ [x,y,z,rx,ry,rz]"""
    fixity = [0, 0, 0, 0, 0, 0]
    valid = [0, 1]
    #
    for x in range(len(fixity)):
        try:
            if isinstance(data[x], (int, float)):
                if not int(data[x]) in valid:
                    raise IOError(f'fixity code {data[x]} not valid')
                fixity[x] = int(data[x])
            elif isinstance(data[x], str):
                if re.match(r"\b(fix(ed)?)\b", data[x], re.IGNORECASE):
                    fixity[x] = 1
                elif re.match(r"\b(free)\b", data[x], re.IGNORECASE):
                    fixity[x] = 0
                else:
                    raise IOError(f"fixity code {data[x]} not valid")            
            else:
                raise IOError('fixity format not valid')
        except IndexError:
            continue
    #
    return fixity
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
