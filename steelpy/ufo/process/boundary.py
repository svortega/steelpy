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
def find_boundary_type(word_in:str) -> str:
    """
    """
    key = {"support": r"\b(node(s)?|support(s)?|constrain(s)?)\b",
           "beam": r"\b(beam(s)?)\b",}
    match = common.find_keyword(word_in, key)
    return match
#
def find_support_item(word_in:str) -> str:
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
def get_support_df(df: DBframework.DataFrame):
    """ [] """
    columns = list(df.columns)
    header = {item:find_support_item(item) for item in columns}
    df.rename(columns=header, inplace=True)
    #
    columns = list(df.columns)
    if 'restrain' in columns:
        df['restrain'] = df['restrain'].apply(get_node_boundary).tolist()
        #fixity = {'x': [], 'y': [], 'z': [], 'rx': [], 'ry': [], 'rz': []}
        #for item in boundary:
        #    fixity['x'].append(item[0])
        #    fixity['y'].append(item[1])
        #    fixity['z'].append(item[2])
        #    fixity['rx'].append(item[3])
        #    fixity['ry'].append(item[4])
        #    fixity['rz'].append(item[5])
        #
        #db = DBframework()
        #fixity = db.DataFrame(fixity)
        #df = db.concat([df, fixity], axis=1)
        #df.drop('restrain')
    else: # here supose to be x,y,z,rx,ry,rz
        1 / 0
    #
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
