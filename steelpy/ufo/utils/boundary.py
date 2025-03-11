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
    boundary_type: str
    def __str__(self) -> str:
        #
        #TODO: include boundary type
        fixity = [self.x, self.y, self.z, self.rx, self.ry, self.rz]
        fixity = [0 if item == None else 1 for item in fixity]
        #
        return "{:>12s} {:12d} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f}\n"\
            .format(str(self.name), self.node, *fixity)
            #.format(str(self.name), self.node, self.x, self.y, self.z, self.rx, self.ry, self.rz)
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
def get_boundary_item(value: list | tuple | dict)->str:
    """
    Returns:
        [force, fx, fy, fz, mx, my, mz, comment]
        [displacement, x, y, z, rx, ry, rz, comment]
        [mass,x, y, z, mx, my, mz, comment]
    """
    #1 / 0
    if isinstance(value, (list, tuple)):
        btype = find_boundary_item(value[0])
    elif isinstance(value, dict):
        columns = list(value.keys())
        columns = [get_support_basic(item) for item in columns]
        if 'beam' in columns:
            btype = 'beam'
        elif 'node' in columns:
            btype = 'node'
        else:
            raise Exception('*** Boundary input format not recognized')
    else:
        raise Exception('*** Boundary input format not recognized')
    return btype   
#
def find_boundary_item(word_in:str) -> str:
    """
    """
    key = {"node": r"\b(node(s)?|support(s)?)\b",
           "beam": r"\b(beam(s)?)\b",}
    match = common.find_keyword(word_in, key)
    return match
#
def get_boundary_beam(values: list|tuple|dict):
    """
    Returns:
        [beam_name, node_end, node_name, values]
    """
    output = [None, None, None]
    if isinstance(values, (list, tuple)):
        steps = len(values)
        idx = []
        for x in range(0, steps, 2):
            item = get_support_basic(values[x])
            try:
                #print(item, values[x + 1])
                if item in ['beam']:
                    output[0] = values[x + 1]
                    idx.extend([x, x + 1])
                elif item in ['node_end']:
                    output[1] = values[x + 1]
                    idx.extend([x, x + 1])
                elif item in ['node']:
                    output[2] = values[x + 1]
                    idx.extend([x, x + 1])
            #try:
            #    # fixity = find_boundary_node([item])
            #    print(item, values[x + 1])
            except IndexError:
                #print(item)
                break
            #    fixity = find_boundary_node([item])
            #    print(item, fixity)
        #
        for x in reversed(idx):
            values.pop(x)

    elif isinstance(values, dict):
        1/0
    else:
        raise Exception('*** Boundary input format not recognized')
    #
    output.extend(values)
    return output
#
def get_boundary_node(values: list|tuple|dict, beams):
    """
    Returns:
        [node_name, force, [fx, fy, fz, mx, my, mz,] [transformation], comment]
        [node_name, displacement, [x, y, z, rx, ry, rz], [transformation], comment]
        [node_name, spring, [Kx, Ky, Kz, Krx, Kry, Krz],  [transformation], comment]
    """
    #
    #output = [None] * 10
    #transformation = None
    #comment = None
    #
    if isinstance(values, (list, tuple)):
        output = get_list_data(values, beams)
    elif isinstance(values, dict):
        output = get_dict_data(values, beams)
    else:
        raise Exception('*** Boundary input format not recognized')
    #
    if not output[0]: # node id
        raise IOError(f'node name missing')
    #
    #output = [node_name, btype, fixity, transformation, comment]
    return output
#
def get_list_data(values: list|tuple, beams) -> list:
    """ """
    transformation = None
    comment = None    
    for x, value in enumerate(values):
        try:
            item = get_support_basic(value)
        except IOError:
            continue
        #
        match item:
            case 'node':
                node_name = values[x+1]
            case 'fixed' | 'pinned' | 'guide' | 'free':
                boundary = find_boundary_node([item])
                btype = boundary[0]
                fixity = boundary[1:]
            case 'constrained':
                fixity = get_node_constrain_list(values[x+1])
                btype = item                   
            case 'displacement':
                fixity = get_node_disp_list(values[x+1])
                btype = item                  
            case 'beam':
                beam_name = values[x+1]
                beam = beams[beam_name]
                transformation = beam.unit_vector
                nodes = beam.connectivity
            case 'node_end':
                node_end = values[x+1]
                if node_end in [0,1]:
                    node_name = nodes[0] #.name
                else:
                    node_name = nodes[1] #.name
                #1 / 0
            case 'transformation':
                transformation = values[x+1]
            case 'angle':
                #comment =
                1 / 0
            
            case _:
                raise IOError(f' {item} not valid')
    #
    return [node_name, btype, fixity, transformation, comment]
#
def get_support_basic(word_in:str):
    """ """
    ''
    try:
        return find_support_basic(word_in)
    except (IOError, TypeError):
        pass   
    # dafault boundary
    try:
        return find_boundary_by_name(word_in)
    except (IOError, TypeError):
        pass
    # user boundary
    try:
        return find_boundary_by_type(word_in)
    except (IOError, TypeError):
        pass
    # not found
    raise IOError(f'*** Boundary {word_in} not recognized') 
#
def find_support_basic(word_in:str) -> str:
    """
    """
    key = {'name': r'\b((boundar(y|ies)(_|-|\s*)?)?(id|name))\b',
           "beam": r"\b(beam(s)?)\b",
           "node": r"\b(node(s)?|support(s)?)\b",
           "node_end": r"\b((node(s)?|support(s)?)(_|-|\s*)?end(s)?)\b",
           'type': r'\b(type)\b',
           'values': r'\b(value(s)?|fixit(y|ies))\b',
           'transformation': r'\b(transformation|direction(_|-|\s*)?cosine(s)?|unit(_|-|\s*)?vector(s)?|translation)\b',
           'angle': r'\b((rotation)?(_|-|\s*)?angle)\b',
           'comment': r'\b(title|comment)\b'}
    match = common.find_keyword(word_in, key)
    return match
#
def find_support_coord(word_in:str) -> str:
    """
    """
    key = {'name': r'\b((boundar(y|ies)(_|-|\s*)?)?(id|name))\b',
           'type': r'\b(type)\b',
           'values': r'\b(value(s)?|fixit(y|ies))\b',
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
def get_dict_data(values: list, beams) -> list:
    """ """
    transformation = None
    comment = None
    for key, value in values.items():
        try:
            item = get_support_basic(key)
        except IOError:
            continue
        #
        #
        match item:
            case 'node':
                node_name = value
                
            case 'fixed' | 'pinned' | 'guide' | 'free':
                boundary = find_boundary_node([item])
                btype = boundary[0]
                fixity = boundary[1:]
                
            case 'type':
                try:
                    btype = find_boundary_by_type(value)
                except IOError:
                    boundary = find_boundary_node([value])
                    btype = boundary[0]
                    fixity = boundary[1:]                    
                
            case 'values':
                if btype in ['constrained']:
                    fixity = get_node_constrain_list(value)
                elif btype in ['displacement']:
                    fixity = get_node_disp_list(value)
                else:
                    raise IOError(f'boundary type {btype} not valid')
            #
            case 'constrained':
                fixity = get_node_constrain_list(value)
                btype = item
            
            case 'displacement':
                fixity = get_node_disp_list(value)
                btype = 'constrained'
                1 / 0
            #
            case 'beam':
                beam_name = value
                beam = beams[beam_name]
                transformation = beam.unit_vector
                nodes = beam.connectivity
            case 'node_end':
                node_end = value
                if node_end in [0,1]:
                    node_name = nodes[0]
                else:
                    node_name = nodes[1]
                #1 / 0
            case 'transformation':
                transformation = value
            case 'angle':
                #comment =
                1 / 0
            #
            case 'comment':
                comment = value
            
            case _:
                pass
    #
    return [node_name, btype, fixity, transformation, comment]    
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
        header = {item:find_support_basic(item) for item in columns}
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
    if not 'type' in columns:
        raise IOError (f'boundary type missing')
    return df
#
# ---------------------------------
#
def find_boundary_by_name(word_in:str) -> str:
    """
    """
    key = {"fixed": r"\b(fix(ed)?|encastre(r)?)\b",
           "pinned": r"\b(pinned|hinged|(on(_|-|\s*)?)?roller(s)?)\b",
           "guide": r"\b(guide)\b",
           "free": r"\b(free)\b",
           "roller": r"\b(roll(er)?)\b"}
    match = common.find_keyword(word_in, key)
    return match
#
def find_boundary_by_type(word_in:str) -> str:
    """
    """
    key = {"constrained": r"\b(constrain(ed|t)?|restrain(ed|t)?|homogeneous)\b",
           "spring": r"\b(k|spring|partially(_|-|\s*)?fixed)\b",
           "displacement": r"\b((prescribed)?(_|-|\s*)?(disp(lacement)?|translation|non(_|-|\s*)?homogeneous))\b"}
    match = common.find_keyword(word_in, key)
    return match
#
def find_boundary_type(word_in:str) -> str:
    """
    """
    try:
        return find_boundary_by_name(word_in)
    except IOError:
        return find_boundary_by_type(word_in)
#
#
#
def find_boundary_node(fixity:list|tuple|dict|str)->list:
    """ [type, x,y,z,rx,ry,rz] """
    if isinstance(fixity, (list, tuple)):
        btype = find_boundary_type(fixity[0])
        match btype:
            case 'fixed'|'pinned'|'guide'|'free'|'roller':
                boundary = get_boundary_by_name(bname=btype)
                boundary.insert(0, 'constrained')
            case 'constrained':
                boundary = get_node_constrain_list(fixity[2])
                boundary.insert(0, btype)
            case 'spring':
                raise NotImplementedError
                #boundary.insert(0, btype)
            case 'displacement':
                boundary = get_node_disp_list(fixity[1])
                boundary.insert(0, btype)
            case _:
                raise NotImplementedError

    elif isinstance(fixity, dict):
        btype = find_boundary_type(fixity['type'])
    else:
        raise Exception('   *** Boundary input format not recognized')
    #
    return boundary
#
def get_node_boundary(fixity:list|tuple|dict|str):
    """ [type, x,y,z,rx,ry,rz] """

    if isinstance(fixity, str):
        boundary = get_boundary_by_name(bname=fixity)
    
    elif isinstance(fixity, (list, tuple)):
        if isinstance(fixity[0], (list, tuple)):
            fixity = fixity[0]
        boundary = get_node_constrain_list(fixity)
    
    elif isinstance(fixity, dict):
        fixity = {find_support_boundary(key): item
                  for key, item in fixity.items()}
        boundary = get_node_constrain_list([fixity['x'], fixity['y'], fixity['z'],
                                        fixity['rx'], fixity['ry'], fixity['rz']])
    
    else:
        raise Exception('   *** Boundary input format not recognized')

    return boundary
#
# -----------------
#
def get_node_constrain_list(data: list) -> list:
    """ [x,y,z,rx,ry,rz]"""
    fixity = [None] * 6
    valid = [0, 1, True, False]
    #
    for x in range(len(fixity)):
        try:
            if isinstance(data[x], (int, float, bool)):
                if not int(data[x]) in valid:
                    raise IOError(f'fixity code {data[x]} not valid')
                elif data[x]:
                    fixity[x] = 0
                #fixity[x] = int(data[x])
            elif isinstance(data[x], str):
                if re.match(r"\b(fix(ed)?|true)\b", data[x], re.IGNORECASE):
                    fixity[x] = 0
                elif re.match(r"\b(free|false)\b", data[x], re.IGNORECASE):
                    fixity[x] = None
                else:
                    raise IOError(f"fixity code {data[x]} not valid")            
            else:
                raise IOError('fixity format not valid')
        except IndexError:
            continue
    #
    return fixity
#
def get_node_disp_list(data:list) -> list:
    """ [x,y,z,rx,ry,rz]"""
    fixity = [None] * 6
    for x in range(len(fixity)):
        try:
            fixity[x] = data[x].value
        except IndexError:
            pass
        except AttributeError:
            if data[x] == None:
                continue
            elif data[x] is False:
                continue
            raise IOError('units missing')
    #
    return fixity
#
def get_nboundary_dict(values: dict):
    """
    return [node, x, y, z, rx, ry, rz]
    """
    1 / 0
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
    """
    Fix = 0
    Free = None
    
    Returns:
        [x,y,z,rx,ry,rz]
    """
    bname = find_boundary_by_name(bname)
    match bname:
        case 'fixed':
            return [0,0,0, 0,0,0]
        case 'pinned':
            return [0,0,0, 0,None,None]
        case 'guide':
            return [0,None,None, 0,None,None]
        case 'free':
            return [None,None,None,None,None,None]
        case 'roller':
            return [None,0,0, 0,None,None]
        case _:
            raise IOError(f"boundary type {bname} not implemented")
    #
    #if re.match(r"\b(fix(ed)?)\b", bname, re.IGNORECASE):
    #    return [0,0,0, 0,0,0]
    
    #elif re.match(r"\b(pinn(ed)?)\b", bname, re.IGNORECASE):
    #    return [0,0,0, 0,None,None]
    #
    #elif re.match(r"\b(roll(er)?)\b", bname, re.IGNORECASE):
    #    return [None,0,0, 0,None,None]

    #elif re.match(r"\b(guide(d)?)\b", bname, re.IGNORECASE):
    #    return [0,None,None, 0,None,None]
    #
    #elif re.match(r"\b(free)\b", bname, re.IGNORECASE):
    #    return [None,None,None,None,None,None]
    #    
    #else:
    #    raise IOError(f"boundary type {bname} not implemented")
    #
    #return value
#
