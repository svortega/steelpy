#
# Copyright (c) 2009-2021 fem2ufo
# 

# Python stdlib imports
from array import array
from collections.abc import Mapping
from collections import defaultdict, Counter
#from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, Iterator, Dict, Iterable, ClassVar, Union
import re

# package imports
#
#
#
def get_beam_line_load(load):
    """
    """
    #print('--->')
    if isinstance(load, dict):
        load = check_beam_dic(load)
    
    elif isinstance(load, (list, tuple)):
        load = check_list_units(load)
    
    else:
        raise Exception('   *** Load input format not recognized')
    return load    
    
#
#
def get_beam_point_load(load):
    """
    """
    #print('--->')
    if isinstance(load, dict):
        load = check_point_dic(load)
    
    elif isinstance(load, (list, tuple)):
        load = check_list_units(load)
    
    else:
        raise Exception('   *** Load input format not recognized')
    return load     
    
#
#
#
def get_beam_concept_load(load):
    """ """
    if isinstance(load, (list, tuple)):
        try:
            load = check_list_units(load)
        except AttributeError:
            load = check_list_number(load, steps=7)
            load.append(load.pop(0))
    elif isinstance(load, dict):
        load = check_point_dic(load)
    else:
        raise Exception('   *** Load input format not recognized')
    return load
#
#
def get_beam_load(data, steps:int=6)->List[float]:
    """ """
    new_data = []
    # q0
    for x in range(3):
        try:
            new_data.append(data[x])
        except IndexError:
            new_data.append(0.0)
    # q1
    for x in range(3,6):
        try:
            new_data.append(data[x])
        except IndexError:
            new_data.append(new_data[x-3])
    # L
    for x in range(6,8):
        try:
            new_data.append(data[x])
        except IndexError:
            new_data.append(0.0)
    return new_data
#
#
#
def check_list_number(data, steps:int=6)->List[float]:
    """ """
    new_data = []
    for x in range(steps):
        try:
            new_data.append(data[x])
        except IndexError:
            new_data.append(0.0)
    return new_data
#
#
def get_value_point(data, label:str, steps:int):
    """ """
    new_data = []
    for x in range(steps):
        try:
            new_data.append(data[label][x])
        except IndexError:
            new_data.append(0)
    return new_data
#
def get_value_line(data, label:str, steps:int):
    """ """
    new_data = [0,0,0, None,None,None]
    for x in range(steps):
        try:
            new_data[x] = data[label][x]
        except IndexError:
            continue
    #
    if new_data[3] == None:
        new_data[3] = new_data[0]

    if new_data[4] == None:
        new_data[4] = new_data[1]

    if new_data[5] == None:
        new_data[5] = new_data[2]
    #
    return new_data
#
#
def check_list_units(data)->List[float]:
    """ """
    load = defaultdict(list)
    L = []
    lendata = len(data)
    for x in range(lendata):
        if data[x].units() == 'metre*second^-2*gram': # N 
            load['N'].append(data[x].convert("newton").value)
        elif data[x].units() == 'metre^2*second^-2*gram': # N*m
            load['N*m'].append(data[x].convert("newton*metre").value)
        elif data[x].units() == 'second^-2*gram': # N/m
            load['N/m'].append(data[x].convert("newton/metre").value)
        elif data[x].units() == 'metre': # m
            L.append(data[x].value)
        else:
            raise IOError("units {:} not compatible"
                          .format(data[x].units()))
    #
    for x in range(2):
        try:
            L[x]
        except IndexError:
            L.append(0)
    #
    new_data = []
    # point and beam point load
    if 'N' in load:
        new_data = get_value_point(load, label='N', steps=3)

    if 'N*m' in load:
        if not new_data:
            new_data = [0,0,0]
        new_data.extend(get_value_point(load, label='N*m', steps=3))

    if new_data:
       for x in range(6):
           try:
               new_data[x]
           except IndexError:
               new_data.append(0)
       new_data.append(L[0])
       return new_data
    #
    # beam line load
    elif 'N/m' in load:
        new_data = get_value_line(load, label='N/m', steps=6)
        new_data.extend(L)
        return new_data
    else:
        raise IOError("units not compatible")
#
#
#
#
def check_point_dic(data)->List[float]:
    """
    new_data: [fx, fy, fz, mx, my, mz, d1]
    """
    new_data = [0,0,0, 0,0,0, 0, 'NULL']
    for key, item in data.items():
        if re.match(r"\b(fx|fa(xial)?)\b", str(key), re.IGNORECASE):
            new_data[0] = item.convert("newton").value
        elif re.match(r"\b(py|fy|in(_)?plane)\b", str(key), re.IGNORECASE):
            new_data[1] = item.convert("newton").value
        elif re.match(r"\b(pz|fz|out(_)?plane)\b", str(key), re.IGNORECASE):
            new_data[2] = item.convert("newton").value
        elif re.match(r"\b(mx|t(orsion)?)\b", str(key), re.IGNORECASE):
            new_data[3] = item.convert("newton*metre").value
        elif re.match(r"\b(my|in(_)?plane)\b", str(key), re.IGNORECASE):
            new_data[4] = item.convert("newton*metre").value
        elif re.match(r"\b(mz|out(_)?plane)\b", str(key), re.IGNORECASE):
            new_data[5] = item.convert("newton*metre").value
        elif re.match(r"\b((l|d(istance)?)(_)?(1|start))\b", str(key), re.IGNORECASE):
            new_data[6] = item.value
        elif re.match(r"\b(title|comment)\b", str(key), re.IGNORECASE):
            new_data[7] = item
    return new_data
#
#
def check_beam_dic(data)->List[float]:
    """
    new_data: [qx1, qy1, qz1, qx2, qy2, qz2, d1, d2]
    """
    new_data = [0,0,0, None,None,None, 0,0, 'NULL']
    for key, item in data.items():
        if re.match(r"\b((qx|t(orsion)?)(_)?(1|start)?)\b", str(key), re.IGNORECASE):
            new_data[0] = item.convert("newton/metre").value
        elif re.match(r"\b((qy|in(_)?plane)(_)?(1|start)?)\b", str(key), re.IGNORECASE):
            new_data[1] = item.convert("newton/metre").value
        elif re.match(r"\b((qz|out(_)?plane)(_)?(1|start)?)\b", str(key), re.IGNORECASE):
            new_data[2] = item.convert("newton/metre").value
        elif re.match(r"\b((qx|t(orsion)?)(_)?(2|end))\b", str(key), re.IGNORECASE):
            new_data[3] = item.convert("newton/metre").value
        elif re.match(r"\b((qy|in(_)?plane)(_)?(2|end))\b", str(key), re.IGNORECASE):
            new_data[4] = item.convert("newton/metre").value
        elif re.match(r"\b((qz|out(_)?plane)(_)?(2|end))\b", str(key), re.IGNORECASE):
            new_data[5] = item.convert("newton/metre").value
        elif re.match(r"\b((l|d(istance)?)(_)?(1|start))\b", str(key), re.IGNORECASE):
            new_data[6] = item.value
        elif re.match(r"\b((l|d(istance)?)(_)?(2|end))\b", str(key), re.IGNORECASE):
            new_data[7] = item.value
        elif re.match(r"\b(title|comment)\b", str(key), re.IGNORECASE):
            new_data[8] = item        
    #
    if new_data[3] == None:
        new_data[3] = new_data[0]

    if new_data[4] == None:
        new_data[4] = new_data[1]

    if new_data[5] == None:
        new_data[5] = new_data[2]
    #
    return new_data
#
#
def duplicates(lst):
    cnt= Counter(lst)
    return [key for key in cnt.keys() if cnt[key]> 1]
#
def indices(lst, items= None):
    items, ind= set(lst) if items is None else items, defaultdict(list)
    for i, v in enumerate(lst):
        if v in items: ind[v].append(i)
    return ind
#
#