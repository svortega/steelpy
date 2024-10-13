#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
from collections import defaultdict, Counter
#from dataclasses import dataclass
#from typing import NamedTuple
#import re

# package imports
#
#
# ---------------------------------
#
def check_list_number(data, steps:int=6)->list[float]:
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
    #
    if not new_data:
        new_data = [0] * steps
    return new_data
#
#
def get_dict_load(items, units:str):
    """ """
    if isinstance(items, (list, tuple)):
        udl = []
        for item in items:
            try:
                udl.append(item.convert(units).value)
            except AttributeError:
                raise IOError(f'units missing')
    else:
        try:
            udl = [items.convert(units).value]
        except AttributeError:
            raise IOError(f'units missing')
    return udl
#
def get_items(data: list, new_data: list,
              rows: int, fill: list|None = None):
    """ """
    if not fill:
        fill = [0] * rows
    
    for x, item in enumerate(data):
        if not item:
            new_data[x] = fill
        else:
            new = []
            for idx in range(rows):
                try:
                    new.append(item[idx])
                except IndexError:
                    new.append(item[idx-1])
            new_data[x] = new
    return new_data
#
# ---------------------------------
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