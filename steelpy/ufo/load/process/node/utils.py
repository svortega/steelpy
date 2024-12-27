#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from collections import defaultdict
import re
#from unittest import case

# package imports
from steelpy.ufo.load.process.utils import (get_value_point,
                                            check_list_number)

import steelpy.utils.io_module.text as common
from steelpy.utils.dataframe.main import DBframework
#
# ---------------------------------
#
def find_NodeLoad_basic(word_in:str) -> str:
    """ """
    key = {"load": r"\b(load(s)?(_|-|\s*)?(id|name)?)\b",
           "type": r"\b((load(_|-|\s*)?)?type)\b",
           "title": r"\b((load(_|-|\s*)?)?title)\b",
           "comment": r"\b((load(_|-|\s*)?)?comment)\b",
           "node": r"\b(node(s)?(_|-|\s*)?(name|id)?)\b",
           "system": r"\b(system)\b",}
    match = common.find_keyword(word_in, key)
    return match
#
def find_NodeLoad_item(word_in:str) -> str:
    """ """
    try:
        return find_NodeLoad_basic(word_in)
    except IOError:
        return find_load_item(word_in)
#
def find_NodeLoad_type(word_in:str) -> str:
    """ """
    key = {"name": r"\b(id|name|load(s)?)\b",
           "type": r"\b((load(_|-|\s*)?)?type)\b",
           "title": r"\b(title|comment)\b",
           "node": r"\b(node(s)?(_|-|\s*)?(name|id)?)\b"}
    
    try:
        match = common.find_keyword(word_in, key)
        return match
    except IOError:
        return find_load_type(word_in)
#
#
#
def find_NodeLoad_force(word_in:str) -> str:
    """ """
    try:
        return find_NodeLoad_basic(word_in)
    except IOError:
        return find_NodeForce_item(word_in)
#
def find_NodeLoad_displacement(word_in:str) -> str:
    """ """
    try:
        return find_NodeLoad_basic(word_in)
    except IOError:
        return find_NodeDisplacement_item(word_in)
#
def find_NodeLoad_mass(word_in:str) -> str:
    """ """
    try:
        return find_NodeLoad_basic(word_in)
    except IOError:
        return find_NodeMass_item(word_in)
#
#
def find_load_type(word_in:str) -> str:
    """ """
    key = {"force": r"\b(force|load)\b",
           "displacement": r"\b((prescribed)?disp(lacement)?)\b",
           "mass": r"\b(mass)\b"}
    match = common.find_keyword(word_in, key)
    return match
#
def find_load_item(word_in:str) -> str:
    """ """
    try:
        return find_NodeForce_item(word_in)
    except IOError:
        pass
    
    try:
        return find_NodeDisplacement_item(word_in)
    except IOError:
        pass
    
    try:
        return find_NodeMass_item(word_in)
    except IOError:
        pass
    #
    return IOError('node load definition')
#
def find_NodeForce_item(word_in:str) -> str:
    """ """
    key = {"fx": r"\b((p|f(orce)?)(_|-|\s*)?(x|a(xial)?))\b",
           "fy": r"\b((p|f(orce)?)(_|-|\s*)?(y|in(_|-|\s*)?plane))\b",
           "fz": r"\b((p|f(orce)?)(_|-|\s*)?(z|out(_|-|\s*)?plane))\b",
           #
           "mx": r"\b(m(oment)?(_|-|\s*)?(x|t(orsion)?))\b",
           "my": r"\b(m(oment)?(_|-|\s*)?(y|out(_|-|\s*)?plane))\b",
           "mz": r"\b(m(oment)?(_|-|\s*)?(z|in(_|-|\s*)?plane))\b"}
    match = common.find_keyword(word_in, key)
    return match
#
def find_NodeDisplacement_item(word_in:str) -> str:
    """ """
    key = {"x": r"\b((delta)?(_|-|\s*)?x)\b",
           "y": r"\b((delta)?(_|-|\s*)?y)\b",
           "z": r"\b((delta)?(_|-|\s*)?z)\b",
           #
           "rx": r"\b(rx|t(orsion)?)\b",
           "ry": r"\b(ry|out(_|-|\s*)?plane)\b",
           "rz": r"\b(rz|in(_|-|\s*)?plane)\b"}
    
    match = common.find_keyword(word_in, key)
    return match
#
def find_NodeMass_item(word_in:str) -> str:
    """ """
    key = {"x": r"\b((m(ass)?)?(_|-|\s*)?x)\b",
           "y": r"\b((m(ass))?(_|-|\s*)?y)\b",
           "z": r"\b(\b((m(ass))?(_|-|\s*)?z)\b",}
           #
           #"rx": r"\b(mx|t(orsion)?)\b",
           #"ry": r"\b(my|out(_|-|\s*)?plane)\b",
           #"rz": r"\b(mz|in(_|-|\s*)?plane)\b"}
    
    match = common.find_keyword(word_in, key)
    return match
#
# ---------------------------------
#
def get_NodaLoad_df(df: DBframework.DataFrame,
                    system: int = 0):
    """
    Returns:
        DataFrame: [force, fx, fy, fz, mx, my, mz, comment]
                   [displacement, x, y, z, rx, ry, rz, comment]
                   [mass,x, y, z, mx, my, mz, comment]
    """
    #
    flag_system = 'global'
    if system == 1:
        flag_system = 'local'     
    #
    columns = list(df.columns)
    #
    # node load input by load item
    try:
        header = {key: find_NodeLoad_item(key)
                  for key in columns}
        df.rename(columns=header, inplace=True)
        df['type'] = df['type'].apply(lambda x: find_load_type(x))
        df['system'] = flag_system
        #
        grpload = df.groupby(['type', 'node'])
        newgrp = defaultdict(list)
        for key, items in grpload:
            nload = items.to_dict(orient='split', index=False)
            #
            if re.match(r"\b(force)\b", key[0], re.IGNORECASE):
                temp = []
                for load in nload['data']:
                    #  [force, fx, fy, fz, mx, my, mz, title]
                    data = dict(zip(nload['columns'], load))
                    force = get_NodeForce_dict(data)
                    temp.append([data['load'], data['system'], key[1], *force, None])
                newgrp['force'].extend(temp)
            
            elif re.match(r"\b(displacement)\b", key[0], re.IGNORECASE):
                temp = []
                for load in nload['data']:
                    #  [displacement, x, y, z, rx, ry, rz, title]
                    data = dict(zip(nload['columns'], load))
                    force = get_NodeDisp_dict(data)
                    temp.append([data['load'], data['system'], key[1], *force, None])
                newgrp['displacement'].extend(temp)
            
            elif re.match(r"\b(mass)\b", key[0], re.IGNORECASE):
                temp = []
                for load in nload['data']:
                    #  [mass, x, y, z, rx, ry, rz, title]
                    data = dict(zip(nload['columns'], load))
                    force = get_NodeMass_dict(data)
                    temp.append([data['load'], data['system'], key[1], *force, None])
                newgrp['mass'].extend(temp)
            
            else:
                raise IOError(f'node load type {key[0]} not available')
    # node input by type
    except re.PatternError:
        header = {key: find_NodeLoad_type(key)
                  for key in columns}
        df.rename(columns=header, inplace=True)
        #
        newgrp = defaultdict(list)
        grpforce = df.groupby(['load', 'node'])
        for key, items in grpforce:
            for item in items.itertuples():
                try: #  [force, fx, fy, fz, mx, my, mz, title]
                    data = item.force
                    force = get_NodeLoad_list_units(['force', *data])
                    newgrp['force'].append([item.load, flag_system, item.node, *force, None])
                except AttributeError:
                    pass
                
                try: #  [mass, x, y, z, rx, ry, rz, title]
                    data = item.mass
                    force = get_NodeLoad_list_units(['mass', *data])
                    newgrp['mass'].append([item.load, flag_system, item.node, *force, None])
                except AttributeError:
                    pass
                
                try: #  [displacement, x, y, z, rx, ry, rz, title]
                    data = item.displacement
                    force = get_NodeLoad_list_units(['displacement', *data])
                    newgrp['displacement'].append([item.load, flag_system, item.node, *force, None])
                except AttributeError:
                    pass
    #
    db = DBframework()
    header = {'force': ['load', 'system',
                        'node', 'type', 
                        'fx', 'fy', 'fz', 'mx', 'my', 'mz',
                        'comment', 'step'],
              'displacement': ['load', 'system',
                               'node', 'type', 
                               'x', 'y', 'z', 'rx', 'ry', 'rz',
                               'comment', 'step'],
              'mass': ['load', 'system',
                       'node', 'type', 
                       'x', 'y', 'z', 'rx', 'ry', 'rz',
                       'comment', 'step']}
    newdf = {}
    for key, item in newgrp.items():
        newdf[key] = db.DataFrame(data=item, columns=header[key])
    #
    return newdf
#
#
def get_nodal_load(load: list|tuple|dict)->list:
    """
    Returns:
        [force, fx, fy, fz, mx, my, mz, comment]
        [displacement, x, y, z, rx, ry, rz, comment]
        [mass,x, y, z, mx, my, mz, comment]
    """
    if isinstance(load, (list, tuple)):
        try:
            load = get_NodeLoad_list_units(load.copy())
        except AttributeError:
            load = get_NodeLoad_list(load)
    elif isinstance(load, dict):
        load = get_NodeLoad_dict(load)
    else:
        raise Exception('   *** Load input format not recognized')
    return load
#
def get_NodeLoad_list(data: list|tuple,
                      steps: int = 6) -> list :
    """
    Returns:
        [force, fx, fy, fz, mx, my, mz, comment]
        [displacement, x, y, z, rx, ry, rz, comment]
        [mass,x, y, z, mx, my, mz, comment]
    """
    loat_type = data.pop(0)
    if isinstance(data[-1], str):
        load_title = data.pop()
    else:
        load_title = None
    #
    new_data = check_list_number(data, steps)
    new_data.append(load_title)
    new_data.insert(0, loat_type)
    return new_data
#
def get_NodeLoad_list_units(data: list|tuple) -> list :
    """
    Returns:
        [force, fx, fy, fz, mx, my, mz, comment]
        [displacement, x, y, z, rx, ry, rz, comment]
        [mass,x, y, z, mx, my, mz, comment]
    """
    loat_type = data.pop(0)
    if isinstance(data[-1], str):
        load_title = data.pop()
    else:
        load_title = None
    #
    new_data = [] # [None, 0,0,0, 0,0,0, 'NULL']
    load = defaultdict(list)
    #
    if re.match(r"\b((prescribed)?disp(lacement)?)\b", loat_type, re.IGNORECASE):
        for item in data:
            if item.units() == 'metre': # m
                load['m'].append(item.value)
            
            elif item.units() == 'radian': # rad
                load['rad'].append(item.value)
        #
        if 'm' in load:
            new_data = get_value_point(load, label='m', steps=3)
        else:
            new_data = [0,0,0]
        
        if 'rad' in load:
            #if not new_data:
            #    new_data = [0,0,0]
            new_data.extend(get_value_point(load, label='rad', steps=3))
        else:
            new_data.extend([0,0,0])
        #
        new_data.insert(0, 'displacement')
    
    elif re.match(r"\b(force|load)\b", loat_type, re.IGNORECASE):
        for item in data:
            if item.units() == 'metre*second^-2*gram': # N 
                load['N'].append(item.convert("newton").value)
            
            elif item.units() == 'metre^2*second^-2*gram': # N*m
                load['N*m'].append(item.convert("newton*metre").value)
        #
        # point and beam point load [Fx, Fy, Fz]
        if 'N' in load:
            new_data = get_value_point(load, label='N', steps=3)
        else:
            new_data = [0,0,0]
        # point and beam point load [Mx, My, Mz]
        if 'N*m' in load:
            #if not new_data:
            #    new_data = [0,0,0]
            new_data.extend(get_value_point(load, label='N*m', steps=3))
        else:
            new_data.extend([0,0,0])
        #
        new_data.insert(0, 'load')
    
    elif re.match(r"\b(mass)\b", loat_type, re.IGNORECASE):
        for item in data:
            # Mass
            if item.units() == 'gram': # m
                load['kg'].append(item.convert("kilogram").value)
        #
        if 'kg' in load: # mass
            new_data = get_value_point(load, label='kg', steps=6)
            #return new_data
        #
        new_data.insert(0, 'mass')
    
    else:
        raise IOError(f'node load type {loat_type} not available')
    #
    new_data.append(load_title)
    return new_data
#
def get_NodeLoad_dict(data: dict)->list:
    """
    Returns:
        [force, fx, fy, fz, mx, my, mz, comment]
        [displacement, x, y, z, rx, ry, rz, comment]
        [mass,x, y, z, mx, my, mz, comment]
    """
    new_data = [None, 0,0,0, 0,0,0, None]
    loat_type = data['type']
    if re.match(r"\b((prescribed)?disp(lacement)?)\b", loat_type, re.IGNORECASE):
        new_data = get_NodeDisp_dict(data)           
    
    elif re.match(r"\b(force|load)\b", loat_type, re.IGNORECASE):
        new_data = get_NodeForce_dict(data)        
    
    elif re.match(r"\b(mass)\b", loat_type, re.IGNORECASE):
        new_data = get_NodeMass_dict(data)
    
    else:
        raise IOError(f'node load type {loat_type} not available')
    #        
    return new_data
#
def get_NodeForce_dict(data: dict)->list:
    """
    Returns:
        [force, fx, fy, fz, mx, my, mz, comment] """
    new_data = ['force', 0,0,0, 0,0,0, None]
    for key, item in data.items():
        load_item = find_NodeLoad_force(key)
        match load_item:
            case 'fx':
                new_data[1] = item.convert("newton").value
            case 'fy':
                new_data[2] = item.convert("newton").value
            case 'fz':
                new_data[3] = item.convert("newton").value
            case 'mx':
                new_data[4] = item.convert("newton*metre").value
            case 'my':
                new_data[5] = item.convert("newton*metre").value
            case 'mz':
                new_data[6] = item.convert("newton*metre").value
            case 'comment':
                new_data[7] = item
    return new_data
#
def get_NodeDisp_dict(data: dict)->list:
    """
    Returns:
         [displacement, x, y, z, rx, ry, rz, comment]
    """
    new_data = ['displacement', 0,0,0, 0,0,0, None]
    for key, item in data.items():
        load_item = find_NodeLoad_displacement(key)
        match load_item:
            case 'x':
                new_data[1] = item.value
            case 'y':
                new_data[2] = item.value
            case 'z':
                new_data[3] = item.value
            case 'rx':
                new_data[4] = item.value
            case 'ry':
                new_data[5] = item.value
            case 'rz':
                new_data[6] = item.value
            case 'comment':
                new_data[7] = item
    return new_data
#
def get_NodeMass_dict(data: dict)->list:
    """
    Returns:
         [mass, x, y, z, rx, ry, rz, comment]
    """
    new_data = ['mass', 0,0,0, 0,0,0, None]
    for key, item in data.items():
        load_item = find_NodeLoad_mass(key)
        match load_item:
            case 'x':
                new_data[1] = item.convert("kilogram").value
            case 'y':
                new_data[2] = item.convert("kilogram").value
            case 'z':
                new_data[3] = item.convert("kilogram").value
            case 'comment':
                new_data[7] = item
    return new_data
#
