#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
from collections import defaultdict #, Counter
#from dataclasses import dataclass
from typing import NamedTuple
import re
#from unittest import case

# package imports
import steelpy.utils.io_module.text as common
from steelpy.ufo.load.process.utils import get_value_point
from steelpy.utils.dataframe.main import DBframework
#   
#
# ---------------------------------
#
def get_beam_line_load(load:list|tuple|dict) -> list:
    """
    Returns:
        [line, qx1,qy1,qz1,qt1,qx2,qy2,qz2,qt2,L0,L1, comment]
    """
    #print('--->')
    if isinstance(load, UDL):
        load = get_BeamLine_load(load)
    elif isinstance(load, dict):
        try:
            load = get_BeamLine_dict(load)
        except AttributeError:
            raise IOError('units missing')
    elif isinstance(load, (list, tuple)):
        try:
            load = get_BeamLine_list(load)
        except AttributeError:
            raise IOError('units missing') 
    else:
        raise IOError('Load input format not recognized')
    #
    return load
#
#
def get_beam_point_load(load:list|tuple|dict) -> list:
    """
    Returns:
        [point,Fx, Fy, Fz, Mx, My, Mz, L0, comment]
        [mass, x, y, z, mx, my, mz, L0, comment]
    """
    #print('--->')
    if isinstance(load, PLoad):
        load = get_BeamPoint_load(load)
    elif isinstance(load, dict):
        try:
            load = get_BeamPoint_dict(load)
        except AttributeError:
            raise IOError('units missing')    
    elif isinstance(load, (list, tuple)):
        try:
            load = get_BeamPoint_list(load)
        except AttributeError:
            raise IOError('units missing')
    else:
        raise IOError('Load input format not recognized')
    #
    return load     
#
#
def get_BeamLine_load(data:list|tuple|UDL)->list:
    """
    Returns:
         [line, q1x,q1y,q1z,q1t,q2x,q2y,q2z,q2t,L0,L1, comment]
    """
    try:
        load_type = data.pop(0)
    except AttributeError:
        load_type = 'line'
    #
    if isinstance(data[-1], str):
        try:
            load_comment = data.pop()
        except AttributeError:
            load_comment = data.comment        
    else:
        load_comment = None
    #
    new_data = []
    # q0 [x,y,z,t]
    for x in range(4):
        try:
            new_data.append(data[x].convert("newton/metre").value)
        except AttributeError:
            new_data.append(data[x])
        except IndexError:
            new_data.append(0.0)
    # q1 [x,y,z,t]
    for x in range(4,8):
        try:
            new_data.append(data[x].convert("newton/metre").value)
        except AttributeError:
            new_data.append(data[x])
        except IndexError:
            new_data.append(new_data[x-3])
    # L0 & L1
    for x in range(8,10):
        try:
            new_data.append(data[x].value)
        except AttributeError:
            new_data.append(data[x])
        except IndexError:
            new_data.append(0.0)
    #
    new_data.insert(0, load_type)
    new_data.append(load_comment)    
    return new_data
#
#
def get_BeamPoint_load(data:list|tuple|PLoad)->list:
    """
    Returns:
        [point, Fx, Fy, Fz, Mx, My, Mz, L0, comment]
    """
    try:
        load_type = data.pop(0)
    except AttributeError:
        load_type = 'point'
    #
    if isinstance(data[-1], str):
        try:
            load_comment = data.pop()
        except AttributeError:
            load_comment = data.comment
    else:
        load_comment = None
    #
    new_data = []
    # Fx,y,z
    for x in range(1, 4):
        try:
            new_data.append(data[x].convert("newton").value)
        except AttributeError:
            new_data.append(data[x])
        except IndexError:
            new_data.append(0.0)
    # Mx,y,z
    for x in range(4, 7):
        try:
            new_data.append(data[x].convert("newton*metre").value)
        except AttributeError:
            new_data.append(data[x])
        except IndexError:
            new_data.append(0.0)
    # L0
    try:
        new_data.append(data[0].value)
    except AttributeError:
        new_data.append(data[0])
    #
    new_data.insert(0, load_type)
    new_data.append(load_comment)
    return new_data
#
# ---------------------------------
#
#
def get_value_line(data, label:str, steps:int = 8)->list[float]:
    """
    Returns:
        [q0x,q0y,q0z,q0t, q1x,q1y,q1z,q1t]
    """
    new_data = [0,0,0,0, None,None,None, None]
    for x in range(steps):
        try:
            new_data[x] = data[label][x]
        except IndexError:
            continue
    #
    start = steps // 2
    for x in range(start, steps):
        #print(x)
        if new_data[x] == None:
            new_data[x] = new_data[x-start]
    #
    return new_data
#
# ---------------------------------
#
#
def get_BeamLoad_df(df: DBframework.DataFrame,
                    system: int):
    """
    Returns:
        Dataframe:
        [line, qx1,qy1,qz1,qt1, qx2,qy2,qz2,qt2,L0,L1, tile]
        [point, fx, fy, fz, mx, my, mz, L0, comment]
        [mass, x, y, z, rx, ry, rz, L0, comment]
    """
    #
    flag_system = 'global'
    if system == 1:
        flag_system = 'local'    
    #
    columns = list(df.columns)
    try:
        header = {item:find_BeamLoad_item(item) for item in columns}
        df.rename(columns=header, inplace=True)
        df['type'] = df['type'].apply(lambda x: find_load_type(x))
        #
        if 'system' not in columns:
            df['system'] = flag_system
        #
        grpload = df.groupby(['type', 'beam'])
        newgrp = defaultdict(list)
        for key, items in grpload:
            bload = items.to_dict(orient='split', index=False)
            #
            if re.match(r"\b(line)\b", key[0], re.IGNORECASE):
                # ['type', 'qx1', 'qy1', 'qz1', 'qt1',
                #  'qx2', 'qy2', 'qz2', 'qt2',
                #  'L1', 'L2']
                temp = []
                #for bname, bload in grpbeam:
                for load in bload['data']:
                    data = dict(zip(bload['columns'], load))
                    line = get_BeamLine_dict(data)
                    temp.append([data['name'], data['system'], key[1], *line, None])
                newgrp['line'].extend(temp)
    
            elif re.match(r"\b(point)\b", key[0], re.IGNORECASE):
                # ['type', Fx, Fy, Fz, Mx, My, Mz, L0, comment]
                temp = []
                for load in bload['data']:
                    data = dict(zip(bload['columns'], load))
                    point = get_BeamPoint_dict(data)
                    temp.append([data['name'], data['system'], key[1], *point, None])
                newgrp['point'].extend(temp)
                
            elif re.match(r"\b(mass)\b", key[0], re.IGNORECASE):
                # ['type', x, y, z, mx, my, mz, L0, comment]
                temp = []
                for load in bload['data']:
                    data = dict(zip(bload['columns'], load))
                    mass = get_BeamPoint_dict(data)
                    temp.append([data['name'], data['system'], key[1], *mass, None])
                newgrp['mass'].extend(temp)
                
            else:
                raise NotImplementedError(f'load type {key} not valid')
    except re.PatternError:
        header = {item:find_BeamLoad_type(item) for item in columns}
        df.rename(columns=header, inplace=True)
        #
        newgrp = defaultdict(list)
        grpforce = df.groupby(['name', 'beam'])        
        for key, items in grpforce:
            for item in items.itertuples():
                try: #  [qx1,qy1,qz1,qt1,qx2,qy2,qz2,qt2,L0,L1, tile]
                    data = item.line
                    force = get_BeamLoad_list_units(['line', *data])
                    newgrp['line'].append([item.name, flag_system, item.beam, *force, None])
                except AttributeError:
                    pass
                
                try: # ['type', Fx, Fy, Fz, Mx, My, Mz, L0, comment]
                    data = item.point
                    force = get_BeamLoad_list_units(['point', *data])
                    newgrp['point'].append([item.name, flag_system, item.beam, *force, None])
                except AttributeError:
                    pass
                
                try: # ['type', x, y, z, mx, my, mz, L0, comment]
                    data = item.mass
                    force = get_BeamLoad_list_units(['mass', *data])
                    newgrp['mass'].append([item.name, flag_system, item.beam, *force, None])
                except AttributeError:
                    pass
    #
    #
    db = DBframework()
    header = {'line': ['name', 'system',
                       'beam', 'type', 
                       'qx0', 'qy0', 'qz0', 'qt0',
                       'qx1', 'qy1', 'qz1', 'qt1',
                       'L0', 'L1', 'comment', 'step'],
              'point': ['name', 'system',
                        'beam', 'type', 
                        'fx', 'fy', 'fz', 'mx', 'my', 'mz',
                        'L0', 'comment', 'step'],
              'mass': ['name', 'system',
                       'beam', 'type', 
                       'x', 'y', 'z', 'rx', 'ry', 'rz',
                       'L0', 'comment', 'step']}
    newdf = {}
    for key, item in newgrp.items():
        newdf[key] = db.DataFrame(data=item, columns=header[key])
    #
    return newdf
#
#
def get_BeamLoad_list_units(data: list|tuple) -> list[float] :
    """
    Returns:
        [line, qx1,qy1,qz1,qt1, qx2,qy2,qz2,qt2,L0,L1, tile]
        [point, fx, fy, fz, mx, my, mz, L0, comment]
        [mass, x, y, z, rx, ry, rz, L0, comment]
    """
    loat_type = data.pop(0)
    if isinstance(data[-1], str):
        load_comment = data.pop()
    else:
        load_comment = None
    #
    #
    load = defaultdict(list)
    L = []
    new_data = []
    if re.match(r"\b(line(_|-|\s*)?d(istributed)?(_|-|\s*)?(l(oad)?)?)\b",
                loat_type, re.IGNORECASE):
        for item in data:
            # Force
            if item.units() == 'second^-2*gram': # N/m
                load['N/m'].append(item.convert("newton/metre").value)
            
            #elif item.units() == 'metre^2*second^-2*gram': # N*m
            #    load['N*m'].append(item.convert("newton*metre").value)
            
            # load distance from end nodes
            elif item.units() == 'metre': # m
                L.append(item.value)
        #
        # L0 & L1 
        for x in range(2):
            try:
                L[x]
            except IndexError:
                L.append(0)
        #
        new_data = get_value_line(load, label='N/m')
        new_data.extend(L)
        new_data.insert(0, 'line')
    
    elif re.match(r"\b(p(oint)?)\b", loat_type, re.IGNORECASE):
        for item in data:
            # Force
            if item.units() == 'metre*second^-2*gram': # N 
                load['N'].append(item.convert("newton").value)
            
            elif item.units() == 'metre^2*second^-2*gram': # N*m
                load['N*m'].append(item.convert("newton*metre").value)
            
            # load distance from end nodes
            elif item.units() == 'metre': # m
                L.append(item.value)
        #
        # point and beam point load [Fx, Fy, Fz]
        if 'N' in load:
            new_data = get_value_point(load, label='N', steps=3)
        # point and beam point load [Mx, My, Mz]
        if 'N*m' in load:
            if not new_data:
                new_data = [0,0,0]
            new_data.extend(get_value_point(load, label='N*m', steps=3))
        else:
            if not new_data:
                raise IOError("units not compatible") 
            new_data.extend([0,0,0])
        #
        new_data.append(L[0])
        new_data.insert(0, 'point')
    
    elif re.match(r"\b(mass)\b", loat_type, re.IGNORECASE):
        for item in data:
            # mass
            if item.units() == 'gram': # m
                load['kg'].append(item.convert("kilogram").value)
                
            # load distance from end nodes
            elif item.units() == 'metre': # m
                L.append(item.value)
        #
        if 'kg' in load: # mass
            new_data = get_value_point(load, label='kg', steps=3)
        
        if not new_data:
            raise IOError("units not compatible")
        #
        new_data.append(L[0])
        new_data.insert(0, 'mass')
        
    else:
        raise IOError(f'beam load type {loat_type} not available')    
    #
    new_data.append(load_comment)
    return new_data
#
# ---------------------------------
#
def get_BeamLoad_dict(data: dict)->list:
    """
    Returns:
        [line, qx1,qy1,qz1,qt1, qx2,qy2,qz2,qt2,L0,L1, tile]
        [point, fx, fy, fz, mx, my, mz, L0, comment]
        [mass, x, y, z, rx, ry, rz, L0, comment]
    """
    load_type = data['type']
    if re.match(r"\b(line|(u|v)?(_|-|\s*)?d(istributed)?(_|-|\s*)?(l(oad)?)?)\b",
                load_type, re.IGNORECASE):
        new_data = get_BeamLine_dict(data)
    
    elif re.match(r"\b(p(oint)?)\b", load_type, re.IGNORECASE):
        new_data = get_BeamPoint_dict(data)

    else:
        raise IOError(f'beam load type {load_type} not available')
    #
    return new_data
#
#
def get_BeamLine_dict(data: dict)->list:
    """
    Returns:
        [type,
        qx1, qy1, qz1, qt1,
        qx2, qy2, qz2, qt2,
        L1, L2,
        comment]
    """
    #
    new_data = [0, 0, 0, 0, None,None,None,None, 0, 0, None]
    for key, item in data.items():
        load_item = find_BeamLoad_line(key)
        #print(load_item)
        match load_item:
            # End 1
            case 'qx1':
                new_data[0] = item.convert("newton/metre").value
            case 'qy1':
                new_data[1] = item.convert("newton/metre").value
            case 'qz1':
                new_data[2] = item.convert("newton/metre").value
            case 'qt1':
                new_data[3] = item.convert("newton/metre").value
           #
           # End 2
            case 'qx2':
                new_data[4] = item.convert("newton/metre").value
            case 'qy2':
                new_data[5] = item.convert("newton/metre").value
            case 'qz2':
                new_data[6] = item.convert("newton/metre").value
            case 'qt2':
                new_data[7] = item.convert("newton/metre").value
           #
           # L0, L2
            case 'L1':
                new_data[8] = item.convert("metre").value
            case 'L2':
                new_data[9] = item.convert("metre").value
           #
           # Comment
            case 'comment':
                new_data[10] = item
    #
    # for uniform load
    # qx2
    if new_data[4] == None:
        new_data[4] = new_data[0]
    # qy2
    if new_data[5] == None:
        new_data[5] = new_data[1]
    # qz2
    if new_data[6] == None:
        new_data[6] = new_data[2]
    # qt2    
    if new_data[7] == None:
        new_data[7] = new_data[3]    
    #
    new_data.insert(0, 'line')
    return new_data
#
#
def get_BeamLine_list(data: list|tuple) -> list :
    """
    Returns:
        line = [qx1,qy1,qz1,qt1,qx2,qy2,qz2,qt2,L0,L1, tile]
    """
    load_type = data.pop(0)
    #
    if isinstance(data[-1], str):
        load_comment = data.pop()
    else:
        load_comment = None
    #
    load = defaultdict(list)
    L = []
    new_data = []
    if re.match(r"\b(line|(u|v)?(_|-|\s*)?d(istributed)?(\_)?(l(oad)?)?)\b",
                load_type, re.IGNORECASE):    
        for item in data:
            # Force
            if item.units() == 'second^-2*gram': # N/m
                load['N/m'].append(item.convert("newton/metre").value)
            
            #elif item.units() == 'metre^2*second^-2*gram': # N*m
            #    load['N*m'].append(item.convert("newton*metre").value)
            
            # load distance from end nodes
            elif item.units() == 'metre': # m
                L.append(item.value)
        #
        # L0 & L1 
        for x in range(2):
            try:
                L[x]
            except IndexError:
                L.append(0)
        #
        new_data = get_value_line(load, label='N/m')
        new_data.extend(L)
        new_data.insert(0, 'line')        
    else:
        raise IOError(f'load type {load_type} not valid')
    #
    new_data.append(load_comment)
    return new_data     
#
#
#
def get_BeamPoint_dict(data: dict)->list:
    """
    Returns:
        [point, fx, fy, fz, mx, my, mz, L0, comment]
        [mass, x, y, z, mx, my, mz, L0, comment]
    """
    load_type = data['type']
    if re.match(r"\b(p(oint)?)\b", load_type, re.IGNORECASE):
        #new_data[0] = 'point'
        new_data = [0,0,0, 0,0,0, 0, None]
        for key, item in data.items():
            load_item = find_BeamLoad_PointForce(key)
            #print(load_item)
            match load_item:
                case 'fx':
                    new_data[0] = item.convert("newton").value
                case 'fy':
                    new_data[1] = item.convert("newton").value
                case 'fz':
                    new_data[2] = item.convert("newton").value
                # moment
                case 'mx':
                    new_data[3] = item.convert("newton*metre").value
                case 'my':
                    new_data[4] = item.convert("newton*metre").value
                case 'mz':
                    new_data[5] = item.convert("newton*metre").value
                # Load distance
                case 'L':
                    new_data[6] = item.convert("metre").value
                # load comment
                case 'comment':
                    new_data[7] = item
        #
        new_data.insert(0, 'point')
    
    elif re.match(r"\b(mass)\b", load_type, re.IGNORECASE):
        #new_data[0] = 'mass'
        new_data = [0,0,0, 0,0,0, 0, None]
        for key, item in data.items():
            load_item = find_BeamLoad_PointMass(key)
            #print(load_item)
            match load_item:
                case 'x':
                    new_data[0] = item.convert("kilogram").value
                case 'y':
                    new_data[1] = item.convert("kilogram").value
                case 'z':
                    new_data[2] = item.convert("kilogram").value
            #
            # Load distance
                case 'L':
                    new_data[6] = item.convert("metre").value
        # load comment
                case 'comment':
                    new_data[7] = item
        #
        new_data.insert(0, 'mass')
    else:
        raise IOError(f'beam load type {load_type} not available')
    #
    return new_data    
#
def get_BeamPoint_list(data: list|tuple) -> list :
    """
    Returns:
        [point, fx, fy, fz, mx, my, mz, L0, comment]
        [mass, x, y, z, mx, my, mz, L0, comment]
    """
    load_type = data.pop(0)
    if isinstance(data[-1], str):
        load_comment = data.pop()
    else:
        load_comment = None
    #
    #
    load = defaultdict(list)
    L = []
    new_data = []    
    if re.match(r"\b(p(oint)?)\b", load_type, re.IGNORECASE):
        for item in data:
            # Force
            if item.units() == 'metre*second^-2*gram': # N 
                load['N'].append(item.convert("newton").value)
            
            elif item.units() == 'metre^2*second^-2*gram': # N*m
                load['N*m'].append(item.convert("newton*metre").value)
            
            # load distance from end nodes
            elif item.units() == 'metre': # m
                L.append(item.value)
        #
        # point and beam point load [Fx, Fy, Fz]
        if 'N' in load:
            new_data = get_value_point(load, label='N', steps=3)
        # point and beam point load [Mx, My, Mz]
        if 'N*m' in load:
            if not new_data:
                new_data = [0,0,0]
            new_data.extend(get_value_point(load, label='N*m', steps=3))
        else:
            if not new_data:
                raise IOError("units not compatible") 
            new_data.extend([0,0,0])
        #
        new_data.append(L[0])
        new_data.insert(0, 'point')
    
    elif re.match(r"\b(mass)\b", load_type, re.IGNORECASE):
        for item in data:
            # mass
            if item.units() == 'gram': # m
                load['kg'].append(item.convert("kilogram").value)
                
            # load distance from end nodes
            elif item.units() == 'metre': # m
                L.append(item.value)
        #
        if 'kg' in load: # mass
            new_data = get_value_point(load, label='kg', steps=3)
        
        if not new_data:
            raise IOError("units not compatible")
        #
        new_data.append(L[0])
        new_data.insert(0, 'mass')
    #
    new_data.append(load_comment)
    return new_data     
#
# ---------------------------------
#
def find_load_type(word_in:str) -> str:
    """ """
    key = {"line": r"\b((line|u|v)(_|-|\s*)?(d(istributed)?)?(_|-|\s*)?(l(oad)?)?)\b",
           "point": r"\b(p(oint)?)\b",
           "mass": r"\b(mass)\b"}
    match = common.find_keyword(word_in, key)
    return match
#
def find_BeamLoad_basic(word_in:str) -> str:
    """ """
    key = {"load": r"\b(load(s)?(_|-|\s*)?(id|name)?)\b",
           "type": r"\b((load(_|-|\s*)?)?type)\b",
           "title": r"\b((load(_|-|\s*)?)?title)\b",
           "comment": r"\b((load(_|-|\s*)?)?comment)\b",
           "system": r"\b((load(_|-|\s*)?)?system)\b",
           "beam": r"\b(beam(s)?(_|-|\s*)?(name|id)?)\b",}
    match = common.find_keyword(word_in, key)
    return match
#
def find_BeamLoad_item(word_in:str) -> str:
    """ """
    try:
        return find_BeamLoad_basic(word_in)
    except IOError:
        return find_load_item(word_in)
#
def find_BeamLoad_type(word_in:str) -> str:
    """ """
    try:
        return find_BeamLoad_basic(word_in)
    except IOError:
        return find_load_type(word_in)
#
def find_BeamLoad_line(word_in:str) -> str:
    """ """
    try:
        return find_BeamLoad_basic(word_in)
    except IOError:
        return find_LineLoad_item(word_in)
#
def find_BeamLoad_PointForce(word_in:str) -> str:
    """ """
    try:
        return find_BeamLoad_basic(word_in)
    except IOError:
        return find_PointLoad_item(word_in)
#
def find_BeamLoad_PointMass(word_in:str) -> str:
    """ """
    try:
        return find_BeamLoad_basic(word_in)
    except IOError:
        return find_PointMass_item(word_in)
#
def find_load_item(word_in:str) -> str:
    """ """
    try:
        return find_LineLoad_item(word_in)
    except IOError:
        pass
    
    try:
        return find_PointLoad_item(word_in)
    except IOError:
        pass
    
    try:
        return find_PointMass_item(word_in)
    except IOError:
        pass
    #
    return IOError('beam load definition')
#
#
#
def find_LineLoad_item(word_in:str) -> str:
    """ """
    key = {"qx1": r"\b((qx|qa(xial)?)(_|-|\s*)?(0|1|start)?)\b",
           "qy1": r"\b((qy|in(_|-|\s*)?plane)(_|-|\s*)?(0|1|start)?)\b",
           "qz1": r"\b((qz|out(_|-|\s*)?plane)(_|-|\s*)?(0|1|start)?)\b",
           "qt1": r"\b((qt(orsion)?|mx)(_|-|\s*)?(0|1|start)?)\b",
           #
           "qx2": r"\b((qx|qa(xial)?)(_|-|\s*)?(2|end))\b",
           "qy2": r"\b((qy|in(_|-|\s*)?plane)(_|-|\s*)?(2|end))\b",
           "qz2": r"\b((qz|out(_|-|\s*)?plane)(_|-|\s*)?(2|end))\b",
           "qt2": r"\b((qt(orsion)?|mx)(_|-|\s*)?(2|end))\b",
           #
           "L1": r"\b((l|d(istance)?)(_|-|\s*)?(0|1|start)?)\b",
           "L2": r"\b((l|d(istance)?)(_|-|\s*)?(2|end))\b"}
    match = common.find_keyword(word_in, key)
    return match
#
def find_PointMass_item(word_in:str) -> str:
    """ """
    key = {"x": r"\b((m(_|-|\s*)?)?x)\b",
           "y": r"\b((m(_|-|\s*)?)?y)\b",
           "z": r"\b(\b((m(_|-|\s*)?)?z)\b",
           #
           "rx": r"\b(mx|t(orsion)?)\b",
           "ry": r"\b(my|out(_|-|\s*)?plane)\b",
           "rz": r"\b(mz|in(_|-|\s*)?plane)\b",
           #
           "L": r"\b((l|d(istance)?)(_|-|\s*)?(0|1|start)?)\b",}
    match = common.find_keyword(word_in, key)
    return match
#
def find_PointLoad_item(word_in:str) -> str:
    """ """
    key = {"fx": r"\b((p|f(orce)?)(_|-|\s*)?(x|a(xial)?))\b",
           "fy": r"\b((p|f(orce)?)(_|-|\s*)?(y|in(_|-|\s*)?plane))\b",
           "fz": r"\b((p|f(orce)?)(_|-|\s*)?(z|out(_|-|\s*)?plane))\b",
           #
           "mx": r"\b(m(oment)?(_|-|\s*)?(x|t(orsion)?))\b",
           "my": r"\b(m(oment)?(_|-|\s*)?(y|out(_|-|\s*)?plane))\b",
           "mz": r"\b(m(oment)?(_|-|\s*)?(z|in(_|-|\s*)?plane))\b",
           #
           "L": r"\b((l|d(istance)?)(_|-|\s*)?(0|1|start)?)\b",}
    
    match = common.find_keyword(word_in, key)
    return match
#
#
# ---------------------------------
#
class UDL(NamedTuple):
    """ Beam UDL load"""
    qx1: float
    qy1: float
    qz1: float
    qt1: float
    #
    qx2: float
    qy2: float
    qz2: float
    qt2: float
    #
    L1: float | None 
    L2: float | None 
    #
    comment: str | int| None
#
class PLoad(NamedTuple):
    """ Beam point load"""
    Fx: float
    Fy: float
    Fz: float
    #
    Mx: float
    My: float
    Mz: float
    #
    L: float | None 
    #
    comment: str | int| None
    #
    @property
    def L1(self):
        """ """
        return self.L