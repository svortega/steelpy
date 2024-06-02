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
import re

# package imports
#
#
# ---------------------------------
#
def get_beam_line_load(load):
    """
    """
    #print('--->')
    if isinstance(load, dict):
        load = get_BeamLine_dic(load)
    
    elif isinstance(load, (list, tuple)):
        load = get_BeamLoad_list_units(load)
    
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
        load = get_BeamLoad_dic(load)
    
    elif isinstance(load, (list, tuple)):
        load = get_BeamLoad_list_units(load)
    else:
        raise Exception('   *** Load input format not recognized')
    return load     
#
#
def get_BeamLine_load(data)->list[float]:
    """
    line = [q1x,q1y,q1z,q1t,q2x,q2y,q2z,q2t,L0,L1, title]
    """
    loat_type = data.pop(0)
    if isinstance(data[-1], str):
        load_title = data.pop()
    else:
        load_title = None
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
    new_data.insert(0, loat_type)
    new_data.append(load_title)    
    return new_data
#
#
def get_BeamNode_load(data)->list[float]:
    """
    point = [Fx, Fy, Fz, Mx, My, Mz, L0, title]
    """
    loat_type = data.pop(0)
    if isinstance(data[-1], str):
        load_title = data.pop()
    else:
        load_title = None
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
    new_data.insert(0, loat_type)
    new_data.append(load_title)
    return new_data
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
    return new_data
#
#
def get_value_line(data, label:str, steps:int = 8):
    """
    new_data = [q0x,q0y,q0z,q0t, q1x,q1y,q1z,q1t]
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
#
def get_BeamLoad_list_units(data) -> list[float] :
    """
    line = [qx1,qy1,qz1,qt1,qx2,qy2,qz2,qt2,L0,L1, tile]
    point = [Fx, Fy, Fz, Mx, My, Mz, L0, title]
    mass = [x, y, z, mx, my, mz, L0, title]
    """
    loat_type = data.pop(0)
    if isinstance(data[-1], str):
        load_title = data.pop()
    else:
        load_title = None
    #
    #
    load = defaultdict(list)
    L = []
    new_data = []
    if re.match(r"\b(line|(u|v)?(\_)?d(istributed)?(\_)?(l(oad)?)?)\b",
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
    #
    new_data.append(load_title)
    return new_data    
    #
    #load = defaultdict(list)
    #L = []
    #lendata = len(data)
    #for x in range(lendata):
    #    # Force
    #    if data[x].units() == 'metre*second^-2*gram': # N 
    #        load['N'].append(data[x].convert("newton").value)
    #    
    #    elif data[x].units() == 'metre^2*second^-2*gram': # N*m
    #        load['N*m'].append(data[x].convert("newton*metre").value)
    #    # line
    #    elif data[x].units() == 'second^-2*gram': # N/m
    #        load['N/m'].append(data[x].convert("newton/metre").value)
    #    # mass
    #    elif data[x].units() == 'gram': # m
    #        load['kg'].append(data[x].convert("kilogram").value)
    #        
    #    # load distance from end nodes
    #    elif data[x].units() == 'metre': # m
    #        L.append(data[x].value)
    #    
    #    else:
    #        raise IOError("units {:} not compatible"
    #                      .format(data[x].units()))
    ##
    ## L0 & L1 
    #for x in range(2):
    #    try:
    #        L[x]
    #    except IndexError:
    #        L.append(0)
    ##
    ##
    ## beam line load 
    #if 'N/m' in load:
    #    # [q1x,q1y,q1z,q1t,q2x,q2y,q2z,q2t,L0,L1]
    #    new_data = get_value_line(load, label='N/m')
    #    new_data.extend(L)
    #    return new_data
    # 
    #elif 'kg' in load: # mass
    #    new_data = get_value_point(load, label='kg', steps=3)
    #    
    #else: # point and beam point load
    #    new_data = []
    #    # point and beam point load [Fx, Fy, Fz]
    #    if 'N' in load:
    #        new_data = get_value_point(load, label='N', steps=3)
    #    # point and beam point load [Mx, My, Mz]
    #    if 'N*m' in load:
    #        if not new_data:
    #            new_data = [0,0,0]
    #        new_data.extend(get_value_point(load, label='N*m', steps=3))
    #    #
    #    #
    #    if new_data:
    #        for x in range(6):
    #            try:
    #                new_data[x]
    #            except IndexError:
    #                new_data.append(0)
    #        new_data.append(L[0])
    #        # [Fx, Fy, Fz, Mx, My, Mz, L0]
    #        return new_data
    #    #
    #    raise IOError("units not compatible")
#
# ---------------------------------
#
def get_BeamLoad_dic(data: dict)->list[float]:
    """
    line = [qx1,qy1,qz1,qt1, qx2,qy2,qz2,qt2,L0,L1, tile]
    point: [fx, fy, fz, mx, my, mz, L0, title]
    mass: [x, y, z, rx, ry, rz, L0, title]
    """
    loat_type = data['type']
    #new_data = [0,0,0, 0,0,0,  0, 'NULL'] # 0,0,0,0,0,0,
    if re.match(r"\b(line|(u|v)?(\_)?d(istributed)?(\_)?(l(oad)?)?)\b",
                loat_type, re.IGNORECASE):
        #new_data[0] = 'line'
        new_data = [0,0,0,0, 0,0,0,0, 0,0, None]
        for key, item in data.items():
            # End 1
            if re.match(r"\b(qx(\_)?(0|1|start)?)\b", key, re.IGNORECASE):
                new_data[0] = item.convert("newton/metre").value
            
            elif re.match(r"\b(qy(\_)?(0|1|start)?)\b", key, re.IGNORECASE):
                new_data[1] = item.convert("newton/metre").value
            
            elif re.match(r"\b(qz(\_)?(0|1|start)?)\b", key, re.IGNORECASE):
                new_data[2] = item.convert("newton/metre").value
            
            elif re.match(r"\b(qt(\_)?(0|1|start)?)\b", key, re.IGNORECASE):
                new_data[3] = item.convert("newton*metre").value
            
            #
            # End 2
            #
            elif re.match(r"\b(qx(\_)?(2|end))\b", key, re.IGNORECASE):
                new_data[4] = item.convert("newton/metre").value
            
            elif re.match(r"\b(qy(\_)?(2|end))\b", key, re.IGNORECASE):
                new_data[5] = item.convert("newton/metre").value
            
            elif re.match(r"\b(qz(\_)?(2|end))\b", key, re.IGNORECASE):
                new_data[6] = item.convert("newton/metre").value
            
            elif re.match(r"\b(qt(\_)?(2|end))\b", key, re.IGNORECASE):
                new_data[7] = item.convert("newton*metre").value
            #
            # L
            elif re.match(r"\b((l|d(istance)?)(\_)?(0|1|start)?)\b", key, re.IGNORECASE):
                new_data[8] = item.value
            
            elif re.match(r"\b((l|d(istance)?)(\_)?(2|end))\b", key, re.IGNORECASE):
                new_data[9] = item.value
            #
            elif re.match(r"\b(title|comment|name|id)\b", key, re.IGNORECASE):
                new_data[10] = item            
        #
        #
        new_data.insert(0, 'line')        
    
    elif re.match(r"\b(p(oint)?)\b", loat_type, re.IGNORECASE):
        #new_data[0] = 'point'
        new_data = [0,0,0, 0,0,0, 0, None]
        
        for key, item in data.items():
            if re.match(r"\b(fx|fa(xial)?)\b", key, re.IGNORECASE):
                new_data[0] = item.convert("newton").value
            
            elif re.match(r"\b(py|fy|in(_)?plane)\b", key, re.IGNORECASE):
                new_data[1] = item.convert("newton").value
            
            elif re.match(r"\b(pz|fz|out(_)?plane)\b", key, re.IGNORECASE):
                new_data[2] = item.convert("newton").value
            #
            elif re.match(r"\b(mx|t(orsion)?)\b", key, re.IGNORECASE):
                new_data[3] = item.convert("newton*metre").value
            
            elif re.match(r"\b(my|out(_)?plane)\b", key, re.IGNORECASE):
                new_data[4] = item.convert("newton*metre").value
            
            elif re.match(r"\b(mz|in(_)?plane)\b", key, re.IGNORECASE):
                new_data[5] = item.convert("newton*metre").value
            #
            elif re.match(r"\b((l|d(istance)?)(_)?(0|1|start)?)\b", key, re.IGNORECASE):
                new_data[6] = item.value
            
            elif re.match(r"\b(title|comment|name|id)\b", key, re.IGNORECASE):
                new_data[7] = item
        #
        new_data.insert(0, 'point')
    
    elif re.match(r"\b(mass)\b", loat_type, re.IGNORECASE):
        #new_data[0] = 'mass'
        new_data = [0,0,0, 0,0,0, 0, None]
        for key, item in data.items():
            
            if re.match(r"\b(x)\b", key, re.IGNORECASE):
                new_data[0] = item.convert("kilogram").value
            
            elif re.match(r"\b(y)\b", key, re.IGNORECASE):
                new_data[1] = item.convert("kilogram").value
            
            elif re.match(r"\b(z)\b", key, re.IGNORECASE):
                new_data[2] = item.convert("kilogram").value
            #
            #elif re.match(r"\b(mx|t(orsion)?)\b", key, re.IGNORECASE):
            #    new_data[3] = item.convert("newton*metre").value
            #
            #elif re.match(r"\b(my|out(_)?plane)\b", key, re.IGNORECASE):
            #    new_data[4] = item.convert("newton*metre").value
            #
            #elif re.match(r"\b(mz|in(_)?plane)\b", key, re.IGNORECASE):
            #    new_data[5] = item.convert("newton*metre").value
            #
            elif re.match(r"\b((l|d(istance)?)(_)?(0|1|start)?)\b", key, re.IGNORECASE):
                new_data[6] = item.value
            
            elif re.match(r"\b(title|comment|name|id)\b", key, re.IGNORECASE):
                new_data[7] = item
        #
        new_data.insert(0, 'mass')
    
    else:
        raise IOError(f'beam load type {loat_type} not available')    
    #
    #        
    #new_data.append(load_title)
    return new_data
#
#
def get_BeamLine_dic(data: dict)->list:
    """
    new_data: [qx1, qy1, qz1, qt1,
               qx2, qy2, qz2, qt2,
               d1, d2,
               title]
    """
    rows = 0
    new_data = [0, 0, 0, 0, None,None,None,None, 0, 0, None]
    for key, item in data.items():
        #
        # End 1
        if re.match(r"\b((qx|qa(xial)?)(_)?(0|1|start)?)\b", key, re.IGNORECASE):
            new_data[0] = item.convert("newton/metre").value
            
        elif re.match(r"\b((qy|in(_)?plane)(_)?(0|1|start)?)\b", key, re.IGNORECASE):
            new_data[1] = item.convert("newton/metre").value
            
        elif re.match(r"\b((qz|out(_)?plane)(_)?(0|1|start)?)\b", key, re.IGNORECASE):
            new_data[2] = item.convert("newton/metre").value
            
        elif re.match(r"\b(qt(orsion)?|mx(_)?(0|1|start)?)\b", key, re.IGNORECASE):
            new_data[3] = item.convert("newton/metre").value
        #
        # End 2
        elif re.match(r"\b((qx|qa(xial)?)(_)?(2|end))\b", key, re.IGNORECASE):
            new_data[4] = item.convert("newton/metre").value           
        
        elif re.match(r"\b((qy|in(_)?plane)(_)?(2|end))\b", key, re.IGNORECASE):
            new_data[5] = item.convert("newton/metre").value
        
        elif re.match(r"\b((qz|out(_)?plane)(_)?(2|end))\b", key, re.IGNORECASE):
            new_data[6] = item.convert("newton/metre").value           
        
        elif re.match(r"\b(qt(orsion)?|mx(_)?(2|end)?)\b", key, re.IGNORECASE):
            new_data[7] = item.convert("newton/metre").value            
        #
        # L0, L2
        elif re.match(r"\b((l|d(istance)?)(_)?(0|1|start)?)\b", key, re.IGNORECASE):
            new_data[8] = item.value
        
        elif re.match(r"\b((l|d(istance)?)(_)?(2|end))\b", key, re.IGNORECASE):
            new_data[9] = item.value
        #
        # Comment
        elif re.match(r"\b(title|comment|name|id)\b", key, re.IGNORECASE):
            new_data[10] = item
    #
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
    #
    new_data.insert(0, 'line')
    #
    return new_data
#
#
def get_BeamLine_dicxx(data: dict, ntiems: int=11)->list:
    """
    new_data: [qx1, qy1, qz1, qt1,
               qx2, qy2, qz2, qt2,
               d1, d2,
               title]
    """
    rows = 0
    #new_data = [0, 0, 0, 0, None,None,None,None, 0, 0, None]
    new_data = [None] * ntiems
    qdata = [None] * 8
    Ldata = [None] * 2
    tdata = []
    for key, item in data.items():
        #
        # End 1
        if re.match(r"\b((qx|qa(xial)?)(_)?(0|1|start)?)\b", key, re.IGNORECASE):
            #new_data[0] = item.convert("newton/metre").value
            qdata[0] = get_dict_load(item, units='newton/metre')
            rows = max(rows, len(qdata[0]))
            
        elif re.match(r"\b((qy|in(_)?plane)(_)?(0|1|start)?)\b", key, re.IGNORECASE):
            #new_data[1] = item.convert("newton/metre").value
            qdata[1] = get_dict_load(item, units='newton/metre')
            rows = max(rows, len(qdata[1]))
            
        elif re.match(r"\b((qz|out(_)?plane)(_)?(0|1|start)?)\b", key, re.IGNORECASE):
            #new_data[2] = item.convert("newton/metre").value
            qdata[2] = get_dict_load(item, units='newton/metre')
            rows = max(rows, len(qdata[2]))
            
        elif re.match(r"\b(qt(orsion)?|mx(_)?(0|1|start)?)\b", key, re.IGNORECASE):
            #new_data[3] = item.convert("newton/metre").value
            qdata[3] = get_dict_load(item, units='newton/metre')
            rows = max(rows, len(qdata[3]))
        #
        # End 2
        elif re.match(r"\b((qx|qa(xial)?)(_)?(2|end))\b", key, re.IGNORECASE):
            #new_data[4] = item.convert("newton/metre").value
            qdata[4] = get_dict_load(item, units='newton/metre')
            rows = max(rows, len(qdata[4]))            
        
        elif re.match(r"\b((qy|in(_)?plane)(_)?(2|end))\b", key, re.IGNORECASE):
            #new_data[5] = item.convert("newton/metre").value
            qdata[5] = get_dict_load(item, units='newton/metre')
            rows = max(rows, len(qdata[5]))
        
        elif re.match(r"\b((qz|out(_)?plane)(_)?(2|end))\b", key, re.IGNORECASE):
            #new_data[6] = item.convert("newton/metre").value
            qdata[6] = get_dict_load(item, units='newton/metre')
            rows = max(rows, len(qdata[6]))            
        
        elif re.match(r"\b(qt(orsion)?|mx(_)?(2|end)?)\b", key, re.IGNORECASE):
            #new_data[7] = item.convert("newton/metre").value
            qdata[7] = get_dict_load(item, units='newton/metre')
            rows = max(rows, len(qdata[7]))            
        #
        # L0, L2
        elif re.match(r"\b((l|d(istance)?)(_)?(0|1|start)?)\b", key, re.IGNORECASE):
            Ldata[0] = get_dict_load(item, units='metre')
            #new_data[8] = item.value
        
        elif re.match(r"\b((l|d(istance)?)(_)?(2|end))\b", key, re.IGNORECASE):
            Ldata[1] = get_dict_load(item, units='metre')
            #new_data[9] = item.value
        #
        # Comment
        elif re.match(r"\b(title|comment|name|id)\b", key, re.IGNORECASE):
            #new_data[10] = item
            if isinstance(item, (list, tuple)):
                tdata.append([x for x in item])
            else:
                tdata.append([item])
    #
    if not rows:
        raise IOError('no line load found')
    #
    # do q1
    #
    new_data[:4] = get_items(data=qdata[:4],
                             new_data=new_data[:4],
                             rows=rows)
    #
    # do L
    new_data[8:10] = get_items(data=Ldata,
                               new_data=new_data[8:10],
                               rows=rows)
    #
    # for uniform load
    # qx2
    new_data[4:5] = get_items(data=qdata[4:5],
                              new_data=new_data[4:5],
                              rows=rows, fill=new_data[0])
    #if new_data[4] == None:
    #    new_data[4] = new_data[0]
    # qy2
    new_data[5:6] = get_items(data=qdata[5:6],
                              new_data=new_data[5:6],
                              rows=rows, fill=new_data[1])
    #if new_data[5] == None:
    #    new_data[5] = new_data[1]
    # qz2
    new_data[6:7] = get_items(data=qdata[6:7],
                              new_data=new_data[6:7],
                              rows=rows, fill=new_data[2])
    #if new_data[6] == None:
    #    new_data[6] = new_data[2]
    # qt2
    new_data[7:8] = get_items(data=qdata[7:8],
                              new_data=new_data[7:8],
                              rows=rows, fill=new_data[3])    
    #if new_data[7] == None:
    #    new_data[7] = new_data[3]    
    #
    new_data[10:] = get_items(data=tdata,
                              new_data=new_data[10:],
                              rows=rows, fill=[None]*rows)
    #
    new_data.insert(0, ['line'] * rows)
    #
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