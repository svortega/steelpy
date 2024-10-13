#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
#from array import array
from typing import NamedTuple
from collections.abc import Mapping
from collections import defaultdict #, Counter
import re

# package imports
from .utils import (get_value_point,
                    check_list_number)

import steelpy.utils.io_module.text as common
from steelpy.utils.dataframe.main import DBframework
#
# ---------------------------------
#
class PointNode(NamedTuple):
    """
    [fx, fy, fz, mx, my, mz, name, title, load_name, system, load_complex, load_type]
    """
    fx: float
    fy: float
    fz: float
    mx: float
    my: float
    mz: float
    name: int|str
    title: str
    load_name: int|str
    system: int
    load_complex: int
    load_type:str = "force"

    #
    def __str__(self, units: str = "si") -> str:
        """ """
        output = (f"{str(self.name):12s} {10 * ' '} "
                  f"{self.fx: 1.3e} {self.fy: 1.3e} {self.fz: 1.3e} "
                  f"{self.coordinate_system.upper():6s} "
                  f"{self.load_complex}\n")
                  
        #
        if (load_title := self.title) == "NULL":
            load_title = ""        
        step = self.load_type
        output += (f"{step:12s} {10 * ' '} "
                   f"{self.mx: 1.3e} {self.my: 1.3e} {self.mz: 1.3e} "
                   f"{str(load_title)}\n")
        return output

    #
    @property
    def coordinate_system(self) -> str:
        if self.system != 0:
            return "local"
        return "global"
    #


#
# ---------------------------------
#
class DispNode(NamedTuple):
    """
    """
    x: float
    y: float
    z: float
    rx: float
    ry: float
    rz: float
    name: int|str
    title: str
    load_name: int|str
    system: int
    load_complex: int
    load_type:str = "Nodal Displacement"
    #
    def __str__(self, units:str="si") -> str:
        """ """
        output = (f"{str(self.name):12s} {10 * ' '} "
                  f"{self.x: 1.3e} {self.y: 1.3e} {self.z: 1.3e} "
                  f"{self.coordinate_system.upper():6s} "
                  f"{self.load_complex}\n")
                  
        #
        if (load_title := self.title) == "NULL":
            load_title = ""        
        step = self.load_type
        output += (f"{step:12s} {10 * ' '} "
                   f"{self.rx: 1.3e} {self.ry: 1.3e} {self.rz: 1.3e} "
                   f"{str(load_title)}\n")
        return output
    #
    #def Cpmt(self):
    #    """ Penalty Method technique"""
    #    factor = self.Cfactor()
    #    
    #    disp = [self.x, self.y, self.z,
    #            self.rx, self.ry, self.rz]
    #    disp = [item *  factor[x]
    #            for x, item in enumerate(disp)]
    #    return disp
    #
    #def Cfactor(self):
    #    """ """
    #    Ktranslation = 10**10
    #    Krotation = 10**12
    #    # x, y, z
    #    translation = [Ktranslation] * 3
    #    # rx, ry, rz
    #    rotation = [Krotation] * 3
    #    factor = translation + rotation
    #    return factor 
    #
    #def Kp(self, Kg):
    #    """ """
    #    kp =  np.zeros((12, 12))
    #    stiffness = Kg.max()
    #    factor = self.Cfactor(stiffness)
    #    
    #
    @property
    def coordinate_system(self):
        if self.system != 0:
            return "local"
        return "global" 
#
# ---------------------------------
#
class NodeLoadBasic(Mapping):
    __slots__ = ['_system_flag']

    def __init__(self) -> None:
        """
        """
        # 0-global/ 1-local
        self._system_flag: int = 0 # Global system
    #
    #
    def __iter__(self):
        """
        """
        items = list(set(dict.fromkeys(self._labels)))
        return iter(items)

    def __contains__(self, value) -> bool:
        return value in self._labels

    def __len__(self) -> int:
        items = list(set(dict.fromkeys(self._labels)))
        return iter(items)

    #
    def __str__(self) -> str:
        """ """
        output = ""
        nodes = list(dict.fromkeys(self._labels))
        #nodes = set(self._labels)
        for node in nodes:
            items = self.__getitem__(node)
            for item in items:
                output += item.__str__()
                # print('---')
        return output
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
# ---------------------------------
#
#
def find_NodeLoad_item(word_in:str) -> str:
    """ """
    key = {"name": r"\b(id|name|load(s)?)\b",
           "type": r"\b((load(_|-|\s*)?)?type)\b",
           "title": r"\b(title|comment)\b",
           "node": r"\b(node(s)?(_|-|\s*)?(name|id)?)\b"}
           #"system": r"\b(system)\b",}
    try:
        match = common.find_keyword(word_in, key)
        return match
    except IOError:
        return find_load_item(word_in)
#
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
        return find_force_item(word_in)
    except IOError:
        pass
    
    try:
        return find_disp_item(word_in)
    except IOError:
        pass
    
    try:
        return find_mass_item(word_in)
    except IOError:
        pass
    #
    return IOError('load definition')
#
def find_force_item(word_in:str) -> str:
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
def find_disp_item(word_in:str) -> str:
    """ """
    key = {"x": r"\b((delta(_|-|\s*)?)?x)\b",
           "y": r"\b((delta(_|-|\s*)?)?y)\b",
           "z": r"\b(\b((delta(_|-|\s*)?)?z)\b",
           #
           "rx": r"\b(rx|t(orsion)?)\b",
           "ry": r"\b(ry|out(_|-|\s*)?plane)\b",
           "rz": r"\b(rz|in(_|-|\s*)?plane)\b"}
    
    match = common.find_keyword(word_in, key)
    return match
#
def find_mass_item(word_in:str) -> str:
    """ """
    key = {"x": r"\b((m(_|-|\s*)?)?x)\b",
           "y": r"\b((m(_|-|\s*)?)?y)\b",
           "z": r"\b(\b((m(_|-|\s*)?)?z)\b",
           #
           "rx": r"\b(mx|t(orsion)?)\b",
           "ry": r"\b(my|out(_|-|\s*)?plane)\b",
           "rz": r"\b(mz|in(_|-|\s*)?plane)\b"}
    
    match = common.find_keyword(word_in, key)
    return match
#
# ---------------------------------
#
#
#
def get_NodaLoad_df(df: DBframework.DataFrame,
                    system: int = 0):
    """ """
    #
    flag_system = 'global'
    if system == 1:
        flag_system = 'local'     
    #
    columns = list(df.columns)
    header = {key: find_NodeLoad_item(key)
              for key in columns}
    #
    #
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
                force = get_NodeForce_dic(data)
                temp.append([data['name'], data['system'], key[1], *force, None])
            newgrp['force'].extend(temp)
        
        elif re.match(r"\b(displacement)\b", key[0], re.IGNORECASE):
            temp = []
            for load in nload['data']:
                #  [displacement, x, y, z, rx, ry, rz, title]
                data = dict(zip(nload['columns'], load))
                force = get_NodeDisp_dic(data)
                temp.append([data['name'], data['system'], key[1], *force, None])
            newgrp['displacement'].extend(temp)
        
        elif re.match(r"\b(mass)\b", key[0], re.IGNORECASE):
            temp = []
            for load in nload['data']:
                #  [mass, x, y, z, rx, ry, rz, title]
                data = dict(zip(nload['columns'], load))
                force = get_NodeMass_dic(data)
                temp.append([data['name'], data['system'], key[1], *force, None])
            newgrp['mass'].extend(temp)
        
        else:
            raise IOError(f'node load type {key[0]} not available')        
    #
    #
    db = DBframework()
    header = {'force': ['name', 'system',
                        'node', 'type', 
                        'fx', 'fy', 'fz', 'mx', 'my', 'mz',
                        'title', 'step'],
              'displacement': ['name', 'system',
                               'node', 'type', 
                               'x', 'y', 'z', 'rx', 'ry', 'rz',
                               'title', 'step'],
              'mass': ['name', 'system',
                       'node', 'type', 
                       'x', 'y', 'z', 'rx', 'ry', 'rz',
                       'title', 'step']}
    newdf = {}
    for key, item in newgrp.items():
        newdf[key] = db.DataFrame(data=item, columns=header[key])
    #
    return newdf
#
#
def get_nodal_load(load: list|tuple|dict)->list:
    """
    force = [Fx, Fy, Fz, Mx, My, Mz, title]
    displacement = [x, y, z, rx, ry, rz, title]
    mass = [x, y, z, mx, my, mz, title]
    """
    if isinstance(load, (list, tuple)):
        try:
            load = get_NodeLoad_list_units(load.copy())
        except AttributeError:
            load = get_NodeLoad_list(load)
    elif isinstance(load, dict):
        load = get_NodeLoad_dic(load)
    else:
        raise Exception('   *** Load input format not recognized')
    return load
#
def get_NodeLoad_list(data: list|tuple,
                      steps: int = 6) -> list :
    """
    force = [Fx, Fy, Fz, Mx, My, Mz, title]
    displacement = [x, y, z, rx, ry, rz, title]
    mass = [x, y, z, mx, my, mz, title]
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
    force = [Fx, Fy, Fz, Mx, My, Mz, title]
    displacement = [x, y, z, rx, ry, rz, title]
    mass = [x, y, z, mx, my, mz, title]
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
def get_NodeLoad_dic(data: dict)->list:
    """
    force : [fx, fy, fz, mx, my, mz, title]
    displacement : [x, y, z, rx, ry, rz, title]
    mass = : [x, y, z, rx, ry, rz, title]
    """
    new_data = [None, 0,0,0, 0,0,0, None]
    loat_type = data['type']
    
    if re.match(r"\b((prescribed)?disp(lacement)?)\b", loat_type, re.IGNORECASE):
        new_data = get_NodeDisp_dic(data)
        #new_data[0] = 'displacement'
        #for key, item in data.items():
        #    
        #    if re.match(r"\b(x)\b", str(key), re.IGNORECASE):
        #        new_data[1] = item.value
        #        
        #    elif re.match(r"\b(y)\b", str(key), re.IGNORECASE):
        #        new_data[2] = item.value
        #        
        #    elif re.match(r"\b(z)\b", str(key), re.IGNORECASE):
        #        new_data[3] = item.value
        #    
        #    elif re.match(r"\b(rx)\b", str(key), re.IGNORECASE):
        #        new_data[4] = item.value
        #        
        #    elif re.match(r"\b(ry)\b", str(key), re.IGNORECASE):
        #        new_data[5] = item.value
        #        
        #    elif re.match(r"\b(rz)\b", str(key), re.IGNORECASE):
        #        new_data[6] = item.value
        #        
        #    elif re.match(r"\b(title|comment|name|id)\b", key, re.IGNORECASE):
        #        new_data[7] = item            
    
    elif re.match(r"\b(force|load)\b", loat_type, re.IGNORECASE):
        new_data = get_NodeForce_dic(data)
        #new_data[0] = 'load'
        #
        #for key, item in data.items():
        #    
        #    if re.match(r"\b(fx|fa(xial)?)\b", key, re.IGNORECASE):
        #        new_data[1] = item.convert("newton").value
        #        
        #    elif re.match(r"\b(py|fy|in(_)?plane)\b", key, re.IGNORECASE):
        #        new_data[2] = item.convert("newton").value
        #        
        #    elif re.match(r"\b(pz|fz|out(_)?plane)\b", key, re.IGNORECASE):
        #        new_data[3] = item.convert("newton").value
        #    # 
        #    elif re.match(r"\b(mx|t(orsion)?)\b", key, re.IGNORECASE):
        #        new_data[4] = item.convert("newton*metre").value
        #        
        #    elif re.match(r"\b(my|out(_)?plane)\b", key, re.IGNORECASE):
        #        new_data[5] = item.convert("newton*metre").value
        #        
        #    elif re.match(r"\b(mz|in(_)?plane)\b", key, re.IGNORECASE):
        #        new_data[6] = item.convert("newton*metre").value
        #    
        #    elif re.match(r"\b(title|comment|name|id)\b", key, re.IGNORECASE):
        #        new_data[7] = item            
    
    elif re.match(r"\b(mass)\b", loat_type, re.IGNORECASE):
        new_data = get_NodeMass_dic(data)
        #new_data[0] = 'mass'
        #for key, item in data.items():
        #    if re.match(r"\b(x)\b", key, re.IGNORECASE):
        #        new_data[1] = item.convert("newton").value
        #        
        #    elif re.match(r"\b(y)\b", key, re.IGNORECASE):
        #        new_data[2] = item.convert("newton").value
        #        
        #    elif re.match(r"\b(z)\b", key, re.IGNORECASE):
        #        new_data[3] = item.convert("newton").value
        #    #
        #    elif re.match(r"\b(rx?)\b", key, re.IGNORECASE):
        #        new_data[4] = item.convert("newton*metre").value
        #        
        #    elif re.match(r"\b(ry)\b", key, re.IGNORECASE):
        #        new_data[5] = item.convert("newton*metre").value
        #        
        #    elif re.match(r"\b(rz)\b", key, re.IGNORECASE):
        #        new_data[6] = item.convert("newton*metre").value            
        #    
        #    elif re.match(r"\b(title|comment|name|id)\b", key, re.IGNORECASE):
        #        new_data[7] = item            
    
    else:
        raise IOError(f'node load type {loat_type} not available')
    #        
    return new_data
#
def get_NodeForce_dic(data: dict)->list:
    """ force : [fx, fy, fz, mx, my, mz, title] """
    new_data = ['force', 0,0,0, 0,0,0, None]
    #
    for key, item in data.items():
        if re.match(r"\b(fx|fa(xial)?)\b", key, re.IGNORECASE):
            new_data[1] = item.convert("newton").value
            
        elif re.match(r"\b(py|fy|in(_)?plane)\b", key, re.IGNORECASE):
            new_data[2] = item.convert("newton").value
            
        elif re.match(r"\b(pz|fz|out(_)?plane)\b", key, re.IGNORECASE):
            new_data[3] = item.convert("newton").value
        # 
        elif re.match(r"\b(mx|t(orsion)?)\b", key, re.IGNORECASE):
            new_data[4] = item.convert("newton*metre").value
            
        elif re.match(r"\b(my|out(_)?plane)\b", key, re.IGNORECASE):
            new_data[5] = item.convert("newton*metre").value
            
        elif re.match(r"\b(mz|in(_)?plane)\b", key, re.IGNORECASE):
            new_data[6] = item.convert("newton*metre").value
        
        elif re.match(r"\b(title|comment|name|id)\b", key, re.IGNORECASE):
            new_data[7] = item
    #
    return new_data
#
def get_NodeDisp_dic(data: dict)->list:
    """ displacement : [x, y, z, rx, ry, rz, title] """
    new_data = ['displacement', 0,0,0, 0,0,0, None]
    #
    for key, item in data.items():
        if re.match(r"\b(x)\b", str(key), re.IGNORECASE):
            new_data[1] = item.value
            
        elif re.match(r"\b(y)\b", str(key), re.IGNORECASE):
            new_data[2] = item.value
            
        elif re.match(r"\b(z)\b", str(key), re.IGNORECASE):
            new_data[3] = item.value
        
        elif re.match(r"\b(rx)\b", str(key), re.IGNORECASE):
            new_data[4] = item.value
            
        elif re.match(r"\b(ry)\b", str(key), re.IGNORECASE):
            new_data[5] = item.value
            
        elif re.match(r"\b(rz)\b", str(key), re.IGNORECASE):
            new_data[6] = item.value
            
        elif re.match(r"\b(title|comment|name|id)\b", key, re.IGNORECASE):
            new_data[7] = item
    #
    return new_data
#
def get_NodeMass_dic(data: dict)->list:
    """ mass = : [x, y, z, rx, ry, rz, title] """
    new_data = ['mass', 0,0,0, 0,0,0, None]
    #
    for key, item in data.items():
        if re.match(r"\b(x)\b", key, re.IGNORECASE):
            new_data[1] = item.convert("newton").value
            
        elif re.match(r"\b(y)\b", key, re.IGNORECASE):
            new_data[2] = item.convert("newton").value
            
        elif re.match(r"\b(z)\b", key, re.IGNORECASE):
            new_data[3] = item.convert("newton").value
        #
        elif re.match(r"\b(rx?)\b", key, re.IGNORECASE):
            new_data[4] = item.convert("newton*metre").value
            
        elif re.match(r"\b(ry)\b", key, re.IGNORECASE):
            new_data[5] = item.convert("newton*metre").value
            
        elif re.match(r"\b(rz)\b", key, re.IGNORECASE):
            new_data[6] = item.convert("newton*metre").value            
        
        elif re.match(r"\b(title|comment|name|id)\b", key, re.IGNORECASE):
            new_data[7] = item
    #
    return new_data
#