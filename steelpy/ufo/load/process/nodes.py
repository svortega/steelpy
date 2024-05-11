#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from array import array
from typing import NamedTuple
from collections.abc import Mapping
from collections import defaultdict #, Counter
import re

# package imports
#import pandas as pd
# steelpy.f2uModel.load.process
from .operations import (get_value_point,
                         check_list_number)
                         #check_point_dic)
 
from steelpy.utils.dataframe.main import DBframework
#import numpy as np
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
    load_type:str = "Nodal Load"

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
    __slots__ = ['_labels', '_title', '_complex',
                 '_load_id', '_system']

    def __init__(self) -> None:
        """
        """
        self._labels: list = []
        self._title: list = []
        self._complex: array = array("I", [])
        self._load_id: list = []
        # 0-global/ 1-local
        self._system: array = array("I", [])

    def __iter__(self):
        """
        """
        items = list(set(dict.fromkeys(self._labels)))
        #items = list(dict.fromkeys(self._labels))
        #items = set(self._labels)
        return iter(items)

    def __contains__(self, value) -> bool:
        return value in self._labels

    def __len__(self) -> int:
        items = list(set(dict.fromkeys(self._labels)))
        #return len(self._labels)
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
#
class NodeLoadMaster(NodeLoadBasic):
    """
    FE Node Load class
    
    NodeLoad
        |_ name
        |_ number
        |_ type
        |_ complex
        |_ point [x, y, z, mx, my, mz]
        |_ acceleration [th[0], th[1], th[2],..., th[n]]
        |_ displacement [th[0], th[1], th[2],..., th[n]]
        |_ mass
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    __slots__ = ['_labels', '_title', '_complex', '_load_id', '_system',
                 '_type', '_fx', '_fy', '_fz', '_mx', '_my', '_mz']

    def __init__(self, load_type:str) -> None:
        """
        """
        super().__init__()
        # real
        self._type = load_type
        self._fx: array = array('f', [])
        self._fy: array = array('f', [])
        self._fz: array = array('f', [])
        self._mx: array = array('f', [])
        self._my: array = array('f', [])
        self._mz: array = array('f', [])
        #self._distance: array = array('f', [])
    #
    #
    def __getitem__(self, node_name: int|str) -> list:
        """
        """
        idx_list: list = [x for x, _item in enumerate(self._labels)
                          if _item == node_name]
        #
        points: list = []
        if self._type in ['load']:
            for idx in idx_list:
                points.append(PointNode(self._fx[idx], self._fy[idx], self._fz[idx],
                                        self._mx[idx], self._my[idx], self._mz[idx],
                                        self._labels[idx], self._title[idx], self._load_id[idx],
                                        self._system[idx], self._complex[idx], self._type))
        
        elif self._type in ['displacement']:
            for idx in idx_list:
                points.append(DispNode(self._fx[idx], self._fy[idx], self._fz[idx],
                                       self._mx[idx], self._my[idx], self._mz[idx],
                                       self._labels[idx], self._title[idx], self._load_id[idx],
                                       self._system[idx], self._complex[idx], self._type))
        
        elif self._type in ['mass']:
            raise NotImplementedError()
        
        else:
            raise IOError(f'load type: {self._type} not valid')
            
        return points
    #
    #
    def __delitem__(self, node_name: int|str) -> None:
        """
        """
        indexes = [i for i, x in enumerate(self._labels)
                   if x == node_name]
        indexes.sort(reverse=True)
        for _index in indexes:
            self._fx.pop(_index)
            self._fy.pop(_index)
            self._fz.pop(_index)
            self._mx.pop(_index)
            self._my.pop(_index)
            self._mz.pop(_index)
            #
            self._labels.pop(_index)
            self._title.pop(_index)
            self._system.pop(_index)
            #self._distance.pop(_index)
            self._complex.pop(_index)
    #
    #
    @property
    def df(self):
        """nodes in dataframe format"""
        db = DBframework()
        # TODO : merge nodes
        #title = []
        data = {'load_name': self._load_id,
                'load_type': ['basic' for item in self._labels],
                'load_id': [idx + 1 for idx, item in enumerate(self._labels)],
                'load_system': self._system, #['global' if item == 0 else 'local'
                           #for item in self._system],
                'load_comment': self._title,
                'node_name':self._labels, 
                'Fx':self._fx, 'Fy':self._fy, 'Fz':self._fz, 
                'Mx':self._mx, 'My':self._my, 'Mz':self._mz}
        #
        # beam point case
        try:
            data.update({'L1': self._L1})
        except AttributeError:
            pass
        #
        # beam node load conversion
        try:
            data.update({'element_name': self._beam})
        except AttributeError:
            pass
        #      
        #
        df_nload = db.DataFrame(data=data, index=None)
        #df = df[['load_name', 'load_type', 'load_id', 'load_system', 'load_comment',
        #         'node_name', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']]          
        return df_nload

    @df.setter
    def df(self, values):
        """nodes in dataframe format"""
        #update
        self._labels.extend(values.node_name.tolist())
        self._load_id.extend(values.load_name.tolist())
        self._title.extend(values.load_title.tolist())
        self._system.extend([1 if item in ['local', 'member', 1] else 0
                             for item in values.system.tolist()])
        self._complex.extend([0 for item in values.system.tolist()])
        #
        self._fx.extend(values.Fx.tolist())
        self._fy.extend(values.Fy.tolist())
        self._fz.extend(values.Fz.tolist())
        self._mx.extend(values.Mx.tolist())
        self._my.extend(values.My.tolist())
        self._mz.extend(values.Mz.tolist())
        #
        # beam point case
        try:
            self._L1.extend(values.L1.tolist())
        except AttributeError:
            pass
        #
        # beam node load conversion
        try:
            self._beam.extend(values.element_name.tolist())
        except AttributeError:
            pass
        #
        #print('nodes df out')
    
#
# ---------------------------------
#
def get_nodal_load(load: list|tuple|dict)->list:
    """ """
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
    froce = [Fx, Fy, Fz, Mx, My, Mz, title]
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
#
def get_NodeLoad_list_units(data: list|tuple) -> list :
    """
    froce = [Fx, Fy, Fz, Mx, My, Mz, title]
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
        
        if 'rad' in load:
            if not new_data:
                new_data = [0,0,0]
            new_data.extend(get_value_point(load, label='rad', steps=3))        
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
        # point and beam point load [Mx, My, Mz]
        if 'N*m' in load:
            if not new_data:
                new_data = [0,0,0]
            new_data.extend(get_value_point(load, label='N*m', steps=3))
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
        new_data[0] = 'displacement'
        
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
    
    elif re.match(r"\b(force|load)\b", loat_type, re.IGNORECASE):
        new_data[0] = 'load'
        
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
    
    elif re.match(r"\b(mass)\b", loat_type, re.IGNORECASE):
        new_data[0] = 'mass'
        
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
    
    else:
        raise IOError(f'node load type {loat_type} not available')
    #        
    return new_data