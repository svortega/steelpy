# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from operator import itemgetter
#from typing import NamedTuple
import re

# package imports
#from steelpy.process.units.units import Units
from steelpy.utils.sqlite.utils import create_connection #, create_table
import numpy as np


class BasicProperty:
    __slots__ = ['name', '_groups', '_elements']
    
    def __init__(self):
        """
        """
        #self._units = Units()
        self._groups:List = []
        self._elements:List = []
        #self._default:bool = False
    
    #
    @property
    def group(self):
        """
        """
        return self._groups
    
    @group.setter
    def group(self, value:list|tuple|str|int):
        """
        """
        if isinstance(value, (list, tuple)):
            self._groups.extend(list(value))
        else:
            self._groups.append(value)
    #
    #
    @property
    def elements(self):
        """
        """
        return self._elements
    
    @elements.setter
    def elements(self, value:list|tuple|str|int):
        """
        """
        if isinstance(value, (list, tuple)):
            self._elements.extend(list(value))
        else:
            self._elements.append(value)
    #
    #def set_default(self):
    #    """
    #    """
    #    self._default = True
#
#
class HydroBasic(Mapping):
    
    __slots__ = ['db_file']
    
    #
    def __init__(self, db_file: str):
        """
        """
        self.db_file = db_file
        #
        # create table
        conn = create_connection(self.db_file)
        with conn:
            self._create_table(conn)         
    #
    #
    def __setitem__(self, name:str|int,
                    value: list|tuple|dict|str) -> None:
        """
        """
        #
        data = (name, 'profile', None)
        conn = create_connection(self.db_file)
        with conn:
            self._push_data(conn, data)
        #
        item = self.__getitem__(name)
        #
        if isinstance(value, list|tuple):
            if isinstance(value[0], list|tuple):
                item.profile = value
            else:
                profile = self.get_profile(value)
                item.profile = profile
        
        elif isinstance(value, dict):
            1 / 0
            
        elif isinstance(value, str):
            1 / 0
            data = self.get_data([value])
        else:
            raise IOError('Input data not valid')
    
    #    
    #
    def __len__(self) -> int:
        return len(self._label)

    def __iter__(self):
        """
        """
        return iter(self._label)
    
    def __contains__(self, value) -> bool:
        return value in self._label
    #
    #
    #
    # ------------------
    #
    def get_profile(self, values, step: int = 10):
        """ [constant/profile, water_depth, thickness, step, title] """
        outval = ['constant', None, None, 10, None]
        #mgprofile = outval[2]
        #
        if isinstance(values, dict):
            1/0
        else:
            #if isinstance(values[-1], str):
            #    mgprofile = values.pop()
            #
            #outval[0] = values[0]
            #
            #try:
            #    outval[1] = values[1].convert('kilogram/metre^3').value
            #except IndexError:
            #    pass
            #
            #try:
            #    outval[3] = values[2].value
            #    #mgprofile = 'uniform'
            #    outval[2] = 'uniform'
            #except IndexError:
            #    #outval[2] = 'profile'
            #    pass
            #
            #
            if re.match(r"\b(constant|linear)\b", values[0], re.IGNORECASE):
                depth = values[1]
                factor = [values[2]]
                #
                try:
                    factor.extend([values[3]])
                except IndexError:
                    pass
                #
                profile = get_profile(elevation=depth, value=factor,
                            steps=step)
                #
                return profile
        #
        #if re.match(r"\b(constant|uniform|mean)\b", mgprofile, re.IGNORECASE):
        #    outval[2] = 'uniform'
        #
        #return outval
    #    
#
#
def get_list(data, steps:int=6)->list[float]:
    """ """
    new_data = []
    for x in range(steps):
        try:
            try:
                new_data.append(data[x].value)
            except AttributeError:
                new_data.append(data[x])
        except IndexError:
            new_data.append(0.0)
    return new_data
#
#
def get_profile(elevation: float, value: float, steps: int):
    """ """
    levels = [item /steps  for item in range(steps, -1, -1)]
    levels.extend([-1 * item for item in reversed(levels[:-1])])
    profile = [[item * elevation, *value] for item in levels]
    return profile
#
#
@dataclass
class HydroItem:
    __slots__ = ['name', '_db_file']
    #
    def __init__(self, name: int|str,
                 db_file: str):
        """ """
        self.name = name
        self._db_file = db_file
    #
    #
    @property
    def profile(self):
        """ """
        #print('-->')
        conn = create_connection(self._db_file)
        with conn:        
            prof = self._pull_profile(conn)
        #1/0
        return prof

    @profile.setter
    def profile(self, value: list[list]):
        """ """
        prof = [[item[0].value, *item[1:]] for item in value]
        #for item in value:
        #    elev = item[0].value
        #    factor = item[1:]
        #    prof.append([elev, factor])
        #
        # Converting velocity to SI units
        try:
            prof = [[item[0], item[1].value] for item in prof]
        except AttributeError:
            pass        
        #        
        #
        prof.sort(key=itemgetter(0), reverse=True)
        #
        conn = create_connection(self._db_file)
        with conn:        
            self._push_profile(conn, profile_data=prof)        
    #
    #
    def get_profile(self, Z):
        """ """
        prof = self.profile
        elev = list(reversed([item[0] for item in prof]))
        value = list(reversed([item[1] for item in prof]))
        items = np.zeros((Z.shape))
        for i, point in enumerate(Z):
            items[i] = np.interp(point, elev, value)
        #
        return items
    #
    #