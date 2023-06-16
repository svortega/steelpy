# 
# Copyright (c) 2019-2023 steelpy
#
# 
# Python stdlib imports
from __future__ import annotations

#from openpyxl.chart import title

#from array import array
from collections.abc import Mapping
#from collections import defaultdict
from dataclasses import dataclass
from typing import NamedTuple
#from math import prod

# package imports
# steelpy.f2uModel.load
#import pandas as pd
from steelpy.process.dataframe.main import DBframework
from steelpy.f2uModel.load.process.combination import LoadCombinationBasic

#
#
#
class MetoceanCombination(Mapping):
    """
    FE Metocean Combination Class
    
    Combination
        |_ name
        |_ number
        |_ surface
        |_ gravity
        |_ water_density
        |_ air_density
        |_ fluid_viscosity
        |_ air_viscosity
        |_ zone
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    __slots__ = ['_combination']

    def __init__(self):
        """
        """
        self._combination:dict = {}
    
    #
    #
    def __getitem__(self, comb_name:str|int):
        """
        """
        return self._combination[comb_name]
    
    
    def __setitem__(self, comb_name: str|int, comb_title: str) -> None:
        """
        """
        self._combination[comb_name] = CombTypes(comb_name, comb_title)
    #
    #
    def __delitem__(self, load_name:str|int):
        """
        """
        del self._combination[load_name]
    #
    def __len__(self) -> int:
        return len(self._combination)
    #
    def __iter__(self):
        """
        """
        return iter(self._combination)
#
@dataclass
class Combination_Metocean:
    """
    """
    #__slots__ = []
    name:str|int
    #number:int
    wave:str|int = 0
    wave_direction:float = 0.0
    wave_kinematics:float = 0.0
    #
    buoyancy : bool = True
    #
    current: str|int = 0
    current_blockage:float  = 0.0
    current_stretching : bool = True
    current_direction:float = 0.0
    #
    wind:str|int = 0
    wind_direction:float = 0.0
    #wind_areas:list = []
#
#
class WaveBasic(NamedTuple):
    #number:int
    wave:str|int
    direction:float
    kinematics:float
    title:str
#
class CurrentBasic(NamedTuple):
    #number:int
    current: str|int
    direction:float
    blockage:float
    stretching : bool
    title:str
#
class WindBasic(NamedTuple):
    #number:int
    wind:str|int
    direction:float
    #wind_areas:list = []
    title:str
#
#
#
class CombTypes:
    """
    """
    __slots__ = ['_wave', '_current', '_wind', 
                 'name', 'number', 'title']
    
    def __init__(self, name:int, title:str):
        """
        """
        self.name = name
        self.title = title
        #
        self._wave = None
        self._current = None
        self._wind = None
    #
    @property
    def wave(self):
        """
        """
        return self._wave
    
    @wave.setter
    def wave(self, values: list|dict):
        """
        [wave_name, Direction(deg), Kinematics, title]
        """
        values, title = get_values(values)
        self._wave = WaveBasic(wave=values[0],
                               direction=values[1],
                               kinematics=values[2],
                               title=title)
    #
    @property
    def current(self):
        """
        """
        return self._current
    
    @current.setter
    def current(self, values: list|dict):
        """
        """
        values, title = get_values(values)
        self._current =  CurrentBasic(current=values[0],
                                      direction=values[1],
                                      blockage=values[2],
                                      stretching=values[3],
                                      title=title)
    #
    @property
    def wind(self):
        """
        """
        return self._wind
    
    @wind.setter
    def wind(self, values: list|dict):
        """
        """
        values, title = get_values(values)
        self._wind = WindBasic(wind=values[0],
                               direction=values[1],
                               title=title)
#
#
def get_values(values):
    """ """
    if isinstance(values, (list|tuple)):
        if isinstance(values[-1], str):
            title = values.pop(-1)
        else:
            title = ""
        
    elif isinstance(values, dict):
        pass
    else:
        raise IOError('Input data not valid')
    #
    return values, title
#
#