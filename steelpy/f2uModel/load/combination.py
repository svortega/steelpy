# 
# Copyright (c) 2019-2020 steelpy
#
# 
# Python stdlib imports
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Union, Dict

# package imports


#
#
class CombTypes:
    """
    """
    __slots__ = ['_basic', '_metocean', 'name', 
                 'number', '_level']
    
    def __init__(self):
        """
        """
        self._basic:Dict = BasicLoad("basic")
        self._metocean = BasicLoad("metocean")
        self._level = 1
    #
    @property
    def metocean(self):
        """
        """
        return self._metocean
    #
    @property
    def basic_load(self):
        """
        """
        return self._basic
    #
    @property
    def level(self):
        """
        """
        return self._level
#
#
class LoadCombination(Mapping):
    """
    FE Load Combination Class
    
    LoadCombination
        |_ name
        |_ number
        |_ design_condition : operating/storm
        |_ functional_load [load number, factor]
        |_ wave_load [load number, factor]
        |_ wind_load [load number, factor]
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    #
    __slots__ = ('number', 'name', '_combination')
    #
    def __init__(self):
        """
        """
        self._combination = {}
    
    def __setitem__(self, load_name:Union[str,int], basic_load:int) -> None:
        """
        """
        #try:
        #    self._combination[load_name].append(basic_load)
        #except KeyError:
        #    self._combination[load_name] = [basic_load]
        self._combination[load_name] = CombTypes()
        self._combination[load_name].name = basic_load
        self._combination[load_name].number = load_name
    
    def __getitem__(self, load_name:Union[str,int]):
        """
        """
        return self._combination[load_name]
    #
    def __delitem__(self, load_name:Union[str,int]):
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
#
class BasicLoad(Mapping):
    """
    FE Metocean Combination Class
    """
    __slots__ = ['_basic', '_type']

    def __init__(self, bl_type:str):
        """
        """
        self._basic:Dict = {}
        self._type = bl_type
    #
    def __getitem__(self, load_name:Union[str,int]):
        """
        """
        return self._basic[load_name]
    
    def __setitem__(self, load_name: Union[str,int], factor: float) -> None:
        """
        """
        self._basic[load_name] = factor
    #
    #
    def __delitem__(self, load_name:Union[str,int]):
        """
        """
        del self._basic[load_name]
    #
    def __len__(self) -> int:
        return len(self._basic)
    #
    def __iter__(self):
        """
        """
        return iter(self._basic)
    #
    #
    @property
    def load_type(self):
        """
        """
        return self._type
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
        self._combination:Dict = {}
    
    #
    #
    def __getitem__(self, comb_name:Union[str,int]):
        """
        """
        return self._combination[comb_name]
    
    
    def __setitem__(self, comb_name: Union[str,int], comb_title: str) -> None:
        """
        """
        self._combination[comb_name] = Combination_Metocean(comb_title)
    #
    #
    def __delitem__(self, load_name:Union[str,int]):
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
    name:Union[str,int]
    #number:int
    wave:Union[str,int] = 0
    wave_direction:float = 0.0
    wave_kinematics:float = 0.0
    #
    buoyancy : bool = True
    #
    current: Union[str,int] = 0
    current_blockage:float  = 0.0
    current_stretching : bool = True
    current_direction:float = 0.0
    #
    wind:Union[str,int] = 0
    wind_direction:float = 0.0
    #wind_areas:list = []
#