# 
# Copyright (c) 2019-2021 steelpy
#
# 
# Python stdlib imports
#from array import array
from collections.abc import Mapping
#from collections import defaultdict
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Union, Dict, List, Union
from math import prod

# package imports
from steelpy.f2uModel.load.operations.combination import LoadCombinationBasic

#
#
class CombTypes:
    """
    """
    __slots__ = ['_basic', '_metocean', 'name', 'number',
                 'title', '_combination']
    
    def __init__(self, name:int, title:str):
        """
        """
        self.name = name
        self.title = title
        #
        self._basic:Dict = BasicLoad("basic")
        self._metocean = BasicLoad("metocean")
        self._combination = BasicLoad("combination")
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
    def load_combination(self):
        """
        """
        return self._combination
#
class LoadCombination(LoadCombinationBasic):
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
    __slots__ = ('number', 'name', '_combination', '_basic')
    #
    def __init__(self, basic_load):
        """
        """
        super().__init__()
        self._combination = {}
        self._basic = basic_load
    
    def __setitem__(self, load_name:int, load_title:str) -> None:
        """
        """
        try:
            self._labels.index(load_name)
            raise Exception('    *** warning load combination name {:} already exist'
                            .format(load_name))
        except ValueError:
            self._labels.append(load_name)
            self._title.append(load_title)
            #
            self._combination[load_name] = CombTypes(name=load_name, title=load_title)
            # TODO: fix numbering sequence
            load_number =  len(self._combination)
            self._combination[load_name].number = load_number
            self._number.append(load_number)
    
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
    def to_basic(self):
        """ """
        # get load combination and convert to basic loads
        for key, item in self._combination.items():
            for comb_name, factor in item._combination.items():
                for precomb, prefactor in self._combination[comb_name]._basic.items():
                    try:
                        item._basic[precomb] += prefactor * factor
                    except KeyError:
                        item._basic[precomb] = prefactor * factor
        # organize basic load
        basic_loads = defaultdict(list)
        for key, item in self._combination.items():
            for bname, factor in item._basic.items():
                try:
                    basic = self._basic[bname]
                    basic_loads[item.name].append([basic.title, factor])
                except KeyError:
                    raise Warning("  warning basic load {:} not found".format(bname))
        #print('-->')
        return basic_loads
    #
    def solve_combinations(self, basic_res, memb_force):
        """
        """
        comb_res = {}
        memb_comb = {}
        bloads = self.to_basic()
        for cname, comb in bloads.items():
            lcomb = self._combination[cname].title
            beam_load = {}
            for bname, factor in comb:
                try:
                    comb_res[lcomb] += basic_res[bname] * factor
                except KeyError:
                    comb_res[lcomb] = basic_res[bname] * factor
                #
                for mname, member in memb_force[bname].items():
                    try:
                        beam_load[mname] += member.end_force * factor
                    except KeyError:
                        beam_load[mname] =  member.end_force * factor
            memb_comb[lcomb] = beam_load
        return comb_res, memb_comb
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
    def __setitem__(self, load_name: Union[str,int], factor: float) -> None:
        """
        """
        self._basic[load_name] = factor
    #
    def __getitem__(self, load_name:Union[str,int]):
        """
        """
        return self._basic[load_name]

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