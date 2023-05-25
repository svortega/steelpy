# 
# Copyright (c) 2009-2023 fem2ufo
#
# Python stdlib imports
from __future__ import annotations
import re
from collections.abc import Mapping
from typing import NamedTuple, Tuple, List, Union

#
# package imports
from steelpy.process.units.main import Units
from steelpy.f2uModel.properties.process.operations import BasicProperty


class MarineGrowth(Mapping):
    """
    """
    __slots__ = ['_marine_growth', 'f2u_units',
                 'sea_water_density', '_default']
    
    def __init__(self):
        """
        """
        global f2u_units #, mg_default
        f2u_units = Units()
        #
        self._default: None|str = None
        self._marine_growth:dict = {}
        self.sea_water_density:float = 1032.0 * f2u_units.kg / f2u_units.m**3
    #
    def __getitem__(self, mg_name: str|int):
        """
        """
        return self._marine_growth[mg_name]
    
    def __setitem__(self, mg_name:str|int, mg_type: str) -> None:
        """
        """
        if re.match(r"\b(constant)\b", mg_type, re.IGNORECASE):
            self._marine_growth[mg_name] = MGconstant()
        elif re.match(r"\b(profile)\b", mg_type, re.IGNORECASE):
            self._marine_growth[mg_name] = MGprofile(self)
        else:
            raise IOError("marine growth type {:} not implemented".format(mg_type))
        #
        self._marine_growth[mg_name].name = mg_name
    #
    #def add_new(self, mg_name, mg_number=None):
    #    """
    #    """
    #    if mg_number:
    #        _number = mg_number
    #    else:
    #        _number = len(self._marine_growth) + self.number        
    #    self.__setitem__(mg_name, _number)
    #    
    #
    def __len__(self) -> int:
        return len(self._marine_growth)

    def __iter__(self):
        """
        """
        return iter(self._marine_growth)
    
    def __contains__(self, value) -> bool:
        return value in self._marine_growth
#
#
class Profile(NamedTuple):
    """
    """
    elevation : List
    thickness : List
#
class MGprofile(BasicProperty):
    """
    FE Marine Growth Class
    
    Marine_Growth
        |_ name
        |_ number
        |_ case : constant_thickness / depth_profile
        |_ profile
        |_ inertia
        |_ factor
        |_ items
        |_ sets
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    #
    __slots__ = ['option', '_units', 'case', #'name',
                 'inertia', 'factor', 'density', 
                 'density_factor', 'absolute_elevations',
                 '_elevation', '_thickness', '_roughness',
                 '_densities', '_label', '_cls']
    #
    def __init__(self, cls):
        """
        """
        BasicProperty.__init__(self)
        #
        self._cls = cls
        self.case = 'profile'
        #
        self.inertia : bool = True
        self.absolute_elevations : bool = False
        self.density_factor : float = 1.36
        self.density : Units = 1400.0 * f2u_units.kg / f2u_units.m**3
        #
        self._label : list[float] = []
        self._elevation : list[float] = []
        self._thickness: list[float] = []
        self._roughness: list[float] = []
        self._densities: list[float] = []
    #
    #
    #def __getitem__(self, level_name:Union[str,int]):
    #    """
    #    """
    #    print('--')
    #
    #def __setitem__(self, level_name:Union[str,int], data:List) -> None:
    @property
    def level(self):
        """
        level_name : Elevation ID
        data = [level, thickness, roughness, density]  
        ----------------------------------------------
        level : elevation above mudline
        thickness : marine growth thickness
        roughness : surface roughness (0m defaul)
        density : dry weight density of marine growth (1400kg/m3 defaul)
        """
        return MGtypes(cls=self)
    #
    #
    def get_density_factor(self, density):
        """
        """
        self.density_factor = density/self.sea_water_density
        return self.density_factor
    #
    def set_default(self):
        """
        """
        #print('--')
        self._cls._default = self.name
    #   
#
class MGconstant(BasicProperty):
    __slots__ = ['case', '_units', '_roughness',
                 'density', 'density_factor', '_thickness']
    
    def __init__(self):
        """
        """
        BasicProperty.__init__(self)
        self.case = 'constant_thickness'
        #
        self._roughness : float = 0.0
        self.density : float = 1400.0 * f2u_units.kg / f2u_units.m**3
        self.density_factor : float = 1.36
    #
    #
    def constant(self, thickness, roughness, 
                 density_factor=None):
        """
        Marine Growth defined with a Constant thickness
        """
        self._thickness = thickness
        self._roughness = roughness
        if density_factor:
            self.density_factor = density_factor
    #
    @property
    def thickness(self):
        """
        """
        return self._thickness
    
    @thickness.setter
    def thickness(self, value):
        """
        """
        self._thickness = value
    #
    @property
    def roughness(self):
        """
        """
        return self._roughness
    
    @roughness.setter
    def roughness(self, value):
        """
        """
        self._roughness = value
#
#
class MGtypes:
    
    __slots__ = ["_cls"]
    
    def __init__(self, cls):
        """
        """
        self._cls = cls
    
    #
    #
    def __getitem__(self, level_name:str|int):
        """
        """
        print('--')
    #
    def __setitem__(self, level_name:str|int,
                    data:list) -> None:
        """
        """
        self._cls._label.append(level_name)
        #
        data = self._get_value(data, steps=4)
        self._cls._elevation.append(data[0])
        self._cls._thickness.append(data[1])
        self._cls._roughness.append(data[2])
        try:
            1/data[3]
            self._cls._densities.append(data[3])
        except ZeroDivisionError:
            self._cls._densities.append(1400.0)
        print('--')
    #
    #
    #
    #
    def _get_value(self, value, steps):
        """
        """
        if isinstance(value, (list, tuple)):
            value = get_list(value, steps)
        elif isinstance(value, dict):
            value = get_dic(value)
        else:
            raise Exception('   *** input format not recognized')
        return value
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
def get_dic(data)->list[float]:
    """ """
    new_data = [0,0,0,0]
    for key, item in data.items():
        if re.match(r"\b(elevation)\b", str(key), re.IGNORECASE):
            new_data[0] = item.value
        elif re.match(r"\b(thickness)\b", str(key), re.IGNORECASE):
            new_data[1] = item.value
        elif re.match(r"\b(roughness)\b", str(key), re.IGNORECASE):
            new_data[2] = item.value
        elif re.match(r"\b(density)\b", str(key), re.IGNORECASE):
            new_data[3] = item.value
    return new_data
#
#