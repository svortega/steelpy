# 
# Copyright (c) 2009 steelpy
#
# 
# Python stdlib imports
from __future__ import annotations
#from array import array
from collections.abc import Mapping
#from collections import defaultdict
#from dataclasses import dataclass
#from math import prod

# package imports
#from steelpy.utils.dataframe.main import DBframework
from ..process.combination import LoadCombinationBasic

#
#
class LoadCombConcept(LoadCombinationBasic):
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
    __slots__ = ('number', 'name', '_component',
                 '_combination', '_basic')
    #
    def __init__(self, basic_load, component:str|int):
        """
        """
        super().__init__()
        self._combination = {}
        self._basic = basic_load
        self._component = component
    
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
            self._combination[load_name] = CombTypes(name=load_name,
                                                     title=load_title,
                                                     component=self._component)
            #
            load_id = next(self.get_number())
            self._combination[load_name].number = load_id
            self._number.append(load_id)
    
    def __getitem__(self, load_name:str|int):
        """
        """
        return self._combination[load_name]
    #
    def __delitem__(self, load_name:str|int):
        """
        """
        del self._combination[load_name]
    #
    #
    #def to_basic(self):
    #    """ """
    #    db = DBframework()
    #    # get combination of combination and convert them to basic loads
    #    # TODO : check this loop works
    #    for key, item in self._combination.items():
    #        for comb_name, factor in item._combination.items():
    #            for precomb, prefactor in self._combination[comb_name]._basic.items():
    #                try:
    #                    item._basic[precomb] += prefactor * factor
    #                except KeyError:
    #                    item._basic[precomb] = prefactor * factor
    #    # organize basic load
    #    #basic_loads = defaultdict(list)
    #    # form combination formed by basic loads only
    #    dftemp = []
    #    for key, item in self._combination.items():
    #        for bl_name, factor in item._basic.items():
    #            #basic_loads[key].append([bl_name, factor])
    #            dftemp.append([key, item.number, 'combination', item.title, bl_name, factor])
    #            #try:
    #            #    #basic = self._basic[bname]
    #            #    #basic_loads[item.name].append([basic.title, factor])
    #            #    basic_loads[key].append([bl_name, factor])
    #            #except KeyError:
    #            #    raise Warning("  warning basic load {:} not found".format(bname))
    #    #
    #    header = ['load_name', 'load_id','load_type', 'load_title', 'basic_load', 'factor']
    #    dfcomb = db.DataFrame(data=dftemp, columns=header, index=None)
    #    return dfcomb #, basic_loads
#
#
class CombTypes:
    """
    """
    __slots__ = ['_basic', '_metocean', 'name', 'number',
                 'title', '_combination','_component']
    
    def __init__(self, name:int, title:str,
                 component:str|int):
        """
        """
        self.name = name
        self.title = title
        self._component = component
        #
        self._basic = BasicLoad("basic")
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
    def basic(self):
        """
        """
        return self._basic
    #
    @property
    def combination(self):
        """
        """
        return self._combination
    #
    @property
    def mesh_name(self):
        """ component name """
        return self._component
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
        self._basic:dict = {}
        self._type = bl_type
    #
    def __setitem__(self, load_name: str|int, factor: float) -> None:
        """
        """
        self._basic[load_name] = factor
    #
    def __getitem__(self, load_name:str|int):
        """
        """
        return self._basic[load_name]

    #
    #
    def __delitem__(self, load_name:str|int):
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
        #comb =  list(dict.fromkeys(self._basic))
        #return iter(comb)
        return iter(self._basic)
    #
    #
    @property
    def load_type(self):
        """
        """
        return self._type
#
