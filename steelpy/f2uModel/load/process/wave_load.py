#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from array import array
from collections.abc import Mapping
#from collections import defaultdict
#from dataclasses import dataclass
#from operator import sub, add
from typing import NamedTuple
# import re
#
# package imports
#
from steelpy.f2uModel.load.process.basic_load import BasicLoadMain
from steelpy.f2uModel.load.concept.beam import BeamIMMaster

#
#
class WaveData(NamedTuple):
    """ """
    name: str | int
    seastate: dict
    design_load: str
    criterion : str
#
#@dataclass
class WaveLoadItem(BeamIMMaster):
    """ """
    __slots__ = ['_seastate', '_name',
                 '_load', '_criterion', '_design_load']

    def __init__(self):
        """ """
        super().__init__()
        #
        #self._name = load_name        
        #
        self._seastate:list = []
        self._design_load:list = []
        self._criterion:list = []
    #
    def __setitem__(self, load_name: int|str,
                    wave_load: list) -> None:
        """
        criterion : 
        """
        try:
            self._labels.index(load_name)
            raise Exception('wave load case {:} already exist'.format(load_name))
        except ValueError:
            self._labels.append(load_name)
            self._seastate.append(wave_load[0])
            self._design_load.append(wave_load[1])
            self._criterion.append(wave_load[2])
    #
    def __getitem__(self, load_name: int | str):
        """
        """
        try:
            index = self._labels.index(load_name)
        except ValueError:
            raise KeyError('Invalid load name : {:}'.format(load_name))
        #
        return WaveData(name=self._labels[index], 
                        seastate=self._seastate[index],
                        design_load=self._design_load[index],
                        criterion=self._criterion[index])        
    #
#
#
#
class HydroLoad(BasicLoadMain):
    
    def __init__(self, plane: NamedTuple) -> None:
        """
        """
        #super().__init__()
        #
        self._labels: list[str|int] = []
        self._title: list[str] = []
        self._number: array = array("I", [])        
        self._hydro: dict = {}        
        #
        #self._plane = plane
        #self._condition = WaveLoadItemSQL(db_file=self.db_file)
        #
    #
    #
    #
    # -----------------------------------------------
    #
    def __setitem__(self, load_name:int|str, load_title:str) -> None:
        """
        """
        try:
            self._labels.index(load_name)
            raise Warning(f'    *** warning load name {load_name} already exist')
        except ValueError:
            self._labels.append(load_name)
            self._title.append(load_title)
            self._number.append(load_number)
            #
            #self._basic[load_name] = LoadTypeSQL(name=load_name,
            #                                     number=load_number,
            #                                     title=load_title,
            #                                     plane=self._plane, 
            #                                     bd_file=self.db_file)
            #self._basic[load_name] = [load_title, load_number]
    #           
    def __getitem__(self, load_name: str|int):
        """
        """
        try:
            #load_name = load_name
            index = self._labels.index(load_name)
            #return self._basic[load_name]
            return LoadTypeSQL(load_name=load_name,
                               plane=self._plane, 
                               bd_file=self.db_file)
        except ValueError:
            raise IOError("load case not defined")    
#
#
#
class MetoceanLoad(Mapping):
    __slots__ = ['_condition', '_labels']
    
    def __init__(self):
        """ """
        self._labels:list = []
        self._condition:list = []
    #
    def __len__(self) -> int:
        return len(self._labels)

    def __iter__(self):
        """
        """
        return iter(self._labels)
    
    def __contains__(self, value):
        return value in self._labels    
    #
    #@property
    #def condition(self):
    #    """ """
    #    return self._condition
    #
    #@condition.setter
    #def condition(self, value):
    #    """ """
    #    self._condition.append(value)
    #
    #
    