# 
# Copyright (c) 2009-2021 fem2ufo
#

# Python stdlib imports
#from typing import NamedTuple, Dict, List, Iterable, Union

# package imports
from steelpy.f2uModel.load.inmemory.basic_load import BasicLoad
from steelpy.f2uModel.load.inmemory.combination import LoadCombination
#
#
class LoadingInmemory:
    
    __slots__ = ['_basic', '_combination']
    
    def __init__(self):
        """
        """
        self._basic = BasicLoad()
        self._combination = LoadCombination(basic_load=self._basic)
    
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
#