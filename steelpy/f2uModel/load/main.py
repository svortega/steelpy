# 
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports

# package imports
from steelpy.f2uModel.load.basic_load import BasicLoad
from steelpy.f2uModel.load.combination import LoadCombination

#
#
class Loading:
    
    __slots__ = ['_basic', '_combination']
    
    def __init__(self):
        """
        """
        self._basic = BasicLoad()
        self._combination = LoadCombination()
    
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