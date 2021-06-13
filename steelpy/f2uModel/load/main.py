# 
# Copyright (c) 2009-2021 fem2ufo
#

# Python stdlib imports
from typing import NamedTuple, Dict, List, Iterable, Union

# package imports
from steelpy.f2uModel.load.inmemory.main import LoadingInmemory
from steelpy.f2uModel.load.sqlite.main import LoadingSQL
#
#
class Load:
    """
    """
    __slots__ = ['_load']
    
    def __init__(self, mesh_type:str,
                 db_file:Union[str,None]):
        """
        """
        if mesh_type != "inmemory":
            self._load = LoadingSQL(db_file=db_file,
                                    db_system=mesh_type)
        else:
            self._load = LoadingInmemory()
    #
    @property
    def basic(self):
        """
        """
        return self._load._basic

    @basic.setter
    def basic(self, values:List):
        """
        """
        number = 1
        for item in values:
            try:
                self._load._basic[number] = item[0]
                number =+ 1
            except Warning:
                pass
            if item[1] in ['node']:
                print('node')
                self._load._basic[number].node = item[2:]
            else:
                print('beam')
    #
    @property
    def combination(self):
        """
        """
        return self._load._combination
    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += self._load._basic.__str__()
        output += self._load._combination.__str__()
        return output
#
#