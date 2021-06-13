#
# Copyright (c) 2009-2021 fem2ufo
# 

# Python stdlib imports
from array import array
from typing import NamedTuple, Tuple, List, Iterator, Dict, Iterable, ClassVar, Union
from collections.abc import Mapping

# package imports
from steelpy.f2uModel.load.operations.operations import (check_list_units, 
                                                         check_list_number, 
                                                         check_point_dic)
#from steelpy.f2uModel.load.sqlite.node import NodeLoadSQL

#
#
class PointNode(NamedTuple):
    """
    """
    fx: float
    fy: float
    fz: float
    mx: float
    my: float
    mz: float
    number: int
    load_name: str
    system:int
    load_complex:int
    #
    def __str__(self,units:str="si") -> str:
        """ """
        output  = (f"{str(self.number):12s} {10*' '} "
                   f"{self.fx: 1.3e} {self.fy: 1.3e} {self.fz: 1.3e} "
                   f"{self.coordinate_system.upper():12s}\n")
                   #f"{0: 1.3e} {0: 1.3e} {0: 1.3e}\n")
        if (step := self.load_name) == "NULL":
            step = 12 * " "
        output += (f"{step:12s} {10*' '} "
                   f"{self.mx: 1.3e} {self.my: 1.3e} {self.mz: 1.3e} {self.load_complex}\n")
        return output
    #
    @property
    def coordinate_system(self) -> str:
        if self.system != 0:
            return "local"
        return "global" 
#
#
class NodeLoadBasic(Mapping):
    """
    """
    #__slots__ = ['_nodes']
    
    def __init__(self) -> None:
        """
        """
        self._labels: List[Union[str, int]] = []
        self._title: List[Union[str, int, None]] = []
        self._complex: array = array("I", [])
        # 0-global/ 1-local
        self._system_flag:int = 0
        self._system: array = array("I", [])        
    
    def __iter__(self)-> Iterable:
        """
        """
        items = list(set(self._labels))
        return iter(items)
    
    def __contains__(self, value) -> bool:
        return value in self._labels
    
    def __len__(self) -> float:
        return len(self._labels)
    #
    def __str__(self) -> str:
        """ """
        print('---')
#
#
class NodeLoadMaster(NodeLoadBasic):
    """
    FE Node Load class
    
    NodeLoad
        |_ name
        |_ number
        |_ type
        |_ complex
        |_ point [x, y, z, mx, my, mz]
        |_ acceleration [th[0], th[1], th[2],..., th[n]]
        |_ displacement [th[0], th[1], th[2],..., th[n]]
        |_ mass
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    __slots__ =  ['_title', '_labels', '_index', '_complex',
                  '_fx', '_fy', '_fz', '_mx', '_my', '_mz',
                  '_system', '_system_flag', '_distance']
    
    def __init__(self) -> None:
        """
        """
        super().__init__()
        # real
        self._fx: array  = array('f',[])
        self._fy: array  = array('f',[])
        self._fz: array  = array('f',[])
        self._mx: array = array('f',[])
        self._my: array = array('f',[])
        self._mz: array = array('f',[])
        self._distance: array = array('f',[])
    #
    #@property
    #def name(self) -> str:
    #    """
    #    """
    #    return self._title[self._index]
    #
    #@name.setter
    #def name(self, load_name:str) -> None:
    #    """
    #    """
    #    try:
    #        self._title[self._index] = load_name
    #    except AttributeError:
    #        #self.load_name = load_name
    #        raise IndexError("load name not found")
    #
    def __delitem__(self, node_number: int) -> None:
        """
        """    
        indexes = [i for i, x in enumerate(self._labels) 
                   if x == node_number]
        indexes.sort(reverse=True)
        for _index in indexes:
            self._fx.pop(_index)
            self._fy.pop(_index)
            self._fz.pop(_index)
            self._mx.pop(_index)
            self._my.pop(_index)
            self._mz.pop(_index)      
            #
            self._labels.pop(_index)
            self._title.pop(_index)
            self._system.pop(_index)
            self._distance.pop(_index)
            self._complex.pop(_index)
#
#
def get_nodal_load(load):
    """ """
    if isinstance(load, (list, tuple)):
        try:
            load = check_list_units(load)
            load = load[:6]
        except AttributeError:
            load = check_list_number(load, steps=6)
    elif isinstance(load, dict):
        load = check_point_dic(load)
    else:
        raise Exception('   *** Load input format not recognized')
    return load
#
#
#

