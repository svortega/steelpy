#
# Copyright (c) 2009-2020 fem2ufo
# 

# Python stdlib imports
from array import array
from collections.abc import Mapping
#from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, Iterator, Dict, Iterable, ClassVar, Union
#import re

# package imports
from steelpy.f2uModel.load.operations import NodeLoadMaster # , check_point_list

# ---------------------------------
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
    name: str
    number: int
    system:str
    load_complex:int
    #
    def __str__(self):
        print('-->')
        return "{: 14.5f} {: 14.5f} {: 14.5f} {: 14.5f} {: 14.5f} {: 14.5f}".format(self.fx, self.fy, self.fz,
                                                                                    self.mx, self.my, self.mz)
    #
    #@property
    #def name(self):
      #  """
      #  """
      #  return self.load_name
    #@name.setter
    #def name(self, load_name:str):
    #    """
    #    """
    #    self.load_name = load_name  
#
#
class NodeLoad(NodeLoadMaster):
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
    #
    def __setitem__(self, node_number:int, 
                    point_load:List) -> None:
        """
        """
        #
        self._labels.append(node_number)
        self._title.append(node_number)
        self._index = self._labels.index(node_number)
        self._system.append(self._system_flag)
        self._distance.append(-1)
        self._complex.append(0)
        #
        point_load = self._get_nodal_load(point_load)
        self._fx.append(point_load[0])
        self._fy.append(point_load[1])
        self._fz.append(point_load[2])
        self._mx.append(point_load[3])
        self._my.append(point_load[4])
        self._mz.append(point_load[5])

    #
    def __getitem__(self, node_number: int)-> List[Tuple]:
        """
        """
        _index_list: List = [x for x, _item in enumerate(self._labels)
                             if _item == node_number]
        #
        _points: List = []
        for _index in _index_list:
            _points.append(PointNode(self._fx[_index], self._fy[_index], self._fz[_index],
                                     self._mx[_index], self._my[_index], self._mz[_index],
                                     self._title[_index], self._labels[_index], 
                                     self._system[_index], self._complex[_index]))
        return _points
    #
    @property
    def get_items_sum(self):
        """
        """
        items = []
        set_items = set(self._labels)
        for item in set_items:
            _index = self._labels.index(item)
            rows = self.__getitem__(item)
            sum_load = [0,0,0,0,0,0]
            for _load in rows:
                sum_load[ 0 ] += _load.fx
                sum_load[ 1 ] += _load.fy
                sum_load[ 2 ] += _load.fz
                sum_load[ 3 ] += _load.mx
                sum_load[ 4 ] += _load.my
                sum_load[ 5 ] += _load.mz
            items.append(PointNode(*sum_load, 
                                   self._title[_index], 
                                   self._labels[_index],
                                   self._system[_index],
                                   self._complex[_index]))
        return items
   
#
#
