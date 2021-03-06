#
# Copyright (c) 2009-2021 fem2ufo
# 

# Python stdlib imports
#from array import array
#from collections.abc import Mapping
#from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, Iterator, Dict, Iterable, ClassVar, Union
#import re

# package imports
from steelpy.f2uModel.load.operations.nodes import (NodeLoadMaster, 
                                                    get_nodal_load, 
                                                    PointNode)

# ---------------------------------
#
class NodeLoadInMemory(NodeLoadMaster):
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
        #
        if isinstance(point_load, dict):
            point_load = get_nodal_load(point_load)
            self._title.append(point_load[-1])
            point_load.pop()            
        elif isinstance(point_load[-1], str):
            title = point_load.pop()
            self._title.append(title)
            point_load = get_nodal_load ( point_load )
        else:
            self._title.append("NULL")
            point_load = get_nodal_load(point_load)
        #
        self._system.append(self._system_flag)
        self._distance.append(-1)
        self._complex.append(0)
        #
        #point_load = get_nodal_load(point_load)
        self._fx.append(point_load[0])
        self._fy.append(point_load[1])
        self._fz.append(point_load[2])
        self._mx.append(point_load[3])
        self._my.append(point_load[4])
        self._mz.append(point_load[5])
        #
        self._index = len(self._labels) - 1
        #print('--')
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
                                     self._labels[_index], self._title[_index],
                                     self._system[_index], self._complex[_index]))
        return _points
    #
    #@property
    #def get_items_sum(self):
    def get_group_nodes(self):
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
                                   self._labels[_index],
                                   self._title[_index], 
                                   self._system[_index],
                                   self._complex[_index]))
        return items
    #
    def update(self, other) -> None:
        """
        """
        #1/0
        pl = other._point
        try:
            #_name = pl.name
            index = pl._index
            # update
            self._labels.extend(pl._labels)
            self._title.extend(pl._title)
            self._system.extend(pl._system)
            self._distance.extend(pl._distance)
            self._complex.extend(pl._complex)
            #
            self._fx.extend(pl._fx)
            self._fy.extend(pl._fy)
            self._fz.extend(pl._fz)
            self._mx.extend(pl._mx)
            self._my.extend(pl._my)
            self._mz.extend(pl._mz)
            # print('?????')
        except AttributeError:
            pass
    #
#
