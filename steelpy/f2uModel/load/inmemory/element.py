#
# Copyright (c) 2009-2021 fem2ufo
# 

# Python stdlib imports
from array import array
from collections.abc import Mapping
#from collections import defaultdict
#from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, Union, Iterable, Dict  


# package imports
from steelpy.f2uModel.load.operations.operations import NodeLoadMaster, get_beam_load
from steelpy.f2uModel.load.operations.element import (LineBeam, PointBeam,  point2node, #line2node,
                                                      get_beam_point_load)

#
# ---------------------------------
#
class BeamDistributed(Mapping):
    """
    """
    __slots__ = ['_type', '_labels', '_index', '_complex',
                 '_L1', '_qx1', '_qy1', '_qz1', 
                 '_L2', '_qx2', '_qy2', '_qz2',
                 '_system', '_system_flag', '_title']

    def __init__(self) -> None:
        """
        """
        # ens 1
        self._L1: array = array("f", [])
        self._qx1: array = array("f", [])
        self._qy1: array = array("f", [])
        self._qz1: array = array("f", [])
        # end 2
        self._L2: array = array("f", [])
        self._qx2: array = array("f", [])
        self._qy2: array = array("f", [])
        self._qz2: array = array("f", [])
        #
        self._labels: List[Union[str, int]] = []
        self._title: List[str] = []
        self._index: int
        self._complex: array = array("I", [])
        # 0-global/ 1-local
        self._system_flag:int = 0
        self._system: array = array("I", [])
    #
    def __setitem__(self, element_name: Union[int, str], 
                    udl: Union[List[float], Dict[str,float]]) -> None:
        """
        """
        self._labels.append(element_name)
        self._title.append(element_name)
        self._index = len(self._labels)-1
        self._system.append(self._system_flag)
        self._complex.append(0)
        #
        # update inputs
        udl = get_beam_load(udl)
        # end 1
        self._qx1.append(udl[0])
        self._qy1.append(udl[1])
        self._qz1.append(udl[2])
        # end 2
        self._qx2.append(udl[3])
        self._qy2.append(udl[4])
        self._qz2.append(udl[5])
        # distance from ends
        self._L1.append(udl[6])
        self._L2.append(udl[7])
    #
    def __getitem__(self, element_name: Union[int, str]) -> List[Tuple]:
        """
        """
        _index_list: List = [x for x, _item in enumerate(self._labels)
                             if _item == element_name]
        
        _udl_list: List[Tuple] = []
        for _index in _index_list:
            _udl_list.append(LineBeam(self._qx1[_index], self._qy1[_index], self._qz1[_index],
                                      self._qx2[_index], self._qy2[_index], self._qz2[_index],
                                      self._L1[_index], self._L2[_index],
                                      self._labels[_index], self._title[_index],
                                      self._system[_index], self._complex[_index]))
        return _udl_list
    #
    @property
    def name(self) -> str:
        """
        """
        return self._title[self._index]
    
    @name.setter
    def name(self, load_name:str) -> None:
        """
        """
        try:
            self._title[self._index] = load_name
        except AttributeError:
            #self.load_name = load_name
            raise IndexError("load name not found")    
    #
    #
    @property
    def coordinate_system(self):
        if self._system_flag != 0:
            return "local"
        return "global"
    
    @coordinate_system.setter
    def coordinate_system(self, system:Union[str,int]):
        """
        Coordinate system for load : global or local (member)
        """
        self._system_flag = 0
        if system in ['local', 'member', 1]:
            self._system_flag = 1
        
    #
    #
    def __len__(self) -> float:
        return len(self._labels)

    def __delitem__(self, element_name: int) -> None:
        """
        """
        _index_list: List = [x for x, _item in enumerate(self._labels)
                             if _item == element_name]

        _index_list.sort(reverse=True)
        for _index in _index_list:
            self._L1.pop(_index)
            self._qx1.pop(_index)
            self._qy1.pop(_index)
            self._qz1.pop(_index)
            #
            self._L2.pop(_index)
            self._qx2.pop(_index)
            self._qy2.pop(_index)
            self._qz2.pop(_index)
            #
            self._labels.pop(_index)
            self._title.pop(_index)
            self._complex.pop(_index)

    def __contains__(self, value) -> bool:
        return value in self._labels

    def __iter__(self) -> Iterable:
        """
        """
        items = list(set(self._labels))
        return iter(items)
    #
    #
    #def get_nodal_load(self, elements, materials, sections) -> List:
    #    """
    #    """
    #    items = [ ]
    #    #set_items = set(self._labels )
    #    #if not set_items:
    #    #    return [ ]
    #    #
    #    for item in self.__iter__():
    #        # _index = self._labels.index(item)
    #        element = elements[item]
    #        material = materials[element.material]
    #        section = sections[element.section].properties
    #        # property = section.properties
    #        #
    #        # nodes = element.nodes
    #        #
    #        nloads = self.__getitem__(item)
    #        for nl in nloads:
    #            nodal_load, member_nload = nl.node_equivalent(element, material, section)
    #            #print ( '-->' )
    #    return items
#
#
class BeamPoint(NodeLoadMaster):
    __slots__ =  ['_title', '_labels', '_index', '_complex',
                  '_fx', '_fy', '_fz', '_mx', '_my', '_mz']
    
    def __init__(self) -> None:
        """
        """
        super().__init__()
    #
    #
    def __setitem__(self, element_name:Union[int, str], 
                    point_load: Union[List[float], Dict[str,float]]) -> None:
        """
        """
        self._labels.append(element_name)
        self._title.append(element_name)
        self._system.append(self._system_flag)
        self._complex.append(0)
        #
        point_load = get_beam_point_load(point_load)
        self._fx.append(point_load[0])
        self._fy.append(point_load[1])
        self._fz.append(point_load[2])
        self._mx.append(point_load[3])
        self._my.append(point_load[4])
        self._mz.append(point_load[5])
        self._distance.append(point_load[6])
        #
        self._index = len(self._labels) - 1
        #print('--')
    
    def __getitem__(self, element_name:Union[int, str])-> List[Tuple]:
        """
        """
        _index_list: List = [x for x, _item in enumerate(self._labels)
                             if _item == element_name]
        #
        _points: List = []
        for _index in _index_list:
            _points.append(PointBeam(self._fx[_index], self._fy[_index], self._fz[_index],
                                     self._mx[_index], self._my[_index], self._mz[_index],
                                     self._distance[_index], self._labels[_index], 
                                     self._title[_index], self._system[_index], 
                                     self._complex[_index]))
        return _points    
    #
    @property
    def coordinate_system(self):
        if self._system_flag != 0:
            return "local"
        return "global"
    
    @coordinate_system.setter
    def coordinate_system(self, system:Union[str,int]):
        """
        Coordinate system for load : global or local (member)
        """
        self._system_flag = 0
        if system in ['local', 'member', 1]:
            self._system_flag = 1
        
    #    
    #
    #def get_nodal_load(self, elements, materials, sections) -> List:
    #    """
    #    """
    #    items = point2node(self, elements, materials, sections)
    #    return items
#
#
