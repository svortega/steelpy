#
# Copyright (c) 2009-2023 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from array import array
#from collections.abc import Mapping
#from collections import defaultdict
#from operator import sub, add
import re
#from typing import NamedTuple, Tuple, List, Union, Iterable, Dict


# package imports
# steelpy
#from ....process.units.buckingham import Number
#from ....process.math.operations import linspace
# steelpy.f2uModel.load
from ..process.nodes import NodeLoadMaster
from ..inmemory.node import NodeLoadIM
from ..process.beam import (LineBeam, PointBeam, BeamDistMaster,
                            BeamLoadItem, BeamLoad)
from ..process.operations import (check_point_dic, check_list_units,
                                  check_beam_dic, #get_beam_point_load, 
                                  get_beam_node_load, get_beam_udl_load)
#
# ---------------------------------
#
class BeamLoadItemIM(BeamLoadItem):
    __slots__ = ['_labels', '_load', '_node_eq',
                 '_beam_id', '_f2u_beams']

    def __init__(self, load_name:str|int, 
                 load_title:str, beams) -> None:
        """
        """
        super().__init__()
        #self._load_name = load_name
        self._load = BeamLoadIM(load_name=load_name,
                                load_title=load_title)
        self._node_eq = BeamToNodeIM(load_name=load_name)
        self._f2u_beams = beams

    #
    def __setitem__(self, beam_name: int|str,
                    beam_load: list) -> None:
        """
        """
        try:
            self._f2u_beams[beam_name]
        except KeyError:
            raise IOError(f"beam {beam_name} not found")
        #
        self._beam_id = beam_name
        self._labels.append(beam_name)
        if re.match(r"\b(point)\b", str(beam_load[0]), re.IGNORECASE):
            self._load._point[beam_name] = beam_load[1:]
            
        elif re.match(r"\b(line|u(niform(ly)?)?d(istributed)?l(oad)?)\b",
                       str(beam_load[0]), re.IGNORECASE):
            self._load._line[beam_name] = beam_load[1:]
        
        else:
            raise IOError(f'Beam lod type {beam_load[0]} not implemented')

    #
    def __getitem__(self, beam_name: int | str):
        """
        """
        #
        try:
            beam = self._f2u_beams[beam_name]
        except KeyError:
            raise IOError(f"beam {beam_name} not found")
        #
        #if not beam_name in self._labels:
        self._labels.append(beam_name)
        return self._load(beam=beam) # ,  beam_name=beam_name
    #
    #
    def fer(self):
        """ Return Fix End Reactions (FER) global system"""
        beams = self._f2u_beams
        for key in set(self._labels):
            beam = beams[key]
            end_nodes = beam.connectivity
            res = self._load(beam=beam).fer()
            for gnload in res:
                self._node_eq[key] = [[end_nodes[0], *gnload[4], gnload[2]],
                                      [end_nodes[1], *gnload[5], gnload[2]]]
        #print('--> get_end_forces')
        #1 / 0
    #    
#
#
#
#
#
class BeamLoadIM(BeamLoad):
    __slots__ = ['_system_flag','_beam_id', '_load_name',
                 '_line', '_point', '_load_title', '_beam']

    def __init__(self, load_name:str|int, load_title:str):
        """
        """
        super().__init__()
        self._load_name = load_name
        self._load_title = load_title
        self._line = BeamDistributedIM(load_name)
        self._point = BeamPointIM(load_name)
        #self._node_eq = BeamToNodeIM(load_name)
    #
    def __call__(self, beam): #beam_name:int|str
        """ """
        #self._beam_id = beam_name
        self._beam = beam
        return self
    #
    #@property
    #def line(self):
    #    """ """
    #    return self._line
    #
    #@property
    #def point(self):
    #    """point load"""
    #    return self._point
    #       
#
#
#
# ---------------------------------
# Beam Loading Calculation
#
class BeamDistributedIM(BeamDistMaster):
    """
    """
    __slots__ = ['_type', '_labels', '_load_name',
                 '_index', '_complex', '_L', 
                 '_L1', '_qx1', '_qy1', '_qz1', 
                 '_L2', '_qx2', '_qy2', '_qz2',
                 '_system', '_title'] #'_system_flag',

    def __init__(self, load_name:str|int) -> None:
        """
        """
        super().__init__()
        #self._L: array = array("f", [])
        self._load_name = load_name
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
    #
    def __setitem__(self, element_name: int|str,
                    line_load: list|dict) -> None:
        """
                        |
             q0         | q1
        o------|        |----------o
        |                          |
        +  L0  +        +    L1    +
        +........... L ............+
        
        [L, qx1, qy2, qz1, qx1, qy1, qz1, L0, L1, comment]
        """
        self._labels.append(element_name)
        self._load_id.append(self._load_name)
        self._complex.append(0)
        #
        #Lbeam =  line_load.pop(0)
        #self._L.append(Lbeam)
        # end 1
        self._qx1.append(line_load[0])
        self._qy1.append(line_load[1])
        self._qz1.append(line_load[2])
        # end 2
        self._qx2.append(line_load[3])
        self._qy2.append(line_load[4])
        self._qz2.append(line_load[5])
        # distance from ends
        self._L1.append(line_load[6])
        self._L2.append(line_load[7])
        #
        self._system.append(line_load[8])
        self._title.append(line_load[9])
    #
    def __getitem__(self, element_name: int|str) :
        """
        """
        index_list: list = [x for x, item in enumerate(self._labels)
                            if item == element_name]
        
        udl_list: list = []
        for index in index_list:
            udl_list.append(LineBeam(self._qx1[index], self._qy1[index], self._qz1[index],
                                      self._qx2[index], self._qy2[index], self._qz2[index],
                                      self._L1[index], self._L2[index],
                                      self._labels[index], self._title[index],  
                                      self._load_id[index], 
                                      self._system[index], self._complex[index]))
        return udl_list
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
    #
    def __delitem__(self, element_name: int) -> None:
        """
        """
        index_list: list = [x for x, _item in enumerate(self._labels)
                            if _item == element_name]
        index_list.sort(reverse=True)

        for index in index_list:
            #self._L.pop(index)
            #
            self._L1.pop(index)
            self._qx1.pop(index)
            self._qy1.pop(index)
            self._qz1.pop(index)
            #
            self._L2.pop(index)
            self._qx2.pop(index)
            self._qy2.pop(index)
            self._qz2.pop(index)
            #
            self._labels.pop(index)
            self._title.pop(index)
            self._complex.pop(index)
    #
    #
    #def to_node(self, beam):
    #    """
    #    return
    #    In plane [V, M, w, theta]
    #    Out plane [V, M, w, theta]
    #    [fx, fy, fz, mx, my, mz]
    #    """
    #    #
    #    items = []
    #    for index, name in enumerate(self._labels):
    #        item = LineBeam(self._qx1[index], self._qy1[index], self._qz1[index],
    #                        self._qx2[index], self._qy2[index], self._qz2[index],
    #                        self._L[index], self._L1[index], self._L2[index],
    #                        self._labels[index], self._title[index],
    #                        self._system[index], self._complex[index])
    #        #
    #        beam.load = item
    #        reactions = beam.reactions()
    #        items.append(reactions)
    #
    #    return items
#
#
class BeamPointIM(NodeLoadMaster):
    __slots__ = ['_title', '_labels', '_load_name',
                 '_index', '_complex',  
                 '_fx', '_fy', '_fz', '_mx', '_my', '_mz',
                 '_L1', '_type']
    
    def __init__(self, load_name:str|int) -> None:
        """
        """
        super().__init__("point")
        self._load_name = load_name
        self._L1: array = array('f', [])
    #
    def __setitem__(self, element_name:str|int,
                    point_load: list|dict) -> None:
        """
        """
        self._labels.append(element_name)
        self._load_id.append(self._load_name)
        self._complex.append(0)
        #
        #Lbeam =  point_load.pop(0)
        #
        #self._L.append(Lbeam)
        #point_load = get_beam_point_load(point_load)
        self._fx.append(point_load[0])
        self._fy.append(point_load[1])
        self._fz.append(point_load[2])
        self._mx.append(point_load[3])
        self._my.append(point_load[4])
        self._mz.append(point_load[5])
        #
        self._L1.append(point_load[6])
        #
        self._system.append(point_load[7])
        self._title.append(point_load[8])
        #print('--')
    
    def __getitem__(self, element_name:str|int):
        """
        """
        index_list: list = [x for x, item in enumerate(self._labels)
                            if item == element_name]
        #
        points: list = []
        for index in index_list:
            points.append(PointBeam(self._fx[index], self._fy[index], self._fz[index],
                                    self._mx[index], self._my[index], self._mz[index],
                                    self._L1[index],
                                    self._labels[index], self._title[index],
                                    self._load_id[index],
                                    self._system[index], self._complex[index]))
        return points
    #
    #
    def __delitem__(self, element_name:int|str) -> None:
        """
        """
        indexes = [i for i, x in enumerate(self._labels)
                   if x == element_name]
        
        indexes.sort(reverse=True)
        
        for _index in indexes:
            #self._L.pop(_index)
            self._L1.pop(_index)
            #
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
            self._complex.pop(_index)    
    #
    #def get_nodal_load(self, elements, materials, sections) -> List:
    #    """
    #    """
    #    items = point2node(self, elements, materials, sections)
    #    return items
#
#
class BeamToNodeIM(NodeLoadIM):
    __slots__ = ['_title', '_labels', '_index', '_complex',  
                 '_fx', '_fy', '_fz', '_mx', '_my', '_mz',
                 '_beam', '_type']
    
    def __init__(self, load_name:str|int) -> None:
        """
        """
        super().__init__(load_name, None, "load")
        self._beam: list = []
    #
    def __setitem__(self, element_name:str|int,
                    node_load: list|dict) -> None:
        """
        """
        for node in node_load:
            self._beam.append(element_name)
            node_number = node[0]
            point_load = node[1:]
            super().__setitem__(node_number, point_load)
    #
    #
    #
    #
#
#
