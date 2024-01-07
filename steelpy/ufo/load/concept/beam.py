#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from array import array
from collections.abc import Mapping
#from collections import defaultdict
#from operator import sub, add
import re
#from typing import NamedTuple, Tuple, List, Union, Iterable, Dict


# package imports
# steelpy
# steelpy.f2uModel.load
from ..process.nodes import NodeLoadMaster
from ..concept.node import NodeLoadIM

from steelpy.ufo.load.process.beam.beam import (LineBeam, PointBeam, 
                                                BeamLoad) 

from steelpy.utils.units.buckingham import Number
from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.math.operations import trnsload

from ..process.operations import (check_point_dic, check_list_units,
                                  check_beam_dic,  
                                  get_beam_node_load, get_beam_udl_load)
#
#
# ---------------------------------
#
class BeamIMMaster(Mapping):
    
    def __init__(self) -> None:
        """
        """
        self._index: int
        self._labels: list[str|int] = []
        self._title: list[str] = []
        self._load_id: list[str|int] = []
        self._complex: array = array("I", [])
        # 0-global/ 1-local
        #self._system_flag: int = 0
        self._system: array = array("I", [])
    #
    def __len__(self) -> int:
        return len(self._labels)
    #
    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    def __iter__(self):
        """
        """
        items = list(set(self._labels))
        return iter(items)
    #
    def __str__(self) -> str:
        """ """
        output = ""
        beams = list(dict.fromkeys(self._labels))
        #beams = list(set(self._labels))
        for beam in beams:
            items = self.__getitem__(beam)
            for item in items:
                output += item.__str__()
        #print('---')
        return output
    #
    #
    #def _get_line_load(self):
    #    """ return line load in correct format"""
    #    print('-->')
    #    1/0
    #
    #
    #
    #
    #@property
    #def coordinate_system(self):
    #    if self._system_flag != 0:
    #        return "local"
    #    return "global"
    #
    #@coordinate_system.setter
    #def coordinate_system(self, system:str|int):
    #    """
    #    Coordinate system for load : global or local (member)
    #    """
    #    self._system_flag = 0
    #    if system in ['local', 'member', 1]:
    #        self._system_flag = 1
    #
    #
    # ------------------
    #   
#
#
#
# ---------------------------------
#
class BeamLoadItemIM(BeamIMMaster):
    __slots__ = ['_labels', '_load', '_node_eq',
                 '_beam_id', '_f2u_beams']

    def __init__(self, load_name:str|int, 
                 load_title:str,
                 component: int, 
                 beams) -> None:
        """
        """
        super().__init__()
        #self._load_name = load_name
        self._load = BeamLoadIM(load_name=load_name,
                                load_title=load_title,
                                component=component)
        #
        #self._node_eq = BeamToNodeIM(load_name=load_name)
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
    def __str__(self, units: str = "si") -> str:
        """ """
        output = ""
        output += self._line.__str__(load_name=self._name)
        output += self._point.__str__(load_name=self._name)
        return output
    #
    # ----------------------------------------
    #
    def _get_line(self, line_load: list|dict):
        """ get line load in beam local system"""
        #
        # update inputs
        if isinstance(line_load, dict):
            udl = check_beam_dic(line_load)
            title = udl.pop()
            
        elif isinstance(line_load[-1], str):
            title = line_load.pop()
            if isinstance(line_load[0], Number):
                udl = check_list_units(line_load)
            else:
                udl = get_beam_udl_load(line_load)
        else:
            title ='NULL'
            udl = get_beam_udl_load(line_load)
        #
        # get system local = 1
        try:
            1 / self._system_flag
            return [*udl, 1, title]
        except ZeroDivisionError:
            # local nodal loading
            nload = [*udl[:3], 0, 0, 0,
                     *udl[3:6], 0, 0, 0,]
            nload = trnsload(nload, self._beam.T3D())
            nload = [*nload[:3], *nload[6:9]] 
            return [*nload, *udl[6:], 1, title]
    #
    def _get_point(self, point_load: list|dict):
        """ get point load in beam local system"""
        # update inputs
        if isinstance(point_load, dict):
            point = check_point_dic(point_load)
            title = point.pop()
        
        elif isinstance(point_load[-1], str):
            title = point_load.pop()
            if isinstance(point_load[0], Number):
                point = check_list_units(point_load)
            else:
                point = get_beam_node_load(point_load)
        
        else:
            title = 'NULL'
            point = get_beam_node_load(point_load)
        #
        # get system local = 1
        try: # Local system
            1 / self._system_flag
            return [*point, 1, title]
        except ZeroDivisionError: # global to local system
            pload = [*point[:6], 0, 0, 0, 0, 0, 0]
            pload = trnsload(pload, self._beam.T3D())
            return [*pload[:6], point[6], 1, title]
    #
    # ----------------------------------------
    #
    def fer(self, beams):
        """ Return Fix End Reactions (FER) global system"""
        #beams = self._f2u_beams
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
    __slots__ = ['_system_flag','_beam_id', 'load_name',
                 '_line', '_point', '_load_title', '_beam']

    def __init__(self, load_name:str|int, load_title:str,
                 component: int):
        """
        """
        super().__init__()
        self.load_name = load_name
        self._load_title = load_title
        self._line = BeamDistributedIM(load_name,
                                       component=component)
        #
        self._point = BeamPointIM(load_name,
                                  component=component)
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
class BeamDistributedIM(BeamIMMaster):
    """
    """
    __slots__ = ['_type', '_labels', '_load_name',
                 '_index', '_complex', '_L', 
                 '_L0', '_qx0', '_qy0', '_qz0', 
                 '_L1', '_qx1', '_qy1', '_qz1',
                 '_system', '_title', '_component'] #'_system_flag',

    def __init__(self, load_name:str|int,
                 component: int) -> None: 
        """
        """
        super().__init__()
        self._component = component
        self._load_name = load_name
        # ens 1
        self._L0: array = array("f", [])
        self._qx0: array = array("f", [])
        self._qy0: array = array("f", [])
        self._qz0: array = array("f", [])
        # end 2
        self._L1: array = array("f", [])
        self._qx1: array = array("f", [])
        self._qy1: array = array("f", [])
        self._qz1: array = array("f", [])
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
        Ltype =  line_load.pop(0)
        #self._L.append(Lbeam)
        # end 1
        self._qx0.append(line_load[0])
        self._qy0.append(line_load[1])
        self._qz0.append(line_load[2])
        # end 2
        self._qx1.append(line_load[3])
        self._qy1.append(line_load[4])
        self._qz1.append(line_load[5])
        # distance from ends
        self._L0.append(line_load[6])
        self._L1.append(line_load[7])
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
            udl_list.append(LineBeam(self._qx0[index], self._qy0[index], self._qz0[index],  # qx0,qy0,qz0
                                      self._qx1[index], self._qy1[index], self._qz1[index], # qx1,qy1,qz1
                                      self._L0[index], self._L1[index],                     # L0, L1
                                      self._labels[index], self._title[index],              # name, title
                                      self._load_id[index],                                 # load_name
                                      self._component,                                      # component
                                      self._system[index], self._complex[index]))           # system, load_complex
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
            self._L0.pop(index)
            self._qx0.pop(index)
            self._qy0.pop(index)
            self._qz0.pop(index)
            #
            self._L1.pop(index)
            self._qx1.pop(index)
            self._qy1.pop(index)
            self._qz1.pop(index)
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
    @property
    def df(self):
        """ """
        #
        db = DBframework()
        #
        data = {'load_name': self._load_id,
                'load_type': ['basic' for item in self._labels],
                'load_id': [idx + 1 for idx, item in enumerate(self._labels)],
                'load_system': self._system, #['global' if item == 0 else 'local'
                                #for item in self._system],
                'load_comment': self._title,
                'element_name':self._labels, 
                'L0':self._L0, 'qx0':self._qx0, 'qy0':self._qy0, 'qz0':self._qz0, 
                'L1':self._L1, 'qx1':self._qx1, 'qy1':self._qy1, 'qz1':self._qz1}        
        #
        cols = ['load_name',  'load_type',
                'load_id', 'load_system',
                'load_comment', 
                'element_name',
                'L0', 'qx0', 'qy0', 'qz0',
                'L1', 'qx1', 'qy1', 'qz1']
        df = db.DataFrame(data=data, columns=cols)
        #
        #df = df[['load_name', 'load_type', 'load_id', 'load_system', 'load_comment',
        #         'element_name',
        #        'L0', 'qx0', 'qy0', 'qz0',
        #        'L1', 'qx1', 'qy1', 'qz1']]
        #       
        #print('--->')
        return df
#
#
class BeamPointIM(NodeLoadMaster):
    __slots__ = ['_title', '_labels', '_load_name',
                 '_index', '_complex', '_component', 
                 '_fx', '_fy', '_fz', '_mx', '_my', '_mz',
                 '_L0', '_type']
    
    def __init__(self,  load_name:str|int,
                 component: int) -> None:
        """
        """
        super().__init__("point")
        self._load_name = load_name
        self._component = component
        self._L0: array = array('f', [])
    #
    def __setitem__(self, element_name:str|int,
                    point_load: list|dict) -> None:
        """
        """
        self._labels.append(element_name)
        self._load_id.append(self._load_name)
        self._complex.append(0)
        #
        Ltype =  point_load.pop(0)
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
        self._L0.append(point_load[6])
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
            points.append(PointBeam(self._fx[index], self._fy[index], self._fz[index],  # fx,fy,fz
                                    self._mx[index], self._my[index], self._mz[index],  # mx,my,mz
                                    self._L0[index],                                    # L0
                                    self._labels[index], self._title[index],            # name,title
                                    self._load_id[index],                               # load_name
                                    self._component,                                    # component_name
                                    self._system[index], self._complex[index]))         # system, load_complex
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
            self._L0.pop(_index)
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
    @property
    def df(self):
        """ """
        #
        db = DBframework()
        #
        data = {'load_name': self._load_id,
                'load_type': ['basic' for item in self._labels],
                'load_id': [idx + 1 for idx, item in enumerate(self._labels)],
                'load_system': self._system, #['global' if item == 0 else 'local'
                                #for item in self._system],
                'load_comment': self._title,
                'element_name':self._labels, 
                'L0':self._L0, 'Fx':self._fx, 'Fy':self._fy, 'Fz':self._fz, 
                'Mx':self._mx, 'My':self._my, 'Mz':self._mz}
        #
        cols = ['load_name',  'load_type',
                'load_id', 'load_system',
                'load_comment', 
                'element_name',
                'L0', 'Fx', 'Fy', 'Fz',
                'Mx', 'My', 'Mz']
        df = db.DataFrame(data=data, columns=cols)
        #
        #df = df[['load_name', 'load_type', 'load_id', 'load_system', 'load_comment',
        #         'element_name',
        #         'L0', 'Fx', 'Fy', 'Fz',
        #         'Mx', 'My', 'Mz']]
        #       
        #print('--->')
        return df    
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
            node_id = node[0]
            point_load = node[1:]
            super().__setitem__(node_id, point_load)
    #
    #
    @property
    def df(self):
        """nodes in dataframe format"""
        #
        db = DBframework()
        #
        data = {'load_name': self._load_id,
                'load_type': ['basic' for item in self._labels],
                #'load_id': [idx + 1 for idx, item in enumerate(self._labels)],
                'load_system': self._system, #['global' if item == 0 else 'local'
                           #for item in self._system],
                'load_comment': self._title,
                'element_name': self._beam, 
                'node_name':self._labels, 
                'Fx':self._fx, 'Fy':self._fy, 'Fz':self._fz, 
                'Mx':self._mx, 'My':self._my, 'Mz':self._mz}      
        #
        cols = ['load_name', 'load_type', 
                'load_system', 'load_comment', 
                'element_name', 'node_name', 
                'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        #
        df = db.DataFrame(data=data, columns=cols)
        #
        df = df[['load_name', 'load_type', 'load_comment', 'load_system',
                 'element_name', 'node_name',
                 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']]
        return df     
    #
    #
#
#

