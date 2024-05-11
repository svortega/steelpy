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
#from typing import NamedTuple


# package imports
# steelpy
# steelpy.f2uModel.load
#from ..process.nodes import NodeLoadMaster
from ..concept.node import NodeLoadIM

from steelpy.ufo.load.process.beam.beam import LineBeam, PointBeam
from steelpy.ufo.load.process.beam.main import BeamTypeBasic, BeamLineBasic, BeamPointBasic

#from steelpy.utils.units.buckingham import Number
from steelpy.utils.dataframe.main import DBframework
#from steelpy.utils.math.operations import trnsload

from steelpy.ufo.load.process.operations import (get_BeamLoad_dic,
                                                 get_BeamLoad_list_units,
                                                 get_BeamLine_dic,  
                                                 get_BeamNode_load,
                                                 get_BeamLine_load)
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
        beams = list(dict.fromkeys(self._labels))
        return len(beams)
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
    # ------------------
    #   
#
#
#
# ---------------------------------
#
class BeamLoadItemIM(BeamIMMaster):
    __slots__ = ['_labels', '_load', '_node_eq',
                 '_beam_id', '_ufo_beam', '_beam']

    def __init__(self, load_name:str|int, 
                 load_title:str,
                 component: int)-> None:
                 #beams) -> None:
        """
        """
        super().__init__()
        #self._ufo_beam = beams
        #self._load_name = load_name
        self._load = BeamLoadTypeIM(load_name=load_name,
                                    load_title=load_title,
                                    #ufo_beam=beams, 
                                    component=component)
        #
        #self._node_eq = BeamToNodeIM(load_name=load_name)

    #
    def __setitem__(self, beam_name: int|str,
                    beam_load: list) -> None:
        """
        """
        #
        if re.match(r"\b(point)\b", str(beam_load[0]), re.IGNORECASE):
            self._load._point[beam_name] = beam_load[1:]
            
        elif re.match(r"\b(line|u(niform(ly)?)?d(istributed)?l(oad)?)\b",
                       str(beam_load[0]), re.IGNORECASE):
            self._load._line[beam_name] = beam_load[1:]
        
        else:
            raise IOError(f'Beam lod type {beam_load[0]} not implemented')
        #
        self._beam_id = beam_name
        self._labels.append(beam_name)
        #
    #
    def __getitem__(self, beam_name: int | str):
        """
        """
        self._labels.append(beam_name)
        return self._load(beam=beam_name)
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
#
#
#
#
#
class BeamLoadTypeIM(BeamTypeBasic):
    __slots__ = ['_system_flag','_beam_id', 'load_name',
                 '_line', '_point', '_load_title',
                 '_beam_id']

    def __init__(self, load_name:str|int,
                 load_title:str,
                 #ufo_beam, 
                 component: int):
        """
        """
        super().__init__()
        self.load_name = load_name
        self._load_title = load_title
        #self._ufo_beam = ufo_beam
        #
        self._line = BeamDistributedIM(load_name,
                                       #ufo_beam=ufo_beam, 
                                       component=component)
        #
        self._point = BeamPointIM(load_name,
                                  #ufo_beam=ufo_beam, 
                                  component=component)
        #self._node_eq = BeamToNodeIM(load_name)
    #
    def __call__(self, beam):
        """ """
        self._beam_id = beam
        #self._line(beam=self._beam)
        #self._point(beam=self._beam)
        return self
    #
    #
    #
    # ------------------
#
#
#
# ---------------------------------
# Beam Loading Calculation
#
class BeamDistributedIM(BeamLineBasic):
    __slots__ = ['_type', '_labels', '_load_name', '_beam', 
                 '_index', '_complex', '_L', '_load_id', 
                 '_L0', '_qx0', '_qy0', '_qz0', '_qt0', 
                 '_L1', '_qx1', '_qy1', '_qz1', '_qt1',
                 '_system', '_title', '_component'] # '_ufo_beam'

    def __init__(self, load_name:str|int,
                 #ufo_beam, 
                 component: int) -> None: 
        """
        """
        super().__init__()
        self._component = component
        self._load_name = load_name
        #self._ufo_beam = ufo_beam
        #
        self._labels: list = []
        self._title: list = []
        #self._complex: array = array("I", [])
        self._load_id: list = []
        self._system: array = array("I", [])        
        # ens 1
        self._L0: array = array("f", [])
        self._qx0: array = array("f", [])
        self._qy0: array = array("f", [])
        self._qz0: array = array("f", [])
        self._qt0: array = array("f", [])
        # end 2
        self._L1: array = array("f", [])
        self._qx1: array = array("f", [])
        self._qy1: array = array("f", [])
        self._qz1: array = array("f", [])
        self._qt1: array = array("f", [])
    #
    #
    def __setitem__(self, beam_name: int|str,
                    line_load: list|dict) -> None:
        """
                        |
             q0         | q1
        o------|        |----------o
        |                          |
        +  L0  +        +    L1    +
        +........... L ............+
        
        [q0x, q0y, q0z, q0t,
         q1x, q1y, q1z, q1t,
         L0, L1,
         system, comment]
        """
        line_load =  self._get_line(line_load)        
        #
        self._labels.append(beam_name)
        self._load_id.append(self._load_name)
        # end 1
        self._qx0.append(line_load[0])
        self._qy0.append(line_load[1])
        self._qz0.append(line_load[2])
        self._qt0.append(line_load[3])
        # end 2
        self._qx1.append(line_load[4])
        self._qy1.append(line_load[5])
        self._qz1.append(line_load[6])
        self._qt1.append(line_load[7])
        # distance from ends
        self._L0.append(line_load[8])
        self._L1.append(line_load[9])
        #
        self._system.append(line_load[10])
        self._title.append(line_load[11])
    #
    def __getitem__(self, beam_name: int|str) :
        """
        """
        idx_list: list = [x for x, item in enumerate(self._labels)
                          if item == beam_name]
        #1/0
        udl_list: list = []
        for index in idx_list:
            udl_list.append(LineBeam(self._labels[index], # name
                                     self._title[index],  # title
                                     self._load_id[index],      # load_name
                                     self._component,           # component
                                     self._system[index],       # system
                                     #
                                     self._qx0[index], self._qy0[index],  # qx0,qy0
                                     self._qz0[index], self._qt0[index],  # qz0,qt0
                                     self._qx1[index], self._qy1[index],  # qx1,qy1
                                     self._qz1[index], self._qt1[index],  # qz1,qt1
                                     self._L0[index], self._L1[index]))   # L0, L1
        return udl_list
    #
    def _get_line(self, line_load: list|dict):
        """ get line load in beam local system"""
        
        if isinstance(line_load, (list, tuple)):
            try:
                udl = get_BeamLoad_list_units(line_load.copy())
            except AttributeError:
                udl = get_BeamLine_load(line_load)
        elif isinstance(line_load, dict):
            udl = get_BeamLine_dic(line_load)
        #
        load_type = udl.pop(0)
        title = udl.pop()
        #
        return [*udl, self._system_flag, title]
    #
    #
    def __delitem__(self, beam_name: int) -> None:
        """
        """
        index_list: list = [x for x, _item in enumerate(self._labels)
                            if _item == beam_name]
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
class BeamPointIM(BeamPointBasic):
    __slots__ = ['_title', '_labels', '_load_name', '_beam', 
                 '_index', '_complex', '_component', 
                 '_fx', '_fy', '_fz', '_mx', '_my', '_mz',
                 '_L0', '_type', '_system', '_load_id'] # '_ufo_beam'
    
    def __init__(self,  load_name:str|int,
                 #ufo_beam, 
                 component: int) -> None:
        """
        """
        super().__init__()
        self._load_name = load_name
        self._component = component
        #self._ufo_beam = ufo_beam        
        #
        self._labels: list = []
        self._title: list = []
        #self._complex: array = array("I", [])
        self._load_id: list = []
        self._system: array = array("I", [])
        #
        self._type = "point"
        self._fx: array = array('f', [])
        self._fy: array = array('f', [])
        self._fz: array = array('f', [])
        self._mx: array = array('f', [])
        self._my: array = array('f', [])
        self._mz: array = array('f', [])        
        #
        self._L0: array = array('f', [])
    #
    def __setitem__(self, beam_name:str|int,
                    point_load: list|dict) -> None:
        """
        """
        point_load =  self._get_point(point_load)        
        #
        self._labels.append(beam_name)
        self._load_id.append(self._load_name)
        #
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
    
    def __getitem__(self, beam_name:str|int):
        """
        """
        index_list: list = [x for x, item in enumerate(self._labels)
                            if item == beam_name]
        #
        points: list = []
        for index in index_list: 
            points.append(PointBeam(self._labels[index], self._title[index],            # name,title
                                    self._load_id[index],                               # load_name
                                    self._component,                                    # component_name
                                    self._system[index], #self._complex[index],         # system, load_complex
                                    #
                                    self._fx[index], self._fy[index], self._fz[index],  # fx,fy,fz
                                    self._mx[index], self._my[index], self._mz[index],  # mx,my,mz
                                    self._L0[index]))                                   # L0       
        return points
    #
    #
    def __delitem__(self, beam_name:int|str) -> None:
        """
        """
        indexes = [i for i, x in enumerate(self._labels)
                   if x == beam_name]
        
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
            #self._complex.pop(_index)    
    #
    # -----------------------------------------------
    #
    # TODO : chekc if works for dict
    def _get_point(self, point_load: list|dict):
        """ get point load in beam local system"""
        # update inputs
        if isinstance(point_load, (list, tuple)):
            try:
                point = get_BeamLoad_list_units(point_load.copy())
            except AttributeError:
                point = get_BeamNode_load(point_load)            
        # update inputs
        elif isinstance(point_load, dict):
            point = get_BeamLoad_dic(point_load)       
        #
        load_type = point.pop(0)
        title = point.pop() 
        #
        return [*point, self._system_flag, title]
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
    @df.setter
    def df(self, values):
        """nodes in dataframe format"""
        #update
        self._labels.extend(values.node_name.tolist())
        self._load_id.extend(values.load_name.tolist())
        self._title.extend(values.load_title.tolist())
        self._system.extend([1 if item in ['local', 'member', 1] else 0
                             for item in values.system.tolist()])
        #self._complex.extend([0 for item in values.system.tolist()])
        #
        self._fx.extend(values.Fx.tolist())
        self._fy.extend(values.Fy.tolist())
        self._fz.extend(values.Fz.tolist())
        self._mx.extend(values.Mx.tolist())
        self._my.extend(values.My.tolist())
        self._mz.extend(values.Mz.tolist())
        #
        # beam point case
        try:
            self._L1.extend(values.L1.tolist())
        except AttributeError:
            pass
        #
        # beam node load conversion
        try:
            self._beam.extend(values.element_name.tolist())
        except AttributeError:
            pass
        #    
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

