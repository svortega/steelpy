#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from array import array
from dataclasses import dataclass
from collections.abc import Mapping
#import re

# package imports
#from steelpy.ufo.load.concept.beam import BeamLoadTypeIM
from steelpy.ufo.load.process.actions import SelfWeight
#from steelpy.f2uModel.load.inmemory.combination import BasicLoad
from steelpy.utils.units.buckingham import Number
from steelpy.utils.math.operations import trnsload, linspace

from steelpy.ufo.load.process.operations import (check_point_dic, check_list_units,
                                                 check_beam_dic,  
                                                 get_beam_node_load, get_beam_udl_load)

#
class BeamBasicLoad:
    """
    FE Load Cases

    LoadType
        |_ name
        |_ number
        |_ basic
        |_ combination_level
        |_ time_series
        |_
        |_ temperature

    **Parameters**:
      :number:  integer internal number
      :name:  string node external name
    """
    __slots__ = ['_cls',  '_selfweight', #'gravity',
                 '_basic', 'combination']

    def __init__(self, beam):
        """
        """
        self._cls = beam
        #self._system_flag = 1 # local
        #self.gravity = 9.80665  # m/s^2
        #
        self._selfweight = SelfWeight()
        #self._basic = BeamLoadTypeIM(load_name="beam_load",
        #                             load_title="basic",
        #                             component=beam.name)
        #self._basic.local_system()
        self.combination = BeamLoadComb("beam_load")
    #
    @property
    def basic(self):
        """return basic beam load"""
        return self._basic(beam=self._cls)
    #
    #@property
    #def combination(self):
    #    """return beam load combinations"""
    #    return self._combinations

    def __str__(self, units: str = "si") -> str:
        """ """
        unit_lenght = " m"
        unit_force = "  N"
        unit_bm = "N*m"
        unit_fl = "N/m"
        output = "\n"
        output += "{:}\n".format(80 * "_")
        output += "\n"
        output += f"{35 * ' '}BASIC LOAD\n"
        output += "\n"
        output += f"--- Beam \n"
        output += f"Element Name{6 * ' '}L1[{unit_lenght}] qx1[{unit_fl}] qy1[{unit_fl}] qz1[{unit_fl}] System Complex\n"
        output += f"Line Load{9 * ' '}L2[{unit_lenght}] qx2[{unit_fl}] qy2[{unit_fl}] qz2[{unit_fl}] Comment\n"
        output += "\n"
        output += f"--- Beam \n"
        output += f"Element Name{6 * ' '}L1[{unit_lenght}] fx [{unit_force}] fy [{unit_force}] fz [{unit_force}] System Complex\n"
        output += f"Point Load{15 * ' '}mx [{unit_bm}] my [{unit_bm}] mz [{unit_bm}] Comment\n"
        output += "\n"
        output += f"--- Gravity/Selfweight\n"
        #output += f"Element Number{4 * ' '}L1[{unit_lenght}] fx [{unit_force}] fy [{unit_force}] fz [{unit_force}] System Complex\n"
        output += "\n"
        output += "{:}\n".format(80 * ".")
        output += "\n"
        output += f"--- Beam \n"
        output += self._basic.__str__()
        output += "\n"
        return output
#
#
@dataclass
class BeamLoadComb:
    
    def __init__(self, b):
        """
        """
        self._labels = []
        self._factor: array = array("f", [])
    #
    def __setitem__(self, load_name: int|str,
                    factor: float) -> None:
        """ """
        
        self._labels.append(load_name)
        self._factor.append(factor)
    
    def __getitem__(self, element_name: int|str) :
        """
        """
        pass
#
#
class BeamBasicItem(Mapping):
    
    def __init__(self, beam_name:str):
        """
        """
        self._labels = {}
    
    
    def __iter__(self):
        """
        """
        items = list(dict.fromkeys(self._labels))
        return iter(items)

    def __contains__(self, value) -> bool:
        return value in self._labels

    def __len__(self) -> int:
        return len(self._labels)
#
#
#
# ---------------------------------
#
#
@dataclass
class BeamTypeBasic:
    """ """
    def __init__(self):
        """
        """
        self._system_flag = 0  # Global system default
    #
    # -----------------------------------------------
    #
    def __call__(self, beam):
        """ """
        self._beam = beam
        self._beam_id = beam.name
        #self._line(beam=self._beam)
        #self._point(beam=self._beam)
        return self
    #
    # -----------------------------------------------
    #
    #
    @property
    def line(self):
        """
        Linear Varying Load (lvl) - Non Uniformly Distributed Load
        
        value : [qx1, qy1, qz1, qx2, qy2, qz2, L1, L2]
    
                        |
             q0         | q1
        o------|        |----------o
        |                          |
        +  L0  +        +    L1    +

        """
        return self._line #[beam_name]

    @line.setter
    def line(self, values: list):
        """
        Linear Varying Load (lvl) - Non Uniformly Distributed Load
                value : [qx1, qy1, qz1, qx2, qy2, qz2, L1, L2]
    
                        |
             q0         | q1
        o------|        |----------o
        |                          |
        +  L0  +        +    L1    +
        """
        #beam_name = self._beam.name
        beam_name = self._beam_id
        #
        #
        if isinstance(values, dict):
            #load = self._get_line(values)
            #load.insert(0, 'load')
            self._line[beam_name] = values
    
        elif isinstance(values[0], (list, tuple)):
            for item in values:
                #load =  self._get_line(item)
                #load.insert(0, 'load')
                self._line[beam_name] = item
        else:
            #load =  self._get_line(values)
            #load.insert(0, 'load')
            self._line[beam_name] = values
    #
    #
    # ------------------
    #
    @property
    def point(self):
        """ Concentrated force """
        return self._point #[beam_name]

    @point.setter
    def point(self, values: list):
        """
        Concentrated force
        """
        #beam_name = self._beam.name
        beam_name = self._beam_id
        #
        if isinstance(values, dict):
            #load = self._get_point(values)
            #load.insert(0, 'force')
            self._point[beam_name] = values
    
        elif isinstance(values[0], (list, tuple)):
            for item in values:
                #value.insert(0, self._Lbeam)
                #load = get_beam_point_load(load=values)
                #load = self._get_point(item)
                #load.insert(0, 'force')
                self._point[beam_name] = item
        else:
            #values.insert(0, self._Lbeam)
            #load = get_beam_point_load(load=values)
            #load =  self._get_point(values)
            #load.insert(0, 'force')
            self._point[beam_name] = values
    #
    #     
    #
    # -----------------------------------------------
    #
    # ---------------------------------
    #
    @property
    def coordinate_system(self):
        if self._system_flag != 0:
            return "local"
        return "global"
    #
    @coordinate_system.setter
    def coordinate_system(self, system:str|int):
        """
        Coordinate system for load : global or local (member)
        """
        self._system_flag = 0
        if system in ['local', 'member', 1]:
            self._system_flag = 1
        #
        self._line.coordinate_system = self._system_flag
        self._point.coordinate_system = self._system_flag        
    #    
    #
    def local_system(self):
        """set load beam local system"""
        self._system_flag = 1
        self._line.coordinate_system = self._system_flag
        self._point.coordinate_system = self._system_flag
        return "local"

    # @property
    def global_system(self):
        """set load beam global system"""
        self._system_flag = 0
        self._line.coordinate_system = self._system_flag
        self._point.coordinate_system = self._system_flag
        return "global"
    #
    # -----------------------------------------------
    #
    def beam_function(self, beams, steps:int = 10):
        """ """
        #
        #beamfun = defaultdict(list)
        loadfun = []
        # line load
        for key, items in self._line.items():
            beam = beams[key]
            mat = beam.material
            sec = beam.section.properties()
            Lsteps = linspace(start=0, stop=beam.L, num=steps+1, endpoint=True)
            for bitem in items:
                lout = bitem.Fx(x=Lsteps, L=beam.L,
                                E=mat.E, G=mat.G, 
                                Iy=sec.Iy, Iz=sec.Iz,
                                J=sec.J, Cw=sec.Cw, Area=sec.area)
                #beamfun[key].extend(lout)
                #beamfun[key].append(lout)
                #
                loadfun.extend(lout)
                #loadfun.append([key, bitem.name, lout])
                #loadfun.extend([[bitem.name, 'local', key, *step]
                #                for step in lout])
                #for step in lout:
                #    # load_name, load_title, [Fx, Fy, Fz, Mx, My, Mz]
                #    loadfun.append([bitem.name, 'local', key, *step])
        # point load
        for key, items in self._point.items():
            beam = beams[key]
            mat = beam.material
            sec = beam.section.properties()
            Lsteps = linspace(start=0, stop=beam.L, num=steps+1, endpoint=True)
            for bitem in items:
                lout = bitem.Fx(x=Lsteps, L=beam.L,
                                E=mat.E, G=mat.G, 
                                Iy=sec.Iy, Iz=sec.Iz,
                                J=sec.J, Cw=sec.Cw, Area=sec.area)
                #beamfun[key].append([bitem.name, lout])
                #beamfun[key].extend(lout)
                #beamfun[key].append(lout)
                #
                loadfun.extend(lout)
                #loadfun.append([key, bitem.name, lout])
                #loadfun.extend([[bitem.name, 'local', key, *step]
                #                for step in lout])                
                #for step in lout:
                #    loadfun.append([bitem.name, 'local', key, *step])
        #print('---')
        return loadfun # beamfun # 
    #
    #
    # -----------------------------------------------
    #
    def __str__(self, units: str = "si") -> str:
        """ """
        #unit_lenght = " m"
        #unit_force = "  N"
        #unit_bm = "N*m"
        #unit_fl = "N/m"
        output = ""
        #if header:
        #    output += "\n"
        #    output += f"--- Beam \n"
        #    output += f"Element Name{6 * ' '}L1[{unit_lenght}] qx1[{unit_fl}] qy1[{unit_fl}] qz1[{unit_fl}] System Complex\n"
        #    output += f"Line Load{9 * ' '}L2[{unit_lenght}] qx2[{unit_fl}] qy2[{unit_fl}] qz2[{unit_fl}] Comment\n"
        #    output += "\n"
        #    output += f"--- Beam \n"
        #    output += f"Element Name{6 * ' '}L1[{unit_lenght}] fx [{unit_force}] fy [{unit_force}] fz [{unit_force}] System Complex\n"
        #    output += f"Point Load{15 * ' '}mx [{unit_bm}] my [{unit_bm}] mz [{unit_bm}] Comment\n"
        #    output += "\n"
        #    output += "{:}\n".format(80 * ".")
        #    output += "\n"
        # 1/0
        # output += "--- Beam Line Load\n"
        output += self._line.__str__()
        # output += "--- Beam Point Load\n"
        output += self._point.__str__()
        # output += self._nodal_displacement.__str__()
        return output
#    
#
# ---------------------------------
#
class BeamLoadBasic(Mapping):
    __slots__ = ['_system_flag', '_beam']
    #
    # ---------------------------------
    #
    #def __call__(self, beam):
    #    """ """
    #    self._beam = beam
    #    return self    
    #
    # ---------------------------------
    #
    def __init__(self):
        """
        """
        self._system_flag = 0  # Global system default       
    
    def __len__(self) -> int:
        return len(self._labels)

    def __iter__(self):
        items = list(dict.fromkeys(self._labels))
        return iter(items)
    
    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    # ---------------------------------
    # TODO : get number from database
    def get_number(self, start:int=0):
        """
        """
        try:
            n = max(self._labels)
        except ValueError:
            n = start
        #
        while True:
            n += 1
            yield n
    #
    # ---------------------------------
    #
    @property
    def coordinate_system(self):
        if self._system_flag != 0:
            return "local"
        return "global"
    #
    @coordinate_system.setter
    def coordinate_system(self, system:str|int):
        """
        Coordinate system for load : global or local (member)
        """
        self._system_flag = 0
        if system in ['local', 'member', 1]:
            self._system_flag = 1
    #
#
#
class BeamLineBasic(BeamLoadBasic):
    __slots__ = ['_system_flag']
    #
    #
    #def __init__(self):
    #    """
    #    """
    #    super().__init__()
    #
    # -----------------------------------------------
    #
    # TODO : chekc if works for dict
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
            nload = [*udl[0:4], 0, 0,
                     *udl[4:8], 0, 0]
            nload = trnsload(nload, self._beam.T3D())
            #nload = [*nload[:4], *nload[6:10]] 
            return [*nload[:4],    # end 1 [qx, qy, qz, qt]
                    *nload[6:10],  # end 2 [qx, qy, qz, qt]
                    *udl[8:],      # [L0, L1]
                    1, title]      # Local system, title
#
#
class BeamPointBasic(BeamLoadBasic):
    __slots__ = ['_system_flag']
    #
    #
    def __init__(self):
        """
        """
        super().__init__()
    #
    # -----------------------------------------------
    #
    # TODO : chekc if works for dict
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
    #
    def __str__(self) -> str:
        """ """
        output = ""
        nodes = list(dict.fromkeys(self._labels))
        for node in nodes:
            items = self.__getitem__(node)
            for item in items:
                output += item.__str__()
                # print('---')
        return output
#