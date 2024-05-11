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
#from steelpy.utils.units.buckingham import Number
from steelpy.utils.math.operations import trnsload, linspace

from steelpy.ufo.load.process.operations import (get_BeamLoad_dic,
                                                 get_BeamLoad_list_units,
                                                 get_BeamLine_dic,  
                                                 get_BeamNode_load,
                                                 get_BeamLine_load)

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
        self._system_flag = 1  # Local system default
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
        beam_name = self._beam_id
        return self._line[beam_name]

    @line.setter
    def line(self, values: list):
        """
        Linear Varying Load (lvl) - Non Uniformly Distributed Load
        
        value : [qx1, qy1, qz1, qt1,  qx2, qy2, qz2, qt2,  L1, L2, title]
    
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
            values.update({'type':'line', })
            self._line[beam_name] = values
    
        elif isinstance(values[0], (list, tuple)):
            for item in values:
                item.insert(0, 'line')
                self._line[beam_name] = item
        else:
            values.insert(0, 'line')
            self._line[beam_name] = values
    #
    #
    # -----------------------------------------------
    #
    @property
    def point(self):
        """ Concentrated force """
        beam_name = self._beam_id
        return self._point[beam_name]

    @point.setter
    def point(self, values: list):
        """
        Concentrated force
        """
        #beam_name = self._beam.name
        beam_name = self._beam_id
        #
        if isinstance(values, dict):
            values.update({'type':'point', })
            self._point[beam_name] = values
    
        elif isinstance(values[0], (list, tuple)):
            for item in values:
                item.insert(0, 'point')
                self._point[beam_name] = item
        else:
            values.insert(0, 'point')
            self._point[beam_name] = values
    #  
    #
    # -----------------------------------------------
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
        self._system_flag = 1
        if system in ['global' , 0]:
        #if system in ['local', 'member', 1]:
            self._system_flag = 0
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
        self._system_flag = 1  # Local system default       
    
    def __len__(self) -> int:
        items = list(set(dict.fromkeys(self._labels)))
        return len(items)

    def __iter__(self):
        items = list(set(dict.fromkeys(self._labels)))
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
        self._system_flag = 1
        if system in ['global' , 0]:
            self._system_flag = 0
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
    def _get_line(self, line_load: list|tuple|dict):
        """ get line load in beam local system"""
        #
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
    #def __init__(self):
    #    """
    #    """
    #    super().__init__()
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