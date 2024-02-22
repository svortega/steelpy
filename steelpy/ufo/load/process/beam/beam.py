#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
from collections import namedtuple
#from collections import defaultdict
from dataclasses import dataclass
from operator import sub, add
from dataclasses import dataclass
from typing import NamedTuple


# package imports
#
from steelpy.utils.units.buckingham import Number
from steelpy.utils.math.operations import trnsload
from steelpy.utils.math.operations import linspace
# steelpy.
from steelpy.ufo.load.process.operations import(check_list_units,
                                                check_beam_dic,
                                                check_point_dic,
                                                get_beam_node_load,
                                                get_beam_udl_load)
#
#
# steelpy.
from steelpy.trave.beam.operation import BeamTorsion, BeamAxial, BeamBending
#from steelpy.trave.beam.roark.chapter10.table103_point import BTWSupPoint
from steelpy.trave.beam.pilkey.chapter14.table_C import BTOpenSupports
from steelpy.trave.beam.operation import BeamBasic
#
import numpy as np
#
# ---------------------------------
#
@dataclass #(kw_only=True)
class BeamLoadClass:
    """ """
    name: str | int
    title: str
    load_name: str | int
    component_name: str | int
    system: int
    #
    # -------------------------
    #
    def local_system(self):
        """
        return q load in local beam system

        retunr:
        [qx0, qy0, qz0, qx1, qy1, qz1]
        """
        try:  # local system
            1 / self.system
        except TypeError:
            if self.system.lower() == 'global':
                raise IOError('line line should be in beam local system')
        except ZeroDivisionError:
            raise IOError('line line should be in beam local system')
    #    
    @property
    def coordinate_system(self):
        """
        Local  : 1
        Global : 0
        """
        try:  # local system
            1 / self.system
            return "local"
        except TypeError:
            return self.system.lower()
        except ZeroDivisionError:
            return 'global'
    #
    # -------------------------
    # Load Process Formulas
    # -------------------------
    #
    def fer_beam(self, beam, system :int) -> list:
        """ """
        L =beam.L
        material = beam.material
        section = beam.section.properties()
        nodes = beam.nodes
        bops = BeamBasic(L=L,
                         area=section.area,
                         Iy=section.Iy, Iz=section.Iz,
                         J=section.J, E=material.E,
                         G=material.G, Cw=section.Cw)
        #
        bops.supports(end1=nodes[0].fixity, end2=nodes[1].fixity)
        #bops.supports(end1='fixed', end2='fixed')
        #
        Fbar = self.Fx(x=L, L=L,
                       E=material.E, G=material.G, 
                       Iy=section.Iy, Iz=section.Iz,
                       J=section.J, Cw=section.Cw,
                       Area=section.area)
        #
        # Axial   [P, blank, blank, u]
        # Torsion [T, Phi, Psi, B, Tw]
        # Bending [V, M, theta, w]
        #
        end1, end2 = bops._reaction(Fbar=Fbar)
        #
        # Force = [fx, fy, fz, mx, my, mz]
        # Node 1
        Fend1 = [end1.x[0], -1 * end1.y[0], -1 *end1.z[0], # Fx, Fy, Fz
                 end1.t[0], 1 * end1.z[1], 1 * end1.y[1]]  # Mx, My, Mz
        # Node 2
        Fend2 = [end2.x[0], -1 * end2.y[0], -1 *end2.z[0], # Fx, Fy, Fz
                 end2.t[0], 1 * end2.z[1], 1 * end2.y[1]]  # Mx, My, Mz
        #
        # Displacement = [x, y, z, rx, ry, rz]
        # Node 1
        Rend1 = [end1.x[3], -1 * end1.y[3], -1 * end1.z[3],  # x, y, z,
                 end1.t[1], -1 * end1.z[2], -1 * end1.y[2]]  # rx, ry, rz
        # Node 2
        Rend2 = [end2.x[3], 1 * end2.y[3], 1 * end2.z[3],  # x, y, z,
                 end2.t[1], 1 * end2.z[2], 1 * end2.y[2]]  # rx, ry, rz
        #
        try: # local
            1 / system
        except ZeroDivisionError:
            pload = [*Fend1, *Fend2]
            pload = trnsload(pload, beam.T3D())
            Fend1 = list(pload[:6])
            Fend2 = list(pload[6:])
            #
            pload = [*Rend1, *Rend2]
            pload = trnsload(pload, beam.T3D())
            Rend1 = list(pload[:6])
            Rend2 = list(pload[6:])
        #
        Pend1 = [*Fend1, *Rend1, *end1.t[2:]]
        Pend2 = [*Fend2, *Rend2, *end2.t[2:]]
        #
        FER_out = [self.load_name, self.component_name,
                   self.title, self.load_type, 
                   'basic', self.system,
                   self.name, Pend1, Pend2]
        #1 / 0
        return FER_out
    #  
#
#
@dataclass
class LineBeam(BeamLoadClass):
    """
    """
    qx0: float
    qy0: float
    qz0: float
    qt0: float
    #
    qx1: float
    qy1: float
    qz1: float
    qt1: float
    #
    L0: float
    L1: float
    #
    load_complex: int = 0
    load_type: str = "Line Load"
    #
    #
    def __str__(self, units:str="si") -> str:
        """ """
        output = (f"{str(self.name):12s} "
                  f"{self.L0: 1.3e} {self.qx0: 1.3e} "
                  f"{self.qy0: 1.3e} {self.qz0: 1.3e} "
                  f"{self.coordinate_system.upper():6s} "
                  f"{self.load_complex}\n")
        
        if (load_title:= self.title) == "NULL":
            load_title = ""        
        step = self.load_type
        output += (f"{step[:12]:12s} {self.L1: 1.3e} {self.qx1: 1.3e} "
                   f"{self.qy1: 1.3e} {self.qz1: 1.3e} "
                   f"{str(load_title)}\n")
        return output
    #
    # -------------------------
    # Load Process Formulas
    # ------------------------
    #
    def Fx(self, x:float|list, L:float,
           E:float, G: float,
           Iy:float, Iz:float,
           J: float, Cw: float, Area: float):
        """
        Beam load local system

        x : distance from end 1
        L : beam length
        E : Elastic module
        Iy : In plane
        Iz : Out plane
        R :

        return:
        [load_name, component,
        load_title, load_type,
        load_level, load_system,
        beam_name, x, Fx, Fy, Fz]
        Fx, Fy, Fz --> [V, M, theta, w] Note positive load is upwards
        """
        # convert load to beam local
        self.local_system()
        #
        # Axial [P, blank, blank, u]
        Fx = BeamAxial(L=L, E=E, A=Area)
        # Torsion [T, B, psi, phi, Tw]
        #Mx_blank = [0, 0, 0, 0, 0]
        Mx = BeamTorsion(E=E, L=L, G=G, J=J, Cw=Cw)
        #Mx.torque(T=self.mx, L1=self.L0)        
        #
        # In plane [V, M, theta, w]
        in_plane =  BeamBending(L=L, E=E, I=Iy)
        # Out plane [V, M, theta, w]
        out_plane = BeamBending(L=L, E=E, I=Iz)
        #
        # [Fx, Fy, Fz, Mx, My, Mz]
        if isinstance(x, (list, tuple)):
            Fx_out = []
            for xstep in x:
                Fx_out.append([self.load_name, self.component_name,
                               self.title, self.load_type, 'basic',
                               self.coordinate_system,
                               self.name, xstep,
                               np.array(Fx.axial(x=xstep,
                                                 P=self.qx0, L1=self.L0,
                                                 P2=self.qx1, L2=self.L1)),
                               # TODO : torsion
                               #np.array(Mx_blank),
                               np.array(Mx.torque(x=xstep,
                                                  T=self.qt0, L1=self.L0,
                                                  T2=self.qt1, L2=self.L1)),
                               #
                               np.array(in_plane.line(x=xstep,
                                                      q1=self.qy0, q2=self.qy1,
                                                      L1=self.L0, L2=self.L1)),
                               np.array(out_plane.line(x=xstep,
                                                       q1=self.qz0, q2=self.qz1,
                                                       L1=self.L0, L2=self.L1))])
        else:
            Fx_out = [self.load_name, self.component_name,
                      self.title, self.load_type, 'basic',
                      self.coordinate_system,
                      self.name, x,
                      np.array(Fx.axial(x=x,
                                        P=self.qx0, L1=self.L0,
                                        P2=self.qx1, L2=self.L1)),
                      # TODO : torsion
                      #np.array(Mx_blank),
                      np.array(Mx.torque(x=x,
                                         T=self.qt0, L1=self.L0,
                                         T2=self.qt1, L2=self.L1)),
                      #
                      np.array(in_plane.line(x=x,
                                             q1=self.qy0, q2=self.qy1,
                                             L1=self.L0, L2=self.L1)),
                      np.array(out_plane.line(x=x,
                                              q1=self.qz0, q2=self.qz1,
                                              L1=self.L0, L2=self.L1))]
        #
        # [load_name, load_title, load_type, load_system, 
        # beam_number, x, Fx, Fy, Fz]
        return Fx_out
    #
    #
    def _FaX(self, L: float):
        """Beam Axial load"""
        # Calculate area
        h = L - self.L0 - self.L1
        m = (self.qx0 + self.qx1) / 2.0
        area = m * h
        try:
            if self.qx0 > self.qx1:
                Lcog = (self.qx0 + 2*self.qx1)/(self.qx0 + self.qx1) * h / 3.0
            else:
                Lcog = h - (2*self.qx0 + self.qx1)/(self.qx0 + self.qx1) * h / 3.0
        except ZeroDivisionError:
            return 0, 0
        #
        h = area / m
        qload = area / h
        L1 = max(self.L0 + Lcog - h / 2.0, 0)
        L2 = max(L - (self.L0 + Lcog + h / 2.0), 0)
        #
        # Axial
        Fa = distributed_axial(w=qload, L=L,
                               l1=L1, l2=L2)
        #       
        return Fa    
    #
    #
#
#
#Pload = namedtuple('PointLoad', ['fx', 'fy', 'fz', 'mx', 'my', 'mz'])
#
@dataclass
class PointBeam(BeamLoadClass):
    """
    """
    fx: float
    fy: float
    fz: float
    mx: float
    my: float
    mz: float
    #
    L0: float
    #
    load_complex: int = 0
    load_type: str = "Point Load"
    #
    #
    def __str__(self, units:str="si") -> str:
        """ """
        output  = (f"{str(self.name):12s} "
                   f"{self.L0: 1.3e} {self.fx: 1.3e} "
                   f"{self.fy: 1.3e} {self.fz: 1.3e} "
                   f"{self.coordinate_system.upper():6s} "
                   f"{self.load_complex}\n")
        
        if (load_title:= self.title) == "NULL":
            load_title = ""        
        step = self.load_type
        output += (f"{step[:12]:12s} {10*' '} {self.mx: 1.3e} "
                   f"{self.my: 1.3e} {self.mz: 1.3e} "
                   f"{str(load_title)}\n")
        return output
    #
    #
    # -------------------------
    # Load Process Formulas
    # ------------------------    
    #
    #
    #def Fa(self, L: float):
    #    """Beam Axial load"""
    #    # Axial 
    #    Fa = axial_load(W=self.fx, L=L, l1=self.L0)
    #    return Fa
    #
    #
    def Fx(self, x:float|list, L: float,
           E:float, G: float, 
           Iy:float, Iz:float,
           J: float, Cw: float, Area: float) -> list:
        """
        Beam load local system

        x : distance from end 1
        L : beam length
        E : Elastic module
        Iy : In plane
        Iz : Out plane
        R :

        Return:
        Fx : [load_name, component,
              load_title, load_type,
              load_level, load_system,
              beam_name, L_step,
              Axial, Torsion, Bending_inplane, Bending_outplane]
        
        Axial, Bending_inplane, Bending_outplane = [V, M, w, theta]
        Torsion = [T, B, psi, phi, Tw]
        """
        
        # convert load to beam local
        self.local_system()
        #
        # Axial [P, blank, blank, u]
        Fx =  BeamAxial(L=L, E=E, A=Area)
        # Torsion [T, B, psi, phi, Tw]
        Mx = BeamTorsion(E=E, L=L, G=G, J=J, Cw=Cw)
        # In plane [V, M, theta, w]
        Finplane =  BeamBending(L=L, E=E, I=Iy)
        # Out plane [V, M, theta, w]
        F_outplane = BeamBending(L=L, E=E, I=Iz)
        #
        # [Fx, Fy, Fz, Mx, My, Mz]
        if isinstance(x, (list,tuple)):
            Fx_out = []
            for xstep in x:
                Axial =  Fx.axial(x=xstep, P=self.fx, L1=self.L0)
                #Axial = Fx.loading_function(x=xstep)
                Torsion =  Mx.torque(x=xstep, T=self.mx, L1=self.L0)
                #Torsion = Mx.loading_function(x=xstep, J=J, Cw=Cw)
                Finp = Finplane.point(x=xstep, P=self.fy, M=self.mz, L1=self.L0)
                #Finp = list(map(sum, zip(F_inplane(xstep, E, Iy), M_outplane(xstep, E, Iy))))
                Foutp = F_outplane.point(x=xstep, P=self.fz, M=self.my, L1=self.L0)
                #Foutp = list(map(sum, zip(F_outplane(xstep, E, Iz), M_inplane(xstep, E, Iz))))
                #
                Fx_out.append([self.load_name, self.component_name,
                               self.title, self.load_type, 'basic',
                               self.coordinate_system,
                               self.name, xstep,
                               np.array(Axial), np.array(Torsion),
                               np.array(Finp), np.array(Foutp)])
        else:
            Axial =  Fx.axial(x=x, P=self.fx, L1=self.L0)
            #Axial = Fx.loading_function(x=x)
            # torsion = [FT, FB, Fpsi, Fphi]
            #Torsion = Mx.loading_function(x=x, J=J, Cw=Cw)
            Torsion =  Mx.torque(x=x, T=self.mx, L1=self.L0)
            #
            Finp = Finplane.point(x=x, P=self.fy, M=self.mz, L1=self.L0)
            #Finp = list(map(sum, zip(F_inplane(x, E, Iy), M_outplane(x, E, Iy))))
            Foutp = F_outplane.point(x=x, P=self.fz, M=self.my, L1=self.L0)
            #Foutp = list(map(sum, zip(F_outplane(x, E, Iz), M_inplane(x, E, Iz))))
            #
            Fx_out = [self.load_name, self.component_name,
                      self.title, self.load_type, 'basic',
                      self.coordinate_system, 
                      self.name, x,
                      np.array(Axial), np.array(Torsion),
                      np.array(Finp), np.array(Foutp)]
        # [axial, torsion, bending_inplane, bending_outplane]
        return Fx_out
#
#
#
#
# ---------------------------------
#
#
@dataclass
class BeamLoadBasic:
    """ """
    def __init__(self):
        """
        """
        self._system_flag = 0  # Global system default
    #
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
    # -----------------------------------------------
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
#
@dataclass
class BeamLoadXXX(BeamLoadBasic):
    __slots__ = ['_line', '_point', '_system_flag']

    def __init__(self): #, beams
        """
        """
        #self._system_flag = 0  # Global system default
        super().__init__()
    #
    #def __call__(self, beam_name: int | str):
    #    """ """
    #    self._beam_id = beam_name
    #    return self
    #
    # ------------------
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
        #beam_name = self._beam.name
        #1 / 0
        #beam_name = self._beam_id
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
        beam_name = self._beam.name
        #beam_name = self._beam_id
        #
        if isinstance(values, dict):
            load = self._get_line(values)
            load.insert(0, 'load')
            self._line[beam_name] = load
    
        elif isinstance(values[0], (list, tuple)):
            for item in values:
                load =  self._get_line(item)
                load.insert(0, 'load')
                self._line[beam_name] = load
        else:
            load =  self._get_line(values)
            load.insert(0, 'load')
            self._line[beam_name] = load
    #
    #def _get_line(self, line_load: list|dict):
    #    """ get line load in beam local system"""
    #    #
    #    # update inputs
    #    if isinstance(line_load, dict):
    #        udl = check_beam_dic(line_load)
    #        title = udl.pop()
    #        
    #    elif isinstance(line_load[-1], str):
    #        title = line_load.pop()
    #        if isinstance(line_load[0], Number):
    #            udl = check_list_units(line_load)
    #        else:
    #            udl = get_beam_udl_load(line_load)
    #    else:
    #        title ='NULL'
    #        udl = get_beam_udl_load(line_load)
    #    #
    #    1 / 0
    #    # get system local = 1
    #    try:
    #        1 / self._system_flag
    #        return [*udl, 1, title]
    #    except ZeroDivisionError:
    #        # local nodal loading
    #        nload = [*udl[:3], 0, 0, 0,
    #                 *udl[3:6], 0, 0, 0,]
    #        nload = trnsload(nload, self._beam.T3D())
    #        nload = [*nload[:3], *nload[6:9]] 
    #        return [*nload, *udl[6:], 1, title]
    #
    # ------------------
    #
    @property
    def point(self):
        """ Concentrated force """
        #beam_name = self._beam.name
        #beam_name = self._beam_id
        #1 / 0
        return self._point #[beam_name]

    @point.setter
    def point(self, values: list):
        """
        Concentrated force
        """
        beam_name = self._beam.name
        #beam_name = self._beam_id
        #
        if isinstance(values, dict):
            load = self._get_point(values)
            load.insert(0, 'force')
            self._point[beam_name] = load
    
        elif isinstance(values[0], (list, tuple)):
            for item in values:
                #value.insert(0, self._Lbeam)
                #load = get_beam_point_load(load=values)
                load = self._get_point(item)
                load.insert(0, 'force')
                self._point[beam_name] = load
        else:
            #values.insert(0, self._Lbeam)
            #load = get_beam_point_load(load=values)
            load =  self._get_point(values)
            load.insert(0, 'force')
            self._point[beam_name] = load
    #
    #def _get_point(self, point_load: list|dict):
    #    """ get point load in beam local system"""
    #    # update inputs
    #    if isinstance(point_load, dict):
    #        point = check_point_dic(point_load)
    #        title = point.pop()
    #    
    #    elif isinstance(point_load[-1], str):
    #        title = point_load.pop()
    #        if isinstance(point_load[0], Number):
    #            point = check_list_units(point_load)
    #        else:
    #            point = get_beam_node_load(point_load)
    #    
    #    else:
    #        title = 'NULL'
    #        point = get_beam_node_load(point_load)
    #    #
    #    # get system local = 1
    #    try: # Local system
    #        1 / self._system_flag
    #        return [*point, 1, title]
    #    except ZeroDivisionError: # global to local system
    #        pload = [*point[:6], 0, 0, 0, 0, 0, 0]
    #        pload = trnsload(pload, self._beam.T3D())
    #        return [*pload[:6], point[6], 1, title]
    #
    # ------------------
    # Load Process
    # ------------------
    #
    def reactions(self):
        """Calculate bean reacition according to boundaries"""
        beam = self._beam
        r1 = 8 * [0]
        r2 = 8 * [0]
        b2n = []
        # line loadreactions
        for key, item in self._line.items():
            #print(f'line load {key}')
            for bload in item:
                res = beam.reactions(load=bload)
                r1 = list(map(add, r1, [*res[0][0], *res[0][1]]))
                r2 = list(map(add, r2, [*res[1][0], *res[1][1]]))
                #b2n.append([bload, item.name, item.title])
        # point load
        for key, item in self._point.items():
            #print(f'point load {key}')
            for bload in item:
                res = beam.reactions(load=bload)
                r1 = list(map(add, r1, [*res[0][0], *res[0][1]]))
                r2 = list(map(add, r2, [*res[1][0], *res[1][1]]))
        return [r1[:4], r1[4:]], [r2[:4], r2[4:]]
    #
    def nodal_equivalent(self):
        """
        Convert beam load [udl and point] to beam's end nodal load local system
        """
        res = self.reactions()
        return self._load2reaction(F0=res[0], F1=res[1])
        #
        #return array('d', nload)
        #return nload
    #
    def _load2reaction(self, F0:list, F1:list):
        """ Change load signs to reacctions"""
        #
        R0y, R0z = F0
        R1y, R1z = F1
        #
        lnload = [0] * 12
        #lnload = array('d', lnload)
        # ---------------------------
        # End 1
        # lnload[ 0 ] # axial
        lnload[1] =  1 * R0y[0]  # y
        lnload[2] =  1 * R0z[0]  # z
        # lnload[3]  # torsion Mx
        lnload[4] = 1 * R0z[1]  # my
        lnload[5] = -1 * R0y[1]  # mz
        # ---------------------------
        # End 2
        # lnload[ 0 ] # axial
        lnload[7] = -1 * R1y[0]  # y
        lnload[8] = -1 * R1z[0]  # z
        # lnload[9]  # torsion Mx
        lnload[10] = -1 * R1z[1]  # my
        lnload[11] = 1 * R1y[1]  # mz
        #
        return lnload
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
    def fer(self, beams):
        """ """
        """Beam reacition global system according to boundaries
        """
        b2n = []
        #beam = self._beam
        global_system = 0
        # line loadreactions
        for key, item in self._line.items():
            beam = beams[key]
            node1, node2 = beam.nodes
            #print(f'line load {key}')
            for bload in item:
                res = bload.fer_beam(L=beam.L)
                # local to global system
                gnload = [*res[4], *res[5]]
                lnload = trnsload(gnload, beam.T3D())
                b2n.append([bload.load_name, bload.title, global_system, 
                            beam.number, node1.number, lnload[:6], node2.number, lnload[6:]])
        # point load
        for key, item in self._point.items():
            beam = beams[key]
            node1, node2 = beam.nodes
            #print(f'point load {key}')
            for bload in item:
                res = bload.fer_beam(L=beam.L)
                gnload = [*res[4], *res[5]]
                lnload = trnsload(gnload, beam.T3D())
                #b2n.append([bload, item.name, item.title])
                #b2n.append(res)
                b2n.append([bload.load_name, bload.title, global_system, 
                            beam.number, node1.number, lnload[:6], node2.number, lnload[6:]])
        #
        #1 / 0
        #return [r1[:4], r1[4:]], [r2[:4], r2[4:]]
        return b2n
    #
    # ------------------
    #
    #def df(self, data):
    #    """ """
    #    print('====')
    #    1/0
    #
    #
#
# ---------------------------------
#
class BeamDistMaster(Mapping):
    
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