#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
from typing import NamedTuple


# package imports
#
from steelpy.utils.units.buckingham import Number
from steelpy.utils.math.operations import trnsload
from steelpy.ufo.load.process.beam.utils import(get_BeamLoad_list_units,
                                                get_BeamLine_dict,
                                                get_BeamLoad_dict,
                                                get_BeamLine_load)
from steelpy.trave.beam.operation import BeamTorsion, BeamAxial, BeamBending
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
    comment: str | None
    load_name: str | int
    mesh_name: str | int
    system: int
    #
    # -------------------------
    #
    def local_system(self):
        """
        return q load in local beam system

        Returns:
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
        # set beam operations
        L =beam.L
        material = beam.material
        section = beam.section.properties(poisson=material.poisson)
        beam_ops = BeamBasic(L=L,
                             area=section.area,
                             Iy=section.Iy, Iz=section.Iz,
                             J=section.J, E=material.E,
                             G=material.G, Cw=section.Cw*0,
                             Asy=section.Asy,
                             Asz=section.Asz,
                             Pdelta=False)
        #
        boundary1 = 'fixed'
        #if nodes[0].boundary:
        #    boundary1 = nodes[0].fixity
        #
        boundary2 = 'fixed'
        #if nodes[1].boundary:
        #    boundary2 = nodes[1].fixity
            
        beam_ops.supports(boundary1, boundary2)
        #bops.supports(end1='fixed', end2='fixed')
        #
        # Get load function along beam length
        Fbar = self.Fx(x=L, L=L,
                       E=material.E, G=material.G, 
                       Iy=section.Iy, Iz=section.Iz,
                       J=section.J, Cw=section.Cw*0,
                       Area=section.area,
                       Asy=section.Asy, Asz=section.Asz,
                       P=0.0)
        #
        # end_1,2 = [axial, torsion, in_plane, out_plane]
        # Axial   [P, blank, blank, u]
        # Torsion [T, Phi, Psi, B, Tw]
        # in_plane [V, M, theta, w]
        # out_plane [V, M, theta, w]
        #
        #end_1, end_2 = beam_ops._reaction(Fbar=Fbar)
        end_1, end_2 = beam_ops._FER_shear_force(Fbar=Fbar)
        #
        # Mapp results to forces
        #
        # Force = [fx, fy, fz, mx, my, mz]
        # Node 1
        Fend1 = [end_1.x[0], 1 * end_1.y[0], 1 * end_1.z[0],  # Fx, Fy, Fz
                 end_1.t[0], 1 * end_1.z[1], 1 * end_1.y[1]]  # Mx, My, Mz
        # Node 2
        Fend2 = [end_2.x[0], 1 * end_2.y[0], 1 * end_2.z[0],  # Fx, Fy, Fz
                 end_2.t[0], 1 * end_2.z[1], 1 * end_2.y[1]]  # Mx, My, Mz
        #
        # Displacement = [x, y, z, rx, ry, rz]
        # Node 1
        Rend1 = [end_1.x[3], 1 * end_1.y[3], 1 * end_1.z[3],  # x, y, z,
                 end_1.t[1], 1 * end_1.z[2], 1 * end_1.y[2]]  # rx, ry, rz
        # Node 2
        Rend2 = [end_2.x[3], 1 * end_2.y[3], 1 * end_2.z[3],  # x, y, z,
                 end_2.t[1], 1 * end_2.z[2], 1 * end_2.y[2]]  # rx, ry, rz
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
        Pend1 = [*Fend1, *Rend1, *end_1.t[2:]]
        Pend2 = [*Fend2, *Rend2, *end_2.t[2:]]
        #
        FER_out = [self.load_name, self.mesh_name,
                   self.comment, self.load_type, 
                   'basic', self.system,
                   self.name, Pend1, Pend2]
        #1 / 0
        #print('--> FER')
        return FER_out
    #  
#
#
class LoadFunction(NamedTuple):
    load_name: str | int
    mesh_name: str | int
    load_comment: str
    load_type: str
    load_level: str
    load_system: str | int
    element_name: str | int
    length: float
    axial: list | np.array
    torsion: list | np.array
    VM_inplane: list | np.array
    VM_outplane: list | np.array
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
    load_step: float | None
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
        
        if (load_comment:= self.comment) == "NULL":
            load_comment = ""        
        step = self.load_type
        output += (f"{step[:12]:12s} {self.L1: 1.3e} {self.qx1: 1.3e} "
                   f"{self.qy1: 1.3e} {self.qz1: 1.3e} "
                   f"{str(load_comment)}\n")
        return output
    #
    # -------------------------
    # Load Process Formulas
    # ------------------------
    #
    def Fx(self, x:float|list, L:float,
           E:float, G:float,
           Iy:float, Iz:float,
           J:float, Cw:float,
           Area:float, Asy:float, Asz:float, 
           P:float, factor:float = 1.0) -> list[NamedTuple]:
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
        load_comment, load_type,
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
        #
        Mx = BeamTorsion(E=E, L=L, G=G, J=J, Cw=Cw)     
        #
        # In plane [V, M, theta, w]
        in_plane =  BeamBending(L=L, E=E, G=G, As=Asy, I=Iy)
        # Out plane [V, M, theta, w]
        out_plane = BeamBending(L=L, E=E, G=G, As=Asz, I=Iz)
        #
        # [Fx, Fy, Fz, Mx, My, Mz]
        if isinstance(x, (list, tuple)):
            Fx_out = []
            for xstep in x:
                #print(f'step --> {xstep}')
                Fx_out.append(LoadFunction(self.load_name, self.mesh_name,
                                           self.comment, self.load_type, 'basic',
                                           self.coordinate_system,
                                           self.name, xstep,
                                           
                                           np.array(Fx.axial(x=xstep,
                                                             P=self.qx0 * factor,
                                                             L1=self.L0,
                                                             P2=self.qx1 * factor,
                                                             L2=self.L1)),
                                           # TODO : torsion
                                           np.array(Mx.torque(x=xstep,
                                                              T=self.qt0 * factor,
                                                              L1=self.L0,
                                                              T2=self.qt1 * factor,
                                                              L2=self.L1)),
                                           
                                           np.array(in_plane.line(x=xstep,
                                                                  axial=P, 
                                                                  q1=self.qy0 * factor,
                                                                  q2=self.qy1 * factor,
                                                                  L1=self.L0, L2=self.L1)),
                                           
                                           np.array(out_plane.line(x=xstep,
                                                                   axial=P,
                                                                   q1=self.qz0 * factor,
                                                                   q2=self.qz1 * factor,
                                                                   L1=self.L0, L2=self.L1))))
        else:
            Fx_out = LoadFunction(self.load_name, self.mesh_name,
                                  self.comment, self.load_type, 'basic',
                                  self.coordinate_system,
                                  self.name, x,
                                  
                                  np.array(Fx.axial(x=x,
                                                    P=self.qx0 * factor,
                                                    L1=self.L0,
                                                    P2=self.qx1 * factor,
                                                    L2=self.L1)),
                                  # TODO : torsion
                                  np.array(Mx.torque(x=x,
                                                     T=self.qt0 * factor,
                                                     L1=self.L0,
                                                     T2=self.qt1 * factor,
                                                     L2=self.L1)),
                                  
                                  np.array(in_plane.line(x=x,
                                                         axial=P, 
                                                         q1=self.qy0 * factor,
                                                         q2=self.qy1 * factor,
                                                         L1=self.L0, L2=self.L1)),
                                  
                                  np.array(out_plane.line(x=x,
                                                          axial=P,
                                                          q1=self.qz0 * factor,
                                                          q2=self.qz1 * factor,
                                                          L1=self.L0, L2=self.L1)))
        #
        # [load_name, load_comment, load_type, load_system, 
        # beam_number, x, Fx, Fy, Fz]
        return Fx_out
    #
    #
    def _Fa(self, L: float):
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
            return 0
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
        return Fa[0]
    #
    #
#
#
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
    load_step: float | None
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
        
        if (load_comment:= self.comment) == "NULL":
            load_comment = ""        
        step = self.load_type
        output += (f"{step[:12]:12s} {10*' '} {self.mx: 1.3e} "
                   f"{self.my: 1.3e} {self.mz: 1.3e} "
                   f"{str(load_comment)}\n")
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
           E:float, G:float, 
           Iy:float, Iz:float,
           J: float, Cw:float,
           Area:float, Asy:float, Asz:float,
           P:float, factor:float = 1.0) -> list[NamedTuple]:
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
              load_comment, load_type,
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
        Finplane =  BeamBending(L=L, E=E, G=G, As=Asy, I=Iy)
                               # alpha_s=alpha_s)
        # Out plane [V, M, theta, w]
        F_outplane = BeamBending(L=L, E=E, G=G, As=Asz, I=Iz)
                                # alpha_s=alpha_s)
        #
        # [Fx, Fy, Fz, Mx, My, Mz]
        if isinstance(x, (list,tuple)):
            Fx_out = []
            for xstep in x:
                Axial =  Fx.axial(x=xstep,
                                  P=self.fx * factor,
                                  L1=self.L0)
                #Axial
                Torsion =  Mx.torque(x=xstep,
                                     T=self.mx * factor,
                                     L1=self.L0)
                #Torsion
                Finp = Finplane.point(x=xstep,
                                      axial=P,
                                      P=self.fy * factor,
                                      M=self.mz * factor,
                                      L1=self.L0)
                #
                Foutp = F_outplane.point(x=xstep,
                                         axial=P,
                                         P=self.fz * factor,
                                         M=self.my * factor,
                                         L1=self.L0)
                #
                Fx_out.append(LoadFunction(self.load_name, self.mesh_name,
                                           self.comment, self.load_type, 'basic',
                                           self.coordinate_system,
                                           self.name, xstep,
                                           np.array(Axial), np.array(Torsion),
                                           np.array(Finp), np.array(Foutp)))
        else:
            Axial =  Fx.axial(x=x,
                              P=self.fx * factor,
                              L1=self.L0)
            # torsion = [FT, FB, Fpsi, Fphi]
            Torsion =  Mx.torque(x=x,
                                 T=self.mx * factor,
                                 L1=self.L0)
            
            Finp = Finplane.point(x=x,
                                  axial=P, 
                                  P=self.fy * factor,
                                  M=self.mz * factor,
                                  L1=self.L0)
            
            Foutp = F_outplane.point(x=x,
                                     axial=P,
                                     P=self.fz * factor,
                                     M=self.my * factor,
                                     L1=self.L0)
            
            Fx_out = LoadFunction(self.load_name, self.mesh_name,
                                  self.comment, self.load_type, 'basic',
                                  self.coordinate_system, 
                                  self.name, x,
                                  np.array(Axial), np.array(Torsion),
                                  np.array(Finp), np.array(Foutp))
        # [axial, torsion, bending_inplane, bending_outplane]
        return Fx_out
#
#
#
# ---------------------------------
#
def distributed_axial(w:float, L:float,
                      l1:float, l2:float) -> list[float]:
    """
    Case 6 from Matrix Analysis of Framed Structures [Aslam Kassimali]
    """
    Fa = w/(2*L) * (L-l1-l2) * (L-l1+l2)
    Fb = w/(2*L) * (L-l1-l2) * (L+l1-l2)
    return [Fa, Fb]
#
# ---------------------------------
#