# 
# Copyright (c) 2019 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
#from collections.abc import Mapping
from typing import NamedTuple



# package imports
#
# Bending 
#
#from steelpy.trave.beam.pilkey.chapter11.table2C import BeamBendingSupports
from steelpy.trave.beam.pilkey.chapter11.table3C import BeamBendingSupports
#
#from steelpy.trave.beam.pilkey.chapter11.table2 import BendingGE
from steelpy.trave.beam.pilkey.chapter11.table3 import BendingGE
#
#from steelpy.trave.beam.pilkey.chapter11.table2B import (Trapezoidal,
#                                                         Point, Moment)
from steelpy.trave.beam.pilkey.chapter11.table3B import (Trapezoidal,
                                                         Point, Moment)
#
# Torsion
from steelpy.trave.beam.pilkey.chapter14.table_C import BTOpenSupports
from steelpy.trave.beam.pilkey.chapter14.table_A import TorsionOpenGE
from steelpy.trave.beam.pilkey.chapter14.table_B import ( TOpenConcentrated,
                                                         TOpenDistributed)
from steelpy.trave.beam.pilkey.chapter12.table_2t import (TBarConcentrated,
                                                          TBarDistributed,
                                                          BTBarSupports,
                                                          TorsionBarGE)
#
# Axial
from steelpy.trave.beam.pilkey.chapter12.table_2axial import (PBarConcentrated,
                                                              PBarDistributed,
                                                              PBarSupports,
                                                              AxialBarGE)
#
from steelpy.utils.math.operations import linspace #, linstep
from steelpy.utils.dataframe.main import DBframework
#
#from steelpy.trave.beam.roark.chapter10.table103 import (TOpenConcentrated,
#                                                         TorsionRoarkTW)
#
#from steelpy.trave.beam.roark.chapter10.table103_point import BTWSupPoint
#
#
#
#
#
@dataclass
class BeamBasic:
    __slots__ = ['L','area', 'Iy', 'Iz', 'J', 'E', 'G', 'Cw',
                 'Asy', 'Asz', '_Pdelta', 
                 '_response', '_support0', '_support1',
                 '_torsion', '_axial', '_bminplane', '_bmoutplane']

    def __init__(self, L: float, area: float,
                 Iy:float, Iz:float, J:float,
                 E:float, G:float, Cw:float,
                 Asy: float, Asz: float,
                 Pdelta: bool) -> None:
        """
        L : beam length
        E : Elastic modulus
        G : Shear modulus
        J : Torsial constant [m^4]
        Iy,z : moment of inertia [m^4]
        Cw : Warping constant [m^6]
        alpha_sy,z: Shear correction coefficient
        Pdelta : True/False
        """
        self.L = L
        self.area = area
        self.Iy = Iy
        self.Iz = Iz
        self.J = J
        self.E = E
        self.G = G
        self.Cw = Cw
        self.Asy = Asy
        self.Asz = Asz
        self._Pdelta = Pdelta
        #
        self._torsion = BeamTorsion(L=L, E=E, G=G,
                                    J=self.J, Cw=self.Cw)
        
        self._axial = BeamAxial(L=L, E=E, A=area)
        #
            
        self._bminplane = BeamBending(L=L, E=E, G=G,
                                      As=self.Asy, I=self.Iy)
        #
        self._bmoutplane = BeamBending(L=L, E=E, G=G,
                                       As=self.Asz, I=self.Iz)
    #
    # ----------------------------------------------
    #
    def supports(self, end1: str|list, end2: str|list,
                 torsion1: str|list|None  = None, torsion2: str|list|None  = None,
                 k1: float | None = None, k2: float | None = None):
        """
        end1 : beam end 1 [fix, free, pinned, guide, spring]
        end2 : beam end 2 [fix, free, pinned, guide, spring]
        """
        # FIXME supports one plane
        self._support0 = get_support(end1)
        self._support1 = get_support(end2)
        #
        if not torsion1:
            torsion1 = end1
            
        if not torsion2:
            torsion2 = end2
        
        self._torsion.supports(end1=torsion1, end2=torsion2)
        #
        self._axial.supports(end1=self._support0[0],
                             end2=self._support1[0])
        #
        self._bminplane.supports(end1=self._support0[0],
                                 end2=self._support1[0])
        #
        self._bmoutplane.supports(end1=self._support0[1],
                                  end2=self._support1[1])        
    #
    # ----------------------------------------------
    # Beam Formulas
    #
    def _reaction(self, Fbar):
        """beam end reactions in local system
        
        Return:
        R0 [Fa, T, Ry, Rz]
        R1 [Fa, T, Ry, Rz]
        """
        Pflag = 0
        if self._Pdelta:
            Pflag = 1        
        # get initial condition beam local system
        # [Axial, Torsion, BM_inplane, BM_outplane]
        R0 = self.R0(Fbar=Fbar)        
        #
        #Fx = [load_name, load_title, load_type, load_system,
        #      beam_name, L_step,
        #      Axial, Torsion, BM_inplane, BM_outplane]
        #
        #Axial, BM_inplane, BM_outplane = [V, M, w, theta]
        #Torsion = [T, B, psi, phi, Tw]
        #
        #(lname, ltitle, ltype, lsystem, bnumber, x,
        # Fx_bar, Mx_bar, Fy_bar, Fz_bar)
        #
        # Fbar =[load_name, component, 
        #        load_title, load_type,
        #        load_level, load_system,
        #        beam_name, x, Fx, Mx, Fy, Fz]
        #
        #
        # Axial [P, blank, blank, u]
        Fa = self._axial.response(x=self.L, R0=R0[0], Fx=Fbar[8])
        Fa = [-1*item for item in Fa]
        #
        # Torsion [T, Phi, Psi, B,Tw]
        T1x = self._torsion.response(x=self.L, R0=R0[1], Fx=Fbar[9])
                                     #Fx=[-1 * item for item in Fbar[9]])
        T1x = [-1*item for item in T1x]
        #
        # Bending in plane [V, M, theta, w, P]
        #Fx = [-1 * item for item in Fbar[10]]
        #Fx.append(Fa[0])
        #Rf = [*R0[2], R0.x[0]]
        R1y = self._bminplane.response(x=self.L,
                                       R0=[*R0[2], R0.x[0] * Pflag],
                                       Fx=[-1 * item for item in Fbar[10]])
        #R1y = [-1*item for item in R1y]
        #R1y[0] *= -1
        #
        # Bending out plane [V, M, theta, w, P]
        #Fx = [-1 * item for item in Fbar[11]]
        #Fx.append(Fa[0])
        R1z = self._bmoutplane.response(x=self.L,
                                        R0=[*R0[3], R0.x[0] * Pflag],
                                        Fx=[-1 * item for item in Fbar[11]])
        #R1z = [-1*item for item in R1z]
        #R1z[0] *= -1
        #
        R1 = R0eq(Fa, T1x, R1y, R1z)
        #
        #
        #K = (1 / 3.0 * (2 * 6.20**3 * 120 + 9.80**3 * 0.240)) / 1000**4
        #rtorsion = TorsionROARK(Cw=self.Cw, K=self.J, E=self.E, G=self.G)
        #rtorsion.moment(Tm=Fbar[9][0], Xe=self.L*0.50, L=self.L)
        #thetas = rtorsion.pinned_ends(coord=[0, 0.25 * self.L, 0.50 * self.L, 0.75 * self.L, self.L])
        #thetas = rtorsion.fixed_ends(coord=[0, 0.25 * self.L, 0.50 * self.L, 0.75 * self.L, self.L])
        #thetas = rtorsion.cantilever(coord=[0, 0.25 * self.L, 0.50 * self.L, 0.75 * self.L, self.L])
        #
        return R0, R1

    #
    def R0(self, Fbar):
        """Reaction end 0 local system
        
        Return: 
        R0 [Fa, T, Ry, Rz]
        """
        #
        Pflag = 0
        if self._Pdelta:
            Pflag = 1         
        #
        # Axial [P, blank, blank, u]
        suppx = self._axial.R0(Fbar=Fbar[8])
        #
        # Torsion [T, Phi, Psi, B, Tw]
        suppt = self._torsion.R0(Fbar=Fbar[9])
        #suppt = [-1*item for item in suppt]
        #
        # Bending in plane [V, M, theta, w, P]
        bload = [-1 * item for item in Fbar[10]]
        bload.append(suppx[0] * Pflag)
        suppy = self._bminplane.R0(Fbar=bload)
        #suppy[1] *= -1
        # Bending out plane [V, M, theta, w, P]
        bload = [-1 * item for item in Fbar[11]]
        bload.append(suppx[0] * Pflag)
        suppz =  self._bmoutplane.R0(Fbar=bload)
        #suppz[1] *= -1
        #
        return R0eq(suppx, suppt, suppy, suppz)
        #return self._sign_inverse(suppx, suppt, suppy, suppz)
    #
    #
    def response(self, x: float, R0:list, Fx:list):
        """
        x : Distance from end 1
        I : Moment of inertia [m^4]
        load_list : list of load to be included in calculation (default use all)

        Return:
        Fx = [Axial, Torsion, BM_inplane, BM_outplane]
        
        ----------------------------- 
        Axial, BM_inplane, BM_outplane = [V, M, w, theta]
        V : Shear force
        M : Bending moment
        theta : Slope
        w : Deflection
        
        -----------------------------
        Torsion =  [T, Phi, Psi, B, Tw]
        T : Twisting moment
        B : Bimoment
        psi : Rate of angle of twist
        phi: Angle of twist
        Tw : Warping torque
        """
        # calculate load at x
        Fax, Ftx, Fyx, Fzx = Fx
        R0x, R0t, R0y, R0z = R0
        #
        Pflag = 0
        if self._Pdelta:
            Pflag = 1         
        #
        # Axial [P, blank, blank, u]
        R1x = self._axial.response(x=x, R0=R0x, Fx=Fax)
        # Torsion [T, Phi, Psi, B, Tw]
        # Torsion [T, theta, theta', theta'', theta''']
        #try:
        #    1 / self.Cw
        #    tload =  Ftx
        #except ZeroDivisionError:
        #tload = [-1 * item for item in Ftx]
        R1t = self._torsion.response(x=x, R0=R0t, Fx=Ftx)
        # Bending in plane [V, M, theta, w, P]
        #Fx = [-1 * item for item in Fyx]
        #Fx.append(R0x[0])
        R1y = self._bminplane.response(x=x,
                                       R0=[*R0y, R0x[0] * Pflag],
                                       Fx=[-1 * item for item in Fyx])
        # Bending out plane [V, M, theta, w, P]
        #Fx = [-1 * item for item in Fzx]
        #Fx.append(R0x[0])
        R1z = self._bmoutplane.response(x=x,
                                        R0=[*R0z, R0x[0] * Pflag],
                                        Fx=[-1 * item for item in Fzx])
        #
        return R1x, R1t, R1y, R1z
    #
    # ----------------------------------------------
    # Beam Calculations
    #
    def reactions(self, bloads) -> dict:
        """Calculate reactions
        
        Return : 
        load_name : [load_title, load_type,
                     load_level, load_system,
                     beam_name,
                     R0[Fa, Tx, Ry, Rz],
                     R1[Fa, Tx, Ry, Rz]]
        """
        #
        # -----------------------------------------------------
        # loop basic load
        # 
        reac = defaultdict(list)
        for item in bloads:
            Fbar = item.Fx(x=self.L, L=self.L,
                           E=self.E, G=self.G, 
                           Iy=self.Iy, Iz=self.Iz,
                           J=self.J, Cw=self.Cw,
                           Area=self.area,
                           Asy=self.Asy, Asz=self.Asz,
                           P=0.0, factor=1.0) # FIXME: comb maybe?

            reac[item.load_name].append([item.load_name, item.component_name, 
                                         item.title, item.load_type, 'basic', 
                                         item.coordinate_system,
                                         item.name, 
                                         *self._reaction(Fbar=Fbar)])
        # TODO: need a dictionary as output?
        return reac
    #
    def solve(self, bloads, steps:int=10) -> list[list]:
        """
        
        Return:
        Beam Forces 
        [load_title, load_type, load_system,
         beam_name, L_steps,
         Axial, Torsion, BM_inplane, BM_outplane]
        
        ----------------------------- 
        Axial, BM_inplane, BM_outplane = [V, M, w, theta]
        V : Shear force
        M : Bending moment
        theta : Slope
        w : Deflection
        
        -----------------------------
        Torsion [T, Phi, Psi, B, Tw]
        T : Twisting moment
        B : Bimoment
        psi : Rate of angle of twist
        phi: Angle of twist
        Tw : Warping torque
        """
        #1 / 0
        header =  ['load_name', 'component', 
                   'load_title', 'load_type',
                   'load_level', 'load_system',
                   'element_name'] 
        #
        #        
        # -----------------------------------------------------
        #
        # {load_name : [load_name, component,
        #               load_title, load_type,
        #               load_level, load_system, 
        #               beam_name,
        #               R0[Fa, Tx, Ry, Rz],
        #               R1[Fa, Tx, Ry, Rz]] }     
        #
        reactions = self.reactions(bloads)
        #reac_df = self.support
        #grpsupp = reac_df.groupby(header)
        #
        #
        #Load function 
        # df [load_name, component, 
        #     load_title, load_type, 
        #     load_level, load_system,
        #     beam_name, 
        #     x, Fx, Mx, Fy, Fz]        
        #
        bfunc_df = self.load_function(bloads, steps)
        beamfun = bfunc_df.groupby(header)        
        #
        # -----------------------------------------------------
        #
        #Fblank = [0, 0, 0, 0]
        lbforce = []
        for items in reactions.values():
            for Rb in items:
                title = tuple(Rb[:7])
                mnload = beamfun.get_group(title)
                # invert load sign 
                mnload = mnload[['x','Fx','Mx','Fy','Fz']]
                #mnload.set_index('x', inplace=True)
                for bstep in mnload.itertuples():
                    lbforce.append([*title, bstep.x,
                                    *self.response(x=bstep.x,
                                                   R0=[*Rb[7]], # Rend 1
                                                   Fx=[bstep.Fx, bstep.Mx,
                                                       bstep.Fy, bstep.Fz])])
                #
        #
        return lbforce
    #
    def forces(self, bloads, steps:int=10):
        """
        Return: 
        Beam force  
        Dataframe [load_name, 'component_name',
                   load_title, load_type,
                   'load_level', load_system,
                   element_name, node_end,
                   Fx, Fy, Fz, Mx, My, Mz,
                   delta_x, delta_y, delta_z,
                   theta_x, theta_y, theta_z,
                   psi, B, Tw]
                   
                   theta_x', theta_x'', theta_x''']
        """
        # solve = [load_title, load_type, load_system,
        #          beam_name, L_steps,
        #          Axial, Torsion, BM_inplane, BM_outplane]
        #
        # Axial, BM_inplane, BM_outplane = [V, M, w, theta]
        # Torsion = [T, phi, psi, B, Tw]
        lbforce = self.solve(bloads, steps)
        lbforce =[[*lbf[:8], *lbf[8], *lbf[9], *lbf[10], *lbf[11]]
                  for lbf in lbforce]
        #
        # Member
        #
        header = ['load_name', 'component_name',
                  'load_title', 'load_type',
                  'load_level', 'load_system',
                  'element_name', 'node_end',
                  'Fx', 'blank1', 'blank2', 'delta_x',
                  'Mx', 'theta_x', 'psi', 'B', 'Tw', 
                  'Fy', 'Mz', 'theta_z', 'delta_y',
                  'Fz', 'My', 'theta_y', 'delta_z']
        #
        db = DBframework()
        df_mload = db.DataFrame(data=lbforce, columns=header, index=None)
        df_mload = df_mload[['load_name', 'component_name',
                             'load_title', 'load_type',
                             'load_level', 'load_system',
                             'element_name', 'node_end',
                             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                             'delta_x', 'delta_y', 'delta_z',
                             'theta_x', 'theta_y', 'theta_z',
                             'psi','B', 'Tw']]
        #
        df_mload.rename({'thetax': 'rx', 'thetay': 'ry', 'thetaz': 'rz'})
        return df_mload
    #
    def load_function(self, bloads, steps:int=10):
        """
        Return:
        Load function 
        Dataframe [load_name, load_title, load_type, load_system,
                   element_name, x, Fx, Mx, Fy, Fz]
        """
        Lsteps = linspace(start=0, stop=self.L,
                          num=steps+1, endpoint=True)
        #Lsteps = linstep(d=d, L=self.L, steps=steps)
        #
        # -----------------------------------------------------
        #
        #Fx = [load_name, component, load_title, 
        #      load_level, load_system,
        #      beam_name, x, Fx, Mx, Fy, Fz]
        #
        #Axial, BM_inplane, BM_outplane = [V, M, w, theta]
        #Torsion =  [T, Phi, Psi, B, Tw]        
        #
        beamfun = []
        for idx, item in enumerate(bloads):
            Fx = item.Fx(x=Lsteps, L=self.L,
                         E=self.E, G=self.G, 
                         Iy=self.Iy, Iz=self.Iz,
                         J=self.J, Cw=self.Cw,
                         Area=self.area,
                         Asy=self.Asy, Asz=self.Asz, 
                         P=0.0, factor=1.0) # FIXME: comb maybe?
            #beamfun[item.load_name].extend(lout)
            beamfun.extend(Fx)
        #
        header = ['load_name', 'component',
                  'load_title', 'load_type',
                  'load_level', 'load_system', 
                  'element_name',
                  'x', 'Fx', 'Mx', 'Fy', 'Fz']
        db = DBframework()
        df_func = db.DataFrame(data=beamfun, columns=header, index=None)
        return df_func
        #return beamfun
#
#
class R0eq(NamedTuple):
    """ """
    x: list
    t: list
    y: list
    z: list
#
# ---------------------------------------------
#
@dataclass
class BeamBending:
    __slots__ = [ 'E', 'L', 'I', 'G', 'As', 
                  '_support', '_gen']
    
    def __init__(self, L: float, E:float, G: float, 
                 As: float, I: float) -> None:
                 #alpha_s: float) -> None:    
        """
        L : Length of beam (L)
        E : Elastic modulus (F/L^2)
        G : Shear modulus of elasticity (F/L^2)
        A : area (L^2)
        I : Moment of inertia about neutral axis (L^4)
        """
        self.L = L
        self.E = E
        self.G = G
        self.As = As
        self.I = I
        #self.alpha_s = alpha_s
        #try:
        #    self.As = self.A / self.alpha_s
        #except ZeroDivisionError:
        #    self.As = 0
    #
    def supports(self, end1: str|list, end2: str|list,
                 k1: float|None = None, k2: float|None = None):
        """
        end1 : beam end 1 [fix, free, pinned, guide, spring]
        end2 : beam end 2 [fix, free, pinned, guide, spring]
        """
        self._support = [end1, end2, k1, k2]
    #
    def R0(self, Fbar:list):
        """Reaction end 1"""
        Mb = BeamBendingSupports(L=self.L, E=self.E, G=self.G,
                                 As=self.As, I=self.I)
        Mb.supports(*self._support)
        return Mb(F_bar=Fbar)
    #
    #
    def line(self, x: float|list,
             q1: float, q2: float,
             axial: float, 
             L1: float, L2: float):
        """line loading
        q1 : line load end 1
        q2 : line load end 2
        L1 : Distant to q1 from end 1
        L2 : Distant to q2 from end 2
        axial : Axial force (compression[-])
        """
        udl = Trapezoidal(q1=q1, q2=q2, P=axial, 
                          L=self.L, L1=L1, L2=L2)
        return udl(x=x, E=self.E, G=self.G, As=self.As, I=self.I)
    #
    def point(self, x: float|list,
              P: float, M:float,
              axial: float, 
              L1: float):
        """
        point loading
        P : Point load
        M : Moment
        axial : Axial force
        L1 : Distant to P,M from end 1
        """
        #
        Fp = Point(W=P, P=axial, L=self.L, L1=L1)
        Fm = Moment(M=M, P=axial, L=self.L, L1=L1)
        #point =  Fp(x=x, E=self.E, G=self.G, A=self.A, I=self.I)
        #moment = Fm(x=x, E=self.E, G=self.G, A=self.A, I=self.I)
        return list(map(sum, zip(Fp(x=x, E=self.E, G=self.G, As=self.As, I=self.I),
                                 Fm(x=x, E=self.E, G=self.G, As=self.As, I=self.I))))

    #
    #def loading_function(self, x: float):
    #    """ Arbitrary Loading functions"""       
    #    return self._gen(x=x, E=self.E, I=self.I)
    #
    def response(self, x: float, R0: list, Fx:list):
        """ General reponse """
        res= BendingGE(E = self.E, G=self.G, As=self.As, I=self.I)
                       #alpha_s=self.alpha_s)
        res.R0(*R0)
        res.load(*Fx)
        return res.response(x=x)
#
#
@dataclass
class BeamTorsion:
    __slots__ = [ 'E', 'L', 'G', 'J', 'Cw', 
                  '_support', '_bar', '_open']
                # 'T1',  'L1', 'T2', 'L2']    
    
    def __init__(self, L: float, E:float, G:float,
                 J: float, Cw: float) -> None:    
        """
        L : beam Lenght 
        E : Elastic modulus
        G : Shear modulus
        J : Torsional constant (length^4)
        Cw : Warping constant (length^6)
        """
        self.E = E
        self.G = G
        self.L = L
        self.J = J
        self.Cw = Cw
    #
    #
    def supports(self, end1: str|list, end2: str|list,
                 k1: float | None = None, k2: float | None = None):
        """
        end1 : beam end 1 [fix, free, pinned, guide, spring]
        end2 : beam end 2 [fix, free, pinned, guide, spring]
        """
        self._support = [end1, end2, k1, k2]
    #  
    #
    def R0(self, Fbar:list):
        """Reaction end 1"""
        #1 / 0
        try:
            # open section
            1 / self.Cw
            #if T2:
            #    # distributed
            #    print('fixme')
            #else:
            #    # concentrated            
            #    supp = BTWSupPoint(L=self.L,
            #                       J=self.J, Cw=self.Cw,
            #                       E=self.E, G=self.G)
            supp = BTOpenSupports(L=self.L,
                                  J=self.J, Cw=self.Cw,
                                  E=self.E, G=self.G)
        except ZeroDivisionError:
            # bar section
            supp = BTBarSupports(L=self.L, J=self.J,
                                 E=self.E, G=self.G)
        #
        supp.supports(*self._support)
        return supp(F_bar=Fbar)
    #
    #
    def torque(self, x: float|list,
               T: float, L1: float,
               T2: float|None = None, L2: float|None = None):
        """
        Torque
        """
        #
        if T2:
            # distributed
            barsec = TBarDistributed(T=T, T2=T2, L=self.L, L1=L1, L2=L2)
            opensec = TOpenDistributed(T=T, T2=T2, L=self.L, L1=L1, L2=L2)
        else:
            # concentrated
            barsec = TBarConcentrated(T=T, L=self.L, L1=L1)
            opensec = TOpenConcentrated(T=T, L=self.L, L1=L1)
        #
        try:
            # open section
            1 / self.Cw
            return opensec(x=x, E=self.E, G=self.G,
                           J=self.J, Cw=self.Cw)
        except ZeroDivisionError:
            # bar section
            return barsec(x=x, E=self.E, G=self.G, J=self.J)        
    #
    #
    #def loading_function(self, x: float):
    #    """ Arbitrary Loading functions"""
    #    try:
    #        # open section
    #        1 / self.Cw
    #        return self._open(x=x, E=self.E, G=self.G,
    #                          J=self.J, Cw=self.Cw)
    #    except ZeroDivisionError:
    #        # bar section
    #        return self._bar(x=x, E=self.E, G=self.G, J=self.J)
    #
    #
    def response(self, x: float, R0: list, Fx:list):
        """ General reponse
        
        [T, Phi, Psi, B, Tw]
        [T, theta, theta1, theta2, theta3] thin walled
        """
        try:
            # open section
            1 / self.Cw
            res = TorsionOpenGE(E=self.E, G=self.G,
                                J=self.J, Cw=self.Cw)
            #res = TorsionRoarkTW(E=self.E, G=self.G,
            #                     J=self.J, Cw=self.Cw)
        except ZeroDivisionError:
            # bar section        
            res = TorsionBarGE(E=self.E, G=self.G, J=self.J)
        #
        res.R0(*R0)
        res.load(*Fx)
        return res.response(x=x)
#
#
@dataclass
class BeamAxial:
    __slots__ = [ 'E', 'L', 'A','_support', '_bar']
    
    def __init__(self, L: float, E:float, A: float) -> None:    
        """
        E : Elastic modulus
        """
        self.E = E
        self.L = L
        self.A = A
    #
    def supports(self, end1: str|list, end2: str|list,
                 k1: float|None = None, k2: float|None = None):
        """
        end1 : beam end 1 [fix, free, pinned, guide, spring]
        end2 : beam end 2 [fix, free, pinned, guide, spring]
        """
        self._support = [end1, end2, k1, k2]
    #
    #
    def R0(self, Fbar:list):
        """Reaction end 1"""
        supp = PBarSupports(L=self.L, E=self.E, A=self.A)
        supp.supports(*self._support)
        return supp(F_bar=Fbar)
    #
    def axial(self, x: float|list,
              P: float, L1: float,
              P2: float|None = None, L2: float|None = None):
        """
        x: 
        """
        #self.P1 = P
        #self.L1 = L1
        if P2:
            bar = PBarDistributed(P=P, P2=P2, L=self.L, L1=L1, L2=L2)
        else:
            bar = PBarConcentrated(P=P, L=self.L, L1=L1)
        if isinstance(x, (list,tuple)):
            1 / 0
        else:
            return bar(x=x, E=self.E, A=self.A)
    #
    #def loading_function(self, x: float):
    #    """ Arbitrary Loading functions"""
    #    return self._bar(x=x, E=self.E, A=self.A)
    #
    def response(self, x: float, R0: list, Fx:list):
        """ General reponse
        
        Results:
        [P, blank, blank, u]
        """
        res =  AxialBarGE(E = self.E, A=self.A)
        res.R0(*R0)
        res.load(*Fx)
        return res.response(x=x)
#
# ---------------------------------------------
#
def list2str(axial: int, vertical: int, moment: int):
    """
    vertical [1:fixed, 0:free]
    moment
    """
    try:
        1 / int(vertical)
        if int(moment) == 1:
            return "fixed"
        else:
            return "pinned"
    except ZeroDivisionError:
        if int(moment) == 1:
            if int(axial) == 1:
                return "guided"
            else:
                return "free"
        else:
            return "free"
#
def get_support(support:list|str):
    """ """
    #
    if isinstance(support, (list, tuple)):
        #suppx = list2str(vertical=support[0], moment=support[3])
        suppy = list2str(axial=support[0], vertical=support[1], moment=support[5])
        suppz = list2str(axial=support[0], vertical=support[2], moment=support[4])
    elif isinstance(support, str):
        suppy = suppz = support
    else:
        raise IOError(f"fixity {support} not valid")
    #
    return suppy, suppz
#
#