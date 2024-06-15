#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
#import math
#from typing import NamedTuple
from dataclasses import dataclass
#
# package imports
#
import numpy as np
#from numpy.matlib import repmat
#import matplotlib.pyplot as plt
#
#from steelpy.process.dataframe.main import DBframework
from steelpy.metocean.process.beamhydro import BeamMorisonWave
#


#
#
@dataclass
class BSOTM:
    """
    Calculate the base shear and overturning moment of single pile

    Parameter
    ----------
    ux : horizontal wave particle velocity
    ax : horizontal wave particle acceleration
    z : z-coordinate
    d : Water depth
    D : diamere of cylinder
    condtion :
        1.Linear wave
        2.Non-linear wave
    rho : Sea water density 
    
    Returns
    -------
    bs  : Base shear
    otm : Overturning moment
    """
    __slots__ = ['wave', 'current', 'properties',
                 'condition', 'rho_w', 'up']
    def __init__(self, #kinematics:list,
                 wave: tuple, 
                 current:tuple,
                 properties:tuple,
                 rho_w:float,
                 condition:int=1):
        """
        """
        self.wave = wave
        #self.kinematics = kinematics
        self.current = current
        self.properties = properties
        self.condition = condition
        self.rho_w = rho_w
        self.up: str = 'y'
    #
    #
    #
    #
    #def base_shear(self, dinertia, ddrag):
    #    """ """
    #    #bdrag = ddrag.sum(dim='z')
    #    #binertia = dinertia.sum(dim='z')
    #    bs = dinertia + ddrag
    #    return bs
    #
    def OTM(self, dinertia, ddrag, Z, d):
        """ """
        oinertia = dinertia * (Z + d)
        #oinertia = oinertia.sum(dim='z')
        #
        odrag = ddrag * (Z + d)
        #odrag = odrag.sum(dim='z')
        otm = oinertia + odrag
        #
        otm = otm.to_dataframe(name='OTM').reset_index()
        #otm.drop('row', axis=1, inplace=True)
        #otm['z'] = Z.flatten()
        return otm
    #
    def BS(self, dinertia, ddrag):
        """ Base Shear"""
        bs = dinertia + ddrag
        #bs = self.base_shear(dinertia, ddrag)
        #
        bs = bs.to_dataframe(name='BS').reset_index()
        #bs.drop('row', axis=1, inplace=True)
        #
        #bs['z'] = Z.flatten()
        return bs
    #
    def solveBSOTM(self, D: float, L: float):
        """
        D : Pile diametre
        L : Pile length
        """
        d = self.kinematics.d
        z = self.kinematics.z
        #ux =  self.kinematics.ux
        #ax =  self.kinematics.ax
        #
        eta = self.kinematics.surface.eta
        #
        crestmax = np.max(eta)        
        #
        # [x, z, time] = value --> irregular
        # [x, z, lenght] = value --> regular
        #
        #UX = ux[:, :-1, :] + ux.diff('z') * 0.50
        #AX = ax[:, :-1, :] + ax.diff('z') * 0.50
        #
        # sea water density
        #rho = self.rho  
        #
        dz = np.diff(z)
        # locating the midle point of each element
        Z = z[:-1] + dz * 0.50
        #
        # -----------------------------------------
        # Kinematis
        #
        kin = self.kinematics.get_kin(Z)
        #
        UX = kin['ux']
        AX = kin['ax']
        #
        # -----------------------------------------
        #
        Dh, At = self.Dh(D, Z)
        #
        # -----------------------------------------
        #
        # get Cd & Cm
        cd, cm = self.CdCm(Z, crestmax)
        #
        # -----------------------------------------
        #
        #eta = self.kinematics.surface.eta
        Vct=1.54
        #zd = (z + d) / d
        Vc = self.Vc(Vct, d, z, eta)
        #
        #        
        # -----------------------------------------
        #
        #
        Z = permute2(Z, (UX.shape[0],UX.shape[2]), 1)
        #
        dz = permute2(dz, (UX.shape[0],UX.shape[2]), 1)
        #
        At = permute2(At, (UX.shape[0],UX.shape[2]), 1)
        #
        cm = permute2(cm, (UX.shape[0],UX.shape[2]), 1)
        #
        cd = permute2(cd, (UX.shape[0],UX.shape[2]), 1)        
        #
        Dh = permute2(Dh, (UX.shape[0],UX.shape[2]), 1)
        #
        # -----------------------------------------
        #        
        dmass = self.mass(At, cm, AX)
        #
        ddrag = self.drag(Dh, cd, UX, np.abs(UX))
        #
        bs = self.BS(dmass, ddrag)
        #
        otm = self.OTM(dmass, ddrag, Z, d)
        #
        return bs, otm
    #
    # -----------------------------------------------
    #
    #
    def Fwave(self, beam):
        """Calculation of wave forces on beam elements

        beam: Beam elements
        wave_angle: Wave angle
        """
        #
        wip = self.properties.WIP
        nelev = wip.nelev(beam, self.up)
        #
        # TODO: Maybe separate beam hydro module
        Bwave = BeamMorisonWave(beam=beam, rho=self.rho_w,
                                nelev=nelev, up=self.up)
        #
        #
        # -----------------------------------------
        # TODO: wtheta
        wtheta = self.wave.direction
        kinematics = self.wave.kinematics()
        #
        eta = kinematics.surface.eta
        #
        #
        #crestmax = np.max(eta)
        #
        # -----------------------------------------
        #
        #section = beam.section
        #D = section.diameter
        #
        # -----------------------------------------
        #
        #dz = np.diff(Elev)
        # locating the middle point of each element
        #Z = Elev[:-1] + dz
        #Z = Bwave.Z(nelev=nelev)
        Elev = Bwave.elevations()
        #coord = Bwave.coordinates()
        #
        # -----------------------------------------
        # Hydro diametre & area
        marine_growth = self.properties.marine_growth
        #mg = marine_growth.MG(Z)
        mg = marine_growth.get_profile(Elev)
        #
        #Dh, At = self.Dh(D, Z)
        #Dh, At = beamhydro.Dh(mg=mg)
        #
        # -----------------------------------------
        # Cd & Cm
        cdcm = self.properties.CdCm
        #cd, cm = cdcm.getCdCm(Z, crestmax, condition=self.condition)
        cd, cm = cdcm.get_profile(Elev)
        #
        # -----------------------------------------
        # Current
        # TODO : ctheta
        #eta = np.hstack((list(reversed(eta[1:])), eta))
        ctheta = self.current.direction
        current = self.current.current
        cbf = self.current.blockage_factor
        cbf = cbf.get_profile(Elev)
        eta2 = np.hstack((list(reversed(eta[1:])), eta))
        Vc = current.get_profile(eta2, Elev, cbf)
        #Vc *= 0
        #
        # -----------------------------------------
        # Kinematis
        #
        wkf = self.wave.kinematic_factor
        wkf = wkf.get_profile(Elev)
        #
        #kin2 = kinematics.get_kin2(beam, nelev)
        #Bwave.Fwave2(Vc=Vc, MG=mg, Cd=cd, Cm=cm,
        #             kinematics=kin2, elev=Elev)
        #
        kin = kinematics.get_kin4(elev=Elev, krf=wkf)
        #
        time = kinematics.surface.time
        return Bwave.Fwave(Vc=Vc, MG=mg, Cd=cd, Cm=cm,
                           kinematics=kin, eta=eta,
                           time=time)
    #
    #
#
#
def permute(A, order):
    """ """
    return np.tile(np.expand_dims(A, axis=(0, 1)), order)
#
#def repmat2(A, n, axis:int=1):
#    """
#    """
#    A1 = repmat(A, 1, 1)
#    A1 = np.transpose(A1)
#    A1 = np.expand_dims(A1, 0)
#    A1 = np.transpose(A1)
#    A1 = np.tile(A1, n)
#    return A1
#
#def permute2(A, order, axis:int=1):
#    """ """
#    A1 = repmat(A, order[1], axis)
#    A1 = np.transpose(A1)
#    A1 = np.expand_dims(A1, axis=0)
#    A1 = np.transpose(A1)
#    A1 = np.tile(A1, order[0])
#    return np.transpose(A1)   
#
#def permute1(A, order, axis:int=1):
#    """ """
#    return np.transpose(repmat(A, order, axis))
#
#