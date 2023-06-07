#
# Copyright (c) 2009-2023 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
#import math
#from typing import NamedTuple, Tuple, Union, List, Dict
from dataclasses import dataclass
#
# package imports
#
import numpy as np
from numpy.matlib import repmat
#
#



#
@dataclass
class BSOTM:
    """
    Calculate the base shear and overturning moment of single pile

    Parameters
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
    kinematics: list
    condition: int = 1
    rho: float = 1025 # 
    #
    #
    def CdCm(self, Z):
        """ """
        cm = np.zeros((Z.shape))
        cd = np.zeros((Z.shape))
        cm[Z > 2] = 1.6
        cm[Z <= 2] = 1.2
    
        # switch condition
        if self.condition == 1:
            cd = cd + 1.15
        else:
            #elif condition == 2:
            cd[Z > 2] = 0.65
            cd[Z <= 2] = 1.05
        #
        return cd, cm
    #
    #
    def MGX(self, Z):
        """ MArine growth"""
        mgr = np.zeros((Z.shape))
        mgr[:, :, :] = 0.90 # m
        mgr[Z<-25.0] = 0.50 # m
        mgr[Z>2] = 0.0      # m
        return mgr
    #
    #
    def MG(self, Z):
        """ MArine growth"""
        mgr = np.zeros((Z.shape))
        mgr[:] = 0.90 # m
        mgr[Z < -25.0] = 0.50 # m
        mgr[Z > 2] = 0.0      # m
        return mgr
    #    
    #
    def Dh(self, D: float, Z):
        """Diamtre hydrodynamic"""
        mg = self.MG(Z)
        Dh = D + 2 * mg * 0
        At = np.pi * np.power(Dh, 2) / 4
        return Dh, At
    #
    #
    def Vc2(self, Vct: float, d: float, Z):
        """Current Velocity"""
        vc =  np.zeros((Z.shape))
        vc[:, :, :] = Vct
        vcr = np.power(np.abs(Vct * (Z + d), 1.0/7.0))
        vidx = Z < 0
        vc[vidx] = vcr[vidx]
        return vc
    #
    def Vc(self, Vct: float, d: float, z: list, eta: list):
        """Current Velocity"""
        zd = (z + d) / d
        vc =  np.zeros((eta.size, zd.size))
        vc[:, :] = Vct
        #vcr1 = np.abs(Vct *  d)
        #
        ze = np.zeros((eta.size, zd.size))
        ze[:, :] = z
        zebool = np.transpose(ze[:, :].T < eta)
        #
        vcr = np.zeros((eta.size, zd.size))
        vcr[:, :] = np.power(np.abs(Vct *  zd), 1.0/7.0)
        #
        vc[zebool] = vcr[zebool]
        return vc    
    #
    #
    def inertia(self, At, cm, AX, dz):
        """ """
        # inertia load per unit length
        pinertia = self.rho * cm * At * AX  
        # figure(4)
        # hold all
        # plot(pinertia(:),Z(:),'.-','LineWidth',1)
    
        # figure(5)
        # hold all
        # plot(pdrag(:),Z(:),'.-','LineWidth',1)
    
        dinertia = pinertia * dz
        #binertia = dinertia.sum(dim='z')
        return dinertia
    #
    def mass(self, D: float, cd, UX, dz):
        """ """
        # drag load per unit length
        pdrag = 0.5 * self.rho * cd * D * UX * np.abs(UX)
        ddrag = pdrag * dz
        #bdrag = np.sum(ddrag, axis=2)
        #bdrag = ddrag.sum(dim='z')
        return ddrag
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
        otm.drop('row', axis=1, inplace=True)
        otm['z'] = Z.flatten()
        return otm
    #
    #
    def BS(self, dinertia, ddrag, Z):
        """ Base Shear"""
        bs = dinertia + ddrag
        #bs = self.base_shear(dinertia, ddrag)
        #
        bs = bs.to_dataframe(name='BS').reset_index()
        bs.drop('row', axis=1, inplace=True)
        #
        bs['z'] = Z.flatten()
        return bs
    #
    def solve(self, D: float, L: float):
        """
        D : Pile diametre
        L : Pile length
        """
        d = self.kinematics.d
        z =  self.kinematics.z
        ux =  self.kinematics.ux
        ax =  self.kinematics.ax
        #
        # [x, z, time] = value --> irregular
        # [x, z, lenght] = value --> regular
        #
        UX = ux[:, :-1, :] + ux.diff('z') * 0.50
        AX = ax[:, :-1, :] + ax.diff('z') * 0.50
        #
        # sea water density
        #rho = self.rho  
        #
        dz = np.diff(z)
        # locating the midle point of each element
        Z = z[:-1] + dz * 0.50
        #
        #
        #
        # -----------------------------------------
        #
        Dh, At = self.Dh(D, Z)
        #
        # -----------------------------------------
        #
        # get Cd & Cm
        cd, cm = self.CdCm(Z)
        #
        # -----------------------------------------
        #
        eta = self.kinematics.surface.eta
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
        dinertia = self.inertia(At, cm, AX, dz)
        #
        ddrag = self.mass(Dh, cd, UX, dz)
        #
        bs = self.BS(dinertia, ddrag, Z)
        #
        otm = self.OTM(dinertia, ddrag, Z, d)
        #
        return bs, otm
    #
    #
    def wave_force(self, mesh):
        """Calculation of wave forces on beam elements"""
        beams = mesh.elements().beams()
        #
        # process
        #
        beam = beams[12]
        uvector = np.array(beam.unit_vector)
        #Duv = uvector / 100
        #DS = beam.L / 100
        #
        n1, n2 = beam.nodes
        #x = n1.x - Duv[0] * 0.5
        #y = n1.y - Duv[1] * 0.5
        #z = n1.z - Duv[2] * 0.5
        #
        # -----------------------------------------
        #
        d = self.kinematics.d
        z =  self.kinematics.z
        ux =  self.kinematics.ux
        ax =  self.kinematics.ax
        #
        #
        # -----------------------------------------
        #
        section = beam.section
        D = section.diameter
        #
        zmax = np.maximum(n1.y, n2.y)
        zmin = np.minimum(n1.y, n2.y)
        Elev = np.linspace(zmin, zmax, 10)
        #
        #b_range = np.abs(n2.y - n1.y)
        #
        #bzidx = (z > zmin) & (z < zmax)
        #
        #
        # -----------------------------------------
        #
        dz = np.diff(Elev)
        # locating the midle point of each element
        Z = Elev[:-1] + dz * 0.50        
        #
        # -----------------------------------------
        # Hydro diametre & area
        Dh, At = self.Dh(D, Z)
        #
        # -----------------------------------------
        # Cd & Cm
        cd, cm = self.CdCm(Z)
        #
        #
        # -----------------------------------------
        # Kinematis
        #
        kin = self.kinematics.get_kin(Z)        
        #        
        #
        #
        # ---------------------------------------
        #
        At = permute2(At, (kin['ax'].shape[0], kin['ax'].shape[2]), 1)
        #
        cm = permute2(cm, (kin['ax'].shape[0], kin['ax'].shape[2]), 1)
        #
        cd = permute2(cd, (kin['ax'].shape[0], kin['ax'].shape[2]), 1)        
        #
        Dh = permute2(Dh, (kin['ax'].shape[0], kin['ax'].shape[2]), 1)        
        #
        dz = permute2(dz, (kin['ax'].shape[0], kin['ax'].shape[2]), 1)
        #
        Z = permute2(Z, (kin['ax'].shape[0], kin['ax'].shape[2]), 1)
        #
        # ---------------------------------------
        #
        dinertia = self.inertia(At, cm, kin['ax'], dz)
        #
        ddrag = self.mass(Dh, cd, kin['ux'], dz)
        #
        #
        bs = self.BS(dinertia, ddrag, Z)
        #
        #
        print('--')
    #
    #
#
def bsotm_reg(kinematics,
              D: float,
              condition: int=1,
              rho:float = 1025):
    """
    Calculate the base shear and overturning moment of single pile

    Parameters
    ----------
    ux : horizontal wave particle velocity
    ax : horizontal wave particle acceleration
    z : z-coordinate
    d : Water depth
    D : diamatre of cylinder
    condtion :
        1.Linear wave
        2.Non-linear wave

    Returns
    -------
    bs   : base shear
    ovtm : overturning moment
    """
    # [Wave surface x coord, time,  Water depth]
    ux =  kinematics.ux
    ax =  kinematics.ax
    z =  kinematics.z
    d = kinematics.d
    #
    # sea water density
    rho = rho  
    #dz = z[1:] - z[:-1]
    dz = np.diff(z)
    # locating the midle point of each element
    Z = z[:-1] + dz * 0.50
    #
    #UX = (ux[:, :, 1:] - ux[:, :, :-1]) / 2 + ux[:, :, :-1]
    #AX = (ax[:, :, 1:] - ax[:, :, :-1]) / 2 + ax[:, :, :-1]
    #
    # [x, z, time] = value --> irregular
    # [dummy, z, x] = value --> regular
    #
    UX = ux[:, :-1, :] + ux.diff('z') * 0.50
    AX = ax[:, :-1, :] + ax.diff('z') * 0.50
    #
    Z = permute2(Z, (UX.shape[0],UX.shape[2]), 1)
    #
    dz = permute2(dz, (UX.shape[0],UX.shape[2]), 1)
    #dz = permute(dz, [1, 1, 1])
    #dz = np.tile(dz, [UX.shape[0], UX.shape[1], 1])
    # dz=repmat(dz, [UX.shape[0],UX.shape[1],1])

    cm = np.zeros((dz.shape))
    cd = np.zeros((dz.shape))
    cm[Z > 2] = 1.6
    cm[Z <= 2] = 1.2

    # switch condition
    if condition == 1:
        cd = cd + 1.15
    elif condition == 2:
        cd[Z > 2] = 0.65
        cd[Z <= 2] = 1.05
    # end
    pinertia = rho * np.pi * D**2 / 4 * cm * AX  # inertia load per unit length
    pdrag = 0.5 * rho * cd * D * UX * np.abs(UX)  # drag load per unit length

    # figure(4)
    # hold all
    # plot(pinertia(:),Z(:),'.-','LineWidth',1)

    # figure(5)
    # hold all
    # plot(pdrag(:),Z(:),'.-','LineWidth',1)

    dinertia = pinertia * dz
    binertia = dinertia.sum(dim='z')
    #binertia = np.sum(dinertia, axis=2)
    # binertia=summ(dinertia)
    ddrag = pdrag * dz
    #bdrag = np.sum(ddrag, axis=2)
    bdrag = ddrag.sum(dim='z')
    bs = binertia + bdrag
    #
    oinertia = dinertia * (Z + d)
    oinertia = oinertia.sum(dim='z')
    #oinertia = np.sum((dinertia * (Z + d)), axis=2)
    #
    odrag = ddrag * (Z + d)
    odrag = odrag.sum(dim='z')
    #odrag = np.sum((ddrag * (Z + d)), axis=2)
    ovtm = oinertia + odrag
    #
    return (bs.to_dataframe(name='BS').reset_index(),
            ovtm.to_dataframe(name='OTM').reset_index())
#
#
def permute(A, order):
    """ """
    return np.tile(np.expand_dims ( A, axis=(0, 1) ), order)
#
def repmat2(A, n, axis:int=1):
    """
    """
    A1 = repmat(A, 1, 1)
    A1 = np.transpose(A1)
    A1 = np.expand_dims(A1, 0)
    A1 = np.transpose(A1)
    A1 = np.tile(A1, n)
    return A1
#
def permute2(A, order, axis:int=1):
    """ """
    A1 = repmat(A, order[1], axis)
    A1 = np.transpose(A1)
    A1 = np.expand_dims(A1, axis=0)
    A1 = np.transpose(A1)
    A1 = np.tile(A1, order[0])
    return np.transpose(A1)   
#
def permute1(A, order, axis:int=1):
    """ """
    return np.transpose(repmat(A, order, axis))