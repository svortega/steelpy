#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
#import math
from typing import NamedTuple
from dataclasses import dataclass
#
# package imports
#
import numpy as np
#import pandas as pd
import xarray as xr
#from numpy.matlib import repmat
import matplotlib.pyplot as plt
#
from steelpy.utils.dataframe.main import DBframework
#from steelpy.metocean.process.beamhydro import BeamMorisonWave
from steelpy.utils.math.operations import linstep
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
                 'condition', 'rho', 'up',
                 'morison', 'Zelev',
                 '_eta', '_x', '_time', 'd']
    def __init__(self, 
                 wave: tuple, 
                 current:tuple,
                 properties:tuple,
                 rho_w:float,
                 condition:int=1):
        """
        """
        self.wave = wave
        #self.current = current
        self.properties = properties
        self.condition = condition
        self.rho = rho_w
        self.up: str = 'y'
        #
        # -----------------------------------------
        # Wave Kinematics
        #
        kinematics = self.wave.kinematics()
        d = kinematics.d
        eta = kinematics.surface.eta
        depth_points = kinematics.depth_points
        #
        # -----------------------------------------
        # Grid elevations
        self.Zelev = self.ZElev(d=d, Hlat=0, eta=eta,
                                points=depth_points)
        #
        # -----------------------------------------
        #
        wkf = self.wave.kinematic_factor
        krf = wkf.get_profile(self.Zelev)        
        kin = kinematics.get_kin5(krf=krf, Zelev=self.Zelev)
        #
        self.morison = BeamMorison(kin)
        #
        # -----------------------------------------
        # Current
        #
        self.current = current.current
        cbf = current.blockage_factor
        #crf = cbf.get_profile(self.Zelev)
        self.current.seastate(grid=self.Zelev, eta=eta, cbf=cbf)
        #
        self._eta = eta
        self._x = kinematics.surface.x
        self._time = kinematics.surface.time
        self.d = d
    #
    # -----------------------------------------------    
    #
    def ZElev(self, d: float, Hlat:float,
              eta: list,points: int,
              Htop: float|None = None):
        """Elevation function"""
        crest = eta.max()
        trough = eta.min()
        if not Htop:
            Htop = 2*crest
        #
        steps = int(np.ceil(points / 4))
        d4 = - d // 4
        return np.hstack([np.linspace(-d, 3*d4, steps, endpoint=False),
                          np.linspace(3*d4, trough, steps, endpoint=False),
                          np.linspace(trough, Hlat, steps, endpoint=False),
                          np.linspace(Hlat, crest, steps, endpoint=False),
                          np.linspace(crest, Htop, steps)])     
    #
    #
    #
    def coordinates(self, beam, steps):
        """get coordinates along beam"""
        n1, n2 = beam.nodes
        #steps = self.steps()
        coord = [n1[:3]]
        for step in steps[1:]:
            coord.append(beam.find_coordinate(node_distance=step))
        coord = list(map(list, zip(*coord)))
        if self.up in ['y']:
            coord = [coord[0], coord[2], coord[1]]        
        return UnitVector(*coord)
    #
    def steps(self, beam, nelev: int):
        """beam steps"""
        sect = beam.section.geometry
        steps = linstep(d=sect.Dh, L=beam.L,
                        steps=nelev)
        return steps    
    #
    #
    # -----------------------------------------------
    #    
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
    def Dh(self, mg, beam):
        """Diametre hydrodynamic"""
        # TODO: section naming, element direction
        section = beam.section.geometry
        Dh = section.Dh
        #mg = self.MG(Z)
        Dh += 2 * mg
        At = np.pi * np.power(Dh, 2) / 4
        return Dh, At    
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
        #
        # TODO: Maybe separate beam hydro module
        #Bwave = BeamMorisonWave(beam=beam, rho=self.rho_w,
        #                        nelev=nelev, up=self.up)
        #
        #
        # -----------------------------------------
        # TODO: wtheta
        wtheta = self.wave.direction
        #kinematics = self.wave.kinematics()
        #
        #d = kinematics.d
        #eta = kinematics.surface.eta
        #depth_points = kinematics.depth_points
        #
        #Zelev = self.ZElev(d=d, Hlat=0, eta=eta,
        #                   points=depth_points)
        #
        # -----------------------------------------
        #
        #section = beam.section
        #D = section.diameter
        #
        # TODO : fix dircos
        dircos = beam.dircosines
        if self.up in ['y']:
            dircos = [dircos[0], dircos[2], dircos[1]]
        uvector = UnitVector(*dircos)        
        #
        # -----------------------------------------
        #
        #dz = np.diff(Elev)
        # locating the middle point of each element
        #Z = Elev[:-1] + dz
        #Z = Bwave.Z(nelev=nelev)
        #Elev = Bwave.elevations()
        #coord = Bwave.coordinates()
        #
        # -----------------------------------------
        # Hydro diametre & area
        marine_growth = self.properties.marine_growth
        #mg = marine_growth.MG(Z)
        mg = marine_growth.get_profile(self.Zelev)
        #
        #Dh, At = self.Dh(D, Z)
        #Dh, At = beamhydro.Dh(mg=mg)
        #
        # -----------------------------------------
        # Cd & Cm
        cdcm = self.properties.CdCm
        #cd, cm = cdcm.getCdCm(Z, crestmax, condition=self.condition)
        Cd, Cm = cdcm.get_profile(self.Zelev)
        #
        # -----------------------------------------
        # Current
        # TODO : ctheta
        #eta = np.hstack((list(reversed(eta[1:])), eta))
        #ctheta = self.current.direction
        #current = self.current.current
        #cbf = self.current.blockage_factor
        #cbf = cbf.get_profile(self.Zelev)
        #eta2 = np.hstack((list(reversed(eta[1:])), eta))
        #Vcp = current.get_profile(eta, Zelev, cbf)
        #Vc = current.Vc(cbf, uvector)
        Vc = self.current.Vc(uvector)
        #Vc *= 0
        #
        #print('----------------------------')
        #print(f'max Vnx = {Vc.Un.max()}')
        #print(f'max Vny = {Vc.Vn.max()}')
        #print(f'max Vnz = {Vc.Wn.max()}')        
        #
        # -----------------------------------------
        # Kinematis
        #
        #wkf = self.wave.kinematic_factor
        #wkf = wkf.get_profile(self.Zelev)
        #
        #kin = kinematics.get_kin5(krf=wkf, Zelev=self.Zelev)
        #kin = kinematics.get_kin4(elev=Elev, krf=wkf)
        #
        #
        #print('----------------------------')
        #print(f'max vel_x = {kin["ux"].max()}')
        #print(f'max vel_z = {kin["uz"].max()}')
        #
        #print(f'max acc_x = {kin["ax"].max()}')
        #print(f'max acc_z = {kin["az"].max()}')
        #        
        #
        #time = kinematics.surface.time
        #theta = kinematics.theta
        #
        # -----------------------------------------
        #
        Dh, At = self.Dh(mg=mg, beam=beam)
        #
        #Ka = self.An(kin, uvector)
        #Kv = self.Un(kin, uvector)
        Ka = self.morison.An(uvector)
        #print('----------------------------')
        #print(f'max Anx = {Ka.Anx.max()}')
        #print(f'max Any = {Ka.Any.max()}')
        #print(f'max Anz = {Ka.Anz.max()}')        
        #
        Kv = self.morison.Un(uvector)
        #print('----------------------------')
        #print(f'max Unx = {Kv.Un.max()}')
        #print(f'max Uny = {Kv.Vn.max()}')
        #print(f'max Unz = {Kv.Wn.max()}')
        #print(f'max vn = {Kv.vn.max()}')
        #
        #
        # --------------------------------
        # Global System
        vnn = Kv.vn + Vc.vn
        
        Fdx = (0.5 * self.rho * Cd * Dh * (Kv.Un + Vc.Un) * vnn
               + self.rho * Cm * At * Ka.Anx)
        
        Fdy = (0.5 * self.rho * Cd * Dh * (Kv.Vn + Vc.Vn) * vnn
               + self.rho * Cm * At * Ka.Any)
        
        Fdz = (0.5 * self.rho * Cd * Dh * (Kv.Wn + Vc.Wn) * vnn
               + self.rho * Cm * At * Ka.Anz)
        #
        #print('----------------------------')
        #print(f'fx = {Fdx.max()}')
        #print(f'fy = {Fdy.max()}')
        #print(f'fz = {Fdz.max()}')
        #
        # --------------------------------
        # member Local System        
        #
        #Tm = beam.unit_vector
        #Fb = Tm @ np.array([Fdx, Fdy, Fdz])
        #print('----------------------------')
        #print(f'fx = {Fb[0].max()}')
        #print(f'fy = {Fb[1].max()}')
        #print(f'fz = {Fb[2].max()}')         
        #
        Fi = xr.Dataset(data_vars={'fx': Fdx, 'fy': Fdy, 'fz': Fdz})
        #
        #      
        #
        df =  self.solve(Fi, beam, nelev)
        #
        #
        qload = self.bload(beam, df, nelev)
        #
        #
        #print('--')
        #1 / 0
        #
        #
        return qload
    #
    #
    # -----------------------------------------------
    #
    def solve(self, Fi, beam, nelev: int):
        """ """
        bsteps = self.steps(beam, nelev)
        coord = self.coordinates(beam, bsteps)
        Tm = beam.unit_vector
        # Fix vertical axis
        dz = np.diff(bsteps) 
        dz = np.array([np.average(dz), *dz]) * 0.50
        elev = self.d + (np.array(coord.z) - dz)
        #
        idx = [x for x in range(len(coord.x))]
        #
        Lstep = np.diff(bsteps)
        Lstep = np.concatenate(([np.average(Lstep)], Lstep))
        #
        #start = self._eta.argmax() + 1
        dftemp = []
        for step, crest in enumerate(self._eta):
            #print(f'step:{step}, crest:{crest}')
            item = Fi.roll(x=step)
            Fstep = item.interp(x=coord.x,
                                y=coord.y,
                                z=coord.z,
                                method="linear",
                                assume_sorted=True,
                                kwargs={"fill_value": 0})
            #print(Fstep)
            fx = Fstep['fx'].data[idx, idx, idx] * Lstep
            fy = Fstep['fy'].data[idx, idx, idx] * Lstep
            fz = Fstep['fz'].data[idx, idx, idx] * Lstep
            #
            Fb = Tm @ np.array([fx, fy, fz])
            #Fb = np.array([fx, fz, fy])
            #
            # TODO : fix torsion
            ft = fz * 0
            #
            BS = np.sqrt(np.power(Fb[1], 2)
                         + np.power(Fb[2], 2))
            #
            OTM = BS * elev
            #
            #crest_step = [crest for _ in range(len(fy))]
            #phase_step = [phase[step] for _ in range(len(fy))]
            length_step = [self._x[step] for _ in range(len(fy))]
            #
            #temp = [length_step, bsteps, fx, fy, fz]
            temp = [length_step, bsteps,
                    Fb[0], Fb[1], Fb[2], ft,
                    BS, OTM]
            temp = tuple(zip(*temp))
            dftemp.extend(temp) # , ft, BS, OTM
            #
            #temp = [fx, fy, fz]
            #dftemp.append(temp)        
        #
        #xdf = xr.DataArray()
        #
        db = DBframework()
        df = db.DataFrame(dftemp,
                          columns=['Lw', 'Lb', 'fx', 'fy', 'fz',
                                   'ft', 'BS', 'OTM'])
        #
        # --------------------------------
        #
        #Fx = df.groupby('Lw')['fx'].sum()
        #Fy = df.groupby('Lw')['fy'].sum()
        #Fz = df.groupby('Lw')['fz'].sum()
        #
        #print('----------------------------')
        #print('Beam Local System')
        #print(f'qx max: {Fx.max()}')
        #print(f'qy max: {Fy.max()}')
        #print(f'qz max: {Fz.max()}')    
        #
        return df
    #
    def bload(self, beam, dftemp, nelev: int):
        """ """
        # distance along beam 
        Lbeam = self.steps(beam, nelev)
        bl = beam.L
        uvector = UnitVector(*beam.dircosines)
        #
        #n1, n2 = beam.nodes
        #magA = np.sqrt(np.dot(n1[:3], n1[:3]))
        #magB = np.sqrt(np.dot(n2[:3], n2[:3]))
        #angle = np.arccos(np.dot(n1[:3], n2[:3])/(magA*magB))
        #deg = np.rad2deg(angle)
        #
        #n1 = getattr(n1, self.up)
        direction = getattr(uvector, self.up)
        #
        #if n1 != elev[0]:
        if direction < 0 :
            Lbeam = Lbeam[::-1]
        #
        grpdf = dftemp.groupby('Lw')
        #
        df = DBframework()
        qload = []
        for idx, (key, step) in enumerate(grpdf):
            q = step[['fx', 'fy', 'fz', 'ft', 'BS', 'OTM']] #.T
            qtemp = [[self._eta[idx], self._time[idx], 
                      *q.iloc[x-1],
                      *q.iloc[x],
                      Lbeam[x-1], bl - Lbeam[x]]
                     for x in range(1, len(q))]
            #
            qload.append(df.DataFrame(data=qtemp,
                                      columns=['eta', 'time', 
                                               'qx0', 'qy0', 'qz0', 'qt0', 'BS0', 'OTM0',
                                               'qx1', 'qy1', 'qz1', 'qt1', 'BS1', 'OTM1',
                                               'L0', 'L1'],
                                      index=None))
        #
        #
        qload = df.concat(qload)
        qload['BS'] = qload['BS0'] + qload['BS1']
        qload['OTM'] = qload['OTM0'] + qload['OTM1']        
        qload['element_type'] = 'beam'
        qload['element_name'] = beam.name
        qload['element_id'] = beam.number
        qload['type'] = 'line'
        #
        qload.drop(columns=['BS0', 'BS1', 'OTM0', 'OTM1'],
                   inplace=True)        
        #
        return qload
#
# TODO : remove permutes
#def permute00(A, order):
#    """ """
#    return np.tile(np.expand_dims(A, axis=(0, 1)), order)
# 
#def permute33(A, order):
#    """ """
#    A1 = np.expand_dims(A, axis=0)
#    A1 = np.transpose(A1)
#    A1 = np.tile(A1, order)
#    A1 = np.transpose(A1)
#    return A1
#
#def permute5(A, order, axis:int=1):
#    """ """
#    A1 = np.transpose(A)
#    A1 = repmat(A1, order[1], axis)
#    A1 = np.expand_dims(A1, axis=0)
#    A1 = np.transpose(A1)
#    A1 = np.tile(A1, order[0])
#    A1 = np.transpose(A1)
#    return A1
#
#
class BeamMorison(NamedTuple):
    """ """
    kin: list
    #
    def Un(self, uvector):
        """Components of velocity Global system"""
        comp0 = uvector.x *  self.kin['ux'] + uvector.z * self.kin['uz']
        #
        Un = self.kin['ux'] - uvector.x * comp0    # x
        Wn = - uvector.y * comp0                   # y
        Vn = self.kin['uz'] - uvector.z * comp0    # z
        # Water velocity normal to the cylinder axis
        vn = self.vn(self.kin, uvector)
        #
        return KinVel(Un, Vn, Wn, vn)
    #
    def vn(self, kin, uvector):
        """ Absolute Water velocity normal to the cylinder axis
        
        Vc : current velocity
        """
        vn = np.sqrt(np.power(self.kin['ux'], 2)
                     + np.power(self.kin['uz'], 2)
                     - np.power(uvector.x * self.kin['ux']
                                + uvector.z * self.kin['uz'], 2))
        return vn    
    #
    def An(self, uvector):
        """
        Instantaneous undisturbed acceleration resolved normal
        to the member
        
        kin : 
        """
        #
        comp0 = (uvector.x * self.kin['ax'] + uvector.z * self.kin['az'])
        # components of acceleration normal to the member 
        # in the x,y and z gloabl directions
        Anx = self.kin['ax'] - uvector.x * comp0
        Any = - uvector.y * comp0
        Anz = self.kin['az'] - uvector.z * comp0
        return KinAcc(Anx, Anz, Any)
    
#
class KinVel(NamedTuple):
    """

    Un : Kinematic components of velocity x
    Vn : Kinematic components of velocity y
    Wn : Kinematic components of velocity z
    vn : Fluid velocity normal to the cylinder axis
    rho : Sea water density (1025)
    """
    Un: list
    Vn: list
    Wn: list
    vn : list
    #rho: float
    #
#
class KinAcc(NamedTuple):
    """
    Anx : Kinematic components of velocity x
    Any : Kinematic components of velocity y
    Anz : Kinematic components of velocity z
    rho : Sea water density (1025)
    """
    Anx: list
    Any: list
    Anz: list
    #rho: float
#
class UnitVector(NamedTuple):
    """Unit vector """
    x: float
    y: float
    z: float
#
#
