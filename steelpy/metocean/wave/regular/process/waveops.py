#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
#from array import array
from dataclasses import dataclass
import re
import math
#from typing import NamedTuple

# package imports
from steelpy.metocean.current.main import MeanCurrent
from steelpy.metocean.wave.regular.process.kinematic import get_kinematic, KinematicResults
from steelpy.metocean.wave.regular.process.inout import get_Height
#
from steelpy.metocean.wave.regular.process.surface import get_surface, SurfaceResults
#
#from steelpy.process.dataframe.main import DBframework
import numpy as np


#
#
#
@dataclass
class WaveItem:
    __slots__ = ['H', 'Tw', 'd', 'title', '_Lw', '_surface',
                 'order', 'nstep', 'niter', 'accuracy',
                 'c_vel', 'c_type', 'method', 'finite_depth',
                 '_wave_length', '_z', '_Y', '_B', '_Tanh', '_Highest']

    def __init__(self, H: float, d: float, title: str,
                 Tw: float | None = None,
                 Lw: float | None = None,
                 infinite_depth: bool = False,
                 current: float = 0.0, c_type: int = 1,
                 order: int = 5, nstep: int = 2,
                 niter: int = 40, accuracy: float = 1e-6):
        """ """
        self.title = title
        self.H = H
        self.d = d
        if Lw:
            self._Lw = Lw / self.d
            self.Tw = None
        else:
            self.Tw = Tw
            self._Lw = None
        # current 
        self.c_vel = current
        self.c_type = c_type
        #
        self.order = order
        self.nstep = nstep
        self.niter = niter
        self.accuracy = accuracy
        #
        self.finite_depth = True
        if infinite_depth:
            self.finite_depth = False
        #
        #self.wave_length:float
        self._surface = None

    #
    #
    def current(self, c_type: int = 1, c_velocity: float = 0.0) -> None:
        """Current data  
            c_type  : Current criterion, 1 - Eularian mean, 2 - Stokes depth integrated mean
            c_velocity : Current magnitude"""
        self.c_vel = c_velocity
        self.c_type = c_type

    #
    #
    def get_parameters(self, g: float = 9.80665):
        """ 
        g: gravity # m/s^2
        
        returns:
        MaxH : H/d
        T : Dimensionless period / None
        L : Dimensionless wavelength / None
        c_type  : Current criterion, 1 for Euler, or 2 for Stokes
        current : Current magnitude
        order :  Number of Fourier components or order of Stokes or cnoidal theory
        nsteps : Number of height steps to reach H/d
        niter  : Maximum number of iterations for each step (20)
        crit   : Criterion for convergence (1e-6)
        Height : wave height
        finite_depth : True/False
        """
        MaxH = self.H / self.d
        current = self.c_vel / math.sqrt(g * self.d)
        if self._Lw:
            case = "wavelength"
            L = self._Lw / self.d
            Height = get_Height(MaxH, case, self.finite_depth, L=L)
            T = None
        else:
            case = 'period'
            T = self.Tw * math.sqrt(g / self.d)
            Height = get_Height(MaxH, case, self.finite_depth, T=T)
            L = None
        #
        return [MaxH, case, T, L, self.c_type, current,
                self.order, self.nstep, self.niter, self.accuracy,
                Height, self.finite_depth]

    #
    def get_surface(self, surface_points: int = 36):
        """ claculate wave surface """
        n = self.order
        kd = self._z[1]
        Y = self._Y
        surface = get_surface(n, kd, Y, self.d,
                              surface_points, self.finite_depth)
        #
        if self.title:
            surface['title'] = self.title
        else:
            surface['title'] = None
        #surface['wave_id'] = self.number
        surface['type'] = 'order_1'
        #
        # TODO : check time
        num = len(surface['length'])
        endtime = self.Tw * 0.50
        time = np.linspace(start=0, stop=endtime,
                           num=num, endpoint=True)
        surface['time'] = time
        #time = np.linspace(start=0, stop=self.Tw,
        #                   num=(2*num)-1, endpoint=True)
        #surface['time'] = time[num-1:]
        #
        self._surface = surface
        #
        return surface[['type', 'length', 'eta', 'phase', 'time']]

    #
    def surface(self, surface_points: int = 36):
        """ wave surface """
        surface = self._surface
        if not surface:
            surface = self.get_surface(surface_points)
        #
        return SurfaceResults(surface=surface,
                              Hw=self.H, Tw=self.Tw, d=self.d,
                              finite_depth=self.finite_depth)

    #
    def get_kinematics(self, depth_points: int = 100,
                       surface_points: int | None = None):
        """get wave kinematics"""
        surface = self._surface
        if surface_points:  #or not surface
            surface = self.get_surface(surface_points)

        kpoints = depth_points
        kindf = get_kinematic(self.order, self._z, self._B, self._Tanh,
                              self.d, surface, kpoints,
                              self.finite_depth)
        #
        #surface['wave_name'] = self.name
        #kindf['wave_id'] = self.number
        kindf['type'] = 'order_1'
        #
        #depth_steps = np.arange(depth_points + 1) / depth_points
        #kindf['surface_id'] = repmat(surface['number'].to_numpy(),
        #                                 depth_steps.size, 1).flatten('F')
        #
        #kindf.drop(columns=['phase', 'x'], inplace=True, axis=1)        
        #
        return kindf

    #
    def kinematics(self, depth_points: int = 100,
                   surface_points: int | None = None):
        """wave kinematics"""
        kin = self.get_kinematics(depth_points, surface_points)
        return KinematicResults(surface=self._surface,
                                kindata=kin,
                                depth_points=depth_points)

    #
    @property
    def Lw(self):
        """ wave_length [m]"""
        if not self._Lw:
            self.solve()
        return self._Lw * self.d

    @Lw.setter
    def Lw(self, L: float):
        """ wave_length [m]"""

        self._Lw = L / self.d

    #
    @property
    def Hw(self):
        """Wave height"""
        return self.H
    #
    @property
    def crest(self):
        """Wave crest elevation"""
        surface = self.surface()
        return surface.eta.max()


#
#
#
#
def surfacePoints(d: float, points: int, eta: list,
                  stickup: float = 1.0):
    """get surface points"""
    crestmax = np.ceil(eta.max()) + stickup
    crestmin = np.floor(eta.min())
    step1 = int(np.ceil(points / 2))
    step2 = int(points - step1)

    return np.hstack([np.linspace(-d, crestmin, step1, endpoint=False),
                      np.linspace(crestmin, 0, step2, endpoint=False),
                      np.linspace(0, crestmax, step2)])


#
#
#
#
#
def get_wave_data(case_data: list | tuple | dict) -> list:
    """ """
    if isinstance(case_data, (list, tuple)):
        data = get_data_list(data=case_data)
    elif isinstance(case_data, dict):
        data = get_data_dic(data=case_data)
    else:
        raise IOError('input data not valid')
    return data


#
def get_data_dic(data: dict) -> list:
    """
    [Hw, Tw, d, wave_theory, title, crest_elevation, Lw]
    """
    new_data = [0, 0, 0, "Fourier", None, 0, 0]
    for key, item in data.items():
        # wave basic
        if re.match(r"\b((wave(\_|\s)?)?h(eight)?(w)?)\b", key, re.IGNORECASE):
            try:
                new_data[0] = item.value
            except AttributeError:
                raise IOError(f'units missing for Hw: wave height')

        elif re.match(r"\b((wave(\_|\s)?)?period|t(w)?|s)\b", key, re.IGNORECASE):
            try:
                new_data[1] = item.value
            except AttributeError:
                raise IOError(f'units missing for Tw: water period')

        elif re.match(r"\b((water(\_|\s)?)?d(epth)?)\b", key, re.IGNORECASE):
            try:
                new_data[2] = item.value
            except AttributeError:
                raise IOError(f'units missing for d: water depth')

        # wave type and title
        elif re.match(r"\b((wave(\_|\s)?)?t(ype|heory))\b", key, re.IGNORECASE):
            new_data[3] = item

        elif re.match(r"\b(title)\b", key, re.IGNORECASE):
            new_data[4] = item

        # user data
        elif re.match(r"\b((wave(\_|\s)?)?l(ength)?)\b", key, re.IGNORECASE):
            try:
                new_data[5] = item.value
            except AttributeError:
                raise IOError(f'units missing for Lw: wave length')

        elif re.match(r"\b((wave(\_|\s)?)?c(rest)?(\_|\s)?e(levation)?)\b", key, re.IGNORECASE):
            try:
                new_data[6] = item.value
            except AttributeError:
                raise IOError(f'units missing for WCe: wave crest elevation')

                #elif re.match(r"\b((wave(\_)?)?order)\b", key, re.IGNORECASE):
        #    new_data[7] = item
    return new_data


#
#
def get_data_list(data: list) -> list:
    """
    [Hw, Tw, d, wave_theory, title, crest_elevation, Lw]
    """
    winput = ['Hw', 'Tw', 'd']
    new_data = [0, 0, 0, "Fourier", None, 0, 0]
    # Wave basic data
    for x, item in enumerate(data[:3]):
        try:
            new_data[x] = item.value
        except AttributeError:
            raise IOError(f'units missing for {winput[x]}')
    #
    uinput = ['WCe', 'Lw']
    start = len(winput)
    stop = start + len(uinput)
    #
    # Wave tytpe, title
    for x, item in enumerate(data[start:stop]):
        try:
            new_data[x + start] = item.value
        except AttributeError:
            new_data[x + start] = item
    #
    # Wave basic data
    for x, item in enumerate(data[stop:]):
        try:
            new_data[x + stop] = item.value
        except AttributeError:
            raise IOError(f'units missing for {uinput[x]}')
    return new_data
#
#
