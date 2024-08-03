#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from array import array
#from dataclasses import dataclass
from typing import NamedTuple
import math

# package imports
import matplotlib.pyplot as plt
import numpy as np
#from numpy.matlib import repmat
#
from steelpy.utils.dataframe.main import DBframework

#
#
#@dataclass
class SurfaceResults(NamedTuple):
    #__slots__ = ['Hw', 'Tw', 'd', 'finite_depth', 'surface']
    #
    #def __init__(self, Hw:float, Tw:float, d:float, 
    #             data:None|dict, finite_depth:bool):
    """
    surface: free surface data frame [x, eta, phase]
    Hw : Wave height [unit length]
    Tw : Wave Period [second]
    d  : Water depth [unit length]
    finite_depth : bool
    """
    surface: list
    Hw : float
    Tw :float
    d : float
    finite_depth: bool
    #
    #df = DBframework()
    #self._data = df.DataFrame(data)
    #
    def __getitem__(self, name:str):
        """
        """
        return self.surface[name].to_numpy()
    #
    #def __getattribute__(self, attr):
    #    """
    #    Getter for myattr
    #    :param attr:
    #    :return:
    #    """
    #    if attr in self.__slots__:
    #        return self[attr]
    #    elif re.search(r"\bd\b", attr, re.IGNORECASE):
    #        return self.get_item(item="diameter")
    #    elif re.search(r"\bt(w)?\b", attr, re.IGNORECASE):
    #        return self.get_item(item="wall_thickess")
    #    else:
    #        raise AttributeError(f"Variable {attr} not found")    
    #
    #
    def __str__(self):
        """ """
        output = "\n"
        output += ("#\n")
        output += ("# Surface of wave - trough-crest-trough, \n")
        if self.finite_depth:
            header1 = "X"
            header2 = "eta"
        else:
            header1 = "kX"
            header2 = "k eta"
        #
        for idx, row in self.surface.iterrows():
            output += ("# {:} = {: 1.3e} Phase = {: 1.3e} {:} = {: 1.3e}"
                       .format(header1, row.x, row.phase, header2, row.eta))
            output += "\n"
        return output
    #
    #
    @property
    def eta(self):
        """eta : Surface elevation [unit length]"""
        return self.surface['eta'].to_numpy()
    #
    @property
    def time(self):
        """eta : Surface elevation [unit length]"""
        return self.surface['time'].to_numpy()
    #
    @property
    def phase(self):
        """phase : Wave phase [degree] """
        return self.surface['phase'].to_numpy()
    #
    @property
    def x(self):
        """Surface of wave coordinates : trough-crest-trough [unit length]"""
        return self.surface['length'].to_numpy()
    #
    #
    #
    def plot(self, phase:bool = False):
        """ """
        if phase:
            x = [-1*item for item in reversed(self.surface['phase'].iloc[::-1])]
            x.extend(self.surface['phase'].iloc[1::])
            x_label = r'$\theta$ (deg)'
            y_label = r'$\eta$ ($\theta$)'
        else:
            x = [-1*item for item in reversed(self.surface['length'].iloc[::-1])]
            x.extend(self.surface['x'].iloc[1::])
            x_label = r'$\lambda$ (m)'
            y_label = r'$\eta$ (m)'
        #
        y = [item for item in reversed(self.surface['eta'].iloc[::-1])]
        y.extend(self.surface['eta'].iloc[1::])
        #
        plt.plot(x, y)
        plt.title('Surface')
        plt.xlabel(x_label)
        #if self.finite_depth:
        #    plt.ylabel('$\zeta$ (m)')
        #else:
        plt.ylabel(y_label)
        plt.show()
    #
#
#
#
def get_surface(n: int, kd: float, Y: list, d:int, Tw: float, Lw: float, 
                nprofiles: int, is_finite: bool) -> DBframework:
    """
    n : order - Number of Fourier components or order of Stokes or cnoidal theory
    kd: wave number
    Y : Fourier components or Stokes/Cnoidal theory
    d : water depth
    Tw : Wave period
    Lw : wave length
    nprofiles : Number of points on free surface
    is_finite: 
    """
    pi = np.pi
    npt = number_steps(nprofiles)
    # position of the horizontal axis (theta : 0-pi)
    #xx = np.arange(npt) * pi / nprofiles
    #phase2 = xx * 180 / pi
    #
    #t = Tw / 2.0 # sec
    #omega = 2 * pi / Tw
    x = np.linspace(0, Lw, npt, endpoint=True)
    Theta =  kd / d * x - np.pi # omega * t
    phase = np.round(Theta * 180 / pi,  decimals=1)
    #
    #eta = np.array([surface(item, Y, n) for item in x])
    Keta = surfacenp(Theta, Y, n)
    #
    if is_finite:
        Keta /= kd
        #xx /= kd
        Theta /= kd
    #
    #xx = d * xx / kd
    #
    data = {'theta': Theta * d,
            'eta': Keta * d,
            'phase': phase,
            #'length': x,
            'time': np.linspace(0, Tw, npt, endpoint=True)}
    #
    df = DBframework()
    return df.DataFrame(data)
#
#def repmat2(A, n, axis:int):
#    """
#    """
#    A1 = np.expand_dims(A, axis)
#    #A1 = np.transpose(A1)
#    A1 = np.tile(A1, n)
#    #A1 = np.transpose(A1)
#    #return np.moveaxis(A1, 0, 1)
#    return A1
#
#  Surface elevation
def surface(x, Y, n):
    """
    Free surface wave profile function
    x : position of the horizontal axis (theta : 0-pi)
    Y : Fourier components or Stokes/Cnoidal theory
    n : order - Number of Fourier components or order of Stokes or cnoidal theory

    Return:
    Surface elevation
    """
    #kEta = 0
    kEta = np.sum([Y[j] * np.cos(j * x) 
                   for j in range(1, n)])
    kEta += 0.5 * Y[n] * np.cos(n * x)
    return kEta
#
def surfacenp(x: list, Y: list, n: int):
    """
    Free surface elevation is calculated from the interpolating
    Fourier series.

    x : position of the horizontal axis (theta : -pi..0..+pi)
    Y : Fourier components or Stokes/Cnoidal theory
    n : order - Number of Fourier components or order of Stokes or cnoidal theory

    Return:
    Keta : Surface elevation
    """
    cs = np.array([np.cos((j + 1) * x) for j in range(n)]) 
    Keta = np.sum(Y[1:n] * cs[:n-1, ].T, axis=1)
    Keta += 0.50 * Y[n] * cs[n-1]
    return Keta
#
def number_steps(StpLgth: int) -> int:
    """
    """
    npt = max(StpLgth, 2)
    if npt % 2 == 0:
        npt += 1
    return int(npt)
#
#
def get_etas(n, z, Y, B, Tanh, nprofiles, is_finite):
    """
    Surface - print out coordinates of points on surface for plotting 
    plus check of pressure on surface.
    """
    pi = math.pi
    surface_points = nprofiles
    kd = z[1]
    c=z[4]/math.sqrt(z[1])
    ce=z[5]/math.sqrt(z[1])
    R=1+z[9]/z[1]
    # Surface - print out coordinates of points on surface for plotting 
    # plus check of pressure on surface
    #print("# %s\n", Title);
    #print("%s\n", Method);
    print("# Surface of wave - trough-crest-trough,")
    print("# note quadratic point spacing clustered around crest")
    if is_finite:
        print("# Non-dimensionalised with respect to depth")
        print("#    X/d   eta/d   check of surface pressure")
        #print("# Dummy point to scale plot")
    else:
        print("# Non-dimensionalised with respect to wavenumber")
        print("#    kX    k eta   check of surface pressure")    
    #
    s_range = surface_points // 2
    for i in range(-s_range, s_range+1):
        #NB Quadratic point spacing, clustered near crest
        X = 4 * pi * (i/surface_points)**2
        X = math.copysign(X, i)
        eta = Surface(X, Y, n)
        (y, Pressure, Bernoulli_check, 
         u, v, dphidt, ut, vt, ux, uy) = Point(X, eta, kd, Tanh, 
                                               B, n, ce, c, R, z, is_finite)
        if is_finite:
            print("{:8.4f} {:7.4f} {:7.0e}".format(X/kd, 1+eta/kd, Pressure))
        else:
            print("{:8.4f} {:7.4f} {:7.0e}".format(X, eta, Pressure))
    print("")
    #
    npt = number_steps(nprofiles)
    xx =  array('f', [0 for i in range(npt)])
    eta = array('f', [0 for i in range(npt)])    
    return xx, eta
#