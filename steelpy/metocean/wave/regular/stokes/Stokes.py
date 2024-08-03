#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
#from array import array
from dataclasses import dataclass
#import math
#from typing import NamedTuple

# package imports
from steelpy.metocean.wave.regular.stokes.subroutines import F, CDE, AB
from steelpy.metocean.wave.regular.process.inout import title_block, output
from steelpy.metocean.wave.regular.process.waveops import WaveItem
import numpy as np
#from scipy.optimize import fsolve

#
#
@dataclass
class WaveStokes(WaveItem):

    def __init__(self, Hw: float, d: float,
                 Tw: float | None = None,
                 Lw: float | None = None,
                 title: str | None = None,
                 infinite_depth: bool = False,
                 current: float = 0.0, c_type: int = 1,
                 order: int = 5, nstep: int = 2,
                 niter: int = 40, accuracy: float = 1e-6) -> None:
        """
        """
        super().__init__(Hw=Hw, Tw=Tw, Lw=Lw, d=d, title=title,
                         order=order, nstep=nstep, niter=niter,
                         accuracy=accuracy,
                         current=current, c_type=c_type,
                         infinite_depth=infinite_depth)

    #
    #
    def solve(self):
        """ Solver """
        try:
            self._z
        except AttributeError:
            self.order = 5
            self.method = f"Stokes method order {self.order}"
            #
            data = self.get_parameters()
            z, Y, B, Tanh = StokesMain(*data)
            self._z = z
            self._Y = Y
            self._B = B
            self._Tanh = Tanh
            self._Lw = 2 * np.pi / z[1]
            #self._Highest = Highest        
        #
#
# 
# Main program
#
def StokesMain(MaxH: float, case: str,
               T: float | None,
               L: float | None,
               c_type: int, current: float,
               norder: int, nstep: int,
               niter: int, accuracy: float,
               Height: float, is_finite: bool):
    """
    Stokes theory calculations
        Input
    -------
    MaxH : H/d
    case : period/wavelength
    T : Dimensionless period / None
    L : Dimensionless wavelength / None
    c_type  : Current criterion, 1 - Eularian mean, 2 - Stokes depth integrated mean
    current : Current magnitude
    order :  Number of Fourier components or order of Stokes or cnoidal theory
    nsteps : Number of height steps to reach H/d
    niter  : Maximum number of iterations for each step (20)
    crit   : Criterion for convergence (1e-6)
    Height : ??
    finite_depth : True/False
    
    Output
    --------
    z : Solution vector
    Y : Discrete Fourier transform of the surface elevations.
    B : Fourier coefficients
    Tanh : 
    """
    # inital values
    #g = 9.80665  # m/s^2
    e = np.zeros(norder + 1)
    H = MaxH
    pi = np.pi
    #
    norder = min(norder, 5)
    print(f"# Solution set by {norder}-order Stokes theory")
    #
    #if Lw:
    if case == 'wavelength':
        if L > 10:
            print(f'The dimensionless wavelength [{L:1.2f}] is greater than 10')
            raise IOError("Stokes theory should not be applied")
        #
        kd = 2.0 * pi / L
        #kH = kd * H
        ckd, skd, ss, t, C, D, E = CDE(kd)
    else:  # Period
        if T > 10:
            print(f'The dimensionless period [{T:1.2f}] is greater than 10')
            raise IOError("Stokes theory should not be applied")
        #
        #def SK5func(*args: list):
        #    """ Stokes 5th Equations"""
        #    kd1, kd2, F22, C2, D2 = args
        #    ckd, skd, ss, t, C, D, E = CDE(kd1)
        #    F1 = F(kd1, H, T, current, c_type, C2, norder, D2)
        #    Fd = (F22 - F1) / (kd2 - kd1)
        #    delta = F1 / Fd
        #    #kd22 = kd1
        #    kd11 = kd1 - delta
        #    #F22 = F11
        #    return kd11, kd1, F22, C, D
        #
        # if period is specified, solve dispersion relation using secant method
        # Until February 2015 the bisection method was used for this.
        # I found that in an extreme case (large current) the bracketting
        # of the solution was not correct, and the program failed,
        # without printing out a proper error message.
        print("# Period has been specified. Now solving for L/d iteratively,")
        print("# Printing to check convergence:")
        omega = 2 * pi / T
        # Fenton & McKee for initial estimate
        kFM = (omega**2
               * np.power(1.0 / np.tanh(np.power(omega, 1.5)), 2/3))
        kd1 = kFM
        kd2 = kFM * 1.01
        ckd, skd, ss, t, C, D, E = CDE(kd2)
        F2 = F(kd2, H, T, current, c_type, C, norder, D)
        #
        #iguess = [kd1, kd2, F2, C, D]
        # System Solution
        #kd12, kd22, F22, C2, D2 = fsolve(SK5func, iguess)
        #
        for _iter in range(niter + 1):
            ckd, skd, ss, t, C, D, E = CDE(kd1)
            F1 = F(kd1, H, T, current, c_type, C, norder, D)
            Fd = (F2 - F1) / (kd2 - kd1)
            delta = F1 / Fd
            kd2 = kd1
            kd1 -= delta
            print("{: 8.4f}".format(2 * pi / kd1))
            if abs(delta / kd1) < accuracy:
                break
            F2 = F1
            if _iter > niter:
                raise RuntimeError("Secant for solution of wavenumber has not converged")
                # print("Contact John Fenton johndfenton@gmail.com")
        kd = kd1
        #kH = kd * H
    #
    kH = kd * H
    #
    z = np.zeros(2 * norder + 10 + 1)
    z[1] = kd
    z[2] = kH
    SU = 0.5 * kH / np.power(kd, 3)
    print("# Stokes-Ursell number (SU) :{: 7.4f}".format(SU))
    if SU > 0.5:
        raise Warning("SU > 1/2. Results are unreliable")
    else:
        print("# SU < 1/2, Stokes theory should be valid")
    #
    e[1] = 0.5 * kH
    for i in range(2, norder + 1):
        e[i] = e[i - 1] * e[1]
    #
    # Calculate coefficients
    Y, z, A, B = AB(skd, ss, t, norder, z, e, C, kd, ckd)
    z[7] = C[0] + e[2] * C[2] + e[4] * C[4]  # ubar
    z[8] = - e[2] * D[2] - e[4] * D[4]
    z[9] = 0.5 * C[0] * C[0] + e[2] * E[2] + e[4] * E[4]
    if c_type == 1:
        z[5] = current * np.sqrt(kd)
        z[4] = z[7] + z[5]
        z[6] = z[4] + z[8] / kd - z[7]
    else:
        z[6] = current * np.sqrt(kd)
        z[4] = z[6] - z[8] / kd + z[7]
        z[5] = z[4] - z[7]
    #
    if case == 'wavelength':
        z[3] = 2 * pi / z[4]
    else:  # Period
        z[3] = T * np.sqrt(kd)
    #
    #  Highest wave - eqn (32) of Fenton (1990)
    Tanh = np.array([np.tanh(i * z[1])
                     for i in range(norder + 1)])
    #
    title_block(is_finite, c_type, current, z)
    output(norder, z, Y, B, Tanh, is_finite)
    return z, Y, B, Tanh
#
