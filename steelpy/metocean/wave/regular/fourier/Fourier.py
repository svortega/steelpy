#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from array import array
from dataclasses import dataclass

# package imports
from steelpy.metocean.wave.regular.fourier.Subroutines import Newton, initial
from steelpy.metocean.wave.regular.process.inout import title_block, output
from steelpy.metocean.wave.regular.process.waveops import WaveItem

import numpy as np


#
#
@dataclass
class WaveFourier(WaveItem):

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
        super().__init__(H=Hw, Tw=Tw, Lw=Lw, d=d, title=title,
                         order=order, nstep=nstep, niter=niter,
                         accuracy=accuracy,
                         current=current, c_type=c_type,
                         infinite_depth=infinite_depth)

    def solve(self):
        """ Solver """
        try:
            self._z
        except AttributeError:
            self.method = "Fourier"
            data = self.get_parameters()
            z, Y, B, Tanh = FourierMain(*data)
            self._z = z
            self._Y = Y
            self._B = B
            self._Tanh = Tanh
            self._Lw = 2 * np.pi / z[1]
        #self._Highest = Highest
        #return z, Y, B, Tanh


#
# 
# Main program
#
#
#def FourierMain(h:float, t:Union[float,None], d:float, 
#                Lw:Union[float,None], is_finite:bool,
#                current:float, c_type:int=1,
#                n:int=20, nstep:int=2, number:int=40, accuracy:float=1e-5):
def FourierMain(MaxH: float, case: str,
                T: float | None,
                L: float | None,
                c_type: int, current: float,
                norder: int, nstep: int,
                niter: int, accuracy: float,
                Height: float, is_finite: bool):
    """
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
    #current=0.31
    crit = accuracy
    #pi = np.pi
    #g = 9.80665  # m/s^2
    #
    print("# Solution by {:}-term Fourier series".format(norder))
    method = ("Fourier method with {:} terms in series".format(norder))
    num = 2 * norder + 10
    dhe = Height / nstep
    dho = MaxH / nstep
    #
    Y = np.zeros(num + 1)
    sol = np.zeros((num + 1, 2 + 1))
    B = np.zeros(norder + 1)
    # Commence stepping through steps in wave height
    for ns in range(1, nstep + 1):
        height = ns * dhe
        hoverd = ns * dho
        print("--- height step {:} of {:}".format(ns, nstep))
        # Calculate initial linear solution
        if ns <= 1:
            sol, z, cosa, sina = initial(height, hoverd, current,
                                         c_type, num, norder,
                                         is_finite, case)
        else:
            # Or, extrapolate for next wave height, if necessary
            z = [2.0 * sol[i][2] - sol[i][1] for i in range(num + 1)]
        #
        # Commence iterative solution
        for itr in range(1, niter + 1):
            print("# Iteration {:}:".format(itr))
            # Calculate right sides of equations and differentiate numerically
            # to obtain Jacobian matrix, : solve matrix equation
            z, error, Tanh = Newton(height, hoverd, z, cosa, sina,
                                    num, norder, current, c_type,
                                    is_finite, case)
            # Convergence criterion satisfied?
            print("# Mean of corrections to free surface: {: 1.4e}".format(error))

            criter = crit
            if ns == nstep:
                criter = 1.e-10
            #
            if error < criter * abs(z[1]) and itr > 1:
                break  # Exit for
            #
            if itr == niter:
                print("# Note that the program still had not converged to the degree specified")
            # Operations for extrapolations if more than one height step used
            if ns == 1:
                sol[1: num + 1, 2] = z[1: num + 1]
            else:
                sol[1: num + 1, 1] = sol[1: num + 1, 2]
                sol[1: num + 1, 2] = z[1: num + 1]
        #
        # Fourier coefficients (for surface elevation by slow Fourier transform)
        Y[0] = 0.0
        for j in range(1, norder + 1):
            B[j] = z[j + norder + 10]
            sumf = 0.5 * (z[10] + z[norder + 10] * np.power(-1.0, float(j)))
            sumf += np.sum([z[10 + m] * cosa[(m * j) % (norder + norder)]
                            for m in range(1, norder)])
            Y[j] = 2.0 * sumf / norder
    #
    # End stepping through wave heights
    method = "Fourier"
    title_block(method, c_type, current, z)
    output(norder, z, Y, B, Tanh, is_finite)
    return z, Y, B, Tanh
# 
# 
#
