#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
#
# Python stdlib imports
from array import array
from dataclasses import dataclass
#import math
#from typing import NamedTuple

# package imports
from steelpy.metocean.wave.regular.process.waveops import WaveItem
from steelpy.metocean.wave.regular.cnoidal.Solution import (Solve, hoverd, Ubar_h, eta_h, u_h, v_h,
                                                            Alpha, Q_h, lambda_d, R_h)


# from steelpy.metocean.regular.cnoidal.Elliptic import K, E
# from steelpy.metocean.regular.cnoidal.Solution import hoverd, Ubar_h
from steelpy.utils.dataframe.main import DBframework
import numpy as np
from numpy.matlib import repmat

@dataclass
class WaveCnoidal(WaveItem):

    def __init__(self, H: float, d: float,
                 T: float|None = None, Lw: float|None = None,
                 title:str|None = None,
                 infinite_depth: bool = False,
                 current: float = 0.0, c_type: int = 1,
                 order: int = 5, nstep: int = 2,
                 niter: int = 40, accuracy: float = 1e-6) -> None:
        """
        """
        super().__init__(H=H, Tw=T, Lw=Lw, d=d, title=title,
                         order=order, nstep=nstep, niter=niter,
                         accuracy=accuracy,
                         current=current, c_type=c_type,
                         infinite_depth=infinite_depth)

    #
    #def surface(self, surface_points: int = 36):
    #    """Free surface (wave) profile
    #    surface_points : Number of points on free surface (the program clusters them near crest)
    #    """
    #    self.__call__()

    def solve(self):
        """ Solver """
        try:
            self._z
        except AttributeError:
            self.order = 5
            self.method = "Cnoidal method order {:}".format(self.order)
            #
            data = self.get_parameters()
            Lw, h, alpha, delta, K, Kd, q1, m1_limit, m1, m, mm, c = CnoidalMain(*data)
            #
            self._Lw = Lw
            self.h = h
            self.alpha = alpha
            self.delta =  delta
            self.epsilon = self.H / h
            self.K = K
            self.Kd = Kd
            self.q1 = q1
            self.m1_limit = m1_limit
            self.m1 = m1
            self.m = m
            self.mm = mm
            self.c = c
        # return z, Y, B, Tanh, L, Highest
        print('------')
    #
    def get_surface(self, surface_points:int = 36):
        """ claculate wave surface """
        x, etas, phase = get_surface(L=self._Lw, h=self.h, alpha=self.alpha, epsilon=self.epsilon ,
                                     m=self.m, mm=self.mm,
                                     m1=self.m1, m1_limit=self.m1_limit,
                                     q1=self.q1, K=self.K, Kd=self.Kd,
                                     norder=self.order, nprofiles=surface_points)
        #
        data = {'x': x ,
                'eta': etas ,
                'phase': phase}
                #'z': [(1+item)*d for item in eta],
                #'time': 0 * x}
        #
        df = DBframework()
        surface = df.DataFrame(data)        
        #
        if self.title:
            surface['title'] = self.title
        else:
            surface['title'] = None
        #
        surface['type'] = 'order_1'
        self._surface = surface
        return surface[['type', 'x', 'eta', 'phase']]
    #
    #
    def get_kinematics(self, depth_points:int = 100,
                       surface_points:int|None = None):
        """get wave kinematics"""
        surface = self._surface
        if surface_points: #or not surface
            surface = self.get_surface(surface_points)
        #
        #
        dfkin = get_kinematic(L=self._Lw, h=self.h,
                              alpha=self.alpha, delta=self.delta,
                              m=self.m, mm=self.mm,
                              m1=self.m1, m1_limit=self.m1_limit,
                              q1=self.q1, K=self.K, Kd=self.Kd,
                              c=self.c, surface=surface, 
                              norder=self.order, nprofiles=depth_points)         
        #
        df = DBframework()
        dfkin = df.DataFrame(dfkin)
        dfkin['type'] = 'order_1'
        return dfkin    
#
#
#
#
# Cnoidal theory 
# Main program
#
# def CnoidalMain(h: float, t: Union[float, None], d: float,
#                Lw: Union[float, None], is_finite: bool,
#                current: float, c_type: int = 1,
#                order: int = 5, nstep: int = 2, number: int = 40, 
#                accuracy: float = 1e-6):
def CnoidalMain(MaxH: float, case: str,
                T: float|None,
                L: float|None,
                c_type: int, current: float,
                norder: int, nstep: int,
                niter: int, accuracy: float,
                Height: float, is_finite: bool):
    """
    Cnoidal theory calculations
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
    # current=0.31
    # inital values
    pi = np.pi
    g = 9.80665  # m/s^2
    H = MaxH
    #
    # This is the limiting value for m1:
    m1_limit = 1.e-8
    # If m is c1oser to 1 than this, a dif erent procedure is used,
    # whereby m is set to 1 and K is calculated iteratively from the data

    # Read data
    #
    if case == "wavelength":
        if L < 10.:
            print("The dimensionless wavelength is less than 10")
            raise Warning("Cnoidal theory should probably not be applied")
        #
        m1 = 16. * np.exp(-1*np.sqrt(0.75 * H * L * L))
        m = 1. - m1
    else:
        # Period
        if T < 10.:
            print("The dimensionless period is less than 10.")
            raise Warning("Cnoidal theory should probably not be applied")
        m1 = 16. * np.exp(-1*np.sqrt(0.75 * H * T * T))
        m = 1. - m1
    #
    norder = min(norder, 6)
    if norder < 6:
        print(f"# Solution by {norder:}-order Cnoidal theory")
    else:
        print("# Solution by 5th order Cnoidal theory with Aitken convergence enhancement")
    #
    if c_type == 1:
        Currentname = "Euler"
        ce = current
    else:
        # if (Case==2):
        Currentname = "Stokes"
        cs = current
    #
    # Input_Title_block(monitor)
    # Solving for parameter m or K iteratively
    # Initial estimate
    # if Case == 'Period':
    #    m1 = 16. * math.exp(-math.sqrt(0.75 * H * T * T))
    #    m = 1. - m1
    # else:
    #    # if Wavelength:
    #    m1 = 16. * math.exp(-math.sqrt(0.75 * H * L * L))
    #    m = 1. - m1
    #
    # Use direct iteration to solve for m, or K in the case of very long waves
    #
    L, T, K, Kd, e, ee, m, mm, q1 = Solve(T, H, norder, current, c_type, case,
                                          m, m1, m1_limit)
    #
    # Highest wave - eqn (32) of Fenton (1990)
    Highest = ((0.0077829 * L * L * L + 0.0095721 * L * L + 0.141063 * L)
               / (0.0093407 * L * L * L + 0.0317567 * L * L + 0.078834 * L + 1))
    #
    # Evaluate all the global quantities
    q = np.exp(-pi * Kd / K)
    h = hoverd(H, e, ee, m, mm, norder)
    epsilon = H / h
    alpha = Alpha(epsilon, m, mm, norder)
    delta = 4. / 3 * alpha * alpha
    Ubar_d = Ubar_h(epsilon, e, m, mm, norder) * np.sqrt(h)
    Q_d = Q_h(epsilon, m, mm, norder) * np.power(h, 1.5)
    #
    if c_type == 1:
        c = ce + Ubar_d
        cs = c - Q_d
    else:
        # if (Case==2):
        c = cs + Q_d
        ce = c - Ubar_d

    if case == "wavelength":
        T = L / c
    else:
        # if (Known,Period):
        L = lambda_d(H, K, e, ee, m, mm, norder)
    #
    R_d = R_h(epsilon, m, mm, norder) * h
    U1 = H * L * L
    U2 = U1 / 8 / pi / pi
    if U2 < 0.5:
        print("# The program is being applied to a wave with Stokes-Ursell number < 1/2.")
        raise RuntimeError("Results are unreliable")
    # Output solution summary
    Method = f"# Solution by {norder:}-order Cnoidal theory"
    print("")
    Title_block(H, T, L, Highest, Currentname, current, Method)
    print("# Solution summary:")
    Output(U2, m, L, H, T, c, ce, cs, Ubar_d, Q_d, R_d)
    #
    #print('--')
    return L, h, alpha, delta, K, Kd, q1, m1_limit, m1, m, mm, c


#
#
def get_etas(L, h, alpha, delta, epsilon, m, mm, m1, m1_limit,
             q1, K, Kd, R_d, Method, surface_points, norder):
    """ """
    # Output free surface
    #
    print(Method)
    print("")
    print("# Surface of wave - trough-crest-trough,")
    print("# note quadratic point spacing clustered around crest")
    print("# Non-dimensionalised with respect to depth")
    print("# X/d, eta/d, & check of surface pressure\n")
    s_range = surface_points // 2
    for i in range(-s_range, s_range + 1):
        x = i / surface_points
        x = L * 2 * x * abs(x)
        eta = eta_h(x / h, alpha, epsilon,
                    m, mm, m1, m1_limit, q1, K, Kd, norder)
        eta_d = eta * h
        u = u_h(x / h, eta, alpha, delta, m, mm, m1, m1_limit, q1, K, Kd, norder)
        v = v_h(x / h, eta, alpha, delta, m, mm, m1, m1_limit, q1, K, Kd, norder)
        pressure = 0.5 * (u * u + v * v) + eta - R_d / h
        print("{:8.4f} {:8.4f} {:8.0e}".format(x, eta_d, pressure))
    #
    #print('--')


#
def get_surface(L, h, alpha, epsilon, m, mm, m1, m1_limit,
                q1, K, Kd, norder, nprofiles):
    """ """
    #
    #pi = np.pi
    npt = number_steps(nprofiles)
    x = np.zeros(npt) #np.array([0.0 for i in range(npt)])
    # xx =  array('f', [0 for i in range(npt)])
    eta = np.zeros(npt) # np.array([0.0 for i in range(npt)])
    phase = np.zeros(npt)
    print("")
    for ii in range(npt):
        # xx[ii] = (pi * (ii / nprofiles) / h)
        x[ii] = L * 0.5 * (ii / nprofiles)
        phase[ii] = x[ii] / L * 360
        eta[ii] = eta_h(x[ii] / h, alpha, epsilon, m, mm, m1, m1_limit, q1, K, Kd, norder)
        eta_d = eta[ii] * h
        print("# x/d ={: 8.4f}, Phase ={: 6.1f} theta/d ={: 8.4f}"
              .format(x[ii], phase[ii], eta_d - 1))
        # print("# x/d = {:8.4f}, Phase = {:6.1f}".format(x[ii], x[ii] * 180 / pi, eta[ii]))
    #print('---')
    return x, eta, phase


#
def summary(L, h, alpha, delta, epsilon, m, mm, m1, m1_limit,
            q1, K, Kd, R_d, c,
            Method, surface_points, nprofiles, norder):
    """ """
    print(Method)
    print("# Velocity profiles\n")
    print("# All quantities are dimensionless with respect to g and/or d\n")
    print("#*********************")
    print("# y        u       v")
    print("# -     -------------")
    print("# d        sqrt(gd)")
    print("#*********************")
    for ii in range(surface_points + 1):
        x = L * 0.5 * ii / surface_points
        eta = eta_h(x / h, alpha, epsilon, m, mm, m1, m1_limit, q1, K, Kd, norder)
        eta_d = eta * h
        print("# x/d = {:8.4f}, Phase = {:6.1f}".format(x, x / L * 360))
        for i in range(nprofiles + 1):
            y = i * eta / nprofiles
            u = u_h(x / h, y, alpha, delta, m, mm, m1, m1_limit, q1, K, Kd, norder) * math.sqrt(h) + c
            v = v_h(x / h, y, alpha, delta, m, mm, m1, m1_limit, q1, K, Kd, norder) * math.sqrt(h)
            print("{:7.4f} {:7.4f} {:7.4f}".format(y * h, u, v))
            #
    # fclose(Flowfield) 
    # fflush(NULL)
    # print("\nTouch key to continue ")
    print("Finished")
    # End main program


#
#
def get_kinematic(L, h, alpha, delta, m, mm, m1, m1_limit,
                  q1, K, Kd, c,
                  surface,
                  norder, nprofiles):
    """ """
    etas = surface['eta'].to_numpy()
    xx = surface['x'].to_numpy()
    phase = surface['phase'].to_numpy()
    #
    depth_steps = np.arange(nprofiles + 1) / nprofiles
    #
    #pi = np.pi
    npt = len(etas)
    x = [xx[j] / h for j in range(npt)]
    # eta = [etas[j]  for j in range(npt)]
    z = [[(i / nprofiles) * etas[j] for i in range(nprofiles + 1)]
         for j in range(npt)]
    #
    u = [[u_h(x[ii], z[ii][i], alpha, delta, m, mm, m1, m1_limit, q1, K, Kd, norder)
          * np.sqrt(h) + c
          for ii in range(npt)] for i in range(nprofiles + 1)]
    #
    v = [[v_h(x[ii], z[ii][i], alpha, delta, m, mm, m1, m1_limit, q1, K, Kd, norder)
          * np.sqrt(h)
          for ii in range(npt)] for i in range(nprofiles + 1)]
    #
    # for ii, eta in enumerate(etas):
    #    print("# x/d = {:8.4f}, Phase = {:6.1f}".format(x[ii]*h, h*x[ii] / L * 360))
    #    for i in range(nprofiles+1):
    #        #y = i * eta / nprofiles
    #        u = u_h(x[ii], y[ii][i], alpha, delta, m, mm, m1, m1_limit, q1, K, Kd, norder) * math.sqrt(h) + c
    #        v = v_h(x[ii], y[ii][i], alpha, delta, m, mm, m1, m1_limit, q1, K, Kd, norder) * math.sqrt(h)
    #        print("{:7.4f} {:7.4f} {:7.4f}".format(y[ii][i] * h, u, v))
    #
    for ii in range(npt):
        print("# X/d = {: 8.4f}, Phase = {: 6.1f}".format(x[ii] * h, h * x[ii] / L * 360))
        # print("# X/d = {: 8.4f}, Phase = {: 6.1f}".format(x[ii]*h, h*x[ii] * 180 / pi))
        for i in range(nprofiles + 1):
            print(f'{z[ii][i] * h: 7.4f} {u[i][ii]: 7.4f} {v[i][ii]: 7.4f}')
    #
    #header = ['z', 'u', 'v']
    #
    dfkin = {'x': repmat(xx, depth_steps.size, 1).flatten('F'),
             'phase': repmat(phase, depth_steps.size, 1).flatten('F')}
    output =  {'z': z, 'u': u, 'v': v}
    dfkin.update(output)
    #print('---')
    return dfkin


#
#
# Two title blocks - at input and output
def Input_Title_block(Ffile):
    """"""
    print("# %s", Heading)
    print("# Printing input data here to check")
    print("# Height/Depth:%6.3f", H)
    if Wavelength:
        print("# Length/Depth:%7.2f", L)
    else:
        # if (Known,Period)
        print("\n# Dimensionless Period T*sqrt(g/d):%7.2f", T)
    print("\n# Current criterion: %s,  Dimensionless value:%6.3lf\n", Currentname, Current)
    print("\n%s\n", Method)


#
#
def number_steps(StpLgth):
    """
    """
    npt = max(StpLgth, 2)
    if npt % 2 == 0:
        npt += 1
    return int(npt)


#
#
##############################
# Output
##############################
#
def Title_block(H, T, L, Highest, Currentname, Current, Method):
    """"""
    # print(f"{Heading:}")
    print(f"{Method:}")
    print("# Height/Depth:{:6.3f}, {:3.0f}% of the maximum of H/d ={:6.3f} for this length:"
          .format(H, H / Highest * 100., Highest))
    print(f"# Length/Depth:{L:7.2f}")
    print(f"# Dimensionless Period T*sqrt(g/d):{T:7.2f}")
    print("# Current criterion: {:},  Dimensionless value:{:6.3f}".format(Currentname, Current))


#
def Output(U2, m, L, H, T, c, ce, cs, Ubar_d, Q_d, R_d):
    """ """
    print("# Stokes-Ursell number {:2.4f} Elliptic parameter m {:2.4f},".format(U2, m))
    # print("# Elliptic parameter m %10.7f\n", m)
    print("# Integral quantities")
    print("# Solution non-dimensionalised by g & mean depth")
    print("# Water depth                        (d) {:8.4f}".format(1.))
    print("# Wave length                   (lambda) {:8.4f}".format(L))
    print("# Wave height                        (H) {:8.4f}".format(H))
    print("# Wave period                      (tau) {:8.4f}".format(T))
    print("# Wave speed                         (c) {:8.4f}".format(c))
    print("# Eulerian current                 (u1_) {:8.4f}".format(ce))
    print("# Stokes current                   (u2_) {:8.4f}".format(cs))
    print("# Mean fluid speed in frame of wave (U_) {:8.4f}".format(Ubar_d))
    print("# Volume flux                        (Q) {:8.4f}".format(Q_d))
    print("# Bernoulli constant                 (R) {:8.4f}".format(R_d))
