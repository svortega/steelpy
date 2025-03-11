# 
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
#from bisect import bisect_right
from dataclasses import dataclass
#from math import cosh, sinh, sqrt, sin, cos
#from typing import NamedTuple
#import re

# package imports
#from steelpy.utils.math.vector import Vector
from steelpy.formulas.pilkey.utils.operations import SingularFunction
from steelpy.formulas.pilkey.chapter11.table3 import TableBasic, e_values

#
#
#
#
# ---------------------------------------------------------------
# Pilkey 2nd ed
# TABLE 11-3
# PART B: BEAMS WITH AXIAL FORCES AND ELASTIC FOUNDATIONS: 
# LOADING FUNCTIONS
# ---------------------------------------------------------------
#
#
@dataclass
class ArbitraryLoading(SingularFunction, TableBasic):
    __slots__ = ['L', 'L1', 'E', 'G', 'As', 'I']

    def __init__(self, L: float, L1: float) -> None:
        """
        L  : length of beam
        L1 : Distance to load from the left end
        """
        if L1 > L:
            raise IOError(f'L1 ({L1}) > beam length ({L})')
        #
        self.L: float = L
        self.L1: float = L1

    #
    def q(self, x: float) -> float:
        """ Loading Function"""
        return 0    
    #
    def V(self, x: float) -> float:
        """ Shear Force"""
        return 0

    #
    def M(self, x: float) -> float:
        """ Bending Moment"""
        return 0

    #
    def T(self, x: float) -> float:
        """ Torque"""
        return 0

    #
    def theta(self, x: float, E: float, G: float, I: float) -> float:
        """ Slope = EIy' """
        return 0

    #
    def w(self, x: float, E: float, G: float, I: float) -> float:
        """ Deflection = EIy"""
        return 0
    #
    #
    def max_steps(self):
        """ """
        return [self.L1]
    

    #
    def __call__(self, x: float, E: float, G: float, As: float, I: float):
        """ 
        Formulas positive (+) is downwards and therefore sign is changed to maintain compatibility
        return: [V, M, theta, w]
        """
        #print('here')
        self.E = E
        self.I = I
        self.G = G
        self.As = As
        #
        step1 = x - self.L1
        xfun = self.function_n(step1, 1)
        ef = self.ei(x=xfun, k=0, load=True)
        #ef = self.ei(x=x)
        return [self.V(x, ef),     # Shear force
                self.M(x, ef),     # Bending moment
                self.theta(x, ef), # Slope
                self.w(x, ef)]     # Deflection
#
#
#
@dataclass
class Trapezoidal(ArbitraryLoading):
    __slots__ = ['q1', 'q2', 'L', 'L1', 'L2',
                 'L3', 'slope', 'P']

    def __init__(self, q1: float, q2: float,
                 L: float, L1: float, L2: float, 
                 P: float) -> None:
        """
        q1 : line load end 1
        q2 : line load end 2
        L1  : Distant to q1 from the left end
        L2  : Distant to q2 from the rigth end
        L  : Length of beam
        P  : axial force
        """
        super().__init__(L=L, L1=L1)
        self.q1: float = q1
        self.q2: float = q2
        self.L2: float = L2
        # Distant to q1 from the rigth end
        self.L3 = self.L - self.L2
        self._slope = (self.q2 - self.q1) / (self.L3 - self.L1)
        # axial force
        self.P = P
    #
    #
    def V(self, x: float, ef: tuple) -> float:
        """ Shear Force"""
        #step1 = x - self.L1
        step2 = x - self.L3
        #ef = self.ei(x=x, k=0)
        xfun = self.function_n(step2, 1)
        ef2 = self.ei(x=xfun, k=0)        
        #
        #func1 = (-1 * self._slope
        #         * (ef.e3 * self.function_n(step1, 1)
        #            + ef.Zeta * ef.e5 * self.function_n(step1, 1)
        #            - ef.e3 * self.function_n(step2, 1)
        #            - ef.Zeta * ef.e5 * self.function_n(step2, 1)))
        func1 = (-1 * self._slope
                 * (ef.e3 + ef.Zeta * ef.e5 - ef2.e3 - ef2.Zeta * ef2.e5))
        
        #func2 = (self.q2 * (ef.e2 * self.function_n(step2, 1)
        #                    + ef.Zeta * ef.e4 * self.function_n(step2, 1))
        #         - self.q1 * (ef.e2 * self.function_n(step1, 1)
        #                      - ef.Zeta * ef.e4 * self.function_n(step1, 1)
        #                      + 2 * ef.Zeta * ef.e4 * self.function_n(step2, 1)))
        func2 = (self.q2 * (ef2.e2 + ef2.Zeta * ef2.e4 )
                 - self.q1 * (ef.e2 - ef.Zeta * ef.e4 
                              + 2 * ef.Zeta * ef.e4))
        
        return func1 + func2

    #
    def M(self, x: float, ef: tuple) -> float:
        """ Bending Moment"""
        #step1 = x - self.L1
        step2 = x - self.L3
        #ef = self.ei(x=x, k=0)
        xfun = self.function_n(step2, 1)
        ef2 = self.ei(x=xfun, k=0)        
        #
        #func1 = (-1 * self._slope
        #         * (ef.e4 * self.function_n(step1, 1)
        #            - ef.e4 * self.function_n(step2, 1)))
        func1 = -1 * self._slope * (ef.e4 - ef2.e4 )      
        
        #func2 = (self.q2 * ef.e3 * self.function_n(step2, 1)
        #         - self.q1 * ef.e3 * self.function_n(step1, 1))
        #
        func2 = self.q2 * ef2.e3 - self.q1 * ef.e3
        
        return func1 + func2

    #
    def theta(self, x: float, ef: tuple) -> float:
        """ Slope = EIy' """
        EI = self.E * self.I
        #step1 = x - self.L1
        step2 = x - self.L3
        #ef = self.ei(x=x, k=0)
        xfun = self.function_n(step2, 1)
        ef2 = self.ei(x=xfun, k=0)        
        #
        #func1 = (-1 * self._slope
        #         * (ef.e5 * self.function_n(step1, 1)
        #            - ef.e5 * self.function_n(step2, 1)))
        func1 = -1 * self._slope * (ef.e5 - ef2.e5 )
        
        #func2 = (self.q2 * ef.e4 * self.function_n(step2, 1)
        #         - self.q1 * ef.e4 * self.function_n(step1, 1))
        func2 = self.q2 * ef2.e4 - self.q1 * ef.e4
        
        return (func1 + func2) / EI

    #
    def w(self, x: float, ef: tuple) -> float:
        """ Deflection = EIy"""
        EI = self.E * self.I
        #step1 = x - self.L1
        step2 = x - self.L3
        #ef = self.ei(x=x, k=0)
        xfun = self.function_n(step2, 1)
        ef2 = self.ei(x=xfun, k=0)
        #
        if self.As == 0:
            func1a = 0
            func2a = 0
        else:
            GAs = self.G * self.As
            #
            func1a = (ef.e4 + ef.Zeta * ef.e6
                      - ef2.e4 - ef2.Zeta * ef2.e6) / GAs
            #
            func2a = (self.q1 * (ef.e3 + ef.Zeta * ef.e5 )
                      - self.q2 * (ef2.e3 + ef2.Zeta * ef2.e5)) / GAs            
        #
        #try:
        #    GAs = self.G * self.As
        #    #func1a = (ef.e4 * self.function_n(step1, 1)
        #    #           + ef.Zeta * ef.e6 * self.function_n(step1, 1)
        #    #           - ef.e4 * self.function_n(step2, 1)
        #    #           - ef.Zeta * ef.e6 * self.function_n(step2, 1)) / GAs
        #    #
        #    func1a = (ef.e4 + ef.Zeta * ef.e6
        #              - ef2.e4 - ef2.Zeta * ef2.e6) / GAs
        #    #
        #    #func2a = (self.q1 * (ef.e3 * self.function_n(step1, 1)
        #    #                   + ef.Zeta * ef.e5 * self.function_n(step1, 1))
        #    #        - self.q2 * (ef.e3 * self.function_n(step2, 1)
        #    #                     + ef.Zeta * ef.e5 * self.function_n(step2, 1))) / GAs
        #    func2a = (self.q1 * (ef.e3 + ef.Zeta * ef.e5 )
        #              - self.q2 * (ef2.e3 + ef2.Zeta * ef2.e5)) / GAs
        #except ZeroDivisionError:
        #    #GAs = 0
        #    func1a = 0
        #    func2a = 0
        #
        #func1 = (self._slope
        #         * ((ef.e6 * self.function_n(step1, 1)
        #             - ef.e6 * self.function_n(step2, 1)) / EI
        #            - func1a))
        func1 = self._slope * ((ef.e6 - ef2.e6 ) / EI - func1a)
        #
        #func2 = ((self.q1 * ef.e5 * self.function_n(step1, 1)
        #          - self.q2 * ef.e5 * self.function_n(step2, 1)) / EI
        #         - func2a)
        func2 = (self.q1 * ef.e5 - self.q2 * ef2.e5 ) / EI - func2a
        
        return func1 + func2

    #
    def max_steps(self):
        """ """
        wl = self.L - self.L1 - self.L2
        try:
            1 / self.q1  # end 1
            try:
                1 / self.q2  # end 2
                if self.q1 == self.q2:  # uniform
                    a = self.L1 + wl / 2.0
                    b = self.L2 + wl / 2.0
                    maxM = (a + wl * (b - a) / (2 * self.L)) / self.L
                else:  # trapezoidal
                    qrad = [0.2, 0.4, 0.6, 0.8, 1.0]
                    xrad = [0.555, 0.536, 0.520, 0.508, 0.50]
                    interp = Interpolate(qrad, xrad)
                    #
                    if abs(self.q1) <= abs(self.q2):
                        rad = interp(self.q1 / self.q2)
                    else:
                        rad = interp(self.q2 / self.q1)
                        rad = 1 - rad
                    maxM = (self.L1 / self.L) + rad
            except ZeroDivisionError:  # triangular
                maxM = (self.L1 / self.L) + (1 - 0.5774)
        except ZeroDivisionError:  # triangular
            maxM = (self.L1 / self.L) + 0.5774
            #
        x_steps = [0, 1 / 4, 3 / 8, 2 / 4, 5 / 8, 3 / 4, 1, maxM]
        x_steps = sorted(list(set(x_steps)))
        x_steps = [item for item in x_steps if item <= 1]
        return x_steps
    #


#
#
@dataclass
class Point(ArbitraryLoading):
    __slots__ = ['W', 'L', 'L1', 'P']

    def __init__(self, W: float, L: float, L1: float,
                 P: float):
        """
        W : Point load
        L1 : Distant to W from the left end
        L : Length of beam
        P : axial force (zero by befault)
        """
        super().__init__(L=L, L1=L1)
        self.W = W
        self.P = P

    #
    #def q(self, x: float) -> float:
    #    """ Loading Function"""
    #    step = x - self.L1
    #    return self.function_n(step, -1) * -self.P

    #
    def V(self, x: float, ef: tuple) -> float:
        """ Shear Force"""
        #step = x - self.L1
        #xfun = self.function_n(step, 1)
        #ef = self.ei(x=xfun, k=0)
        #return (-1 * self.W
        #        * (ef.e1 * self.function_n(step, 1)
        #           + ef.Zeta * ef.e3 * self.function_n(step, 1)))
        #
        return -1 * self.W * (ef.e1 + ef.Zeta * ef.e3 )
    #
    def M(self, x: float, ef: tuple) -> float:
        """ Bending Moment"""
        #step = x - self.L1
        #xfun = self.function_n(step, 1)
        #ef = self.ei(x=xfun, k=0)
        return -1 * self.W * ef.e2 #* self.function_n(step, 1)

    #
    def theta(self, x: float, ef: tuple) -> float:
        """ Slope = EIy' """
        #step = x - self.L1
        EI = self.E * self.I
        #xfun = self.function_n(step, 1)
        #ef = self.ei(x=xfun, k=0)
        # -self.W * ef.e3 * self.function_n(step, 1) / EI
        return -self.W * ef.e3 / EI

    #
    def w(self, x: float, ef: tuple) -> float:
        """ Deflection = EIy"""
        #step = x - self.L1
        EI = self.E * self.I
        #xfun = self.function_n(step, 1)
        #ef = self.ei(x=xfun, k=0)
        #
        if self.As == 0:
            func1 = 0
        else:
            GAs = self.G * self.As
            func1 = (ef.e2 - ef.Zeta * ef.e4 ) / GAs            
        #
        #try:
        #    GAs = self.G * self.As
        #    func1 = (ef.e2 - ef.Zeta * ef.e4 ) / GAs
        #except ZeroDivisionError:
        #    func1 = 0
        #
        #return (self.W * (ef.e4 * self.function_n(step, 1) / EI
        #                  - (ef.e2 * self.function_n(step, 1)
        #                     - ef.Zeta * ef.e4 * self.function_n(step, 1)) / (G * As)))
        #
        return self.W * (ef.e4 / EI - func1)

#
#
@dataclass
class Moment(ArbitraryLoading):
    __slots__ = ['C', 'L', 'L1', 'P']

    def __init__(self, M: float, L: float, L1: float,
                 P: float) -> None:
        """
        M : Moment load
        L1 : Distant to M from the left end
        L : Length of beam
        P : axial force
        """
        super().__init__(L=L, L1=L1)
        self.C = M
        self.P = P
    #
    def V(self, x: float, ef: tuple) -> float:
        """ Shear Force"""
        #step = x - self.L1
        #ef = self.ei(x=x, k=0)
        return self.C * ef.Lambda * ef.e4 #* self.function_n(step, 1)
    #
    def M(self, x: float, ef: tuple) -> float:
        """ Bending Moment"""
        #step = x - self.L1
        #ef = self.ei(x=x, k=0)
        # (-1 * self.C
        #        * (ef.e1 * self.function_n(step, 1)
        #           - ef.Eta * ef.e3 * self.function_n(step, 1)))
        return -1 * self.C * (ef.e1 - ef.Eta * ef.e3 ) 

    #
    def theta(self, x: float, ef: tuple) -> float:
        """ Slope = EIy' """
        #step = x - self.L1
        EI = self.E * self.I
        #ef = self.ei(x=x, k=0)
        #(-1 * self.C
        # * (ef.e3 * self.function_n(step, 1)
        #    - ef.Eta * ef.e4 * self.function_n(step, 1)) / EI)        
        return -1 * self.C * (ef.e3 - ef.Eta * ef.e4 ) / EI

    #
    def w(self, x: float, ef: tuple) -> float:
        """ Deflection = EIy"""
        #step = x - self.L1
        EI = self.E * self.I
        #ef = self.ei(x=x, k=0)
        # self.C * ef.e3 * self.function_n(step, 1) / EI
        return self.C * ef.e3 / EI
#
#


