#
# Copyright (c) 2019-2023 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
from dataclasses import dataclass
from math import factorial, sqrt, cos, sin
#from typing import NamedTuple
#import re
#
# package imports

#
#
#
#
@dataclass
class ArbitraryLoading:
    __slots__ = ['L', 'L1']

    def __init__(self, L: float, L1: float) -> None:
        """
        """
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
    def theta(self, x: float, E: float, I: float) -> float:
        """ Slope = EIy' """
        return 0

    #
    def phi(self, x: float, E: float, G: float, Cw: float, K: float) -> float:
        """ angle of rotation at a distance x from the left end (radians) """
        return 0

    #
    def w(self, x: float, E: float, I: float) -> float:
        """ Deflection = EIy"""
        return 0

    #
    def function_n(self, step: float, n: int) -> float:
        """ <x-L>^n """
        if n < 0:
            return 0
        elif step < 0:
            return 0
        elif n == 0:
            return 1
        else:
            return step ** n

    #
    def max_steps(self):
        """ """
        return [self.L1]

    #
    def __call__(self, x: float, E: float, I: float):
        """

        return: [V, M, w, theta]
        """
        # return Vector([self.V(x), self.M(x),
        #               self.theta(x, E, I), self.w(x, E, I)])
        return [self.V(x), self.M(x), self.theta(x, E, I), self.w(x, E, I)]
    #
    #
    def k(self, I:float, E:float):
        """ """
        return sqrt(self.P/(E*I))
    #
    def Fa(self, step: float, E: float, I: float):
        """ """
        k = self.k(I=I, E=E)
        #
        F1 = cos(k * step)
        F2 = k * self.function_n(step=step, n=1)
        F3 = 1 - cos(k * step)
        #
        Fa1 = self.function_n(step=step, n=0) * F1
        Fa2 = sin(F2)
        Fa3 = self.function_n(step=step, n=0) * F3
        Fa4 = F2 - Fa2
        Fa5 = k**2/2 * self.function_n(step=step, n=2) - Fa3
        Fa6 = k**3/6 * self.function_n(step=step, n=3) - Fa4
        return [Fa1, Fa2, Fa3, Fa4, Fa5, Fa6]
#
#
@dataclass
class Trapezoidal(ArbitraryLoading):
    __slots__ = ['wa', 'wl', 'L', 'L1', 'L2',
                 '_L3', '_slope', 'P']

    def __init__(self, P:float, q1: float, q2: float,
                 L: float, L1: float, L2: float) -> None:
        """
        """
        super().__init__(L, L1)
        self.wa: float = q1
        self.wl: float = q2
        self.L2: float = L2
        self.P: float = P
        #
        self._L3 = self.L - self.L2
        self._slope = (self.q2 - self.q1) / (self._L3 - self.L1)
    #
    #
    def V(self, x: float, E: float, I: float) -> float:
        """ Transversal Shear Force"""
        k = self.k(I=I, E=E)
        step1 = x - self.L1
        step2 = x - self._L3
        FaI = self.Fa(step=step1, I=I, E=E)
        FaII = self.Fa(step=step2, I=I, E=E)
        #
        func1 = (-1 * self.wa * FaI[1] + self.wl * FaII[1]) / k
        func2 = -1 * self._slope / k**2 * (FaI[2] - FaII[2])
        return func1 + func2
    #
    def M(self, x: float, E: float, I: float) -> float:
        """ Bending Moment"""
        k = self.k(I=I, E=E)
        step1 = x - self.L1
        step2 = x - self._L3
        FaI = self.Fa(step=step1, I=I, E=E)
        FaII = self.Fa(step=step2, I=I, E=E)
        #
        func1 = -1 * (self.wa * FaI[2] - self.wl * FaII[2]) / k**2
        func2 = -1 * self._slope / k**3 * (FaI[3] - FaII[3])
        return func1 + func2
    #
    def theta(self, x: float, E: float, I: float) -> float:
        """ Slope = EIy' """
        k = self.k(I=I, E=E)
        step1 = x - self.L1
        step2 = x - self._L3
        FaI = self.Fa(step=step1, I=I, E=E)
        FaII = self.Fa(step=step2, I=I, E=E)
        #
        func1 = -1 * (self.wa * FaI[3] - self.wl * FaII[3]) / (k * self.P)
        func2 = -1 * self._slope / (k**2 * self.P) * (FaI[4] - FaII[4])
        return func1 + func2
    #
    def w(self, x: float, E: float, I: float) -> float:
        """ Deflection = EIy"""
        k = self.k(I=I, E=E)
        step1 = x - self.L1
        step2 = x - self._L3
        FaI = self.Fa(step=step1, I=I, E=E)
        FaII = self.Fa(step=step2, I=I, E=E)
        #
        func1 = -1 * (self.wa * FaI[4] - self.wl * FaII[4]) / (k**2 * self.P)
        func2 = -1 * self._slope / (k**3 * self.P) * (FaI[5] - FaII[5])
        return func1 + func2