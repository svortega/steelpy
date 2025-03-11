# 
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
#from bisect import bisect_right
from dataclasses import dataclass
from math import factorial #, cosh, sinh
#import re

# package imports
from steelpy.trave.beam.pilkey.utils.operations import SingularFunction

#
#
#
# ---------------------------------------------------------------
# Pilkey 2nd ed
# TABLE 11-2
# PART B: SIMPLE BEAMS WITH ARBITRARY LOADINGS: 
# LOADING FUNCTIONS
# ---------------------------------------------------------------
#
#
@dataclass
class ArbitraryLoading(SingularFunction):
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
    #def function_n(self, step: float, n: int) -> float:
    #    """ <x-L>^n """
    #    if n < 0:
    #        return 0
    #    elif step < 0:
    #        return 0
    #    elif n == 0:
    #        return 1
    #    else:
    #        return step ** n

    #
    def max_steps(self):
        """ """
        return [self.L1]
    

    #
    def __call__(self, x: float, E: float, I: float):
        """ 
        Formulas positive (+) is downwards and therefore sign is changed to maintain compatibility
        return: [V, M, theta, w]
        """
        #
        return [self.V(x),           # Shear force
                self.M(x),           # Bending moment
                self.theta(x, E, I), # Slope
                self.w(x, E, I)]     # Deflection
#
#
@dataclass
class Trapezoidal(ArbitraryLoading):
    __slots__ = ['q1', 'q2', 'L', 'L1', 'L2',
                 '_L3', '_slope']

    def __init__(self, q1: float, q2: float,
                 L: float, L1: float, L2: float) -> None:
        """
        """
        super().__init__(L, L1)
        self.q1: float = q1
        self.q2: float = q2
        self.L2: float = L2
        #
        self._L3 = self.L - self.L2
        self._slope = (self.q2 - self.q1) / (self._L3 - self.L1)

    #
    def q(self, x: float) -> float:
        """ Loading Function"""
        step1 = x - self.L1
        step2 = x - self._L3
        func1 = (self.function_n(step1, 1) - self.function_n(step2, 1)) * -self._slope
        func2 = -self.function_n(step1, 0) * self.q1 + self.function_n(step2, 0) * self.q2
        return func1 + func2

    #
    def V(self, x: float) -> float:
        """ Shear Force"""
        step1 = x - self.L1
        step2 = x - self._L3
        func1 = -self._slope / 2 * (self.function_n(step1, 2) - self.function_n(step2, 2))
        func2 = -self.function_n(step1, 1) * self.q1 + self.function_n(step2, 1) * self.q2
        return func1 + func2

    #
    def M(self, x: float) -> float:
        """ Bending Moment"""
        step1 = x - self.L1
        step2 = x - self._L3
        func1 = - (self._slope / factorial(3) * (self.function_n(step1, 3) - self.function_n(step2, 3)))
        func2 = -0.50 * (self.q1 * self.function_n(step1, 2) - self.q2 * self.function_n(step2, 2))
        return func1 + func2

    #
    def theta(self, x: float, E: float, I: float) -> float:
        """ Slope = EIy' """
        step1 = x - self.L1
        step2 = x - self._L3
        func1 = (self.function_n(step1, 4) - self.function_n(step2, 4)) * -self._slope / (factorial(4) * E * I)
        func2 = -1 / (factorial(3) * E * I) * (
                    self.function_n(step1, 3) * self.q1 - self.function_n(step2, 3) * self.q2)
        return func1 + func2

    #
    def w(self, x: float, E: float, I: float) -> float:
        """ Deflection = EIy"""
        step1 = x - self.L1
        step2 = x - self._L3
        func1 = (self.function_n(step1, 5) - self.function_n(step2, 5)) * self._slope / (factorial(5) * E * I)
        func2 = 1 / (factorial(4) * E * I) * (self.function_n(step1, 4) * self.q1 - self.function_n(step2, 4) * self.q2)
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
    __slots__ = ['P', 'L', 'L1']

    def __init__(self, P: float, L: float, L1: float):
        """
        """
        super().__init__(L, L1)
        self.P: float = P

    #
    def q(self, x: float) -> float:
        """ Loading Function"""
        step = x - self.L1
        return self.function_n(step, -1) * -self.P

    #
    def V(self, x: float) -> float:
        """ Shear Force"""
        step = x - self.L1
        return -1 * self.P * self.function_n(step, 0)
        #

    def M(self, x: float) -> float:
        """ Bending Moment"""
        step = x - self.L1
        return -1 * self.P * self.function_n(step, 1)

    #
    def theta(self, x: float, E: float, I: float) -> float:
        """ Slope = EIy' """
        step = x - self.L1
        return -1 * self.P * self.function_n(step, 2) / (2 * E * I)

    #
    def w(self, x: float, E: float, I: float) -> float:
        """ Deflection = EIy"""
        step = x - self.L1
        return self.P * self.function_n(step, 3)  / (factorial(3) * E * I)

#
#
@dataclass
class Moment(ArbitraryLoading):
    __slots__ = ['m', 'L', 'L1']

    def __init__(self, M: float, L: float, L1: float) -> None:
        """
        """
        super().__init__(L, L1)
        self.m: float = M

    #
    def M(self, x: float) -> float:
        """ Bending Moment"""
        step = x - self.L1
        return -self.m * self.function_n(step, 0)

    #
    def theta(self, x: float, E: float, I: float) -> float:
        """ Slope = EIy' """
        step = x - self.L1
        return -self.m / (E * I) * self.function_n(step, 1)

    #
    def w(self, x: float, E: float, I: float) -> float:
        """ Deflection = EIy"""
        step = x - self.L1
        return self.m / (2 * E * I) * self.function_n(step, 2)


#
#
