# 
# Copyright (c) 2019-2021 steelpy
# 

# Python stdlib imports
#from collections.abc import Mapping
#from typing import NamedTuple, Dict, List, Tuple
from dataclasses import dataclass
from math import factorial

# package imports

#
#
@dataclass
class SingFunction:
    __slots__ = [ 'L', 'L1' ]

    def __init__(self, L: float, L1: float) -> None:
        """
        """
        self.L: float = L
        self.L1: float = L1
    #
    def q(self, x:float) -> float:
        """ Loading Function"""
        return 0
    #
    def V(self, x:float)->float:
        """ Shear Force"""
        return 0
    #
    def M(self, x:float)->float:
        """ Bending Moment"""
        return 0
    #
    def theta(self, x:float, E:float, I:float)->float:
        """ Slope = EIy' """
        return 0
    #
    def w(self, x:float, E:float, I:float)->float:
        """ Deflection = EIy"""
        return 0
    #
    def _step(self, step:float,  n:int) -> float:
        """ <x-L>^n """
        if n < 0:
            return 0
        elif step < 0:
            return 0
        elif n == 0:
            return 1
        else:
            return step**n
#
#
@dataclass
class Trapezoidal(SingFunction):
    
    __slots__ = [ 'q1', 'q2', 'L', 'L1', 'L2',
                  '_L3', '_slope']

    def __init__(self, q1:float, q2:float,
                 L:float, L1:float, L2:float) -> None:
        """
        """
        SingFunction.__init__ (self, L, L1)
        self.q1: float = q1
        self.q2: float = q2
        self.L2: float = L2
        #
        self._L3 = self.L - self.L2
        self._slope = (self.q2 - self.q1) / (self._L3 - self.L1)
    #
    def q(self, x:float) -> float:
        """ Loading Function"""
        step1 = x - self.L1
        step2 = x - self._L3
        func1 = (self._step(step1, 1) - self._step(step2, 1)) * -self._slope
        func2 = -self._step(step1, 0)*self.q1 + self._step(step2, 0)*self.q2
        return func1 + func2
    #
    def V(self, x:float)->float:
        """ Shear Force"""
        step1 = x - self.L1
        step2 = x - self._L3
        func1 = -self._slope/2 * (self._step(step1, 2) - self._step(step2, 2))
        func2 = -self._step(step1, 1)*self.q1 + self._step(step2, 1)*self.q2
        return func1 + func2
    #
    def M(self, x:float)->float:
        """ Bending Moment"""
        step1 = x - self.L1
        step2 = x - self._L3
        func1 = - (self._slope/factorial(3) * (self._step(step1, 3) - self._step(step2, 3)))
        func2 = -0.50*(self.q1*self._step(step1, 2) - self.q2 * self._step(step2, 2))
        return func1 + func2
    #
    def theta(self, x:float, E:float, I:float)->float:
        """ Slope = EIy' """
        step1 = x - self.L1
        step2 = x - self._L3
        func1 = (self._step(step1, 4) - self._step(step2, 4))* -self._slope/(factorial(4)*E*I)
        func2 = -1/(factorial(3)*E*I)*(self._step(step1, 3)*self.q1 - self._step(step2, 3)*self.q2)
        return func1 + func2
    #
    def w(self, x:float, E:float, I:float)->float:
        """ Deflection = EIy"""
        step1 = x - self.L1
        step2 = x - self._L3
        func1 = (self._step(step1, 5) - self._step(step2, 5))*self._slope/(factorial(5)*E*I)
        func2 = 1/(factorial(4)*E*I) * (self._step(step1, 4)*self.q1 - self._step(step2, 4)*self.q2)
        return func1 + func2
    #
#
#
@dataclass
class Point(SingFunction):
    
    __slots__ = ['P', 'L', 'L1']

    def __init__(self, P: float, L: float, L1: float):
        """
        """
        SingFunction.__init__(self, L, L1)
        self.P: float = P
    #
    def q(self, x:float) -> float:
        """ Loading Function"""
        step = x - self.L1
        return self._step(step, -1) * -self.P
    #
    def V(self, x:float) -> float:
        """ Shear Force"""
        step = x-self.L1
        return -self.P * self._step(step, 0) 
    #
    def M(self, x:float) -> float:
        """ Bending Moment"""
        step = x - self.L1
        return -self.P * self._step(step, 1)
    #
    def theta(self, x: float, E:float, I:float) -> float:
        """ Slope = EIy' """
        step = x - self.L1
        return self._step(step, 2) * -self.P/(2*E*I)
    #
    def w(self, x:float, E:float, I:float) -> float:
        """ Deflection = EIy"""
        step = x - self.L1
        return self._step(step, 3) * self.P/(factorial(3)*E*I)
#
#
@dataclass
class Moment(SingFunction):
    
    __slots__ = [ 'm', 'L', 'L1' ]

    def __init__(self, m: float, L: float, L1: float) -> None:
        """
        """
        SingFunction.__init__(self, L, L1)
        self.m: float = m
    #
    def M(self, x:float) -> float:
        """ Bending Moment"""
        step = x - self.L1
        return -self.m * self._step(step, 0)
    #
    def theta(self, x: float, E:float, I:float) -> float:
        """ Slope = EIy' """
        step = x-self.L1
        return -self.m/(E*I) * self._step(step, 1)
    #
    def w(self, x:float, E:float, I:float) -> float:
        """ Deflection = EIy"""
        step = x-self.L1
        return self.m/(2*E*I) * self._step(step, 2)
    #
#
#
class SingularFunction:
    __slots__ = ['_function', '_n', '_divisor', '_C']

    def __init__(self, function:float, power:int) -> None:
        """ """
        self._function = function
        self._n = power
        self._divisor = 1.0
    #
    #def __pow__(self, power:int, modulo=None):
    #    """ """
    #    self._power = power
    #    if power < 0:
    #        return 0
    #    elif self._function < 0:
    #        return 0
    #    elif power == 0:
    #        return 1
    #    else:
    #        return self._function**power
    #
    def __repr__(self):
        """Return a string representation of self."""
        return '<{:}>^{:}'.format(self._function, self._n)
    #
    def integration(self):
        """ """
        if self._n < 0:
            self._n += 1
        else:
            self._n += 1
            self._divisor = self._n
    #
    def differentiation(self):
        """ """
        if self._n < 0:
            self._n -= 1
        else:
            self._divisor /= self._n
            self._n -= 1