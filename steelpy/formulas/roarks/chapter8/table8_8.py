#
# Copyright (c) 2019-2023 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
from dataclasses import dataclass
import math
from typing import NamedTuple
#import re
#
# package imports

#
#
#



#
# ---------------------------------------------------------------
#    Roark's Formulas Stress and Strain 7ed
#    TABLE 8.8
#    Shear, moment, slope, and deflection formulas for beams
#    under simultaneous axial compression and transverse loading
# ---------------------------------------------------------------
#
#
#
@dataclass
class GeneralExpressions:
    """
    Roark's Formulas Stress and Strain 7ed
    TABLE 8.8
    Shear, moment, slope, and deflection formulas for beams
    under simultaneous axial compression and transverse loading
    """

    __slots__ = ['FR', 'FM', 'Ftheta', 'Fy', 'FA',
                 'RA', 'MA', 'thetaA', 'yA', 'E']

    def __init__(self, E: float) -> None:
        """
        """
        self.E = E
        #pass
    #
    def load(self, FV: float, FM: float,
             Ftheta: float, Fw: float, FA:float) -> None:
        """
        Load @ x
        """
        self.FR = FV
        self.FM = FM
        self.Ftheta = Ftheta
        self.Fy = Fw
        self.FA = FA

    #
    def R0(self, V0: float, M0: float,
           theta0: float, w0: float) -> None:
        """
        Initial Parameters
        """
        self.RA = V0
        self.MA = M0
        self.thetaA = theta0
        self.yA = w0

    #
    def V(self) -> float:
        """ Transverse Shear Force"""
        return self.RA * F1 - self.MA*k*F2 - self.thetaA*self.P*F1 - self.FV*Fa1

    #
    def M(self, x: float) -> float:
        """ Bending moment"""
        return (self.MA*F1
                + (self.RA * F2 + self.thetaA*self.P*F2 - self.FV*Fa2) / k)

    #
    def theta(self, x: float, E: float, I: float) -> float:
        """ Slope"""
        return (self.thetaA*F1
                + (self.MA * k * F2  + self.RA * F3 - self.FV * Fa3 ) / self.P)

    #
    def w(self, x: float, E: float, I: float) -> float:
        """ Deflection"""
        return (self.yA + self.tethaA*F2/k + self.MA*F3/self.P
                + (self.RA*F4 - self.FA * Fa4) / (k*self.P))
    #
    def response(self, x:float, I: float, E: float=2.05e11) -> list[float]:
        """
        x : distance from end 1

        results:
        [V, M, theta, w]
        """
        return [self.V(), self.M(x),
                self.theta(x, E=self.E, I=I), self.w(x, E=self.E, I=I)]
    #
    def k(self, P:float, I:float):
        """ """
        return math.sqrt(P/(self.E*I))
    #
    def Fn(self, x:float, I=float):
        """ """
        k = self.k(P=self.FA, I=I)
        F1 = math.cos(k*x)
        F2 = math.sin(k * x)
        F3 = 1 - math.cos(k * x)
        F4 = k*x - math.sin(k * x)
        return [F1, F2, F3, F4]
#
#
# ---------------------------------------------------------------
#
#
@dataclass
class ArbitrarySupport:
    __slots__ = ['L', 'E', 'I',
                 'FV_bar', 'FM_bar', 'Ftheta_bar', 'Fw_bar']

    def __init__(self, L:float, E: float) -> None:
        """
        Roark's Formulas Stress and Strain 7ed
        TABLE 8.8
        Shear, moment, slope, and deflection formulas for beams
        under simultaneous axial compression and transverse loading
        """
        self.L:float = L
        self.E:float = E

    #
    @property
    def V0(self) -> float:
        """ Shear """
        return 0
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return 0
    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        return 0
    #
    @property
    def w0(self) -> float:
        """ Deflection"""
        return 0
    #
    def k(self, P:float, I:float):
        """ """
        return math.sqrt(P/(self.E*I))
    #
    def Ca(self, I:float):
        """ """
        k = self.k(P=self.P, I=I)
        alpha = self.L1
        kl_alpha = k*(self.L-alpha)
        Ca1 = math.cos(kl_alpha)
        Ca2 = math.sin(kl_alpha)
        Ca3 = 1 - Ca1
        Ca4 = kl_alpha - Ca2
        Ca5 = k**2 / 2 * (self.L - alpha)**2 - Ca3
        Ca6 = k**3/6 * (self.L-alpha)**3 -Ca5
        return [Ca1, Ca2, Ca3, Ca4, Ca5, Ca6]
    #
    def Cn(self, I:float):
        """ """
        k = self.k(P=self.P, I=I)
        C1 = math.cos(k*self.L)
        C2 = math.sin(k*self.L)
        C3 = 1 - C1
        C4 = k*self.L - C2
        return [C1, C2, C3, C4]
#
#
#
# ---------------------------------------------------------------
# Pinned
# ---------------------------------------------------------------
#
@dataclass
class PinnedFixed(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw']

    def __init__(self, L: float, E: float,
                 k1: float | None = None, k2: float | None = None) -> None:
        """
        """
        super().__init__(L, E)
    #
    @property
    def V0(self) -> float:
        """ Transverse Shear """
        return self.W * (C2* - C1*Ca4)/(C2*C3 - C1*C4)
    #
    @property
    def theta0(self) -> float:
        """ Slope """
        return self.W/self.P * (C4*Ca3 - C3*Ca4)/(C2*C3 - C1*C4)
    #
    def V1(self):
        """ Reaction end 2 """
        RB = self.W - self.V0
        MB = (-1*self.W/k
              * ((k*self.L*math.sin(k*alpha) - k*alpha*math.sin(k*self.L))
                /(math.sin(k*self.L - k*self.L * math.cos(k*self.L)))))
        thetaB = 0
        yB = 0
        return [RB, MB, thetaB, yB]
#
#
@dataclass
class PinnedPinned(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw']

    def __init__(self, L: float, E: float,
                 k1: float | None = None, k2: float | None = None) -> None:
        """
        """
        super().__init__(L, E)
    #
    @property
    def V0(self) -> float:
        """ Transverse Shear """
        alpha = self.L1
        return self.W/self.L * (self.L * - alpha)
    #
    @property
    def theta0(self) -> float:
        """ Slope """
        alpha = self.L1
        return self.W/self.P * (math.sin(k*(self.L-alpha))/math.sin(k*self.L)
                                - (self.L-alpha)/self.L)
#
# ---------------------------------------------------------------
# free
# ---------------------------------------------------------------
#
#
@dataclass
class FreeFixed(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw']

    def __init__(self, L: float, E: float,
                 k1: float | None = None, k2: float | None = None) -> None:
        """
        """
        super().__init__(L, E)

    #
    @property
    def theta0(self) -> float:
        """ Slope """
        return self.W/self.P * Ca[2]/C[0]

    #
    @property
    def w0(self) -> float:
        """ Deflection """
        return -1 * self.W/(k * self.P) * (C2 * Ca3 - C1 * Ca4)/C1

    #
    def V1(self):
        """ Reaction end 2 """
        RB = self.W
        MB = -1*self.W/k * (C2*Ca3 + C2*Ca4)/C2
        thetaB = 0
        yB = 0
        return [RB, MB, thetaB, yB]
#
#
# ---------------------------------------------------------------
# guide
# ---------------------------------------------------------------
#
@dataclass
class GuidedFixed(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw']

    def __init__(self, L: float, E: float,
                 k1: float | None = None, k2: float | None = None) -> None:
        """
        """
        super().__init__(L, E)
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return self.W/k * Ca3/C2
    #
    @property
    def w0(self) -> float:
        """ Deflection """
        return -1 * self.W/(k * self.P) * (C3 * Ca3 - C2 * Ca4)/C2
    #
    def V1(self):
        """ Reaction end 2 """
        alpha = self.L1
        RB = self.W
        MB = -1*self.W/k * (math.cos(k*alpha) - math.cos(k*self.L))/math.sin(k*self.L)
        thetaB = 0
        yB = 0
        return [RB, MB, thetaB, yB]
#
@dataclass
class GuidedPinned(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw']

    def __init__(self, L: float, E: float,
                 k1: float | None = None, k2: float | None = None) -> None:
        """
        """
        super().__init__(L, E)
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        alpha = self.L1
        return self.W/k * math.sin(k*(self.L-alpha))/math.cos(k*self.L)
    #
    @property
    def w0(self) -> float:
        """ Deflection """
        alpha = self.L1
        return -1 * self.W/(k * self.P) * (math.sin(k*(self.L-alpha))/math.cos(k*self.L)
                                           - k*(self.L-alpha))
    #
    def V1(self):
        """ Reaction end 2 """
        alpha = self.L1
        RB = self.W
        MB = 0
        thetaB = self.W/self.P * (math.cos(k*alpha)/math.cos(k*self.L) - 1.0)
        yB = 0
        return [RB, MB, thetaB, yB]
#
# ---------------------------------------------------------------
# Fixed
# ---------------------------------------------------------------
#
@dataclass
class FixedFixed(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw']

    def __init__(self, L: float, E: float,
                 k1: float | None = None, k2: float | None = None) -> None:
        """
        """
        super().__init__(L, E)
    #
    @property
    def V0(self) -> float:
        """ Transverse Shear """
        return self.W * (C3 * Ca3 - C2*Ca4) / (C3**2 - C2*C4)
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return -1.0*self.W/k * (C3 * Ca3 - C3*Ca4) / (C3**2 - C2*C4)
    #
    def V1(self):
        """ Reaction end 2 """
        alpha = self.L1
        RB = self.W - self.V0
        MB = self.M0 + self.V0*self.L - self.W * (self.L - alpha)
        thetaB = 0
        yB = 0
        return [RB, MB, thetaB, yB]