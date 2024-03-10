#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
from dataclasses import dataclass
from math import pi, sin, cos, sinh, cosh, sqrt
from typing import NamedTuple
#import re
#
# package imports
#from .table2B import BeamLoad
#from .table2C import BeamBendingSupports
#
#
#
# ---------------------------------------------------------------
#
class e_values(NamedTuple):
    """ ei values """
    e0: float
    e1: float
    e2: float
    e3: float
    e4: float
    e5: float
    e6: float
    Zeta: float
    Lambda: float
    Eta: float
#
#
@dataclass
class TableBasic:
    """
    Pilkey 2nd ed
    TABLE 11-3 (continued)
    PART A: BEAMS WITH AXIAL FORCES AND ELASTIC FOUNDATIONS:
            GENERAL RESPONSE EXPRESSIONS
    """
    def ei(self, x: float, k: float):
        """Values of ei
        
        x : is measured from the left end
        k : Elastic foundation modulus 
        """
        try:
            # Ordinary Beam with or without
            # Shear Deformation
            1 / float(self.P)
            try:
                # Beam on Elastic Foundation, k, with Shear Deformation Effects
                1 / k
                print('elastic foundation')
                raise NotImplementedError('soil k to be included')
            
            except ZeroDivisionError:
                alpha2 = abs(self.P / (self.E * self.I))
                alpha = sqrt(alpha2)
                # Axial force 
                if self.P > 0: # Tension
                    # Beam with Tensile Axial Force P
                    print('tension')
                    return e_values(e0= alpha*sinh(alpha*x),
                                    e1= cosh(alpha*x),
                                    e2= sinh(alpha*x)/alpha,
                                    e3= (cosh(alpha*x) - 1.0) / alpha**2,
                                    e4= (sinh(alpha*x) - alpha*x) / alpha**3,
                                    e5= (-1*alpha**2 * x**2 / 2.0 + cosh(alpha*x) - 1.0) / alpha**4,
                                    e6= (-1*alpha**3 * x**3 / 6.0 + sinh(alpha*x) - alpha*x) / alpha**5,
                                    Zeta=alpha2, Lambda=0, Eta=0) 
                
                else: # compression
                    # Beam with Compressive Axial Force P
                    print('compression')
                    return e_values(e0= -1*alpha*sin(alpha*x),
                                    e1= cos(alpha*x),
                                    e2= sin(alpha*x)/alpha,
                                    e3= (1.0 - cos(alpha*x)) / alpha**2,
                                    e4= (alpha*x - sin(alpha*x)) / alpha**3,
                                    e5= (alpha**2 * x**2 / 2.0 + cos(alpha*x) - 1.0) / alpha**4,
                                    e6= (alpha**3 * x**3 / 6.0 + sin(alpha*x) - alpha*x) / alpha**5,
                                    Zeta= alpha2, Lambda=0, Eta=0)
                 
        except ZeroDivisionError:
            # set values
            return e_values(e0=0,
                            e1=1,
                            e2=x,
                            e3=x**2 / 2.0,
                            e4=x**3 / 6.0,
                            e5=x**4 / 24.0,
                            e6=x**5 / 120.0,
                            Zeta=0,
                            Lambda=0,
                            Eta=0)
    #
    def Pcr(self):
        """P critical"""
        Pcr = K1 * pi**2 * self.E * self.I / self.L ** 2
#
#
#
# ---------------------------------------------------------------
#    Pilkey 2nd ed
#    TABLE 11-3
#    PART A: BEAMS WITH AXIAL FORCES AND ELASTIC FOUNDATIONS: 
#    GENERAL RESPONSE EXPRESSIONS
# ---------------------------------------------------------------
#
#
#
@dataclass
class BendingGE(TableBasic):
    """
    Pilkey 2nd ed
    TABLE 11-3
    PART A: BEAMS WITH AXIAL FORCES AND ELASTIC FOUNDATIONS:
            GENERAL RESPONSE EXPRESSIONS
    """

    __slots__ = ['E', 'G', 'I', 'P', 
                 'FV', 'FM', 'Ftheta', 'Fw',
                 'V0', 'M0', 'theta0', 'w0']

    def __init__(self, E: float, G: float,
                 A: float, I: float,
                 k: float = 0) -> None:
        """
        E : modulus of elasticity
        I : moment of inertia
        G : Shear modulus
        k : elastic foundation modulus
        """
        self.E = E
        self.G = G
        self.A = A
        self.I = I
        self.k = k
    #
    def load(self, FV: float, FM: float,
             Ftheta: float, Fw: float)-> None:
             #P: float) -> None:
        """
        Load @ x
        """
        self.FV = FV
        self.FM = FM
        self.Ftheta = Ftheta
        self.Fw = Fw
        #self.P = P
    #
    def R0(self, V0: float, M0: float,
           theta0: float, w0: float,
           P: float) -> None:
        """ 
        Initial Parameters
        """
        self.V0 = V0
        self.M0 = M0
        self.theta0 = theta0
        self.w0 = w0
        #
        self.P = P
    #
    def V(self, x: float, ef: tuple) -> float:
        """ Shear force"""
        EI = self.E * self.I
        #func = self.ei(x=x, k=0)
        return (self.w0 * ef.Lambda * EI * (ef.e2 + ef.Zeta * ef.e4)
                - self.theta0 * ef.Lambda * EI * ef.e3
                + self.V0 * (ef.e1 + ef.Zeta * ef.e3)
                + self.FV)

    #
    def M(self, x: float, ef: tuple) -> float:
        """ Bending moment"""
        EI = self.E * self.I
        #func = self.ei(x=x, k=0)
        return (self.w0 * ef.Lambda * EI * ef.e3
                + self.theta0 * EI * (ef.e0 - ef.Eta * ef.e3)
                + self.V0 * ef.e2
                + self.M0 * (ef.e1 - ef.Eta * ef.e3)
                + self.FM)

    #
    def theta(self, x: float, ef: tuple) -> float:
        """ Slope"""
        EI = self.E * self.I
        #func = self.ei(x=x, k=0)
        return (self.w0 * ef.Lambda * ef.e4
                + self.theta0 * (ef.e1 - ef.Eta * ef.e3)
                + self.V0 * ef.e3 / EI
                + self.M0 * (ef.e2 - ef.Eta * ef.e4) / EI
                + self.Ftheta)

    #
    def w(self, x: float, ef: tuple) -> float:
        """ Deflection"""
        EI = self.E * self.I
        #func = self.ei(x=x, k=0)
        try:
            GAs = self.G * self.A / ef.Zeta
            func1 = (ef.e2 + ef.Zeta * ef.e4) / GAs
        except ZeroDivisionError:
            func1 = 0
        
        return (self.w0 * (ef.e1 + ef.Zeta * ef.e3)
                - self.theta0 * ef.e2
                - self.V0 * (ef.e4 / EI - func1)
                - self.M0 * ef.e3 / EI
                + self.Fw)
    #
    def response(self, x:float) -> list[float]:
        """
        x is measured from the left end
        
        results: 
        [V, M, theta, w]
        """
        # shear & deflection signed change
        ef = self.ei(x=x, k=0)
        return [-1 * self.V(x, ef),     # Shear force
                self.M(x, ef),          # Bending moment
                self.theta(x, ef),      # Slope
                -1 * self.w(x, ef)]     # Deflection
#
#
