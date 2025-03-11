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
    alpha: float
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
    def ei(self, x: float, k: float, load: bool=False):
        """Values of ei
        
        x : is measured from the left end
        k : Elastic foundation modulus
        load : load flag
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
                #alpha2 = abs(self.P / (self.E * self.I))
                #alpha = sqrt(alpha2)
                # Axial force 
                if self.P > 0: # Tension
                    # Beam with Tensile Axial Force P
                    #print(f'tension : {self.P}')
                    return self._beam_tension(x, load=load)
                    #return self._simple_beam(x, load=load)
                
                else: # compression
                    # Beam with Compressive Axial Force P
                    #print(f'compression : {self.P}')
                    return self._beam_compression(x, load=load)
                 
        except ZeroDivisionError:
            #print('no axial')
            return self._simple_beam(x, load=load)
    #
    def _simple_beam(self, x: float, load: bool):
        """ Ordinary Beam with or without Shear Deformation"""
        e1 = 1
        if load and x == 0:
            e1 = 0
        
        return e_values(e0=0,
                        e1= e1, # 1,
                        e2=x,
                        e3=x**2 / 2.0,
                        e4=x**3 / 6.0,
                        e5=x**4 / 24.0,
                        e6=x**5 / 120.0,
                        Zeta=0.0,
                        alpha=0.0,
                        Lambda=0.0,
                        Eta=0.0)        
    #
    def _beam_compression(self, x: float, load: bool):
        """ Beam with Compressive Axial Force P"""
        Zeta = abs(self.P / (self.E * self.I))
        alpha = sqrt(Zeta)
        
        e1 = round(cos(alpha*x), 6)
        if load and x == 0:
            e1 = 0
        
        return e_values(e0= -1*round(alpha*sin(alpha*x), 6),
                        e1= e1, #cos(alpha*x),
                        e2= sin(alpha*x)/alpha,
                        e3= round((1.0 - cos(alpha*x)) / alpha**2, 6),
                        e4= round((alpha*x - sin(alpha*x)) / alpha**3, 6),
                        e5= round((alpha**2 * x**2) / 2.0 + (cos(alpha*x) - 1.0), 6) / alpha**4,
                        e6= round((alpha**3 * x**3) / 6.0 + (sin(alpha*x) - alpha*x), 6) / alpha**5,
                        Zeta= Zeta,
                        alpha=alpha, 
                        Lambda=0.0,
                        Eta=0.0)
    #
    def _beam_tension(self, x: float, load: bool):
        """ Beam with Tensile Axial Force P"""
        Zeta = abs(self.P / (self.E * self.I))
        alpha = sqrt(Zeta)
        
        e1 = round(cosh(alpha*x), 6)
        if load and x == 0:
            e1 = 0
        
        return e_values(e0= round(alpha*sinh(alpha*x), 6),
                        e1= e1, #cosh(alpha*x),
                        e2= round(sinh(alpha*x)/alpha, 6),
                        e3= round((cosh(alpha*x) - 1.0) / alpha**2, 6),
                        e4= round((sinh(alpha*x) - alpha*x) / alpha**3, 6),
                        e5= round((-1*(alpha**2 * x**2) / 2.0 + (cosh(alpha*x) - 1.0)) / alpha**4, 6),
                        e6= round((-1*(alpha**3 * x**3) / 6.0 + (sinh(alpha*x) - alpha*x)) / alpha**5, 6),
                        Zeta= Zeta,
                        alpha=alpha,
                        Lambda=0.0,
                        Eta=0.0)
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

    __slots__ = ['E', 'G', 'As', 'I', #'alpha_s', 'k', 
                 'FV', 'FM', 'Ftheta', 'Fw',
                 'V0', 'M0', 'theta0', 'w0']

    def __init__(self, E: float, G: float,
                 As: float, I: float,
                 #alpha_s: float, 
                 k: float = 0) -> None:
        """
        E : modulus of elasticity
        G : Shear modulus
        As : Equivalent shear Area
        I : moment of inertia
        alpha_s : Shear correction factor
        k : elastic foundation modulus
        """
        self.E = E
        self.G = G
        self.As = As
        self.I = I
        self.k = k
        #self.alpha_s = alpha_s
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
                - self.M0 * ef.Lambda * ef.e4
                + self.FV)

    #
    def M(self, x: float, ef: tuple) -> float:
        """ Bending moment"""
        EI = self.E * self.I
        #func = self.ei(x=x, k=0)
        return (self.w0 * ef.Lambda * EI * ef.e3
                + self.theta0 * EI * (ef.e0 - ef.Eta * ef.e2)
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
        #
        if self.As == 0:
            func1 = 0
        else:
            GAs = self.G * self.As
            func1 = (ef.e2 + ef.Zeta * ef.e4) / GAs
        #
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
        return [1 * self.V(x, ef),     # Shear force
                1 * self.M(x, ef),     # Bending moment
                1 * self.theta(x, ef), # Slope
                1 * self.w(x, ef)]     # Deflection
#
#
