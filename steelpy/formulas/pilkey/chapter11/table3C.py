#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
from dataclasses import dataclass
#from math import factorial
#from typing import NamedTuple
import re
#
# package imports
#
from steelpy.formulas.pilkey.chapter11.table3 import TableBasic
#
#
# ---------------------------------------------------------------
#
# 
#
@dataclass
class BeamBendingSupports:
    """Beam Bending Initial Parameters """

    __slots__ = ['L', 'E', 'I', 'G', 'As', '_support0']

    def __init__(self, L:float, E:float, G: float, As: float, I: float): #
        """
        L : beam lenght [m]
        E : Elastic modulus [Pa] (default steel)
        """
        self.L:float = L
        self.E:float = E
        self.G:float = G
        self.As:float = As
        self.I:float = I
    #
    def supports(self, supp1:str, supp2:str,
                 k1:float|None=None, k2:float|None=None):
        """
        support1 : suppor 1 (left) [pinned/fixed/guide/spring]
        support2 : suppor 2 (right) [pinned/fixed/guide/spring]
        if spring :
        k1 : spring end 1 [N/m]
        k2 : spring end 2 [N/m]
        """
        self._support0 = self.support_func(supp1, supp2, k1, k2)
    #
    def support_func(self, supp1:str, supp2:str, 
                     k1:float|None=None, k2:float|None=None):
        """
        support1 : suppor 1 (left) [pinned/fixed/guide/spring]
        support2 : suppor 2 (right) [pinned/fixed/guide/spring]
        if spring :
        k1 : spring 1 [N/m]
        k2 : spring 2 [N/m]
        """
        if re.match(r"\b(pinn(ed)?|hinged|roller(s)?)\b", supp1, re.IGNORECASE):
            if re.match(r"\b(pinn(ed)?)\b", supp2, re.IGNORECASE):
                return PinnedPinned(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I)
            elif re.match(r"\b(fix(ed)?|encastre(r)?|infinite)\b", supp2, re.IGNORECASE):
                return PinnedFixed(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I)
            elif re.match(r"\b(guide(d)?)\b", supp2, re.IGNORECASE):
                return PinnedGuided(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I)
            elif re.match(r"\b(spring)\b", supp2, re.IGNORECASE):
                if not k2:
                    raise IOError("Spring k2 value missing")
                return PinnedSpring(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I,
                                    k1=k1, k2=k2)
            else:
                raise IOError("unstable")
        #
        elif re.match(r"\b(fix(ed)?|encastre(r)?|infinite)\b", supp1, re.IGNORECASE):
            if re.match(r"\b(pinn(ed)?|hinged|roller(s)?)\b", supp2, re.IGNORECASE):
                return FixedPinned(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I)
            elif re.match(r"\b(fix(ed)?|encastre(r)?|infinite)\b", supp2, re.IGNORECASE):
                return FixedFixed(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I)
            elif re.match(r"\b(guide(d)?)\b", supp2, re.IGNORECASE):
                return FixedGuided(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I)
            elif re.match(r"\b(free)\b", supp2, re.IGNORECASE):
                return FixedFree(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I)
            elif re.match(r"\b(spring)\b", supp2, re.IGNORECASE):
                if not k2:
                    raise IOError("Spring k2 value missing")
                return FixedSpring(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I,
                                   k1=k1, k2=k2)
            else:
                raise IOError(f"boundary {supp2} not supported")
        #
        elif re.match(r"\b(free)\b", supp1, re.IGNORECASE):
            if re.match(r"\b(fix(ed)?|encastre(r)?|infinite)\b", supp2, re.IGNORECASE):
                return FreeFixed(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I)
            elif re.match(r"\b(spring)\b", supp2, re.IGNORECASE):
                if not k2:
                    raise IOError("Spring k2 value missing")
                return FreeSpring(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I,
                                   k1=k1, k2=k2)
            else:
                raise IOError("unstable")
        #
        elif re.match(r"\b(guide(d)?)\b", supp1, re.IGNORECASE):
            if re.match(r"\b(pinn(ed)?|hinged|roller(s)?)\b", supp2, re.IGNORECASE):
                return GuidedPinned(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I)
            elif re.match(r"\b(fix(ed)?|encastre(r)?|infinite)\b", supp2, re.IGNORECASE):
                return GuidedFixed(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I)
            elif re.match(r"\b(spring)\b", supp2, re.IGNORECASE):
                if not k2:
                    raise IOError("Spring k2 value missing")
                return GuidedSpring(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I,
                                   k1=k1, k2=k2)
            else:
                raise IOError("unstable")
        #
        elif re.match(r"\b(spring)\b", supp1, re.IGNORECASE):
            if not k1:
                raise IOError("Spring k1 value missing")
            #
            if re.match(r"\b(pinn(ed)?|hinged|roller(s)?)\b", supp2, re.IGNORECASE):
                return SpringPinned(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I,
                                   k1=k1, k2=k2)
            elif re.match(r"\b(fix(ed)?|encastre(r)?|infinite)\b", supp2, re.IGNORECASE):
                return SpringFixed(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I,
                                   k1=k1, k2=k2)
            elif re.match(r"\b(guide(d)?)\b", supp2, re.IGNORECASE):
                return SpringGuided(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I,
                                    k1=k1, k2=k2)
            elif re.match(r"\b(free)\b", supp2, re.IGNORECASE):
                return SpringFree(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I,
                                  k1=k1, k2=k2)
            elif re.match(r"\b(spring)\b", supp2, re.IGNORECASE):
                if not k2:
                    raise IOError("Spring k2 value missing")
                return SpringSpring(L=self.L, E=self.E, G=self.G, As=self.As, I=self.I,
                                    k1=k1, k2=k2)
            else:
                raise IOError(f"boundary {supp2} not supported")
        #
        else:
            raise IOError(f"boundary {supp1} not supported")
    #
    def __call__(self, F_bar:list[float]):
        """
        I : moment of inertia [m^4]
        F_bar : [FV_bar, FM_bar, Ftheta_bar, Fw_bar, P]

        return :
        Initial Parameters x = 0 (support 1)
        [V0, M0, theta0, w0]
        """
        return self._support0.initial_parameters(F_bar=F_bar)
#
# ---------------------------------------------------------------
# Pilkey 2nd ed
# TABLE 11-2
# PART C: SIMPLE BEAMS WITH ARBITRARY LOADING: 
# INITIAL PARAMETERS
# ---------------------------------------------------------------
#
#
@dataclass
class ArbitrarySupport(TableBasic):
    __slots__ = ['L', 'E', 'G', 'I', 'ef', 'P', 'As', 
                 'FV_bar', 'FM_bar', 'Ftheta_bar', 'Fw_bar']

    def __init__(self, L:float, E: float, G: float,
                 As: float, I: float) -> None:
        """
        Pilkey 2nd ed
        TABLE 11-3
        PART C: BEAMS WITH AXIAL FORCES AND ELASTIC FOUNDATIONS:
                INITIAL PARAMETER
        
        Parameters:
        L : Length of beam
        E : modulus of elasticity
        G : Shear modulus
        A : Area
        I : Moment of inertia
        """
        self.L = L
        self.E = E
        self.G = G
        self.As = As
        self.I = I
    #
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
    #
    def _spring_2(self, k2:float):
        """ """
        ef = self.e_bar
        A5 = k2 * ef.e3 - self.E * self.I * ef.e2
        A6 = self.E * self.I * (ef.e0 - ef.Eta * ef.e2) - k2 * (ef.e1 - ef.Eta * ef.e3)
        A7 = k2 * (ef.e2 - ef.Eta * ef.e4) - self.E * self.I * (ef.e1 - ef.Eta * ef.e3)
        A8 = self.E * self.I *ef.e3 - k2 * ef.e4
        A9 = k2 * self.Ftheta_bar - self.FM_bar
        return A5, A6, A7, A8, A9
    #
    def _spring_1(self, k1:float):
        """ """
        ef = self.e_bar
        A1 = self.E*self.I * (ef.e0 - ef.Eta * ef.e2) + k1 * (ef.e1 - ef.Eta * ef.e3)
        A2 = ef.e2 + k1 * ef.e3 /  (self.E * self.I)
        A3 = ef.e1 - ef.Eta * ef.e3 + k1 * (ef.e2 - ef.Eta * ef.e4) / (self.E * self.I)
        A4 = self.E * self.I * ef.e3 + k1 * ef.e4
        return A1, A2, A3, A4
    #
    #
    def Fbar(self, FV_bar:float, FM_bar:float, 
             Ftheta_bar:float, Fw_bar:float, P: float):
        """
        FV_bar:     Shear force @ x=L
        FM_bar:     Bending moment @ x=L
        Ftheta_bar: Slope @ x=L
        Fw_bar:     Deflection @ x=L
        P: Axial load
        """
        self.FV_bar = FV_bar
        self.FM_bar = FM_bar
        self.Ftheta_bar = Ftheta_bar
        self.Fw_bar = Fw_bar
        self.P = P
    #
    def initial_parameters(self, F_bar:list[float]):
        """
        F_bar : [FV_bar, FM_bar, Ftheta_bar, Fw_bar]
        """
        self.Fbar(*F_bar)
        self.e_bar = self.ei(x=self.L, k=0)
        return [self.V0,     # Shear force
                self.M0,     # Bending moment
                self.theta0, # Slope
                self.w0]     # Deflection
#
#
# ---------------------------------------------------------------
# Pinned
# ---------------------------------------------------------------
#
@dataclass
class PinnedPinned(ArbitrarySupport):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
    #
    def operator(self):
        """ """
        eb = self.e_bar
        EI = self.E * self.I
        #
        if self.As == 0:
            func1 = 0
        else:
            GAs = self.G * self.As
            func1 =  (eb.e2 + eb.Zeta * eb.e4) / GAs            
        #
        #try:
        #    GAs = self.G * self.As
        #    func1 =  (eb.e2 + eb.Zeta * eb.e4) / GAs
        #except ZeroDivisionError:
        #    func1 = 0
        
        return eb.e0 * (eb.e4 / EI - func1) - eb.e2**2     
    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()
        return ((EI * eb.e0 * self.Fw_bar
                 + eb.e2 * self.FM_bar) / op)
    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()
        #
        if self.As == 0:
            func1 = 0
        else:
            GAs = self.G * self.As
            func1 =  (eb.e2 + eb.Zeta * eb.e4) / GAs            
        #
        #try:
        #    GAs = self.G * self.As
        #    func1 =  (eb.e2 + eb.Zeta * eb.e4) / GAs
        #except ZeroDivisionError:
        #    func1 = 0
        #
        return ((-1 * eb.e2 * self.Fw_bar
                - (eb.e4 / EI -  func1) * self.FM_bar) / op)
#
@dataclass
class PinnedFixed(ArbitrarySupport):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
    #
    def operator(self):
        """ """
        eb = self.e_bar
        EI = self.E * self.I
        #
        if self.As == 0:
            func1 = 0
        else:
            GAs = self.G * self.As
            func1 =  (eb.e2 + eb.Zeta * eb.e4) / GAs            
        #
        #try:
        #    GAs = self.G * self.As
        #    func1 =  (eb.e2 + eb.Zeta * eb.e4) / GAs
        #except ZeroDivisionError:
        #    func1 = 0
        #
        return ((eb.e1 - eb.Eta * eb.e3)
                * (eb.e4 / EI - func1) - eb.e2 * eb.e3 / EI)
    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        eb = self.e_bar
        op = self.operator()
        return ((eb.e2 * self.Ftheta_bar
                 + (eb.e1 - eb.Eta * eb.e3) * self.Fw_bar) / op)
    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        eb = self.e_bar
        op = self.operator()
        EI = self.E * self.I
        #
        if self.As == 0:
            func1 = 0
        else:
            GAs = self.G * self.As
            func1 =  (eb.e2 + eb.Zeta * eb.e4) / GAs            
        #
        #try:
        #    GAs = self.G * self.As
        #    func1 =  (eb.e2 + eb.Zeta * eb.e4) / GAs
        #except ZeroDivisionError:
        #    func1 = 0
        #
        return ((-1 * self.Fw_bar * eb.e3 / EI
                 - self.Ftheta_bar * (eb.e4 / EI - func1)) / op)
#
@dataclass
class PinnedGuided(ArbitrarySupport):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
    #
    def operator(self):
        """ """
        eb = self.e_bar
        #EI = self.E * self.I
        return ((eb.e1 - eb.Eta * eb.e3)
                * (eb.e1 + eb.Zeta * eb.e3)
                +  eb.Lambda * eb.e3**2)
    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()
        return (((-1 * eb.Lambda * EI * eb.e3 * self.Ftheta_bar)
                - (eb.e1 - eb.Eta * eb.e3) * self.FV_bar)
                / op)
    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()        
        return (((-1 * (eb.e1 + eb.Zeta * eb.e3) * self.Ftheta_bar)
                + eb.e3 / EI * self.FV_bar)
                / op)
#
@dataclass
class PinnedSpring(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1']

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
        # spring stiffness
        self.k1 = k1
        self.k2 = k2
    #
    def operator(self):
        """ """
        eb = self.e_bar
        EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        return ((A5 * eb.e2
                 + A6 * (eb.e4 - (eb.e2 + eb.Zeta * eb.e4) * EI
                         / (self.G * self.As))))
    #
    @property
    def theta0(self) -> float:
        """ Bending moment"""
        eb = self.e_bar
        EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        op = self.operator()
        return ((A5 * self.Fw_bar
                 + A9 * (eb.e4 - (eb.e2 + eb.Zeta * eb.e4) * EI
                         / (self.G * self.As)))
                / op)
    
    #
    @property
    def V0(self) -> float:
        """ Shear """
        eb = self.e_bar
        EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        op = self.operator()
        return ((EI * A6 * self.Fw_bar
                 + EI * A9 * (eb.e2
                              + (eb.e2 + eb.Zeta * eb.e4) * EI
                              / (self.G * self.As)))
                / op)
#
# ---------------------------------------------------------------
# Fixed
# ---------------------------------------------------------------
#
@dataclass
class FixedPinned(ArbitrarySupport):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
    #
    def operator(self):
        """ """
        eb = self.e_bar
        EI = self.E * self.I
        #
        if self.As == 0:
            func1 = 0
        else:
            GAs = self.G * self.As
            func1 = (eb.e2 + eb.Zeta * eb.e4) / GAs
        #
        return ((eb.e2 * eb.e3 / EI
                 - (eb.e4 / EI - func1) * (eb.e1 - eb.Eta * eb.e3)))
    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()        
        return ((- self.FM_bar * eb.e3 / EI
                 - self.Fw_bar * (eb.e1 - eb.Eta * eb.e3)) / op)
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()
        #
        if self.As == 0:
            func1 = 0
        else:
            GAs = self.G * self.As
            func1 = (eb.e2 + eb.Zeta * eb.e4) / GAs
        #
        return ((self.Fw_bar * eb.e2
                 + self.FM_bar * (eb.e4 / EI - func1)) / op)
#
@dataclass
class FixedFixed(ArbitrarySupport):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
    #
    def operator(self):
        """ """
        eb = self.e_bar
        EI = self.E * self.I
        #
        if self.As == 0:
            func1 = 0
        else:
            GAs = self.G * self.As
            func1 = (eb.e2 + eb.Zeta * eb.e4) / GAs
            
        return (eb.e3**2 / EI - (eb.e2 - eb.Eta * eb.e4) * (eb.e4 / EI - func1))
    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        eb = self.e_bar
        #EI = self.E * self.I
        op = self.operator()         
        return ((- eb.e3 * self.Ftheta_bar
                 - (eb.e2 - eb.Eta * eb.e4) * self.Fw_bar) / op)

    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()
        #
        if self.As == 0:
            func1 = 0
        else:
            GAs = self.G * self.As
            func1 = (eb.e2 + eb.Zeta * eb.e4) / GAs
        #try:
        #    GAs = self.G * self.As
        #    func1 = (eb.e2 + eb.Zeta * eb.e4) / GAs
        #except ZeroDivisionError:
        #    func1 = 0
        
        return ((eb.e3 * self.Fw_bar
                 + (eb.e4 / EI - func1) * self.Ftheta_bar * EI) / op)
#
@dataclass
class FixedFree(ArbitrarySupport):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)

    def operator(self):
        """ """
        eb = self.e_bar
        #EI = self.E * self.I
        return ((eb.e1 - eb.Eta * eb.e3)
                * (eb.e1 + eb.Zeta * eb.e3)
                +eb.Lambda * eb.e2 * eb.e4)
    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        eb = self.e_bar
        #EI = self.E * self.I
        op = self.operator() 
        return ((-1 * eb.Lambda * eb.e4 * self.FM_bar
                 - (eb.e1 - eb.Eta * eb.e3) * self.FV_bar) / op)
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        eb = self.e_bar
        #EI = self.E * self.I
        op = self.operator()         
        return ((-1 * (eb.e1 + eb.Zeta * eb.e3) * self.FM_bar
                 + eb.e2 * self.FV_bar) / op)
#
@dataclass
class FixedGuided(ArbitrarySupport):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
    #
    def operator(self):
        """ """
        eb = self.e_bar
        #EI = self.E * self.I
        return ((eb.e2 - eb.Eta * eb.e4)
                * (eb.e1 + eb.Zeta * eb.e3)
                + eb.Lambda * eb.e4 * eb.e3)
    #
    @property
    def V0(self) -> float:
        """ Shear force"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()
        return ((-1 * eb.Lambda * EI * eb.e4 * self.Ftheta_bar
                - (eb.e2 - eb.Eta * eb.e4) * self.FV_bar)
                / op)
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()        
        return ((-1 * (eb.e1 + eb.Zeta * eb.e3) * EI * self.Ftheta_bar
                + eb.e3 * self.FV_bar)
                / op)
#
@dataclass
class FixedSpring(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1']

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
        # spring stiffness
        self.k1 = k1
        self.k2 = k2
    #
    def operator(self):
        """ """
        eb = self.e_bar
        EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        return (A5 * eb.e3
                - A7 * (eb.e4 - (eb.e2 + eb.Zeta * eb.e4) * EI / (self.G * self.As)))
    #
    @property
    def V0(self) -> float:
        """ Shear """
        eb = self.e_bar
        EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        op = self.operator()
        return (-1 * EI * A7 * self.Fw_bar - EI * A9 * eb.e3) / op
    
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        eb = self.e_bar
        EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        op = self.operator()
        return (EI * A5 * self.Fw_bar + EI * A9 * eb.e4) / op
#
#
# ---------------------------------------------------------------
# free
# ---------------------------------------------------------------
#
@dataclass
class FreeFixed(ArbitrarySupport):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
    #
    def operator(self):
        """ """
        eb = self.e_bar
        #EI = self.E * self.I
        return (eb.Lambda * eb.e2 * eb.e4
                + (eb.e1 + eb.Zeta * eb.e3)
                * (eb.e1 - eb.Eta * eb.e3))
    #
    @property
    def theta0(self) -> float:
        """ Bending moment"""
        eb = self.e_bar
        #EI = self.E * self.I
        op = self.operator()        
        return ((eb.Lambda * eb.e4 * self.Fw_bar
                - (eb.e1 + eb.Zeta * eb.e3) * self.Ftheta_bar)
                / op)
    
    #
    @property
    def w0(self) -> float:
        """ Shear force"""
        eb = self.e_bar
        #EI = self.E * self.I
        op = self.operator()         
        return (-1 * ((eb.e1 - eb.Eta * eb.e3) * self.Fw_bar
                      + eb.e2 * self.Ftheta_bar)
                / op)
#
@dataclass
class FreeSpring(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1']

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
        # spring stiffness
        self.k1 = k1
        self.k2 = k2
    #
    def operator(self):
        """ """
        eb = self.e_bar
        #EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        return (A6 * (eb.e1 + eb.Zeta * eb.e3)
                + A8 * eb.e2)
    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        eb = self.e_bar
        #EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        op = self.operator()
        return ((eb.Lambda * A8 * self.Fw_bar
                - A9 * (eb.e1 + eb.Zeta * eb.e3))
                / op)
    
    #
    @property
    def w0(self) -> float:
        """ Deflection"""
        eb = self.e_bar
        #EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        op = self.operator()
        return ((-1 * A6 * self.Fw_bar
                 + A9 * eb.e2) / op)
#
#
# ---------------------------------------------------------------
# guide
# ---------------------------------------------------------------
#
@dataclass
class GuidedPinned(ArbitrarySupport):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
    #
    def operator(self):
        """ """
        eb = self.e_bar
        #EI = self.E * self.I
        return ((eb.e1 - eb.Eta * eb.e3)
                * (eb.e1 * eb.Zeta * eb.e3)
                + eb.Lambda * eb.e3**2)
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator() 
        return ((eb.Lambda * EI * eb.e3 * self.Fw_bar
                - (eb.e1 + eb.Zeta * eb.e3) * self.FM_bar)
                / op)
    
    #
    @property
    def w0(self) -> float:
        """ Shear force"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()        
        return ((-1 * (eb.e1 - eb.Eta * eb.e3) * self.Fw_bar
                - (eb.e3 / EI) * self.FM_bar)
                / op)
#
@dataclass
class GuidedFixed(ArbitrarySupport):
    """ """
    __slots__ = [ 'E', 'I', 'FV', 'FM', 'Ftheta', 'Fw' ]

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
    #
    def operator(self):
        """ """
        eb = self.e_bar
        #EI = self.E * self.I
        return ((eb.e2 - eb.Eta * eb.e4)
                * (eb.e1 + eb.Zeta * eb.e3)
                + eb.Lambda * eb.e3 * eb.e4)
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()
        return ((eb.Lambda * EI * eb.e4 * self.Fw_bar
                 - EI * (eb.e1 + eb.Zeta * eb.e3) * self.Ftheta_bar)
                / op)
    
    #
    @property
    def w0(self) -> float:
        """ Shear force"""
        eb = self.e_bar
        EI = self.E * self.I
        op = self.operator()        
        return (-1 * (eb.e2 - eb.Eta * eb.e4) * self.Fw_bar
                - eb.e3 * self.Ftheta_bar) / op
#
@dataclass
class GuidedSpring(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k2', 'k1']

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
        # spring stiffness
        self.k1 = k1
        self.k2 = k2
    #
    def operator(self):
        """ """
        eb = self.e_bar
        #EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        return (eb.Lambda * A8 * eb.e3
                - A7 * (eb.e1 + eb.Zeta * eb.e3))
    #
    @property
    def w0(self) -> float:
        """ Slope"""
        eb = self.e_bar
        #EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        op = self.operator()
        return (A7 * self.Fw_bar + A9 * eb.e3) / op
    
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        eb = self.e_bar
        EI = self.E * self.I
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        op = self.operator()
        return ((eb.Lambda * EI * A8 * self.Fw_bar
                 + EI * A9 * (eb.e1 + eb.Zeta * eb.e3))
                / op)
#
# TODO : update all cases below
# ---------------------------------------------------------------
# Spring
# ---------------------------------------------------------------
#
@dataclass
class SpringPinned(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1']

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
        # spring stiffness
        self.k1 = k1
        #self.k2 = k2
    #
    def operator(self):
        """ """
        eb = self.e_bar
        EI = self.E * self.I
        A1, A2, A3, A4 = self._spring_1(self.k1)
        return ((A1 * (eb.e4 - (eb.e2 + eb.Zeta * eb.e4)
                       / (self.G * self.As)) / EI)
                - A2 * eb.e2)
    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        eb = self.e_bar
        EI = self.E * self.I
        #A1, A2, A3, A4 = self._spring_1(self.k1)
        op = self.operator()
        return ((- eb.e2 * self.Fw_bar
                - self.FM_bar * (eb.e4 - (eb.e2 + eb.Zeta * eb.e4)
                                 / (self.G * self.As)) / EI)
                / op)
    
    #
    @property
    def V0(self) -> float:
        """ Shear """
        #eb = self.e_bar
        #EI = self.E * self.I
        A1, A2, A3, A4 = self._spring_1(self.k1)
        op = self.operator()
        return (A1 * self.Fw_bar + A2 * self.FM_bar) / op
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return self.k1 * self.theta0
#
@dataclass
class SpringFixed(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1']

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
        # spring stiffness
        self.k1 = k1
        #self.k2 = k2
    #
    def operator(self):
        """ """
        eb = self.e_bar
        EI = self.E * self.I
        A1, A2, A3, A4 = self._spring_1(self.k1)
        return (A3 * (eb.e4 - (eb.e2 + eb.Zeta * eb.e4)
                      / (self.G * self.As)) / EI
                - A2 * eb.e3 / EI)
    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        eb = self.e_bar
        EI = self.E * self.I
        #A1, A2, A3, A4 = self._spring_1(self.k1)
        op = self.operator()
        return ((-1 * (eb.e3 / EI) * self.Fw_bar
                 - self.Ftheta_bar * (eb.e4 - (eb.e2 + eb.Zeta * eb.e4)
                                      / (self.G * self.As)) / EI)
                / op)
    
    #
    @property
    def V0(self) -> float:
        """ Shear """
        #eb = self.e_bar
        #EI = self.E * self.I
        A1, A2, A3, A4 = self._spring_1(self.k1)
        op = self.operator()
        return (A3 * self.Fw_bar + A2 * self.Ftheta_bar) / op
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return self.k1 * self.theta0
#
@dataclass
class SpringFree(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1']

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
        # spring stiffness
        self.k1 = k1
        #self.k2 = k2
    #
    def operator(self):
        """ """
        eb = self.e_bar
        EI = self.E * self.I
        A1, A2, A3, A4 = self._spring_1(self.k1)
        return (((self.k1 * eb.e1 + EI * eb.e0)
                 * (eb.e1 + eb.Zeta * eb.e3))
                + eb.Lambda * A4 * eb.e2)
    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        eb = self.e_bar
        #EI = self.E * self.I
        #A1, A2, A3, A4 = self._spring_1(self.k1)
        op = self.operator()
        return (-(eb.e1 + eb.Zeta * eb.e3) * self.FM_bar
                + eb.e2 * self.FV_bar) / op
    
    #
    @property
    def V0(self) -> float:
        """ Shear """
        eb = self.e_bar
        EI = self.E * self.I
        A1, A2, A3, A4 = self._spring_1(self.k1)
        op = self.operator()
        return ((- eb.Lambda * A4 * self.FM_bar
                - (self.k1 * eb.e1 + EI * eb.e0) * self.FV_bar)
                / op)
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return self.k1 * self.theta0
#
@dataclass
class SpringGuided(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw', 'k1']

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
        # spring stiffness
        self.k1 = k1
        #self.k2 = k2
    #
    def operator(self):
        """ """
        eb = self.e_bar
        EI = self.E * self.I
        A1, A2, A3, A4 = self._spring_1(self.k1)
        return (A3 * (eb.e1 + eb.Zeta * eb.e3)
                + eb.Lambda * A4 * eb.e3 / EI)
    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        eb = self.e_bar
        EI = self.E * self.I
        #A1, A2, A3, A4 = self._spring_1(self.k1)
        op = self.operator()
        return ((-1 * (eb.e1 + eb.Zeta * eb.e3) * self.Ftheta_bar
                + (eb.e3 / EI) * self.FV_bar) / op)
    
    #
    @property
    def V0(self) -> float:
        """ Shear """
        eb = self.e_bar
        #EI = self.E * self.I
        A1, A2, A3, A4 = self._spring_1(self.k1)
        op = self.operator()
        return (- eb.Lambda * A4 * self.Ftheta_bar
                - A3 * self.FV_bar) / op
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return self.k1 * self.theta0
#
@dataclass
class SpringSpring(ArbitrarySupport):
    """ """
    __slots__ = ['E', 'I', 'FV', 'FM', 'Ftheta', 'Fw',
                 'k2', 'k1']

    def __init__(self, L:float, E: float,
                 G: float, As: float, I: float, 
                 k1:float|None=None, k2:float|None=None) -> None:
        """
        """
        super().__init__(L=L, E=E, G=G, As=As, I=I)
        # spring stiffness
        self.k1 = k1
        self.k2 = k2
    #
    def operator(self):
        """ """
        eb = self.e_bar
        EI = self.E * self.I
        A1, A2, A3, A4 = self._spring_1(self.k1)
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        return (A2 * A5 + (A1 - self.k2 * A3)
                * (eb.e4 - (eb.e1 + eb.Zeta * eb.e4) * EI / (self.G * self.As)))
    #
    @property
    def theta0(self) -> float:
        """ Slope"""
        eb = self.e_bar
        EI = self.E * self.I        
        #A1, A2, A3, A4 = self._spring_1(self.k1)
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        op = self.operator()
        return ((A5 * self.Fw_bar
                + A9 * (eb.e4 - (eb.e1 + eb.Zeta * eb.e4) * EI / (self.G * self.As)))
                / op)
    
    #
    @property
    def V0(self) -> float:
        """ Shear """
        #eb = self.e_bar
        EI = self.E * self.I        
        A1, A2, A3, A4 = self._spring_1(self.k1)
        A5, A6, A7, A8, A9 = self._spring_2(self.k2)
        op = self.operator()
        return ((EI * (A1 - self.k2 * A3) * self.Fw_bar
                 - EI * A2 * A9) / op)
    #
    @property
    def M0(self) -> float:
        """ Bending moment"""
        return self.k1 * self.theta0
    #
  
#
# ---------------------------------------------------------------
#
#