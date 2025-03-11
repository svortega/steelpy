#
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations

# from bisect import bisect_right
from dataclasses import dataclass
from math import cosh, sinh, sqrt
# import re

# package imports
from steelpy.formulas.pilkey.utils.operations import SingularFunction

#
#
# ---------------------------------------------------------------
# Pilkey 2nd ed
# TABLE 14-1
# PART B: TWISTING OF THIN-WALLED BEAMS WITH ARBITRARY LOADING:
# LOADING FUNCTIONS
# ---------------------------------------------------------------
#
#
@dataclass
class ArbitraryLoading(SingularFunction):
    __slots__ = ["L", "L1"]

    def __init__(self, L: float, L1: float) -> None:
        """ """
        self.L: float = L
        self.L1: float = L1

    #
    def Fphi(self, x: float, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        return 0

    def Fpsi(self, x: float, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        return 0

    def FT(self, x: float, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        return 0

    def FB(self, x: float, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        return 0

    def FTw(self, x: float, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        return 0

    #
    def C(self, E: float, G: float, J: float, Cw: float):
        """ """
        # try:
        return sqrt(G * J / (E * Cw))
        # except ZeroDivisionError:
        #    return 0
    #
    #def function_n(self, step: float, n: int) -> float:
    #    """<x-L>^n"""
    #    if n < 0:
    #        return 0
    #    elif step < 0:
    #        return 0
    #    elif n == 0:
    #        return 1
    #    else:
    #        return step**n
    #
    #def fun_sincos(self, fun, step: float, C: float) -> float:
    #    """<x-L>^n"""
    #    if step < 0:
    #        return 0
    #    else:
    #        return fun(C * step)
    #
    #
    def __call__(self, x: float, E: float, G: float, J: float, Cw: float):
        """
        Formulas positive (+) is downwards

        Return:
        FT : Twisting moment
        Fphi: Angle of twist
        Fpsi : Rate of angle of twist
        FB : Bimoment
        Ftw : Warping torque

        [FT, Fphi, Fpsi, FB, Tw]
        """
        return [
            1 * self.FT(x, E, G, J, Cw),
            1 * self.Fphi(x, E, G, J, Cw),
            1 * self.Fpsi(x, E, G, J, Cw),
            1 * self.FB(x, E, G, J, Cw),
            1 * self.FTw(x, E, G, J, Cw),
        ]


#
#
#
@dataclass
class TOpenConcentrated(ArbitraryLoading):
    __slots__ = ["T", "L", "L1"]

    def __init__(self, T: float, L: float, L1: float):
        """
        Concentrated Torque
        """
        super().__init__(L, L1)
        self.T: float = T

    #
    def Fphi(self, x: float, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        step = x - self.L1
        C = self.C(E=E, G=G, J=J, Cw=Cw)
        Cfun = C * self.function_n(step, 1)
        return (self.T / (C * G * J)
                * (-1 * Cfun
                   + self.fun_sincos(fun=sinh, step=step, C=C)))

    #
    def Fpsi(self, x, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        step = x - self.L1
        C = self.C(E=E, G=G, J=J, Cw=Cw)
        return (self.T / (G * J)
                * (self.function_n(step, 0)
                   - self.fun_sincos(fun=cosh, step=step, C=C)))

    #
    def FT(self, x, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        step = x - self.L1
        return -1 * self.T * self.function_n(step, 0)

    #
    def FB(self, x, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        step = x - self.L1
        C = self.C(E=E, G=G, J=J, Cw=Cw)
        return -1 * self.T / C * self.fun_sincos(fun=sinh, step=step, C=C)


#
#
@dataclass
class TOpenDistributed(ArbitraryLoading):
    __slots__ = ["T1", "T2", "L", "L1", "L2", "L3", "_slope"]

    def __init__(self, T: float, T2: float, L: float, L1: float, L2: float) -> None:
        """
        Distributed Torque
        """
        super().__init__(L, L1)
        self.T1: float = T
        self.T2: float = T2
        self.L2: float = L2
        #
        self.L3 = self.L - self.L2
        self._slope = (self.T2 - self.T1) / (self.L3 - self.L1)

    #
    def Fphi(self, x, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        step1 = x - self.L1
        step2 = x - self.L3
        C = self.C(E=E, G=G, J=J, Cw=Cw)
        #Cfun1 = C * self.function_n(step1, 1)
        #Cfun2 = C * self.function_n(step2, 1)
        #
        t1 = (self.T1 / (C**2 * G * J)
            * (
                self.fun_sincos(fun=cosh, step=step1, C=C) # cosh(Cfun1)
                - self.function_n(step1, 0)
                - C**2 * self.function_n(step1, 2) / 2.0
                - self.fun_sincos(fun=cosh, step=step2, C=C) # cosh(Cfun2)
                + self.function_n(step2, 0)
                + C**2 * self.function_n(step2, 2) / 2.0
            )
        )
        #
        t2 = (self._slope / (E * Cw)
            * (
                self.fun_sincos(fun=sinh, step=step1, C=C) / C #sinh(Cfun1) / C
                - self.function_n(step1, 1)
                - C**2 * self.function_n(step1, 2) / 6.0
                - self.fun_sincos(fun=sinh, step=step2, C=C) / C # sinh(Cfun2) / C
                + self.function_n(step2, 1)
                + C**2 * self.function_n(step2, 2) / 6.0
            )
        )

        return t1 + t2

    #
    def Fpsi(self, x, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        step1 = x - self.L1
        step2 = x - self.L3
        C = self.C(E=E, G=G, J=J, Cw=Cw)
        #Cfun1 = C * self.function_n(step1, 1)
        #Cfun2 = C * self.function_n(step2, 1)
        #
        t1 = ( self.T1 / (C * G * J)
            * (
                C * self.function_n(step1, 1)
                - self.fun_sincos(fun=sinh, step=step1, C=C) #sinh(Cfun1)
                - C * self.function_n(step2, 1)
                + self.fun_sincos(fun=sinh, step=step2, C=C) # sinh(Cfun2)
            )
        )
        #
        t2 = (-1 * self._slope / (C * G * J)
            * (
                self.fun_sincos(fun=cosh, step=step1, C=C) / C # cosh(Cfun1) / C
                - self.function_n(step1, 0) / C
                - C / 2.0 * self.function_n(step1, 2)
                - self.fun_sincos(fun=cosh, step=step2, C=C) / C # cosh(Cfun2) / C
                + self.function_n(step2, 0) / C
                + C / 2.0 * self.function_n(step2, 2)
            )
        )

        return t1 + t2

    #
    def FT(self, x, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        step1 = x - self.L1
        step2 = x - self.L3
        #
        t1 = -1 * self.T1 * (self.function_n(step1, 1) - self.function_n(step2, 1))
        #
        t2 = (-1 * self._slope / 2
              * (self.function_n(step1, 2) - self.function_n(step2, 2))
        )

        return t1 + t2

    #
    def FB(self, x, E: float, G: float, J: float, Cw: float) -> float:
        """ """
        step1 = x - self.L1
        step2 = x - self.L3
        C = self.C(E=E, G=G, J=J, Cw=Cw)
        #Cfun1 = C * self.function_n(step1, 1)
        #Cfun2 = C * self.function_n(step2, 1)
        #
        t1 = (-1 * self.T1 / C**2
              * (
                self.fun_sincos(fun=cosh, step=step1, C=C)  # cosh(Cfun1)
                - self.function_n(step1, 0)
                - self.fun_sincos(fun=cosh, step=step2, C=C) # cosh(Cfun2)
                + self.function_n(step2, 0)
            )
        )
        #
        t2 = self._slope * (
            self.function_n(step1, 1)
            - self.fun_sincos(fun=sinh, step=step1, C=C) / C # sinh(Cfun1) / C
            - self.function_n(step2, 1)
            + self.fun_sincos(fun=sinh, step=step2, C=C) / C # sinh(Cfun2) / C
        )

        return t1 + t2


#
#
