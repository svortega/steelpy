#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
#
# Python stdlib imports
#from array import array
#from dataclasses import dataclass
#import cmath
#import math
#from typing import NamedTuple

# package imports
import numpy as np

#
#
#
###############################################################################################################
# Evaluation of elliptic integrals and functions, using approximations of Fenton and Gardener-Garden for m > 1/2
###############################################################################################################
#
#
def Elliptic_integrals(m: float, m1_limit: float):
    """ """
    pi = np.pi
    m1 = 1. - m
    m4 = np.pow(m, 0.25)

    if m1 > m1_limit:
        K = 2. / np.pow(1 + m4, 2) * np.log(2 * (1 + m4) / (1 - m4))

    Kd = 2 * pi / np.pow(1. + m4, 2)

    q1 = np.exp(-pi * K / Kd)
    e = ((2. - m) / 3. + pi / 2 / K / Kd
         + 2. * np.pow(pi / Kd, 2) * (-1. / 24. + q1 * q1 / np.pow(1 - q1 * q1, 2)))
    #
    ee = np.zeros(6)
    mm = np.zeros(6)    
    ee[1] = e
    mm[1] = m
    for i in range(2, 5 + 1):
        ee[i] = ee[i - 1] * e
        mm[i] = mm[i - 1] * m

    return e, ee, mm, Kd

#
# **********************************************
# Elliptic functions
# **********************************************
#
# Bring back to range [-2K,+2K]
def Shift(u: float, K: float) -> float:
    """ """
    N = np.trunc(u / (4 * K))
    uu = u - N * 4 * K
    
    if uu < -2 * K:
        uu += 4 * K

    if uu > 2 * K:
        uu -= 4 * K
    
    return uu


#
def sn(u: float, m: float, m1: float, m1_limit: float,
       q1: float, K: float, Kd: float) -> float:
    """"""
    factor = 1.
    if m1 > m1_limit:
        factor = np.pow(m, -0.25)
    w = 0.5 * np.pi * Shift(u, K) / Kd
    term = (factor * (np.sinh(w) - q1 * q1 * np.sinh(3 * w))
            / (np.cosh(w) + q1 * q1 * np.cosh(3 * w)))
    return term


#
def cn(u: float, m: float, m1: float, m1_limit: float,
       q1: float, K: float, Kd: float) -> float:
    """ """
    factor = 1.
    if m1 > m1_limit:
        factor = 0.5 * np.pow(m1 / m / q1, 0.25)
    w = 0.5 * np.pi * Shift(u, K) / Kd
    term = (factor * (1 - 2 * q1 * np.cosh(2 * w))
            / (np.cosh(w) + q1 * q1 * np.cosh(3 * w)))
    return term
#
def dn(u: float, m: float, m1: float, m1_limit: float,
       q1: float, K: float, Kd: float) -> float:
    """ """
    factor = 1.
    if m1 > m1_limit:
        factor = 0.5 * np.pow(m1 / q1, 0.25)
    w = 0.5 * np.pi * Shift(u, K) / Kd
    term = (factor * (1 + 2 * q1 * np.cosh(2 * w))
            / (np.cosh(w) + q1 * q1 * np.cosh(3 * w)))
    return term
#
#
