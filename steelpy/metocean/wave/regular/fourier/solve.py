# Copyright (c) 2009 steelpy
#
from __future__ import annotations
#
# Python stdlib imports
from array import array
#from dataclasses import dataclass
#import math

# package imports
from steelpy.metocean.wave.regular.fourier.Dsvdcmp import dsvdcmp
#from steelpy.metocean.wave.regular.fourier.Dsvbksb import dsvbksb
import numpy as np
#
#
def solver(a: list, b: list, m: int,
           n: int, NP: int):
    """
    """
    # Perform decomposition
    a, w, v = dsvdcmp(a, m, n, NP)
    # Set up: see p65 of Press et al.
    wmax = np.max([w[i] if w[i] >= 0 else 0
                   for i in range(n+1)])
    
    wmin = wmax * 1.e-12
    w = [w[i] if w[i] > wmin else 0
         for i in range(n+1)]
    # Back substitute
    solution = dsvbksb(a, w, v, m, n, b)
    return solution, a
#
#
def dsvbksb(u: array, w: array, v: array, m: int, n: int, b: array):
    """
    """
    #1 / 0
    tmp = np.array([np.sum([u[i][j] * b[i] for i in range(1, m + 1)]) / w[j]
                    if w[j] != 0 else 0 for j in range(n + 1)])
    #
    x = np.sum(v[:n + 1, :n + 1] * tmp[:n + 1], axis=1)
    return np.array(x)
#