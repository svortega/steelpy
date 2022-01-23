#
# Copyright (c) 2009-2022 steelpy
#
# Python stdlib imports
from array import array
from dataclasses import dataclass
import math
from typing import NamedTuple, Tuple, Union, List, Dict

# package imports
from steelpy.metocean.regular.fourier.Dsvdcmp import dsvdcmp
from steelpy.metocean.regular.fourier.Dsvbksb import dsvbksb
#
#
def solver(a: List[array], b: array, m: int, n: int, NP: int):
    """
    """
    # Perform decomposition
    a, w, v = dsvdcmp(a, m, n, NP)
    # Set up: see p65 of Press et al.
    wmax = max([w[i] if w[i] >= 0 else 0
                for i in range(n+1)])
    
    wmin = wmax * 1.e-12
    w = [w[i] if w[i] > wmin else 0
         for i in range(n+1)]
    # Back substitute
    solution = dsvbksb(a, w, v, m, n, b)
    return solution, a
#
#
