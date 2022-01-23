#
# Copyright (c) 2009-2022 steelpy
#
# Python stdlib imports
from array import array
#from dataclasses import dataclass
#import math
#from typing import NamedTuple, Tuple, Union, List, Dict


# package imports

#
#
def dsvbksb(u: array, w: array, v: array, m: int, n: int, b: array):
    """
    """
    tmp = [sum([u[i][j] * b[i] for i in range(1, m + 1)]) / w[j]
           if w[j] != 0 else 0 for j in range(n + 1)]

    x = [sum([v[j][jj] * tmp[jj] for jj in range(n + 1)])
         for j in range(n + 1)]
    return x
#
#
#
