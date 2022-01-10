#
# Copyright (c) 2009-2022 steelpy
#
# Python stdlib imports
from array import array
from dataclasses import dataclass
import math
from typing import NamedTuple, Tuple, Union, List, Dict

# package imports

#
def dpythag(a, b):
    """
    """
    absa = abs(a)
    absb = abs(b)
    
    try:
        1/absb
        if absa > absb:
            return absa * math.sqrt(1.0 + DSQR(absb / absa))
        return absb * math.sqrt(1.0 + DSQR(absa / absb))
    except ZeroDivisionError:
    #if absb == 0.0:
        return 0.0
    #
    # return absb * math.sqrt(1.0 + DSQR(absa / absb))
# 
#
def DSQR(a):
    """
    """
    dsqrarg = float(a)
    try:
        1/dsqrarg
        return dsqrarg * dsqrarg
    except ZeroDivisionError:
    #if dsqrarg == 0.0:
        return 0.0
    #return dsqrarg * dsqrarg
#
#