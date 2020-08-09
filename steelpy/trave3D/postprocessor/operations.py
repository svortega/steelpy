# 
# Copyright (c) 2009-2020 fem2ufo
#
# Python stdlib imports
from typing import List, ClassVar, Dict, NamedTuple, Union
#from array import array
#from itertools import chain
#from copy import  copy
#import math as math
#import pickle
#
# package imports
#from steelpy.frame3D.processor.operations import zeros, to_matrix
#
#


class Results(NamedTuple):
    """ Basic load transfer"""
    name: Union[str, int]
    title: str
    load_type: str
    items: List[List[float]]