# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from collections import defaultdict
#import multiprocessing
#import pickle
#from dataclasses import dataclass
#from datetime import datetime as dt

# package imports
#from steelpy.trave.process.dynamic import eigen, trnsient
#from steelpy.trave.utils.solution import UnSolver
#from steelpy.trave.process.static import StaticSolver
#from steelpy.trave.postprocess.main import PostProcess
#
from steelpy.trave.preprocess.main import TraveItemBasic
#
#
#
#
class Trave3D(TraveItemBasic):
    """
    A program for static & dynamic analysis of 3D framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u', '_results', '_plane']

    def __init__(self, mesh = None,
                 sql_file:str|None = None) -> None:
        """
        """
        self._plane2D:bool = False
        super().__init__(mesh=mesh, sql_file=sql_file)
#
#
#
#
class Trave2D(TraveItemBasic):
    """
    A program for static & dynamic analysis of 2D framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u', '_results', '_plane']
    
    def __init__(self, mesh = None,
                 sql_file:str|None = None) -> None:
        """
        """
        self._plane2D:bool = True
        super().__init__(mesh=mesh, sql_file=sql_file)
#
#    
#
#
#

    
    