#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from array import array
#from copy import copy
#from math import fsum
#import pickle
from dataclasses import dataclass
#from typing import NamedTuple
#from itertools import chain
#import time
#
# package imports
#from steelpy.utils.math.operations import to_matrix
#from steelpy.utils.dataframe.main import DBframework
#from steelpy.utils.math.operations import remove_column_row
#
from steelpy.trave.process.static.operations import LinearStatic, PDelta
from steelpy.trave.postprocess.main import PostProcess
#


#
#
#
# ---------------------------------------------
#
@dataclass
class StaticSolver:
    """ This analysis option is used for general linear analysis for framed structures and finite elements"""
    __slots__ = ['_mesh', '_method', '_postprocess',
                 '_log', '_result_name', 'db_file',
                 'second_order', 'nonlinear',
                 '_static', '_PDelta']

    def __init__(self, mesh,
                 result_name: int | str,
                 db_file: str, log: bool,
                 second_order: bool = False,
                 nonlinear: bool = False,
                 #method: str = 'linear',
                 ) -> None:
        """
        plane : Plane system (3D/2D)
        """
        self._mesh = mesh
        self._result_name = result_name
        self.db_file = db_file
        self._postprocess = PostProcess(mesh=self._mesh,
                                        result_name=result_name,
                                        db_file=db_file)
        self._log = log
        self.second_order = second_order
        self.nonlinear = nonlinear
        #
        #
        self._static = LinearStatic(mesh=self._mesh,
                                    nonlinear=nonlinear,
                                    result_name=result_name)

        self._PDelta = PDelta(mesh=self._mesh,
                              nonlinear=nonlinear,
                              result_name=result_name)
        #
        #
        #
    #
    # ------------------------------------------
    #
    def solve(self,  
              max_iter: int = 30,
              beam_steps: int = 10,
              method: str = 'simple_step'):
        """
        Solves the static system by the Direct Stiffness Method (DSM)

        method : banded, frontal
        """
        if self._mesh:
            if self.second_order:
                order = "2nd"
                Pdelta = True
                #self._postprocess._Pdelta = True
                Un, Rn = self._PDelta.solve(max_iter=max_iter,
                                        method=method)
            else:
                order = "1st"
                #self._postprocess._Pdelta = False
                Pdelta = False
                Un, Rn = self._static.solve(max_iter=max_iter)
            #
            print(f"** Solving Linear Static [{order} order] ")
            #
            #Un = run(sparse=sparse, max_iter=max_iter)
            #
            self._postprocess.Un = Un
            self._postprocess.Rn = Rn
        else:
            raise IOError('** error: mesh missing')
        # run postprocessing
        self._postprocess.run(Pdelta=Pdelta,
                              beam_steps=beam_steps)
        #
        return self._postprocess._results
#
#
# ---------------------------------------------
#