#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass


#
# package imports
#from steelpy.trave.postprocess.operations import MainProcess
from steelpy.trave.postprocess.sql.main import PostProcessSQL
#from steelpy.trave.postprocess.sql.Un import UnSQL
#
#
# --------------------
# Results
# --------------------
#
#
class PostProcess(PostProcessSQL):
    __slots__ = ['_mesh', '_process', '_Un',
                 'db_file', '_result_name']
    
    def __init__(self, mesh,
                 result_name:int|str,
                 db_file: str) -> None:
        """
        """
        self._mesh = mesh
        self._result_name = result_name
        #if not name:
        #    self._name = self._mesh._name
        #
        super().__init__(mesh=self._mesh,
                         result_name=self._result_name,
                         db_file=db_file)
        #
    #
    #def mesh(self, mesh):
    #    """ """
    #    self._mesh = mesh
    #    self._plane = self._mesh._plane
    #
    # --------------------
    # Results
    # --------------------
    #
    @property
    def Un(self):
        """ Nodal displacement solution"""
        return self._Un
    #
    #@Un.setter
    #def Un(self, df):
    #    """Nodal displacement solution"""
    #    self._Un = df
    #
    #
    # --------------------
    #
    #
    def results(self, beam_steps: int,
                Pdelta: bool):
        """ """
        res = self._process.results(Un=self.Un.df,
                                    beam_steps=beam_steps,
                                    Pdelta=Pdelta)
        #print('-->')
        return res
    #

#
#
