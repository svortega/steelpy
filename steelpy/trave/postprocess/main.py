#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass


#
# package imports
from steelpy.trave.postprocess.sql.main import PostProcessSQL
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
    # ---------------------------------
    # Results
    # ---------------------------------
    #
    def results(self, beam_steps: int,
                Pdelta: bool):
        """ """
        return self._process.results(Un=self.Un.df,
                                     beam_steps=beam_steps,
                                     Pdelta=Pdelta)
    #
    # ---------------------------------
    #
    def run(self, beam_steps: int= 10):
        """ """
        print("** Postprocessing")
        Un = self.Un
        Pdelta: bool = self._Pdelta
        beam_force= self._process.solve(Un=Un,
                                        steps=beam_steps,
                                        Pdelta=Pdelta)
        # load comb update
        #if not Pdelta:
        #    # combination
        #    combination = self._mesh._load.combination()
        #    comb2basic = combination.to_basic()
       #    beam_force, Qn = self._process._add_comb(beam_force, Qn, comb2basic)
        # ---------------------------------
        # Process
        # element force
        #self._process._push_beam_result(beam_force)
        # Node
        #self._process._push_node_reaction(Qn)
        # element stress
        self._process.solve_stress(beam_force=beam_force)
        print('-->')
    #
    # ---------------------------------
    #
    
#
#
#
