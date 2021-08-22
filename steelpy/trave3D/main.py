# 
# Copyright (c) 2009-2021 steelpy
#

# Python stdlib imports
from array import array
#import datetime
from typing import List, ClassVar, Dict, NamedTuple, Union
from itertools import chain


#import multiprocessing

# package imports
#from scipy.linalg import solve_banded
#
from steelpy.trave3D.processor.static_solver import UDUt, solve_forces, solve_displacement
from steelpy.trave3D.preprocessor.assemble import assemble_banded_matrix




#
#
class Trave3D:
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u']
    
    def __init__(self) -> None:
        """
        """       
        #self._units = Units()
        #self._mesh = Mesh()
        print ("-- module : trave3D version 2.10")
        print ('{:}'.format(52 * '-'))
    
    
    @property
    def f2u_model(self):
        """
        """
        return self._f2u
    
    @f2u_model.setter
    def f2u_model(self, model):
        """
        """
        self._f2u = model
        #self._mesh = model.mesh
        #self._load = model.load
    
    #
    @property
    def load(self):
        """
        """
        return self._load
    
    @load.setter
    def load(self, value):
        """
        """
        self._load = value
    #
    def run_static(self):
        """
        Solves the static system 
        """
        materials = self._f2u.materials
        sections = self._f2u.sections
        elements = self._f2u.mesh.elements
        nodes = self._f2u.mesh.nodes
        boundaries = self._f2u.mesh.boundaries
        #
        basic_load = self._f2u.load.basic
        load_combination = self._f2u.load.combination
        #
        #
        assemble_banded_matrix(elements, nodes, boundaries)
        ndisp = solve_displacement(elements, nodes, boundaries,
                                   materials, sections,
                                   basic_load, load_combination)
        self._f2u._results.node.deflection = ndisp
        #
        membf, nreacs = solve_forces(elements, self._f2u._results.node.deflection)
        #                             basic_load, load_combination)
        #
        nbound = {key:nodes[key].index for key, item in boundaries.node.items()
                  if sum(item[:6])!= 0}
        reactions = []
        for lname, nreac in nreacs:
            reactions.extend([[key, lname, 'global', *nreac[index]] 
                              for key, index in nbound.items()])
        self._f2u._results.node.reaction = reactions
        #
        self._f2u._results.element.force = membf
        #print("-->")
    #
    #
    def print_results(self):
        """
        """
        self._f2u._results.node.print_deflections()
        self._f2u._results.node.print_reactions()
        self._f2u._results.element.print_forces()
    #
    #
#
#    
    
    
    
    