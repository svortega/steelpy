# 
# Copyright (c) 2009-2023 steelpy
#
from __future__ import annotations
# Python stdlib imports
#from collections import defaultdict
#import multiprocessing
#import pickle
#from dataclasses import dataclass

# package imports
#from steelpy.f2uModel.mesh.main import MeshPlane
# steelpy.trave3D
#from .preprocessor.mass import form_mass
from .processor.dynamic_solver import eigen, trnsient
from .processor.static_solver import StaticSolver
#from .postprocessor.main import Results
#
from .processor.operations import ElementProcess
#
from .beam.main import Beam
#
#
#
#
#
#
class TraveItem:
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u', '_results', 
                 '_process', '_plane2D']
    
    def __init__(self, log: bool = False) -> None:
        """
        """
        plane = "3D"
        if self._plane2D:
            plane = "2D"
        #
        print (f"-- module : Trave{plane} version 2.50")
        print ('{:}'.format(52 * '-'))
    #
    @property
    def mesh(self):
        """
        """
        return self._mesh
    
    @mesh.setter
    def mesh(self, value):
        """
        """
        self._mesh = value
        self._mesh.plane(self._plane2D)
        #
        #
        
    #
    def static(self, method:str|None=None,
               second_order: bool = False):
        """
        Solves the static system
        method : banded, frontal
        second_order : Second order (True/False)
        """
        # get mesh
        mesh = self._mesh
        elements = self._mesh.elements()
        boundaries = self._mesh.boundaries()        
        # loading
        load =  mesh.load()
        #basic_load = load.basic()
        #df_nload = basic_load.node_df()
        #
        # ------------------------------
        #
        static = StaticSolver(plane=mesh._plane)
        #
        self._process = ElementProcess(elements=elements,
                                       boundaries=boundaries, 
                                       #load=load, 
                                       plane=mesh._plane)        
        #
        # ------------------------------
        # Static solution
        # ------------------------------
        # get K
        K = mesh.K(solver=method)
        # get displacements
        #df_ndisp = solve_deflections(df_nload, method, m2D=self._m2D)
        #
        jbc = mesh.jbc()
        static.Kglobal(jbc=jbc, Ka=K)
        static.load(load)
        df_ndisp = static.deflection(method=method)
        # ------------------------------
        # get beam end node forces
        # ------------------------------
        #
        #elements = mesh.elements()
        #
        #df_nforce = beam_end_force(elements, df_ndisp, m2D=self._m2D)
        #
        #self._elements.elements(elements)
        #df_nforce = self._elements.end_force(df_ndisp)
        #
        # -----------------------------------
        # get beam force along lenght
        # -----------------------------------
        #df_membf = beam_int_force(elements, 
        #                          basic_load,
        #                          df_ndisp=df_ndisp, 
        #                          df_nforce=df_nforce,
        #                          df_nload=df_nload, 
        #                          m2D=self._m2D)
        #
        self._process.input(load, df_ndisp)
        #results = self._process.internal_force(df_ndisp, beam_steps= 10)
        #
        # -----------------------------------
        # Updating load combinations
        # -----------------------------------
        # get load combinations
        #load_combination = load.combination()
        #df_comb = load_combination.to_basic()
        #
        #(df_nload, df_ndisp,
        # df_nforce, df_membf) = comb_update(df_comb, df_nload, df_ndisp,
        #                                    df_nforce, df_membf, m2D=self._m2D)
        #
        # -------------------------------------
        # node reactions
        # -------------------------------------
        #
        #boundaries = mesh.boundaries()
        #
        #df_reactions = get_reactions(boundaries, df_nforce, m2D=self._m2D)
        #
        # -------------------------------------
        # Results
        # -------------------------------------        
        #
        #results = Results(mesh=mesh, postprocess=df_process)
        #results.postprocess(df_ndisp, df_nforce,
        #                    df_membf, df_reactions)
        #
        # -----------------------------------
        #print("-->")
        #return results
    #
    def dynamic(self):
        """ """
        pass
    #
    #
    def solve(self, beam_steps: int= 10):
        """ """
        results = self._process.internal_force(beam_steps=beam_steps)
        return results
#
#
class Trave3D(TraveItem):
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u', '_results', '_plane']
    
    def __init__(self) -> None:
        """
        """
        self._plane2D=False
        super().__init__()
        #
        #print ("-- module : Trave3D version 2.50")
        #print ('{:}'.format(52 * '-'))
        #self._results = Results()        
    #
    #
    def dynamic(self, end_time: float, delta_t: float,
                    type:str|None=None):
        """
        Solves the dynamic system
        
        end_time: simulation time [seconds]
        deltat : time increment [seconds]
        type : modal, time history
        mass : lumped, consistent
        damping : 
        """
        file = open( "stfmx.f2u", "rb" )
        jbc = pickle.load( file )
        stf = pickle.load( file )
        mass =  pickle.load( file )
        file.close()
        #
        ibandm = 1
        #load = []
        npt =  int(end_time // delta_t)
        #
        trnsient(stf, mass, jbc, npt, 
                 #load,
                 #disp, vel, acc, fmag, olddis,
                 #wk, damp, loadin, maxnode,
                 ibandm)        
        #
        print('--')
    #
    def nfreq(self):
        """ """
        # geometry
        elements = self._mesh.elements()
        nodes = self._mesh.nodes()
        boundaries = self._mesh.boundaries()
        # pure python solution
        assemble_banded_Kmatrix(elements, nodes, boundaries)
        #
        mss, ibandm = form_mass(elements, nodes, boundaries)
        eigen(ibandm=ibandm, ivib=2)
        print('--')
    #
    #
#
#
class Trave2D(TraveItem):
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u', '_results', '_plane']
    
    def __init__(self) -> None:
        """
        """
        self._plane2D=True
        super().__init__()
        #print ("-- module : Trave2D version 2.50")
        #print ('{:}'.format(52 * '-'))
        #self._results = Results()
    #
#    
#
#
#
    
    
    