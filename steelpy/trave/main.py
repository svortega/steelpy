# 
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
#from collections import defaultdict
#import multiprocessing
#import pickle
#from dataclasses import dataclass
#from datetime import datetime as dt

# package imports
from steelpy.trave.process.dynamic import eigen, trnsient
from steelpy.trave.process.solution import UnSolver
from steelpy.trave.postprocess.main import PostProcess
#
from steelpy.trave.beam.main import Beam
#
#
#
#
class TraveItem:
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_plane2D', '_postprocess', '_solver', 
                 'db_file', '_build', '_name']
    #'_mesh', '_load', '_f2u','_results', 
    
    def __init__(self, mesh,
                 name: str|None = None, 
                 sql_file: str|None = None, 
                 log: bool = False) -> None:
        """
        """
        plane = "3D"
        if self._plane2D:
            plane = "2D"
        #
        print (f"-- module : Trave{plane} version 2.50")
        print ('{:}'.format(52 * '-'))
        #
        # -----------
        # TODO : if not mesh, get file and open mesh class
        self._mesh = mesh
        self._mesh.plane(self._plane2D)
        #
        #try:
        #    self._name = f'{self._mesh._name}_res'
        #except AttributeError:
        #    self._name = None
        
        #super().__init__(name=name, sql_file=sql_file)
        #
        #self._postprocess = PostProcess(mesh=self._mesh,
        #                                sql_file=self.db_file)
        #
        self._solver = UnSolver(mesh=self._mesh)
        #
        #if not name:
        #    name = self._mesh._name
        #
        self._postprocess = PostProcess(mesh=self._mesh,
                                        name=name)
                                        #sql_file=sql_file)
        
    #
    #@property
    #def mesh(self):
    #    """
    #    """
    #    return self._mesh
    
    #@mesh.setter
    #def mesh(self, mesh):
    #    """
    #    """
    #    self._mesh = mesh
    #    self._mesh.plane(self._plane2D)
    #
    #
    #
    def static(self, #mesh= None, 
               method:str|None=None,
               second_order: bool = False):
        """
        Solves the static system by the Direct Stiffness Method (DSM)
        
        method : banded, frontal
        second_order : Second order (True/False)
        """
        #
        #static = StaticSolver(plane2D=self._plane2D)
        static = self._solver.static()
        #
        if self._mesh:
            #self._mesh = mesh
            #self._mesh.plane(self._plane2D)
            #
            # ------------------------------
            # Get K matrix
            # ------------------------------
            #mesh = self._mesh
            Kg = self._mesh.K(solver=method) 
            jbc = self._mesh.jbc()
            #
            # ------------------------------
            # Get load vector
            # ------------------------------        
            # 
            #load =  mesh.load()
            #basic_load = load.basic()
            #Fn_df = self._mesh.load().case().Fn()
            Fn = self._mesh._load._basic.Fn()
            #
            # ------------------------------
            # Static solution
            # ------------------------------
            #
            # plane=mesh._plane, method=method
            #      
            Udf = static.solve(Kg=Kg, Fn=Fn, jbc=jbc)
            #
            # -----------------------------------
            # Postprocess
            # -----------------------------------
            #
            #self._mesh.Un = Udf
            #
            self._postprocess.Un.df = Udf
            #self._postprocess = PostProcess(mesh=self._mesh)
        else:
            return static
        #
        #print('-->')
    #
    def dynamic(self):
        """ """
        1 / 0
    #
    #
    def results(self,
                sql_file:str|None = None,
                name:str|None = None,
                beam_steps: int= 10):
        """ """
        #1 / 0
        #
        #if mesh:
        #    mesh.plane(self._plane2D)
        #    self._postprocess = PostProcess(mesh=mesh)
        #
        # -------------------------------------
        # Results
        # -------------------------------------        
        #
        #if not name:
        #    name = f'{self._mesh._name}_res'
        #
        #
        #postprocess = PostProcess(mesh=self._mesh,
        #                          name=name,
        #                          sql_file=sql_file)
        #
        #self._postprocess = PostProcess(mesh=self._mesh)
        results = self._postprocess.results(beam_steps)
        return results
    #   
#
#
class Trave3D(TraveItem):
    """
    A program for static & dynamic analysis of 3D framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u', '_results', '_plane']
    
    def __init__(self, mesh = None,
                 sql_file:str|None = None) -> None:
        """
        """
        self._plane2D=False
        super().__init__(mesh=mesh, sql_file=sql_file)
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
    A program for static & dynamic analysis of 2D framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u', '_results', '_plane']
    
    def __init__(self, mesh = None,
                 sql_file:str|None = None) -> None:
        """
        """
        self._plane2D=True
        super().__init__(mesh=mesh, sql_file=sql_file)
        #print ("-- module : Trave2D version 2.50")
        #print ('{:}'.format(52 * '-'))
        #self._results = Results()
    #
#    
#
#
#
    
    
    