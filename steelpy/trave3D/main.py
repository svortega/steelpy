# 
# Copyright (c) 2009-2023 steelpy
#
from __future__ import annotations
# Python stdlib imports
#from collections import defaultdict
#import multiprocessing
#import pickle

# package imports
#from steelpy.f2uModel.results.main import Results
# steelpy.trave3D
from .preprocessor.mass import form_mass
from .processor.dynamic_solver import eigen, trnsient
from .processor.static_solver import solve_deflections
from .postprocessor.main import Results
#
#
class Trave3D:
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u', '_results']
    
    def __init__(self) -> None:
        """
        """       
        print ("-- module : trave3D version 2.50")
        print ('{:}'.format(52 * '-'))
    #
    #
    @property
    def mesh(self):
        """
        """
        return self._results.mesh
    
    @mesh.setter
    def mesh(self, value):
        """
        """
        #self._mesh = value
        self._results = Results(mesh=value)
    #
    def run_static(self, method:str|None=None):
        """
        Solves the static system 
        method : banded, frontal
        """
        # -----------------------------------
        df_ndisp = self.deflection(method=method)
        # -----------------------------------
        #
        self._results.postprocess(df_ndisp=df_ndisp)
        #
        #self._results.node.deflection = df_ndisp
        #
        # -----------------------------------
        #print("-->")
        return self._results
    #
    def run_dynamic(self, end_time: float, delta_t: float,
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
    def deflection(self, method: str, log: bool = False):
        """Calculate node deflection
        
        method: numpy/banded/sparse
        """
        # geometry
        mesh = self._results.mesh
        # loading
        load = self._results.load
        basic_load = load.basic()
        #df_nload = basic_load.get_node_load2(elements, nodes, boundaries)
        df_nload = basic_load.node_df()
        #
        # solution
        #
        jbc, K = mesh.K(solver=method, log=log)
        #
        #
        df_ndisp, df_nload = solve_deflections(df_nload, method)
        #
        load._df_nodal = df_nload
        return df_ndisp
#
#
#    
    
    
    
    