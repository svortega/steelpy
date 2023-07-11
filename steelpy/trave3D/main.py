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
#from .preprocessor.mass import form_mass
from .processor.dynamic_solver import eigen, trnsient
from .processor.static_solver import solve_deflections
from .postprocessor.main import Results
#
from .processor.operations import beam_end_force, beam_int_force, comb_update, get_reactions
#
#
class TraveMain:
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u', '_results', '_m2D']
    
    def __init__(self) -> None:
        """
        """       
        print ("-- module : Trave3D version 2.50")
        print ('{:}'.format(52 * '-'))
        #self._results = Results()
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
        self._mesh = value
    #
    def run_static(self, method:str|None=None,
                   log: bool = False):
        """
        Solves the static system 
        method : banded, frontal
        """
        # get mesh
        mesh = self._mesh
        # loading
        load =  mesh.load()
        basic_load = load.basic()        
        df_nload = basic_load.node_df()
        #
        # ------------------------------
        # Static solution
        # ------------------------------
        # get K
        jbc, K = mesh.K(solver=method, log=log, m2D=self._m2D)
        # get displacements
        df_ndisp = solve_deflections(df_nload, method, m2D=self._m2D)
        #
        # ------------------------------
        # get beam end node forces
        # ------------------------------
        #
        elements = mesh.elements()
        #
        df_nforce = beam_end_force(elements, df_ndisp, m2D=self._m2D)
        #
        # -----------------------------------
        # get beam force along lenght
        # -----------------------------------
        df_membf = beam_int_force(elements, 
                                  basic_load,
                                  df_ndisp=df_ndisp, 
                                  df_nforce=df_nforce,
                                  df_nload=df_nload, 
                                  m2D=self._m2D)
        #
        # -----------------------------------
        # Updating load combinations
        # -----------------------------------
        # get load combinations
        load_combination = load.combination()
        df_comb = load_combination.to_basic()
        #
        (df_nload, df_ndisp,
         df_nforce, df_membf) = comb_update(df_comb, df_nload, df_ndisp,
                                            df_nforce, df_membf, m2D=self._m2D)
        #
        # -------------------------------------
        # node reactions
        # -------------------------------------
        #
        boundaries = mesh.boundaries()
        #
        df_reactions = get_reactions(boundaries, df_nforce, m2D=self._m2D)
        #
        # -------------------------------------
        # Results
        # -------------------------------------        
        #
        results = Results(mesh=mesh, m2D=self._m2D)
        results.postprocess(df_ndisp, df_nforce,
                            df_membf, df_reactions)
        #
        # -----------------------------------
        #print("-->")
        return results
    #
    #
#
#
#
class Trave3D(TraveMain):
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u', '_results', '_m2D']
    
    def __init__(self) -> None:
        """
        """       
        super().__init__()
        self._m2D:bool = False
    #
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
#
#
#
class Trave2D(TraveMain):
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_mesh', '_load', '_f2u', '_results', '_m2D']
    
    def __init__(self) -> None:
        """
        """       
        super().__init__()
        self._m2D:bool = True
    #
#    
    
    
    
    