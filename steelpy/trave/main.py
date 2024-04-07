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
from steelpy.utils.dataframe.main import DBframework
#
import numpy as np
#
#
#
class TraveItem:
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_plane2D', '_postprocess', '_solver', 
                 'db_file', '_build', '_name', '_Pdelta']
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
        self._solver = UnSolver(mesh=self._mesh)
        #
        self._postprocess = PostProcess(mesh=self._mesh,
                                        name=name)
                                        #sql_file=sql_file)
        
    #
    #
    def staticX(self, #mesh= None, 
               sparse:bool=True,
               second_order: bool = False):
        """
        Solves the static system by the Direct Stiffness Method (DSM)
        
        method : banded, frontal
        second_order : Second order (True/False)
        """
        #
        #static = StaticSolver(plane2D=self._plane2D)
        static = self._solver.static()
        self._Pdelta = False
        #
        if (mesh := self._mesh):
            #self._mesh = mesh
            #self._mesh.plane(self._plane2D)
            #
            # ------------------------------
            # Get K matrix
            # ------------------------------
            #mesh = self._mesh
            Ks = mesh.Ke(sparse = sparse)
                              # condensed = False) 
            jbc = mesh.jbc()
            #Ddof = self._mesh._nodes.DOF_unreleased()
            #
            # ------------------------------
            # Get load vector
            # ------------------------------        
            #
            #self._mesh._load._basic._nodes.displacement
            #
            #load =  self._mesh.load()
            #basic_load = load.basic()
            #Fn_df = self._mesh.load().case().Fn()
            Fn = mesh._load._basic.Fn()
            #
            # ------------------------------
            # Static solution
            # ------------------------------
            #
            #      
            Udf = static.solve(Ks=Ks,
                               Fn=Fn,
                               jbc=jbc,
                               sparse = sparse)
            #
            #static.PMT(K=Kg, basic=basic_load, jbc=jbc)
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
    def static(self,
               sparse:bool=True,
               second_order: bool = False,
               max_iter: int = 30):
        """
        Solves the static system by the Direct Stiffness Method (DSM)
        
        method : banded, frontal
        second_order : Second order (True/False)
        """
        #
        static = self._solver.static()
        #
        self._Pdelta = False
        solve = static.solveLinear
        if second_order:
            self._Pdelta = True
            solve = static.solvePdelta
        #
        if (mesh := self._mesh):
            #
            # ------------------------------
            # Get Ke matrix
            # ------------------------------
            #
            Ke = mesh.Ke
            Kg = mesh.Kg
            Kt = mesh.Kt
            jbc = mesh.jbc
            #
            # ------------------------------
            # Get load vector
            # ------------------------------        
            #
            Fn = mesh._load._basic.Fn()
            #
            
            colgrp = ['load_name', 'load_id', 
                      'load_level', 'load_title',
                      'load_system', 'component_name',
                      'node_name', 'node_index']
            #
            hforce =  mesh._plane.hforce
            hdisp = mesh._plane.hdisp
            #
            Fi = Fn[colgrp+hforce]
            Di = Fn[colgrp+hdisp]
            #
            # ------------------------------
            # Linear solution
            # ------------------------------                
            Un = solve(Ke=Ke, Kg=Kg, Kt=Kt,
                       Fn=Fi, Dn=Di, 
                       jbc=jbc,
                       sparse=sparse,
                       max_iter=max_iter)
            #
            self._postprocess.Un.df = Un
        
        else:
            return static        
    #
    def Pdelta(self, sparse:bool=True,
               max_iter: int = 30):
        """ Performs second order (P-Delta) analysis """
        #
        pdelta = self._solver._Pdelta()
        self._Pdelta = True
        #
        # ------------------------------
        if (mesh := self._mesh):
            # ------------------------------
            # Get K stiffness matrix
            # ------------------------------
            #
            Ke = mesh.Ke(sparse = sparse)
            jbc = mesh.jbc()
            #
            # ------------------------------
            # Get load vector
            # ------------------------------        
            #
            Fn = mesh._load._basic.Fn()
            #
            # ------------------------------
            # Linear solution
            # ------------------------------
            #
            colgrp = ['load_name', 'load_id', 
                      'load_level', 'load_title',
                      'load_system', 'component_name',
                      'node_name', 'node_index']
            #
            hforce =  mesh._plane.hforce
            hdisp = mesh._plane.hdisp
            #
            Fi = Fn[colgrp+hforce]
            Di = Fn[colgrp+hdisp]
            #
            Un = pdelta.solveLinear(Ke=Ke,
                                    Fn=Fi,
                                    Dn=Di, 
                                    jbc=jbc,
                                    sparse = sparse)
            #
            # ------------------------------
            # Pdelta solution
            # ------------------------------            
            #
            colgrp = ['load_name', 'load_id', 'load_level',
                      'load_system', 'component_name',]
            #
            Ugrp = Un.groupby(colgrp)
            Dgrp = Di.groupby(colgrp)
            #
            hdisp = ['node_name', *mesh._plane.hdisp]
            #
            Utemp = []
            for key, noded in Ugrp:
                # TODO: select basic only
                if key[2] != 'basic':
                    continue
                Uii = noded[hdisp].set_index('node_name')
                #ndisp.set_index('node_name', inplace=True)
                #
                # ------------------------------
                # Assambly matrix start            
                Kg = mesh.Kg(D=Uii)
                Kt = Ke + Kg
                #
                Ustep = Dgrp.get_group(key)
                #Ui = Ustep.copy()
                #
                Ui = pdelta.solveLinear(Ks=Kt,
                                        Fn=Fi,
                                        Dn=Ustep,
                                        jbc=jbc,
                                        sparse = sparse)
                #
                #
                Us = Ui[hdisp].set_index('node_name')
                Utemp.extend([[*key,  noded['load_title'].iloc[0],
                              nname, *Us.loc[nname]]
                              for x, nname in enumerate(jbc.index)])
                #
                #
                #for x in range(max_iter):
                #    #
                #    # ------------------------------
                #    # solve displacements
                #    #
                #    Ui = pdelta.solveLinear(Ks=Kt,
                #                            Fn=Fi,
                #                            Dn=Ustep,
                #                            jbc=jbc,
                #                            sparse = sparse)
                #    #
                #    #
                #    # ------------------------------
                #    #
                #    Ui = Ui[hdisp].set_index('node_name')
                #    print(Ui)
                #    #
                #    #try:
                #    #Ud = Uii + Ui
                #    #print(Ud)
                #    #
                #    #Ui.loc[:,'x'] = Ui.loc[:,'x'].add(Ud.loc[:,'x'])
                #    #Ui = ndisp + Ud
                #    Uii.loc[:,'x'] =  Ui.loc[:,'x'] - Uii.loc[:,'x']
                #    print('------------')
                #    print(Uii)
                #    #except UnboundLocalError:
                #    #    pass
                #    #
                #    # ------------------------------
                #    # Assambly matrix start i
                #    #
                #    #Uii = Ui[hdisp]
                #    #Uii.set_index('node_name', inplace=True)
                #    Kg = self._mesh.Kg(D=Uii)
                #    # Tangent stiffness matrix Kt
                #    Kt = Ke + Kg
                #    #
                #    # M = scipy.linalg.det(Kt)
                #    # M = np.linalg.det(Kt)
                #    #
                #    Uii = Ui.copy()
                #    #1 / 0
        #
        #1 / 0
        Us = self.df(Utemp)
        self._postprocess.Un.df = Us
    #
    #
    def df(self, dftemp):
        """displacement dataframe"""
        db = DBframework()
        header = ['load_name', 'load_id', 'load_level',
                  'load_system', 'component_name', 'load_title', 
                  'node_name', *self._mesh._plane.hdisp]
        return db.DataFrame(data=dftemp, columns=header, index=None)
    # 
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
        results = self._postprocess.results(beam_steps,
                                            Pdelta=self._Pdelta)
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
    
    
    