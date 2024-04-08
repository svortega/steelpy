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
#from steelpy.trave.process.solution import UnSolver
from steelpy.trave.process.static import StaticSolver
from steelpy.trave.postprocess.main import PostProcess
#
from steelpy.trave.beam.main import Beam
from steelpy.utils.dataframe.main import DBframework
#
#import numpy as np
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
        #self._solver = UnSolver(mesh=self._mesh)
        #
        self._postprocess = PostProcess(mesh=self._mesh,
                                        name=name)
                                        #sql_file=sql_file)
        
    #
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
        if (mesh := self._mesh):
            static =  StaticSolver(mesh=mesh)
            #
            if second_order:
                order = "2nd"
                self._Pdelta = True
                solve = static.solvePdelta
                #Fn = mesh._load._combination.Fn()
            else:
                order = "1st"
                self._Pdelta = False
                solve = static.solveLinear
                #Fn = mesh._load._basic.Fn()
            #
            print(f"** Solving Linear Static [{order} order] ")
            #start_time = time.time()            
            Un = solve(sparse=sparse,
                       max_iter=max_iter)
            #
            self._postprocess.Un.df = Un
        else:
            raise IOError('** error: mesh missing')        
    #
    #
    def staticXX(self,
                 sparse:bool=True,
                 second_order: bool = False,
                 max_iter: int = 30):
        """
        Solves the static system by the Direct Stiffness Method (DSM)
        
        method : banded, frontal
        second_order : Second order (True/False)
        """
        #
        #static = self._solver.static()
        static =  StaticSolver(plane=self._mesh)
        #          
        #
        if (mesh := self._mesh):
            # prepare inputs
            if second_order:
                self._Pdelta = True
                solve = static.solvePdelta
                Fn = mesh._load._combination.Fn()
            else:
                self._Pdelta = False
                solve = static.solveLinear
                Fn = mesh._load._basic.Fn()
            #
            # ------------------------------
            # Get K matrix
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
            #
            Un = solve(Ke=Ke, Kg=Kg, Kt=Kt,
                       Fn=Fi, Dn=Di, 
                       jbc=jbc,
                       sparse=sparse,
                       max_iter=max_iter)
            #
            # ------------------------------
            # Load Combination
            # ------------------------------            
            #
            if not second_order:
                load_comb = mesh._load.combination()
                df_comb = load_comb.to_basic()
                #
                # displacments
                Un = self._update_ndf(dfnode=Un, dfcomb=df_comb)
            #            
            #
            self._postprocess.Un.df = Un
        
        else:
            #return static
            raise IOError('** error: mesh missing')
    #
    #
    def _update_ndf(self, dfnode, dfcomb, 
                   values:list[str] = ['x', 'y', 'z', 'rx', 'ry', 'rz']):
        """
        Update node displacements to include lcomb
        """
        db = DBframework()
        # group basic load by name
        ndgrp = dfnode.groupby('load_name')
        # get combinations with basic loads 
        #
        combgrp = dfcomb.groupby('load_name')
        for key, combfactors in combgrp:
            for row in combfactors.itertuples():
                comb = ndgrp.get_group(row.basic_load).copy()
                comb.loc[:, values] *= row.factor
                comb['load_level'] = 'combination'
                comb['load_name'] = row.load_name
                comb['load_id'] = row.load_id
                comb['load_title'] = row.load_title
                comb['component_name'] = row.component_name
                #
                try:
                    dftemp = db.concat([dftemp, comb], ignore_index=True)
                except UnboundLocalError:
                    dftemp = comb
        #
        try:
            #check = dftemp.groupby(['node_name', 'c']).sum().reset_index()
            dftemp = dftemp.groupby(['load_name', 'load_id','load_level',
                                     'load_title', 'load_system',
                                     'component_name', 'node_name'],
                                      as_index=False)[values].sum()
            #test
            dfnode = db.concat([dfnode, dftemp], ignore_index=True)
        except UnboundLocalError:
            pass
        #
        return dfnode #, memb_comb
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
    
    
    