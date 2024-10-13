# 
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
#from collections import defaultdict
#import multiprocessing
#import pickle
#from dataclasses import dataclass
from datetime import datetime as dt

# package imports
from steelpy.trave.process.dynamic import eigen, trnsient
#from steelpy.trave.utils.solution import UnSolver
from steelpy.trave.process.static import StaticSolver
from steelpy.trave.postprocess.main import PostProcess
#
#from steelpy.trave.beam.main import Beam
#from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.sqlite.utils import create_connection, create_table
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
    __slots__ = ['_plane2D', '_postprocess', '_solver', '_name', 
                 'db_file', '_build', '_name', '_Pdelta']
    
    def __init__(self, mesh,
                 name: str|None = None, 
                 sql_file: str|None = None, 
                 log: bool = False) -> None:
        """
        """
        self._name = name
        #
        plane = "3D"
        if self._plane2D:
            plane = "2D"
        #
        print('{:}'.format(52 * '-'))
        print (f"-- module : Trave{plane} version 2.50")
        print ('{:}'.format(52 * '-'))
        #
        # -----------
        # TODO : if not mesh, get file and open mesh class
        self._mesh = mesh
        self._mesh.plane(self._plane2D)
        #
        #self._solver = UnSolver(mesh=self._mesh)
        self.db_file = self._mesh.db_file
        #
        self._postprocess = PostProcess(mesh=self._mesh,
                                        name=self._name)
                                        #sql_file=self.db_file)
        #
        # create table
        conn = create_connection(self.db_file)
        with conn:
            self._new_table(conn)
    #
    # --------------------------------------------
    # SQL ops
    #
    def _new_table(self, conn) -> None:
        """ """
        table = "CREATE TABLE IF NOT EXISTS Solution (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    component_id INTEGER NOT NULL REFERENCES Component(number), \
                    type TEXT NOT NULL, \
                    Pdelta BOOL NOT NULL, \
                    plane TEXT NOT NULL, \
                    date TEXT NOT NULL,\
                    units TEXT NOT NULL,\
                    title TEXT);"
        create_table(conn, table)    
    #
    #
    def _push_data(self, conn,
                   name: int|str,
                   component: int,
                   analysis_type: str,
                   Pdelta: int, 
                   plane: str, 
                   title: str|None = None):
        """ """
        table = 'INSERT INTO Solution(name, component_id, type, Pdelta,\
                                      plane, date, units, title)\
                            VALUES(?,?,?,?,?,?,?,?)'
        #
        date = dt.now().strftime('%Y-%m-%d')
        data = (name, component, analysis_type, Pdelta, 
                plane, 'si', date, title)
        # push
        cur = conn.cursor()
        out = cur.execute(table, data)
        return out.lastrowid    
    #
    # --------------------------------------------
    #
    def static(self,
               second_order: bool = False,
               ineleastic: bool = False,
               sparse:bool=True,
               max_iter: int = 30):
        """
        Solves the static system by the Direct Stiffness Method (DSM)
        
        method : banded, frontal
        second_order : Second order (True/False)
        """
        #
        # update table        
        conn = create_connection(self.db_file)
        with conn:
            name = 'test' #self._name
            component = self._mesh._component
            analysis_type = 'static'
            Pdelta = second_order
            #
            plane = "3D"
            if self._plane2D:
                plane = "2D"
            #
            self._push_data(conn, name, component, 
                            analysis_type, Pdelta, 
                            plane)
        #
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
    def modal(self):
        """ Natural Period"""
        pass
    #
    #
    def buckling(self,
                 ineleastic: bool = False):
        """ Eigen Buckling"""
        pass
    #
    #
    def time_history(self,
                     second_order: bool = False,
                     ineleastic: bool = False,):
        """ """
        pass
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
    
    
    