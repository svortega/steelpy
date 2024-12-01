#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
#from array import array
#from copy import copy
#from math import fsum
#import pickle
from dataclasses import dataclass
#from typing import NamedTuple
#from itertools import chain
import time
#
# package imports
from steelpy.utils.math.operations import to_matrix
from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.math.operations import remove_column_row
#
from steelpy.trave.postprocess.main import PostProcess
#
import numpy as np
#from scipy.linalg import cholesky_banded, cho_solve_banded
from scipy.sparse.linalg import spsolve
#
#
# ---------------------------------------------
# solver using numpy
# ---------------------------------------------
#
def solver_np(a, b):
    """
    a : Global stiffness matrix
    b : Global force vector
    """
    #nloads = df2.stack() #.values
    ndisp = np.linalg.solve(a, b)
    #if not np.allclose(np.dot(stf, ndisp), nload): #, rtol=1e-05, atol=1e-08
    #    raise RuntimeError('Solution fail')
    #
    return ndisp
#
def solver_sparse(a, b):
    """
    """
    ndisp = spsolve(a.tocsr(), b)
    return ndisp
#
#
# ---------------------------------------------
#
# ---------------------------------------------
#
@dataclass
class StaticSolverX:
    """ Linear static solver class"""
    __slots__ = ['_plane',  '_method']
    
    def __init__(self, plane) -> None:
        """
        plane : Plane system (3D/2D)
        """
        self._plane = plane
    #
    #
    def solve(self, Ks, Fn, jbc, #Ddof, 
              sparse: bool=False, 
              method:str|None = None):
        """
        Linear Static Analysis (1st Order)
        
        Input: 
        Ks  : Global stiffness matrix
        Fn  : FER Node load dataframe
        jbc : Node boundary condition dataframe
        
        Return: 
        Udf : Node displacement global system dataframe
        """
        #
        order = "1st"
        print(f"** Solving U = K^-1 F [{order} order] ")
        start_time = time.time()
        #
        maxstiff = Ks.max()
        # Select solver
        if sparse:
            solver = solver_sparse
            Ks = Ks.tolil() #copy=True
        else:
            solver = solver_np
        #if self._method == 'banded':
        #    solver = solver_Mbanded
        # Solution f = Ku
        #Us = self.DSM(Kg, Fn, jbc, solver)
        Us = self.PMT(Ks, Fn, jbc,
                      solver, maxstiff)
        #
        #
        uptime = time.time() - start_time
        print(f"** Finish Time: {uptime:1.4e} sec")
        return Us
    #
    def DSM(self, Ks, Fn, jbc, solver):
        """
        Direct Stiffness Method
        
        Ks : Global stiffness matrix
        Fn : Global force & displacement matrix
        jbc = Nodes with boundary
        
        Return:
        U : Global node displacement
        """
        # TODO : must be better ways to manipulate dfs
        db = DBframework()
        # group FER node load
        colgrp = ['load_name', 'load_id', 
                  'load_level','load_system',
                  'component_name']        
        #
        # Select nodes' free DOF
        DOFbool, DOFzeros = self._mask_DOF(jbc)
        # jbc Matrix to vector 
        jbcflat = jbc.stack()
        # group FER node load
        fergrp = Fn.groupby(colgrp)
        # Matrix condensed
        #Ks = self._Mcond(Ks, jbcflat)
        #
        Utemp = []
        for key, item in fergrp:
            # set FER df by nodes
            fernodes = item.set_index(['node_name'])
            # Select load nodes with free DOF
            DOFindex = DOFzeros.index.isin(item['node_name'])
            #
            # Select load vector comprising nodes with free DOF 
            Fs = DOFzeros.copy() #.astype('float64')
            Fs.loc[DOFindex] = fernodes[self._plane.hforce] 
            Fs = Fs.stack() # get load matrix as vector
            # Select invidual free nodes DOF
            Fs = Fs[DOFbool]
            #
            #Fs2 = fernodes[self._plane.hforce].stack()
            # Select invidual free nodes DOF 
            #Fs2 = Fs2[DOFbool]
            #
            # Solve displacements U
            Us = iter(solver(Ks, Fs))
            #
            # reshape vector in matrix form [row, col]
            Us = [next(Us) if ieqnum != 0 else ieqnum
                  for ieqnum in jbcflat]
            # bak to matrix form
            Us = to_matrix(Us, self._plane.ndof)
            # pack basic load data for df's dumping 
            Utemp.extend([[*key,  item['load_title'].iloc[0],
                           nname, *Us[x]]
                           for x, nname in enumerate(jbc.index)])
        #
        #
        Us = self.df(Utemp)
        #
        # add displacements from beam solution
        #fergrp = Fn.groupby(colgrp)
        # group Us
        Usgrp = Us.groupby(colgrp)
        cols = self._plane.hdisp
        Utemp = []
        for key, item in fergrp:
            #Usitem = Usgrp.get_group(key)
            Usitem = Usgrp.get_group(key).set_index(['node_name'])
            fernodes = item.set_index(['node_name'])
            # Select load nodes with free DOF
            DOFindex = DOFzeros.index.isin(item['node_name'])            
            #
            #Usitem.loc[DOFindex][cols] = Usitem.loc[DOFindex][cols].sub(fernodes[cols])
            Usitem.loc[item['node_name']][cols] = Usitem.loc[item['node_name']][cols].sub(fernodes[cols])
            #
            Utemp.append(Usitem.reset_index())
        # Update U results
        Us = db.concat(Utemp, ignore_index=True)        
        #
        return Us
    #
    #
    def PMT(self, Ks, Fn,
            jbc, solver,
            maxstiff: float):
        """
        Penalty Method
        
        Ks : Global stiffness matrix
        Fn : Global force & displacement matrix
        jbc = Nodes with boundary
        
        Return:
        U : Global node displacement
        """
        # TODO : must be better ways to manipulate dfs
        #db = DBframework()
        # group FER node load
        colgrp = ['load_name', 'load_id', 
                  'load_level','load_system',
                  'component_name']        
        #
        ndof = self._plane.ndof
        hforce =  self._plane.hforce
        hdisp = self._plane.hdisp
        colrename = self._plane.colrename
        #
        Ktrans = maxstiff * 10_000
        Krot = maxstiff * 100_000
        #
        kb = np.zeros((ndof, ndof), dtype=np.float64)
        #
        ditem = {'x':Ktrans, 'y':Ktrans, 'z':Ktrans,
                 'rx': Krot, 'ry': Krot, 'rz': Krot,}
        #
        hstiff = {key : ditem[key] for key in hdisp}
        kp = np.array([ditem[key] for key in hdisp])
        #
        #
        # Select nodes' free DOF
        Fbool, Fzeros = self._mask_DOF(jbc)
        Dbool, Dzeros = self._mask_DOF(jbc, rename=False)
        # jbc Matrix to vector 
        jbcflat = jbc.stack()
        # group FER node load
        fergrp = Fn.groupby(colgrp)
        #
        Utemp = []
        for key, item in fergrp:
            #Kitem = Ks.tolil(copy=True)
            Kitem = Ks.copy()
            # set FER df by nodes
            fernodes = item.set_index(['node_name'])
            # Select load nodes with free DOF
            Findex = Fzeros.index.isin(item['node_name'])
            #
            # Force Section
            #
            # Select load vector comprising nodes with free DOF 
            Fs = Fzeros.copy() #.astype('float64')
            #Fs = Fs[Findex]
            #Fs = Fs.add(fernodes[hforce].astype('float64'))
            #
            #Fs = fernodes[hforce].astype('float64').copy()
            #
            #Fs.loc[Findex] = fernodes[hforce].astype('float64')
            Fs.loc[Findex] = fernodes[hforce].astype('float64')
            #Fs = Fs.stack() # get load matrix as vector
            # Select invidual free nodes DOF
            #Fs = Fs[Fbool]
            #
            #
            Ds = Dzeros.copy()
            Dindex = Dzeros.index.isin(item['node_name'])
            Ds.loc[Dindex] = fernodes[hdisp].astype('float64')
            #
            #
            # update Ks + Kg if axial load
            #1 / 0
            #
            # Displacement Section
            if Ds.any(axis=None):
            #try:
                #1 / len(Dzeros)
                #Ds = Dzeros.copy()
                #Dindex = Dzeros.index.isin(item['node_name'])
                #Ds.loc[Dindex] = fernodes[hdisp].astype('float64')
                #print(f'displacement all: {Ds.all(axis=None)} empty :{Ds.empty}')
                Ds = Ds.mul(hstiff)
                Ds.rename(columns=colrename, inplace=True)
                # Update gloabl load vextor with displacement load
                #Fs = Fs.add(Ds, axis='columns')
                #print(f'---> {Ds} {Fs}')
                #
                for nodeid, idof in zip(item['node_name'], item['node_index']):
                    # check if node with displacement
                    try:
                        nload = Ds.loc[nodeid] #.to_numpy()
                        if not nload.any():
                            continue
                        Fs.loc[nodeid] = Fs.loc[nodeid].add(nload, axis='rows')
                        #Ditem = [nodeid + x
                        #         for x, step in enumerate(nload.tolist())
                        #         if step != 0]
                    except KeyError:
                        continue
                    # DOF
                    niqi = idof * ndof
                    niqj = niqi + ndof                    
                    #
                    kc = kp.copy()
                    kc[nload == 0] = 0
                    ksup = kb.copy()
                    ksup.flat[0::ndof+1] = kc
                    # update global K matrix with dummy stiffness
                    Kitem[niqi:niqj, niqi:niqj] += ksup
                    #1 / 0
                    #Ddof =  list(set(Ddof + Ditem))
                    #print('--->', niqi, niqj)
            #except ZeroDivisionError:
            #    pass
            #
            # FIXME: Matrix condensed
            K11, K12 = self._Msplit(Kitem, jbcflat)
            #Kitem1 = self._Mcond(Kitem, jbcflat)
            #
            # get load matrix as vector
            #Fs = Fs[Findex].stack()
            Fs = Fs.stack()
            Fs = Fs.loc[Fbool]
            #
            # Solve displacements U
            Us = self._Us(K11, Fs, solver,
                          jbcflat, ndof)
            #Us = iter(solver(K11, Fs))
            ##
            ## reshape vector in matrix form [row, col]
            #Us = [next(Us) if ieqnum != 0 else ieqnum
            #      for ieqnum in jbcflat]
            # bak to matrix form
            #Us = to_matrix(Us, ndof)
            # pack basic load data for df's dumping 
            Utemp.extend([[*key,  item['load_title'].iloc[0],
                           nname, *Us[x]]
                           for x, nname in enumerate(jbc.index)])
        #
        #
        Us = self.df(Utemp)
        return Us
    #
    def _Mcond(self, Kitem: list, jbcflat: list):
        """ """
        jbcc = jbcflat.values
        index = list(reversed([i for i, item in enumerate(jbcc)
                               if item == 0]))
        for i in index:
            Kitem = remove_column_row(Kitem, i, i)
        #
        return Kitem
    #
    def _Msplit(self, Km, jbcflat: list):
        """
        Km: The unpartitioned matrix (or vector) to be partitioned.
        D1: A list of the indices for degrees of freedom that have unknown displacements.
        D2: A list of the indices for degrees of freedom that have known displacements.
        """
        jbcc = jbcflat.values
        D2 = [i for i, item in enumerate(jbcc)
              if item == 0]
        
        D1 = [i for i, item in enumerate(jbcc)
              if item != 0]        
        # 1D vectors
        if Km.shape[1] == 1:
            # Partition the vector into 2 subvectors
            m1 = Km[D1, :]
            m2 = Km[D2, :]
            return m1, m2
        # 2D matrices
        else:
            # Partition the matrix into 4 submatrices
            m11 = Km[D1, :][:, D1]
            m12 = Km[D1, :][:, D2]
            #m21 = Km[D2, :][:, D1]
            #m22 = Km[D2, :][:, D2]
            return m11, m12 #, m21, m22
    #
    #
    def _Us(self, Ks, Fs, solver, 
            jbcflat: list, ndof: float):
        """ """
        # Solve displacements U
        Us = iter(solver(Ks, Fs))
        # reshape vector in matrix form [row, col]
        Us = [next(Us) if ieqnum != 0 else ieqnum
              for ieqnum in jbcflat]
        # bak to matrix form
        Us = to_matrix(Us, ndof)
        return Us
    #
    #
    #@property
    def df(self, dftemp):
        """displacement dataframe"""
        db = DBframework()
        header = ['load_name', 'load_id', 'load_level',
                  'load_system', 'component_name', 'load_title', 
                  'node_name', *self._plane.hdisp]
        return db.DataFrame(data=dftemp, columns=header, index=None)
    #
    #
    def _mask_DOF(self, jbc, rename: bool = True):
        """ """
        # remove rows with zeros
        dfjbc = jbc.copy()
        if rename :
            dfjbc.rename(columns=self._plane.colrename, inplace=True)
        dfjbc = dfjbc[jbc.any(axis=1)]
        dfjbc = dfjbc.replace(float(0.0), np.nan)
        #
        dfbool = dfjbc.copy()
        dfbool = dfbool.notnull()
        dfbool = dfbool.stack(future_stack=True)
        # Copy dataframe 
        #dfzeros = dfjbc.copy()
        #dfzeros.iloc[:] = float(0.0)
        #
        dfjbc.iloc[:] = float(0.0)
        #for col in dfjbc.columns:
        #    dfjbc[col].values[:] = float(0.0)
        #
        return dfbool.astype('bool'), dfjbc.astype('float64')
    #
    # -----------------------------------------------------------
    # Post-utils
    # -----------------------------------------------------------
    #    
#
# ---------------------------------------------
#
@dataclass
class StaticSolver:
    """ Linear static solver class"""
    __slots__ = ['_mesh',  '_method', '_postprocess',
                 '_log', '_result_name', 'db_file',
                 'second_order', 'ineleastic']
    
    def __init__(self, mesh,
                 result_name: int|str,
                 db_file: str, log: bool,
                 second_order: bool = False,
                 ineleastic: bool = False,
                 ) -> None:
        """
        plane : Plane system (3D/2D)
        """
        self._mesh = mesh
        self._result_name = result_name
        self.db_file = db_file
        self._postprocess = PostProcess(mesh=self._mesh,
                                        result_name=result_name,
                                        db_file=db_file)
        self._log = log
        self.second_order = second_order
        self.ineleastic = ineleastic
    #
    def Ke(self, sparse: bool):
        """Stiffness matrix"""
        return self._mesh.Ke(sparse=sparse)
    #
    # ------------------------------------------
    #
    def FD_basic(self):
        """
        Fn & Dn Nodal Basic Load

        Return
        Fi : Nodal force
        Di : Nodal displacement
        """
        FDn = self._mesh._load._basic.NFD_global(plane=self._mesh._plane)
        #
        if len(FDn) == 0:
            raise IOError('Load data is required')
        #
        col_grp = ['load_name', 'load_id',
                   'load_level', 'load_title',
                   'load_system', 'mesh_name',
                   'node_name', 'node_index']
        #
        head_force =  self._mesh._plane.hforce
        head_disp = self._mesh._plane.hdisp
        #
        Fi = FDn[col_grp + head_force]
        Di = FDn[col_grp + head_disp]
        return Fi, Di
    #
    def FD_combination(self):
        """
        Fn & Dn Nodal Combination Load

        Return
        Fi : Nodal force
        Di : Nodal displacement
        """
        FDn = self._mesh._load._combination.NFD_global(plane=self._mesh._plane)
        col_grp = ['load_name', 'load_id',
                   'load_level', 'load_title', 'load_system',
                   'mesh_name', 'node_name', 'node_index']
        #
        head_force =  self._mesh._plane.hforce
        head_disp = self._mesh._plane.hdisp
        #
        Fi = FDn[col_grp + head_force]
        Di = FDn[col_grp + head_disp]
        return Fi, Di
    #
    # ------------------------------------------
    #
    def _Msplit(self, Km, jbcflat: list):
        """
        Km: The un-partitioned matrix (or vector) to be partitioned.
        D1: A list of the indices for degrees of freedom that have unknown displacements.
        D2: A list of the indices for degrees of freedom that have known displacements.
        """
        jbcc = jbcflat.values
        D2 = [i for i, item in enumerate(jbcc)
              if item == 0]

        D1 = [i for i, item in enumerate(jbcc)
              if item != 0]
        # 1D vectors
        if Km.shape[1] == 1:
            # Partition the vector into 2 subvectors
            m1 = Km[D1, :]
            m2 = Km[D2, :]
            return m1, m2
        # 2D matrices
        else:
            # Partition the matrix into 4 submatrices
            m11 = Km[D1, :][:, D1]
            m12 = Km[D1, :][:, D2]
            # m21 = Km[D2, :][:, D1]
            # m22 = Km[D2, :][:, D2]
            return m11, m12  # , m21, m22
    #
    #
    def _Un_update(self, Un, lcomb,
                   values:list[str]): #  = ['x', 'y', 'z', 'rx', 'ry', 'rz']
        """
        Update node displacements to include lcomb

        Un : node displacement
        lcomb : load combination
        values:
        """
        db = DBframework()
        # group basic load by name
        ndgrp = Un.groupby('load_name')
        # get combinations with basic loads
        #
        combgrp = lcomb.groupby('load_name')
        for key, combfactors in combgrp:
            for row in combfactors.itertuples():
                comb = ndgrp.get_group(row.basic_load).copy()
                comb.loc[:, values] *= row.factor
                comb['load_level'] = 'combination'
                comb['load_name'] = row.load_name
                comb['load_id'] = row.load_id
                comb['load_title'] = row.load_title
                comb['mesh_name'] = row.mesh_name
                #comb['result_id'] = row.result_id
                #
                try:
                    dftemp = db.concat([dftemp, comb], ignore_index=True)
                except UnboundLocalError:
                    dftemp = comb
        #
        try:
            dftemp = dftemp.groupby(['load_name', 'load_id','load_level',
                                     'load_title', 'load_system',
                                     'mesh_name', 'result_name',
                                     'node_name'],
                                      as_index=False)[values].sum()
            #test
            Un = db.concat([Un, dftemp], ignore_index=True)
        except UnboundLocalError:
            pass
        #
        return Un
    #
    # ------------------------------------------
    #
    def solve_Linear(self,
                     sparse: bool = False,
                     max_iter: int | None = None):
        """
        Linear Static Analysis

        Input:
        Ke  : Global stiffness matrix
        Fn  : FER Node load dataframe
        jbc : Node boundary condition dataframe

        Return:
        Udf : Node displacement global system dataframe
        """
        start_time = time.time()
        #
        Ke = self._mesh.Ke(sparse=sparse)
        jbc = self._mesh.jbc()
        # Ke, Kfactor, solver = self._get_solver(sparse)
        #
        # Get global nodal load and displacement
        Fn, Dn = self.FD_basic()
        #
        Un = self.PMT(Ke, Fn, Dn,
                      jbc=jbc, sparse=sparse)
        #
        # Post-processing load combinations
        #
        load_comb = self._mesh._load.combination()
        df_comb = load_comb.to_basic()
        #
        # Update load comb displacements
        Un = self._Un_update(Un=Un,
                             lcomb=df_comb,
                             values=self._mesh._plane.hdisp)
        #
        uptime = time.time() - start_time
        print(f"** {{F}} = [Ke] {{U}} Solution: {uptime:1.4e} sec")
        return Un
    #
    def solve_Pdelta(self,
                    sparse: bool=False, 
                    max_iter: int = 30):
        """
        Linear Static Analysis (2nd Order)
        Approximate method: Two cycles iterative method.
        
        Input: 
        Ks  : Global stiffness matrix
        Fn  : FER Node load dataframe
        jbc : Node boundary condition dataframe
        
        Return: 
        Udf : Node displacement global system dataframe
        """
        # start process
        start_time = time.time()
        # ------------------------------
        # Get basic data
        Ke = self._mesh.Ke(sparse=sparse)
        jbc = self._mesh.jbc()
        # ------------------------------
        # Get global nodal load
        # and displacement (combination)
        Fn, Dn = self.FD_combination()
        # ------------------------------
        # Step 1 : linear solution
        Un = self.PMT(Ke, Fn, Dn,
                      jbc=jbc, sparse=sparse)
        #
        # ------------------------------
        # Step 2 : Pdelta solution
        # ------------------------------            
        #
        col_grp = ['load_name', 'load_id', 'load_level',
                   'load_system', 'mesh_name']
        # grouping df
        Un_grp = Un.groupby(col_grp)
        Fn_grp = Fn.groupby(col_grp)
        Dn_grp = Dn.groupby(col_grp)
        #
        Un_temp = []
        for key, Un_step in Un_grp:
            Fn_step = Fn_grp.get_group(key)
            Dn_step = Dn_grp.get_group(key)
            # ------------------------------
            # Assembly matrix including
            # node displacement to force
            kl = self._mesh.Kt(Dn=Un_step)
            # ------------------------------
            Un_temp.append(self.PMT(Ke=kl,
                                    Fn=Fn_step,
                                    Dn=Dn_step,
                                    jbc=jbc,
                                    sparse=sparse))
        #
        # Results to dataframe
        db = DBframework()
        Us = db.concat(Un_temp, ignore_index=True)
        # end process
        uptime = time.time() - start_time
        print(f"** {{F}} = [Kt] {{U}} Solution: {uptime:1.4e} sec")
        return Us
    #
    # ------------------------------------------
    # TODO : remove redundant code
    #
    def solveLinearX(self, Ke, Kg, Kt,
                     Fn, Dn, jbc,
                     sparse: bool = False,
                     max_iter: int | None = None):
        """
        Linear Static Analysis

        Input:
        Ke  : Global stiffness matrix
        Fn  : FER Node load dataframe
        jbc : Node boundary condition dataframe

        Return:
        Udf : Node displacement global system dataframe
        """
        #
        order = "1st"
        print(f"** Solving U = K^-1 F [{order} order] ")
        # start_time = time.time()
        #
        #
        (Ke, Kfactor,
         jbc, solver) = self._get_solver(sparse)
        #
        Us = self.PMT(Ke, Fn, Dn, jbc,
                      solver, Kfactor)
        #
        return Us

    #
    def solvePdeltaX(self,
                    Ke, Kg, Kt,
                    Fn, Dn, jbc,
                    sparse: bool=False, 
                    max_iter: int = 30):
        """
        Linear Static Analysis (2nd Order)
        Aproximate method: Two cycles iterative method.
        
        Input: 
        Ks  : Global stiffness matrix
        Fn  : FER Node load dataframe
        jbc : Node boundary condition dataframe
        
        Return: 
        Udf : Node displacement global system dataframe
        """
        #
        order = "2nd"
        print(f"** Solving U = K^-1 F [{order} order] ")
        #start_time = time.time()
        #
        # ------------------------------
        #
        (Ke, Kfactor,
         jbc, solver) = self._get_solver(Ke, jbc, sparse)
        #
        # ------------------------------
        # Step 1 
        Un = self.PMT(Ke, Fn, Dn, jbc,
                      solver, Kfactor)
        #
        # ------------------------------
        # Pdelta solution
        # ------------------------------            
        #
        colgrp = ['load_name', 'load_id', 'load_level',
                  'load_system', 'mesh_name',]
        #
        Ugrp = Un.groupby(colgrp)
        Dgrp = Dn.groupby(colgrp)
        #
        hdisp = ['node_name', *self._plane.hdisp]
        #
        Utemp = []
        for key, noded in Ugrp:
            #
            # ------------------------------
            # TODO: select basic only
            #if key[2] != 'basic':
            #    continue
            #
            # ------------------------------
            #
            Dstep = Dgrp.get_group(key)          
            #            
            # ------------------------------
            #
            Uii = noded[hdisp].set_index('node_name')
            #
            # ------------------------------
            #
            # Assembly matrix start
            #kg = Kg(D=Uii)
            #kt = Ke + kg
            #
            kl = Kt(D=Uii)
            kl = kl.tolil()
            #
            # ------------------------------
            #
            #Ui = self.PMT(Ke=kt, F=Fn,
            #              D=Dstep,
            #              jbc=jbc,
            #              solver=solver,
            #              maxstiff=Kfactor)
            #
            Uni = self.PMT(Ke=kl, F=Fn,
                           D=Dstep,
                           jbc=jbc,
                           solver=solver,
                           maxstiff=Kfactor)
            # ------------------------------
            Utemp.extend(Uni)
        #
        Us = self.df(Utemp)
        return Us
    #
    def solvePdeltaXX(self, elements,
                    Ke, Fn, Un, jbc, 
                    sparse: bool=False, 
                    max_iter: int = 30):
        """
        Linear Static Analysis (2nd Order)
        
        Input: 
        Ks  : Global stiffness matrix
        Fn  : FER Node load dataframe
        jbc : Node boundary condition dataframe
        
        Return: 
        Udf : Node displacement global system dataframe
        """
        #
        order = "2nd"
        print(f"** Solving U = K^-1 F [{order} order] ")
        #start_time = time.time()
        #
        Ke, Kfactor, solver = self._get_solver(Ke, sparse)
        #
        colgrp = ['load_name', 'load_id', 
                  'load_level', 'load_title',
                  'load_system', 'mesh_name',
                  'node_name', 'node_index']
        #
        hforce =  self._plane.hforce
        hdisp = self._plane.hdisp
        #
        F = Fn[colgrp+hforce]
        D = Fn[colgrp+hdisp]
        #
        #
        colgrp = ['load_name', 'mesh_name',
                  'load_level',  'load_system']      
        #
        #
        # group FER node load
        Fgrp = F.groupby(colgrp)
        Dgrp = D.groupby(colgrp)
        Ugrp = Un.groupby(colgrp)
        #
        Utemp = []
        for key, item in Fgrp:
            # set FER df by nodes
            Fitem = item.set_index(['node_name'])
            #Fitem
            ndisp = Ugrp.get_group(key).set_index(['node_name'])
            ndisp
            for mname, element in elements.items():
                nodes = element.connectivity
                Tb = element.T
                # ---------------------------------------------
                # displacement
                nd_global = np.concatenate((ndisp.loc[nodes[0]],
                                            ndisp.loc[nodes[1]]), axis=None)
                #
                # ---------------------------------------------
                # convert global end-node disp in beam's local system
                # [x,y,z,rx,ry,rz]
                nd_local = Tb @ nd_global
                #
                # --------------------------------------------
                # get kg matrix            
    #    
    # ------------------------------------------
    #
    def _get_solver(self, Ke, sparse: bool):
        """ """

        # jbc = self._mesh.jbc()
        # Ke = self._mesh.Ke(sparse=sparse)

        if sparse:
            solver = solver_sparse
            # Ke = Ke.tolil()
            Kfactor = np.max(Ke.data.max())
        else:
            solver = solver_np
            Kfactor = Ke.max()

        return Kfactor, solver
    #
    def _Us(self, Ke, Fs, solver,
            jbc: list, ndof: float):
        """ """
        # jbc Matrix to vector
        jbcflat = jbc.stack()
        # FIXME: Matrix condensed
        K11, K12 = self._Msplit(Ke, jbcflat)
        # Solve displacements U
        #try:
        Us = iter(solver(K11, Fs))
        #except Warning:
        #    print('-->')
        # reshape vector in matrix form [row, col]
        Us = [next(Us) if ieqnum != 0 else ieqnum
              for ieqnum in jbcflat]
        # bak to matrix form
        Us = to_matrix(Us, ndof)
        return Us
    #
    def _mask_DOF(self, jbc, rename: bool = True):
        """ """
        # remove rows with zeros
        dfjbc = jbc.copy()
        if rename:
            dfjbc.rename(columns=self._mesh._plane.colrename,
                         inplace=True)
        dfjbc = dfjbc[jbc.any(axis=1)]
        dfjbc = dfjbc.replace(float(0.0), np.nan)
        #
        dfbool = dfjbc.copy()
        dfbool = dfbool.notnull()
        dfbool = dfbool.stack(future_stack=True)
        # Copy dataframe
        # dfzeros = dfjbc.copy()
        # dfzeros.iloc[:] = float(0.0)
        #
        dfjbc.iloc[:] = float(0.0)
        # for col in dfjbc.columns:
        #    dfjbc[col].values[:] = float(0.0)
        #
        return dfbool.astype('bool'), dfjbc.astype('float64')
    #
    def _to_df(self, Us):
        """displacement dataframe"""
        db = DBframework()
        header = ['load_name', 'load_id', 'load_level',
                  'load_system', 'mesh_name', 'result_name',
                  'load_title', 'node_name', *self._mesh._plane.hdisp]
        dftemp = db.DataFrame(data=Us, columns=header, index=None)
        #dftemp.rename(columns={'load_system':'system'}, inplace=True)
        #header = ['load_name', 'load_id', 'load_level',
        #          'load_title', 'result_name',
        #          'node_name', 'system',  *self._mesh._plane.hdisp]
        return dftemp #[header]
    #
    # ------------------------------------------
    #
    def PMT(self, Ke, Fn, Dn,
            jbc, sparse: bool=False):
        """
        Penalty Method
        
        Ke : Global stiffness matrix
        Fn : Global force  df
        Dn : Global displacement df
        jbc = Nodes with boundary
        
        Return:
        U : Global node displacement
        """
        Kfactor, solver = self._get_solver(Ke, sparse)
        #
        # TODO : must be better ways to manipulate dfs
        # grouping node items
        colgrp = ['load_name', 'load_id', 
                  'load_level','load_system',
                  'mesh_name']
        #
        ndof = self._mesh._plane.ndof
        head_force = ['node_name', *self._mesh._plane.hforce]
        head_disp = ['node_name', *self._mesh._plane.hdisp]
        #
        # Select nodes' free DOF
        Fn_bool, Fn_zeros = self._mask_DOF(jbc)
        Dn_bool, Dn_zeros = self._mask_DOF(jbc, rename=False)
        # group global node load and displacement
        Fn_grp = Fn.groupby(colgrp)
        Dn_grp = Dn.groupby(colgrp)
        #
        Utemp = []
        for key, Fn_item in Fn_grp:
            # set df by nodes
            Fn_set = Fn_item[head_force].set_index(['node_name'])
            Dn_set = Dn_grp.get_group(key) # .set_index(['node_name'])
            Dn_set = Dn_set[head_disp].set_index(['node_name'])
            # Select load nodes with free DOF
            Fn_index = Fn_zeros.index.isin(Fn_item['node_name'])
            # Select load vector comprising nodes with free DOF 
            Fs = Fn_zeros.copy()
            Fs.loc[Fn_index] = Fn_set.astype('float64')
            #
            Kitem = Ke.copy()
            # Displacement Section
            if Dn_set.any(axis=None):
                Kitem, Fs = self._PTM_Dn(Kitem=Kitem, Fs=Fs,
                                         Fn_item=Fn_item,
                                         Dn_set=Dn_set,
                                         Dn_zeros=Dn_zeros,
                                         ndof=ndof,
                                         K_factor=Kfactor)
            #
            # get load matrix as vector
            Fs = Fs.stack()
            Fs = Fs.loc[Fn_bool]
            # Solve displacements U
            Us = self._Us(Kitem, Fs, solver, jbc, ndof)
            # pack basic load data for df's dumping 
            Utemp.extend([[*key,  self._result_name,
                           Fn_item['load_title'].iloc[0],
                           nname, *Us[x]]
                           for x, nname in enumerate(jbc.index)])
        #
        Us = self._to_df(Utemp)
        return Us
    #
    def _PTM_Dn(self, Kitem, Fs,
                Fn_item, Dn_set, Dn_zeros,
                ndof:int, K_factor:float):
        """
        """
        # set up
        Dn_index = Dn_zeros.index.isin(Fn_item['node_name'])
        Ds = Dn_zeros.copy()
        Ds.loc[Dn_index] = Dn_set.astype('float64')
        #
        Ke_sup = np.zeros((ndof, ndof), dtype=np.float64)
        #
        K_trans = K_factor * 10_000
        K_rot = K_factor * 100_000
        ditem = {'x':K_trans, 'y':K_trans, 'z':K_trans,
                 'rx': K_rot, 'ry': K_rot, 'rz': K_rot,}
        #
        col_rename = self._mesh._plane.colrename
        head_disp = self._mesh._plane.hdisp
        head_stiff = {key : ditem[key] for key in head_disp}
        kp = np.array([ditem[key] for key in head_disp])
        #
        # Displacement check and postprocess
        #
        # calculate force {F} = {D} x Kmod
        FDs = Ds.mul(head_stiff)
        FDs.rename(columns=col_rename, inplace=True)
        # Update global load vector with displacement load
        # print(f'---> {Ds} {Fs}')
        for nodeid, idof in zip(Fn_item['node_name'], Fn_item['node_index']):
            # check if node with displacement
            try:
                nload = FDs.loc[nodeid]  # .to_numpy()
                if not nload.any():
                    continue
                Fs.loc[nodeid] = Fs.loc[nodeid].add(nload, axis='rows')
            except KeyError:
                continue
            #
            # DOF
            niqi = idof * ndof
            niqj = niqi + ndof
            #
            kc = kp.copy()
            kc[nload == 0] = 0
            #
            Ke_sup.flat[0::ndof + 1] = kc
            # update global K matrix with dummy stiffness
            Kitem[niqi:niqj, niqi:niqj] += Ke_sup
            # print('--->', niqi, niqj)
        #
        return Kitem, Fs
    #
    #
    # ------------------------------------------
    #
    def solve(self,
              sparse:bool =True,
              max_iter:int = 30,
              beam_steps:int = 10):
        """
        Solves the static system by the Direct Stiffness Method (DSM)

        method : banded, frontal
        """
        if self._mesh:
            if self.second_order:
                order = "2nd"
                self._postprocess._Pdelta = True
                run = self.solve_Pdelta
                # Fn = mesh._load._combination.Fn()
            else:
                order = "1st"
                self._postprocess._Pdelta = False
                run = self.solve_Linear
                # Fn = mesh._load._basic.Fn()
            #
            print(f"** Solving Linear Static [{order} order] ")
            #
            Un = run(sparse=sparse,
                     max_iter=max_iter)
            #
            self._postprocess.Un.df = Un
        else:
            raise IOError('** error: mesh missing')
        # run postprocessing
        self._postprocess.run(beam_steps=beam_steps)
        #
        return self._postprocess._results
#