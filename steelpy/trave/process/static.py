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
#from steelpy.f2uModel.mesh.main import MeshPlane
from steelpy.utils.math.operations import remove_column_row
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
    # Post-process
    # -----------------------------------------------------------
    #    
#
# ---------------------------------------------
#
@dataclass
class StaticSolver:
    """ Linear static solver class"""
    __slots__ = ['_mesh',  '_method']
    
    def __init__(self, mesh) -> None:
        """
        plane : Plane system (3D/2D)
        """
        self._mesh = mesh
    #
    #
    def _get_solver(self, sparse: bool):
        """ """
        
        jbc = self._mesh.jbc()
        Ke = self._mesh.Ke(sparse=sparse)
        Kfactor = Ke.max()
        
        if sparse:
            solver = solver_sparse
            Ke = Ke.tolil()
        else:
            solver = solver_np
        
        return Ke, Kfactor, jbc, solver
    #
    def _basicload(self):
        """ """
        Fn = self._mesh._load._basic.Fn()
        #
        if len(Fn) == 0:
            raise IOError('Load Combination is required')
        #
        colgrp = ['load_name', 'load_id', 
                  'load_level', 'load_title',
                  'load_system', 'component_name',
                  'node_name', 'node_index']
        #
        hforce =  self._mesh._plane.hforce
        hdisp = self._mesh._plane.hdisp
        #
        Fi = Fn[colgrp+hforce]
        Di = Fn[colgrp+hdisp]
        #
        return Fi, Di
    #
    def solveLinear(self,
                    sparse: bool=False,
                    max_iter: int|None = None):
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
        (Ke, Kfactor,
         jbc, solver) = self._get_solver(sparse)
        #
        # Get basic nodal load and displacement
        Fn, Dn = self._basicload()
        #
        Un = self.PMT(Ke, Fn, Dn, jbc,
                      solver, Kfactor)
        #
        # Post processing load combinations
        # 
        load_comb = self._mesh._load.combination()
        df_comb = load_comb.to_basic()
        #
        # Update load comb displacements
        Un = self._update_ndf(dfnode=Un,
                              dfcomb=df_comb,
                              values=self._mesh._plane.hdisp)
        #
        uptime = time.time() - start_time
        print(f"** {{F}} = [Ke] {{U}} Solution: {uptime:1.4e} sec")
        return Un
    #
    def _update_ndf(self, dfnode, dfcomb, 
                   values:list[str]): #  = ['x', 'y', 'z', 'rx', 'ry', 'rz']
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
    def solveLinearX(self, Ke, Kg, Kt,
                    Fn, Dn, jbc,
                    sparse: bool=False,
                    max_iter: int|None = None):
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
        #start_time = time.time()
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
    #
    #
    def _combload(self):
        """ """
        Fn = self._mesh._load._combination.Fn()
        colgrp = ['load_name', 'load_id', 
                  'load_level', 'load_title',
                  'load_system', 'component_name',
                  'node_name', 'node_index']
        #
        hforce =  self._mesh._plane.hforce
        hdisp = self._mesh._plane.hdisp
        #
        Fi = Fn[colgrp+hforce]
        Di = Fn[colgrp+hdisp]
        #
        return Fi, Di    
    #
    def solvePdelta(self,
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
        start_time = time.time()
        #
        # ------------------------------
        # Get basic data
        (Ke, Kfactor,
         jbc, solver) = self._get_solver(sparse)
        #
        # ------------------------------
        # Get basic nodal load and displacement
        Fn, Dn = self._combload()
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
                  'load_system', 'component_name',]
        #
        Fgrp = Fn.groupby(colgrp)
        Ugrp = Un.groupby(colgrp)
        Dgrp = Dn.groupby(colgrp)
        #
        hdisp = ['node_name', *self._mesh._plane.hdisp]
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
            Fstep = Fgrp.get_group(key)  
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
            kl = self._mesh.Kt(D=Uii)
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
            Uni = self.PMT(Ke=kl,
                           F=Fstep,
                           D=Dstep,
                           jbc=jbc,
                           solver=solver,
                           maxstiff=Kfactor)
            # ------------------------------
            Utemp.append(Uni)
        #
        db = DBframework()
        Us = db.concat(Utemp, ignore_index=True)
        #Us = self.df(Utemp)
        uptime = time.time() - start_time
        print(f"** {{F}} = [Kt] {{U}} Solution: {uptime:1.4e} sec")
        return Us
    #    
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
                  'load_system', 'component_name',]
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
                  'load_system', 'component_name',
                  'node_name', 'node_index']
        #
        hforce =  self._plane.hforce
        hdisp = self._plane.hdisp
        #
        F = Fn[colgrp+hforce]
        D = Fn[colgrp+hdisp]
        #
        #
        colgrp = ['load_name', 'component_name',
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
    #
    #
    def PMT(self, Ke, F, D, 
            jbc, solver,
            maxstiff: float):
        """
        Penalty Method
        
        Ke : Global stiffness matrix
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
        ndof = self._mesh._plane.ndof
        hforce =  self._mesh._plane.hforce
        hdisp = self._mesh._plane.hdisp
        colrename = self._mesh._plane.colrename
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
        Fgrp = F.groupby(colgrp)
        Dgrp = D.groupby(colgrp)
        #
        Utemp = []
        for key, item in Fgrp:
            Kitem = Ke.copy()
            # set FER df by nodes
            Fn = item.set_index(['node_name'])
            # Select load nodes with free DOF
            Findex = Fzeros.index.isin(item['node_name'])
            #
            # Force Section
            #
            # Select load vector comprising nodes with free DOF 
            Fs = Fzeros.copy()
            Fs.loc[Findex] = Fn[hforce].astype('float64')
            #
            # Displacement Section
            #
            Dn = Dgrp.get_group(key).set_index(['node_name'])
            Dindex = Dzeros.index.isin(item['node_name'])
            #
            Ds = Dzeros.copy()
            Ds.loc[Dindex] = Dn[hdisp].astype('float64')
            #
            #
            # Displacement check and postprocess
            if Ds.any(axis=None):
                Ds = Ds.mul(hstiff)
                Ds.rename(columns=colrename, inplace=True)
                # Update gloabl load vextor with displacement load
                #print(f'---> {Ds} {Fs}')
                #
                for nodeid, idof in zip(item['node_name'], item['node_index']):
                    # check if node with displacement
                    try:
                        nload = Ds.loc[nodeid] #.to_numpy()
                        if not nload.any():
                            continue
                        Fs.loc[nodeid] = Fs.loc[nodeid].add(nload, axis='rows')
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
                    #print('--->', niqi, niqj)
            #
            # FIXME: Matrix condensed
            K11, K12 = self._Msplit(Kitem, jbcflat)
            #
            # get load matrix as vector
            Fs = Fs.stack()
            Fs = Fs.loc[Fbool]
            #
            # Solve displacements U
            Us = self._Us(K11, Fs, solver,
                          jbcflat, ndof)
            #
            # pack basic load data for df's dumping 
            Utemp.extend([[*key,  item['load_title'].iloc[0],
                           nname, *Us[x]]
                           for x, nname in enumerate(jbc.index)])
        #
        Us = self.df(Utemp)
        return Us
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
    def _mask_DOF(self, jbc, rename: bool = True):
        """ """
        # remove rows with zeros
        dfjbc = jbc.copy()
        if rename :
            dfjbc.rename(columns=self._mesh._plane.colrename,
                         inplace=True)
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
    #
    def _Us(self, Ks, Fs, solver, 
            jbcflat: list, ndof: float):
        """ """
        # Solve displacements U
        #try:
        Us = iter(solver(Ks, Fs))
        #except Warning:
        #    print('-->')
        # reshape vector in matrix form [row, col]
        Us = [next(Us) if ieqnum != 0 else ieqnum
              for ieqnum in jbcflat]
        # bak to matrix form
        Us = to_matrix(Us, ndof)
        return Us
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