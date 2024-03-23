#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
from array import array
from copy import copy
from math import fsum
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
#
#
# --------------------
# solver pure python
# --------------------
# 
#
def BAK(a: list[float], b: list[float]) -> list:
    """
    back substitution
    """
    #neq, iband = np.shape(a)
    neq = len(a)
    iband = len(a[0])
    wk = copy(b)
    # forward substitution
    for i in range(neq):
        j = max(i - iband + 1, 0)
        wk[i] -= fsum([a[k][i - k] * wk[k]
                       for k in range(j, i)])
        #
        #wk[i] -= np.sum(a[j: i, i - j: 0] * wk[j: i])
    # middle terms
    wk = array('d', [wk[i] / a[i][0] for i in range(neq)])
    #wk = wk[:neq] / a[:neq, 0]
    #wkk = copy(wk)
    # backward substitution
    for i in range(neq - 1, -1, -1):
        j = min(i + iband, neq)
        wk[i] -= fsum([a[i][k - i] * wk[k]
                       for k in range(i+1, j)])
        #
        #wk[i] -= np.sum(a[i, 1: j - i] * wk[i+1: j])
    #
    return wk
#
#
def solver_Mbanded(stf, nloads):
    """ """
    nloads = nloads.values        
    #
    # get displacement in global system
    #x = cho_solve_banded(stf, nloads)
    ndisp = BAK(stf, nloads)
    #
    return ndisp    
#
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
#
#
# ---------------------------------------------
#
# ---------------------------------------------
#
@dataclass
class StaticSolver:
    """ Linear static solver class"""
    __slots__ = ['_plane',  '_method']
    
    def __init__(self, plane) -> None:
        """
        plane : Plane system (3D/2D)
        """
        self._plane = plane
    #
    #
    #
    #
    #
    def solve(self, Ks, Fn, jbc,
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
        # Select solver
        solver = solver_np
        #if self._method == 'banded':
        #    solver = solver_Mbanded
        # Solution f = Ku
        #Us = self.DSM(Kg, Fn, jbc, solver)
        Us = self.PMT(Ks, Fn, jbc, solver)
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
    def PMT(self, Ks, Fn, jbc, solver):
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
        maxstiff = Ks.max()
        Ktrans = maxstiff * 10
        Krot = maxstiff * 100
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
            Fs[Findex] = fernodes[hforce].astype('float64')
            #Fs = Fs.stack() # get load matrix as vector
            # Select invidual free nodes DOF
            #Fs = Fs[Fbool]
            #
            # update Ks + Kg if axial load
            #1 / 0
            #
            # Displacement Section
            try:
                1 / len(Dzeros)
                Ds = Dzeros.copy()
                Dindex = Dzeros.index.isin(item['node_name'])
                Ds.loc[Dindex] = fernodes[hdisp].astype('float64')
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
                            #print('-->')
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
            except ZeroDivisionError:
                pass
            #
            # FIXME: Matrix condensed
            #
            jbcc = jbcflat.values
            index = list(reversed([i for i, item in enumerate(jbcc)
                                   if item == 0]))
            for i in index:
                Kitem = remove_column_row(Kitem, i, i)            
            #
            # get load matrix as vector
            #Fs = Fs[Findex].stack()
            Fs = Fs.stack()
            Fs = Fs.loc[Fbool]
            #
            # Solve displacements U
            Us = iter(solver(Kitem, Fs))
            #
            # reshape vector in matrix form [row, col]
            Us = [next(Us) if ieqnum != 0 else ieqnum
                  for ieqnum in jbcflat]
            # bak to matrix form
            Us = to_matrix(Us, ndof)
            # pack basic load data for df's dumping 
            Utemp.extend([[*key,  item['load_title'].iloc[0],
                           nname, *Us[x]]
                           for x, nname in enumerate(jbc.index)])
        #
        #
        Us = self.df(Utemp)
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
#