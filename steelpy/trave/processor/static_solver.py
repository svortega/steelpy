#
# Copyright (c) 2009-2023 steelpy
#
from __future__ import annotations
# Python stdlib imports
from array import array
from copy import copy
from math import fsum
import pickle
from dataclasses import dataclass
from typing import NamedTuple
#from itertools import chain
#import time
#
# package imports
from steelpy.process.math.operations import to_matrix
#from steelpy.trave.processor.operations import ElementProcess
from steelpy.process.dataframe.main import DBframework
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
def solver_np(stf, nloads):
    """ """
    #nloads = df2.stack() #.values
    ndisp = np.linalg.solve(stf, nloads)
    #if not np.allclose(np.dot(stf, ndisp), nload): #, rtol=1e-05, atol=1e-08
    #    raise RuntimeError('Solution fail')
    #
    return ndisp
#
#
# ---------------------------------------------
#
# ---------------------------------------------
#
#
def solve_deflections(df_nload, method: str,
                      stf=None, jbc=None, 
                      m2D:bool = False): 
    """
    m2D: 2D Marix (False default)
    """
    #
    print("** Calculating Joint Displacements")
    print("** reloaded [k] & {p}")    
    #
    if not stf:
        file = open("stfmx.f2u", "rb")
        df_jbc = pickle.load( file )
        stf = pickle.load( file )
        file.close()
    else:
        stf = stf
        df_jbc = jbc
    #
    #
    ndof:int = 6
    headforce = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
    headdisp = ['x', 'y', 'z', 'rx', 'ry', 'rz']
    if m2D:
        ndof:int = 3
        #df_jbc = df_jbc[['x', 'y', 'rz']]
        headforce = ['Fx', 'Fy', 'Mz']
        headdisp = ['x', 'y', 'rz']
    #
    solver = solver_np
    if method == 'banded':
        solver = solver_Mbanded
    #
    jbcc = df_jbc.stack()
    dfbool, dfzeros = get_mask(df_jbc, m2D=m2D)
    #dfbool = dfjbc.stack()
    #
    dfnload = update_load(df_nload, headforce)
    blgrp = dfnload.groupby(['load_name', 'load_number', 
                             'load_type','load_system'])
    #
    dftemp = []
    for key, litem in blgrp:
        # map loading 
        df1 = litem.set_index(litem['node_name'])
        df2 = dfzeros.copy()
        df2.loc[df2.index.isin(df1['node_name'])] = df1[headforce]
        # get load vector flatted
        df2 = df2.stack()
        nloads = df2[dfbool]
        #df2.mask(dfjbc, other=df2, inplace=True)
        #nloads = df2.stack()
        # Solve displcements
        ndisp = iter(solver(stf, nloads))
        # reshape vector in matrix form [row, col]
        #xxx = df2.mask(dfjbc, other=df2)
        #ndisp = [ndisp[ieqnum - 1] if ieqnum != 0 else ieqnum
        #         for ieqnum in jbcc]
        ndisp = [next(ndisp) if ieqnum != 0 else ieqnum
                 for ieqnum in jbcc]
        #1/0
        ndisp = to_matrix(ndisp, ndof)
        filldata = [[*key,  litem['load_title'].iloc[0],
                    nname, *ndisp[x]]
                    for x, nname in enumerate(df_jbc.index)]
        #
        dftemp.extend(filldata)
    #
    df_ndisp = get_dfdisp(dftemp, headdisp)
    print("** Finished Calculating Joint Displacements")
    #
    return df_ndisp #, df_nload
#
#
def get_mask(df_jbc, plane2D:bool = False):
    """ """
    cols =  {'x':'Fx', 'y':'Fy', 'z':'Fz',
             'rx':'Mx', 'ry':'My', 'rz':'Mz'}
    if plane2D:
        cols = {'x':'Fx', 'y':'Fy','rz':'Mz'}
    # remove rows with zeros
    dfjbc = df_jbc.rename(columns=cols) # , inplace=True
    dfjbc = dfjbc[df_jbc.any(axis=1)]
    dfjbc = dfjbc.replace(0, np.nan)
    dfjbc = dfjbc.notnull()
    #
    dfbool = dfjbc.stack()
    # Copy dataframe 
    #dfzeros = dfjbc.copy()
    #dfzeros.iloc[:] = 0
    #
    dfjbc.iloc[:] = 0
    #
    return dfbool, dfjbc
#
def update_load(df_nload, headforce):
    """ """
    dfnload = (df_nload.groupby(['load_name', 'load_number', 'load_type',
                                 'load_title','load_system', 'node_name'])
               [headforce].sum())
    #
    dfnload.reset_index(inplace=True)
    return dfnload
#
def get_dfdisp(dftemp, headdisp):
    """ """
    db = DBframework()
    header = ['load_name', 'load_number', 'load_type',
              'load_system', 'load_title',
              'node_name', *headdisp]
    return db.DataFrame(data=dftemp, columns=header, index=None)
#
#
#
# ---------------------------------------------
#
# ---------------------------------------------
#
@dataclass
class StaticSolver:
    __slots__ = ['_plane', '_load']
    
    def __init__(self, plane: NamedTuple) -> None:
        """
        """
        self._plane = plane
    #
    #
    def load(self, load):
        """ """
        self._load = load
    #
    #
    def Kglobal(self, jbc, Ka, Kg: list|bool = None):
        """ """
        with open("stfmx.f2u", "wb") as f:
            pickle.dump(jbc, f)
            pickle.dump(Ka, f)
    #
    def deflection(self, method: str):
        """ """
        #
        file = open("stfmx.f2u", "rb")
        df_jbc = pickle.load( file )
        stf = pickle.load( file )
        file.close()        
        #
        basic_load = self._load.basic()
        df_nload  = basic_load.node_df()        
        dfnload = self._load_update(df_nload)
        #
        solver = solver_np
        if method == 'banded':
            solver = solver_Mbanded
        #
        df_ndisp = self.solve(stf, dfnload, df_jbc, solver)
        return df_ndisp
    #
    def solve(self, stf, dfnload, df_jbc, solver):
        """ """
        dfbool, dfzeros = self._mask(df_jbc)
        jbcc = df_jbc.stack()
        blgrp = dfnload.groupby(['load_name', 'load_number', 
                                 'load_type','load_system'])
        #       
        dftemp = []
        for key, litem in blgrp:
            # map loading 
            df1 = litem.set_index(litem['node_name'])
            df2 = dfzeros.copy()
            df2.loc[df2.index.isin(df1['node_name'])] = df1[self._plane.hforce]
            # get load vector flatted
            df2 = df2.stack()
            nloads = df2[dfbool]
            # Solve displcements
            ndisp = iter(solver(stf, nloads))
            # reshape vector in matrix form [row, col]
            ndisp = [next(ndisp) if ieqnum != 0 else ieqnum
                     for ieqnum in jbcc]
            #
            ndisp = to_matrix(ndisp, self._plane.ndof)
            filldata = [[*key,  litem['load_title'].iloc[0],
                        nname, *ndisp[x]]
                        for x, nname in enumerate(df_jbc.index)]
            #
            dftemp.extend(filldata)
        #
        df_ndisp = self.df(dftemp)
        return df_ndisp
    #
    #@property
    def df(self, dftemp):
        """displacement dataframe"""
        db = DBframework()
        header = ['load_name', 'load_number', 'load_type',
                  'load_system', 'load_title',
                  'node_name', *self._plane.hdisp]
        return db.DataFrame(data=dftemp, columns=header, index=None)
    #
    def _load_update(self, df_nload):
        """ """
        dfnload = (df_nload.groupby(['load_name', 'load_number', 'load_type',
                                     'load_title','load_system', 'node_name'])
                   [self._plane.hforce].sum())
        #
        dfnload.reset_index(inplace=True)
        return dfnload
    #
    def _mask(self, df_jbc):
        """ """
        #cols =  {'x':'Fx', 'y':'Fy', 'z':'Fz',
        #         'rx':'Mx', 'ry':'My', 'rz':'Mz'}
        #if self._plane.m2D:
        #    cols = {'x':'Fx', 'y':'Fy','rz':'Mz'}
        # remove rows with zeros
        dfjbc = df_jbc.rename(columns=self._plane.colrename) # , inplace=True
        dfjbc = dfjbc[df_jbc.any(axis=1)]
        dfjbc = dfjbc.replace(0, np.nan)
        dfjbc = dfjbc.notnull()
        #
        dfbool = dfjbc.stack()
        # Copy dataframe 
        #dfzeros = dfjbc.copy()
        #dfzeros.iloc[:] = 0
        #
        dfjbc.iloc[:] = 0
        #
        return dfbool, dfjbc
    #
    # -----------------------------------------------------------
    # Post-process
    # -----------------------------------------------------------
    #    
#
#