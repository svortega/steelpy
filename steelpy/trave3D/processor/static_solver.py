#
# Copyright (c) 2009-2023 steelpy
#
from __future__ import annotations
# Python stdlib imports
from array import array
from copy import copy
from math import fsum
import pickle
#from typing import List #, NamedTuple, Union
from itertools import chain
import time
#
# package imports
from steelpy.process.math.operations import to_matrix #, zeros_vector, matAbd, trns_3Dv
#from steelpy.process.math.vector import Vector
#from steelpy.trave3D.processor.operations import get_deflection
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
    neq = len(a)
    iband = len(a[0])
    wk = copy(b)
    # forward substitution
    for i in range(neq):
        j = max(i - iband + 1, 0)
        wk[i] -= fsum([a[k][i - k] * wk[k]
                       for k in range(j, i)])
    # middle terms
    wk = array('d', [wk[i] / a[i][0] for i in range(neq)])
    #wkk = copy(wk)
    # backward substitution
    for i in range(neq - 1, -1, -1):
        j = min(i + iband, neq)
        wk[i] -= fsum([a[i][k - i] * wk[k]
                       for k in range(i+1, j)])
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
def solve_deflections(df_nload, method: str): 
    """
    """
    print("** Calculating Joint Displacements")
    print("** reloaded [k] & {p}")    
    #
    file = open("stfmx.f2u", "rb")
    df_jbc = pickle.load( file )
    stf = pickle.load( file )
    file.close()
    #
    solver = solver_np
    if method == 'banded':
        solver = solver_Mbanded
    #
    jbcc = df_jbc.stack()
    dfbool, dfzeros = get_mask(df_jbc)
    #dfbool = dfjbc.stack()
    #
    dfnload = update_load(df_nload)
    blgrp = dfnload.groupby(['load_name', 'load_number', 
                             'load_type','load_system'])
    #
    dftemp = []
    for key, litem in blgrp:
        # map loading 
        df1 = litem.set_index(litem['node_name'])
        df2 = dfzeros.copy()
        df2.loc[df2.index.isin(df1['node_name'])] = df1[['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']]
        # get load vector flatted
        df2 = df2.stack()
        nloads = df2[dfbool]
        #df2.mask(dfjbc, other=df2, inplace=True)
        #nloads = df2.stack()
        # Solve displcements
        ndisp = solver(stf, nloads)
        # reshape vector in matrix form [row, col]
        #xxx = df2.mask(dfjbc, other=df2)
        ndisp = [ndisp[ieqnum - 1] if ieqnum != 0 else ieqnum
                 for ieqnum in jbcc]
        #1/0
        ndisp = to_matrix(ndisp, 6)
        filldata = [[*key,  litem['load_title'].iloc[0],
                    nname, *ndisp[x]]
                    for x, nname in enumerate(df_jbc.index)]
        #
        dftemp.extend(filldata)
    #
    df_ndisp = get_dfdisp(dftemp)
    print("** Finished Calculating Joint Displacements")
    return df_ndisp, df_nload
#
#
def get_mask(df_jbc):
    """ """
    # remove rows with zeros
    dfjbc = df_jbc.rename(columns={'x':'Fx', 'y':'Fy', 'z':'Fz',
                           'rx':'Mx', 'ry':'My', 'rz':'Mz'}) # , inplace=True
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
def update_load(df_nload):
    """ """
    dfnload = (df_nload.groupby(['load_name', 'load_number', 'load_type',
                                 'load_title','load_system', 'node_name'])
             [['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']].sum())
    #
    dfnload.reset_index(inplace=True)
    return dfnload
#
def get_dfdisp(dftemp):
    """ """
    db = DBframework()
    header = ['load_name', 'load_number', 'load_type',
              'load_system', 'load_title',
              'node_name', 'x', 'y', 'z', 'rx', 'ry', 'rz']
    return db.DataFrame(data=dftemp, columns=header, index=None)
        