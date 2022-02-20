#
# Copyright (c) 2009-2021 steelpy
#
# Python stdlib imports
from array import array
from copy import copy
from math import fsum
import pickle
from typing import List #, NamedTuple, Union
from itertools import chain
import time
#
# package imports
from steelpy.process.math.operations import to_matrix, zeros_vector, matAbd, trns_3Dv
from steelpy.process.math.vector import Vector
from steelpy.trave3D.processor.operations import get_deflection
#
#
# --------------------
# solver pure python
# --------------------
#
def UDUtx(a: List[float]) -> List:
    """
    Solution of banded system asing UDU decomposition 
    based on Weaver & Gare pp 469-472.  
    """
    neq = len(a)
    iband = len(a[0])
    print("** Processing UDU")
    if a[0][0] == 0.0:
        raise RuntimeError("error:  zero diagonal term")
    if neq == 1:
        return a
    #
    imult = 0
    for j in range(1, neq):
        j2 = max(j - iband + 1, 0)
        # off-diagonal terms
        for i in range(j2+1, j):
            a[i][j - i] -= fsum([a[k][i - k] * a[k][j - k]
                                 for k in range(j2, i)])
            imult += abs(i - j2)
        # diagonal  terms
        for k in range(j2, j):
            temp = a[k][j - k] / a[k][0]
            a[j][0] -= a[k][j - k] * temp
            a[k][j - k] = temp
            imult += 2
        # check diagonal zero terms
        try:
            1.0 / a[j][0]
        except ZeroDivisionError:
            raise RuntimeError("error:  zero diagonal term")
    #
    print("** Finished Processing UDU: mults & divs {:}".format(imult))
    return a
#
def BAKx(a: List[float], b: List[float]) -> List:
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
                       for k in range(i + 1, j)])
    #
    # for i in range(neq-1 , -1, -1):
    #    wkk[i] -= sum([a[i][k - i] * wkk[k]
    #                  for k in range(iband-1, i+1, -1)])
    return wk
#
# --------------------
#
def UDUt(a: List[float]) -> List:
    """
    Solution of banded system asing UDU decomposition 
    based on Weaver & Gare pp 469-472.  
    """
    start_time = time.time()
    neq = len(a)
    iband = len(a[0])
    #print("** Processing UDU")
    if a[0][0] == 0.0:
        raise RuntimeError("error:  zero diagonal term")
    if neq == 1:
        return a
    #
    for j in range(1, neq):
        j2 = max(j - iband + 1, 0)
        # off-diagonal terms
        for i in range(j2+1, j):
            a[i][j - i] -= fsum([a[k][i - k] * a[k][j - k]
                                 for k in range(j2, i)])
        # diagonal  terms
        for k in range(j2, j):
            temp = a[k][j - k] / a[k][0]
            a[j][0] -= a[k][j - k] * temp
            a[k][j - k] = temp
        # check diagonal zero terms
        try:
            1.0 / a[j][0]
        except ZeroDivisionError:
            raise RuntimeError("error:  zero diagonal term")
    #
    uptime = time.time() - start_time
    print("** [UDU] Finish Process Time: {:1.4e} sec".format(uptime))
    return a
#
def BAK(a: List[float], b: List[float]) -> List:
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
#
# ---------------------------------------------
# solver pure python
#
#
def solve_displacement(elements, nodes, boundaries,
                       basic_load, load_combination):
    """
    """
    print("** Calculating Joint Displacements")
    print("** reloaded [k] & {p}")    
    #
    file = open( "stfmx.f2u", "rb" )
    jbc = pickle.load( file )
    stf = pickle.load( file )
    file.close()
    #
    #from scipy.linalg import cholesky_banded, cho_solve_banded
    #
    jbcc = list(chain.from_iterable(jbc))
    #fmemb = {}
    memf_basic = {}
    basic_res = {}
    for key, item in basic_load.items():
        load_res = basic_load.get_basic_load(key, elements, nodes, boundaries)
        # member end load in local system
        memf_basic[item.title] = load_res.member_load
        # node load in global system
        nodal_load = load_res.nodal_load
        # reduce degree of fredooms
        nloads = [nodal_load[i][j] 
                  for i in range(len(jbc)) for j in range(6) 
                  if jbc[i][j] != 0]
        # get displacement in global system
        #x = cho_solve_banded(stf, nloads)
        basic = BAK(stf, nloads)
        # expand degree of fredooms
        basic = [basic[ieqnum - 1] if ieqnum != 0 else ieqnum
                 for ieqnum in jbcc]
        basic_res[item.title] = Vector(basic)
    #
    comb_res, memf_comb = load_combination.solve_combinations(basic_res, memf_basic)
    #
    with open("elemfout.f2u", "wb") as f:
        pickle.dump(memf_basic, f)
        pickle.dump(memf_comb, f)
    #
    #
    deflection = get_deflection(nodes, basic_res, comb_res)
    print("** Finished Calculating Joint Displacements")
    return deflection
#
#
# ---------------------------------------------
