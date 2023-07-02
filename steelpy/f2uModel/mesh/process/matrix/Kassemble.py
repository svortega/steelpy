# 
# Copyright (c) 2009-2023 steelpy
#
from __future__ import annotations
#
# Python stdlib imports
import time
#from multiprocessing import Process, Manager

# package imports
import numpy as np

#
#
#
# -------------------- 
#     assemble
# --------------------
#
#
# ---------------------------------------------
# Solution using numpy
#
#
def assemble_Kmatrix_np(elements, jbc, neq,
                        m2D: bool=False):
    """
    Asseable the element matrices
    -------------------------------------------------
    Ka : stiffness matrix
    jbc : nodes freedom
    m2D : 2D matrix
    """
    print("** Processing Global [K] Matrix")
    start_time = time.time()
    #
    #if m2D:
    #    jbc = jbc[['x', 'y', 'rz']]
    #
    nn = len(jbc)
    Ka = form_Kmatrix_np(elements, nn, m2D=m2D)
    #
    #jbcc = list(it.chain.from_iterable(jbc))
    jbcc = jbc.stack().values
    #
    index = list(reversed([i for i, item in enumerate(jbcc)
                           if item == 0]))
    #
    for i in index:
        Ka = remove_column_row(Ka, i, i) 
    #
    #
    #
    #steps = len(jbcc)
    #for i, bn in enumerate(jbcc):
    #    try:
    #        1/bn
    #    except ZeroDivisionError:
    #        aa[i][:steps] = 0  # Row to zeros
    #        aa[:steps][i] = 0  # Column to zeros
    #        aa[i][i] = 1       # Diagonal to 1
    #
    #from scipy.linalg import cholesky_banded
    #c = cholesky_banded(aa)
    #aa = np.linalg.cholesky(aa)
    #
    #with open("stfmx.f2u", "wb") as f:
    #    pickle.dump(jbc, f)
    #    pickle.dump(aa, f)
    #
    uptime = time.time() - start_time
    print("** [K] assembly Finish Process Time: {:1.4e} sec".format(uptime))
    return Ka
#
#
def remove_column_row(a: list, row: float, col: int):
    """ """
    without_row = np.delete(a, row, axis=0)
    return np.delete(without_row, col, axis=1)
#
#
def form_Kmatrix_np(elements, nn:int,
                    m2D:bool = False):
    """
    Global system stiffness matrix 
    
    elements : 
    nn  : node number
    m2D : Matrix 2D (False default --> 3D)
    :return
    Ka : global stiffness matrix
    """
    #nodal degrees of freedom per node 
    ndof:int = 6
    if m2D:
        ndof:int = 3
    Ka = np.zeros((nn*ndof, nn*ndof), dtype=np.float64)
    for key, element in elements.items():
        # TODO : check applicable to all element type
        keg = np.array(element.K(m2D=m2D))
        idof, jdof = element.DoF
        # node and corresponding dof (start, finish), used to define the
        # elements of the system stiffness and force matrices
        niqi, niqj = idof*ndof, idof*ndof + ndof
        njqi, njqj = jdof*ndof, jdof*ndof + ndof
        # assemble global stiffness matrix, quadrant 1 to 4
        Ka[niqi:niqj, niqi:niqj] += keg[:ndof, :ndof]             # 2nd
        Ka[niqi:niqj, njqi:njqj] += keg[:ndof, ndof:2*ndof]       # 1st
        Ka[njqi:njqj, niqi:niqj] += keg[ndof:2*ndof, :ndof]       # 3rd
        Ka[njqi:njqj, njqi:njqj] += keg[ndof:2*ndof, ndof:2*ndof] # 4th
    return Ka
#
#
def remove_column_of_zeros_and_shift_row(a, row, col):
    """ """
    without_row = np.delete(a, row, axis=0)
    without_row_and_col = np.delete(without_row, col, axis=1)
    z = np.zeros((1, len(without_row_and_col[0])))
    without_col_shifted_row = np.append(z, without_row_and_col, axis=0)
    return without_col_shifted_row
#
#
def swap_col(arr, start_index, last_index):
    """ """
    arr[:, [start_index, last_index]] = arr[:, [last_index, start_index]]
#
#
def swap_row(arr, start_index, last_index):
    """ """
    arr[[start_index, last_index],:] = arr[[last_index, start_index],:]
#
#
#
#
#
#
# ---------------------------------------------
#  Banded Matrix
# ---------------------------------------------
#
# --------------------
#
def UDUt(a: list[float]) -> list:
    """
    Solution of banded system asing UDU decomposition 
    based on Weaver & Gare pp 469-472.  
    """
    start_time = time.time()
    #neq, iband = np.shape(a)
    neq = len(a)
    iband = len(a[0])
    #print("** Processing UDU")
    if a[0, 0] == 0.0:
        raise RuntimeError("error:  zero diagonal term")
    if neq == 1:
        return a
    #
    for j in range(1, neq):
        j2 = max(j - iband + 1, 0)
        # off-diagonal terms
        for i in range(j2+1, j):
            a[i, j - i] -= np.sum([a[k, i - k] * a[k, j - k]
                                   for k in range(j2, i)])
            #a[i, j - i] -= np.sum(a[j2: i, i - j2: 0]
            #                      * a[j2: i, j - j2: j - i])
        # diagonal  terms
        for k in range(j2, j):
            temp = a[k, j - k] / a[k, 0]
            a[j, 0] -= a[k, j - k] * temp
            a[k, j - k] = temp
        # check diagonal zero terms
        try:
            1.0 / a[j, 0]
        except ZeroDivisionError:
            raise RuntimeError("error:  zero diagonal term")
    #
    uptime = time.time() - start_time
    print("** [UDU] Finish Process Time: {:1.4e} sec".format(uptime))
    return a
#
# --------------------
#
def assemble(ipv:list, a:list, aa:list):
    """
    ipv : member ends equations
    a  : member stiffness matrix
    aa : global stiffness matrix
    """
    for i in range(12):
        try:
            1.0 / (ieqn1 := ipv[i])
            for j in range(i, 12):
                try:
                    1.0 / (ieqn2 := ipv[j])
                    iband = np.minimum(ieqn1, ieqn2) - 1
                    jband = np.abs(ieqn1 - ieqn2)
                    aa[iband, jband] += a[i, j]
                except ZeroDivisionError:
                    continue
        except ZeroDivisionError:
            continue
#
# --------------------
#
def max_bandwidth(elements,  jbc):
    """
    calculate max bandwidth
    ------------------------  
    npi : connectivity end 1
    npj : connectivity end 2
    jbc : nodes freedom
    nel: number of elements
    if we
    npi ,npj, jbc, nel
    """
    #TODO : plates 
    ibndm4 = [0]
    for key, element in elements.items():        
        nodes = element.connectivity
        bc1 = jbc.loc[nodes[0]].tolist()
        bc2 = jbc.loc[nodes[1]].tolist()
        ieqn = bc1 + bc2
        try:
            ibndm4.append(max([abs(ieqn1 - ieqn2)
                               for x, ieqn1 in enumerate(ieqn) if ieqn1 > 0
                               for ieqn2 in ieqn[x+1:] if ieqn2 > 0]))
        except ValueError:
            continue
    #
    return max(ibndm4) + 1
#
def form_Kmatrix(elements, jbc:list, neq:int, iband:int):
    """
    elements
    jbc : node equations
    neq : number of equations
    iband : bandwidth

    :return
    a : global banded stiffness matrix (neq X iband)
    """
    aa = np.zeros((neq, iband))
    for key, element in elements.items():
        #idof, jdof = element.DoF
        nodes = element.connectivity
        # TODO : check applicable to all element type
        a = element.K()
        #ipv = jbc[idof] + jbc[jdof]
        bc1 = jbc.loc[nodes[0]].tolist()
        bc2 = jbc.loc[nodes[1]].tolist()        
        ipv = bc1 + bc2
        assemble(ipv, a, aa)
    #
    return aa
#
# ---------------------
#
def assemble_banded_Kmatrix(elements, jbc, neq, iband):
    """
    Asseable the element matrices in upper band form;
    call separatly from formstif, formmass, formgeom
    -------------------------------------------------
    Ka : stiffness matrix
    jbc : nodes freedom
    """
    print("** Processing Global [K] Banded Matrix")
    start_time = time.time()
    #
    Ka = form_Kmatrix(elements, jbc=jbc, neq=neq, iband=iband)
    # set UDUt matrix
    Ka = UDUt(Ka)
    #
    uptime = time.time() - start_time
    print("** [Kb] assembly Finish Process Time: {:1.4e} sec".format(uptime))
    return Ka
#