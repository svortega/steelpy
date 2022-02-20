# 
# Copyright (c) 2009-2022 steelpy
#
# Python stdlib imports
import itertools as it
#from itertools import chain
#import math as math
import pickle
import time
from typing import Dict, List, Tuple # , ClassVar, Iterable, Union
#from multiprocessing import Process, Manager

# package imports
from steelpy.process.math.operations import zeros, to_matrix #, transposeM
from steelpy.trave3D.processor.operations import shape_cond
import numpy as np

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
def assemble_Kmatrix(elements, nodes, boundaries):
    """
    Asseable the element matrices
    -------------------------------------------------
    aa : stiffness matrix
    jbc : nodes freedom
    """
    print("** Processing Global [K] Matrix")
    start_time = time.time()
    aa, jbc = get_Kmatrix_np(elements, nodes, boundaries)
    #
    jbcc = list(it.chain.from_iterable(jbc))
    #
    index = list(reversed([i for i, item in enumerate(jbcc)
                           if item == 0]))
    #
    for i in index:
        aa = remove_column_row(aa, i, i) 
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
    with open("stfmx.f2u", "wb") as f:
        pickle.dump(jbc, f)
        pickle.dump(aa, f)
    #
    uptime = time.time() - start_time
    print("** [K] assembly Finish Process Time: {:1.4e} sec".format(uptime))    
    #
    #return aa
#
#
def get_Kmatrix_np(elements, nodes, boundaries):
    """ """
    free_nodes = elements.get_free_nodes
    jbc, neq = shape_cond(elements=elements, 
                           nodes=nodes, 
                           boundaries=boundaries, 
                           free_nodes=free_nodes)
    aa = form_Kmatrix_np(elements=elements, jbc=jbc)
    return aa, jbc
#
#
#
def remove_column_row(a, row, col):
    """ """
    without_row = np.delete(a, row, axis=0)
    return np.delete(without_row, col, axis=1)
#
#
def form_Kmatrix_np(elements, jbc:List):
    """
    elements
    jbc : node equations

    :return
    a : global stiffness matrix
    """
    ndof:int = 6    # nodal degrees of freedom
    # global system stiffness matrix
    nn = len(jbc)
    aa = np.zeros((nn*ndof, nn*ndof), dtype=np.float64)
    for key, element in elements.items():
        idof, jdof = element.DoF
        keg = np.array(element.Kglobal)
        # node and corresponding dof (start, finish), used to define the
        # elements of the system stiffness and force matrices
        niqi, niqj = idof*ndof, idof*ndof + ndof
        njqi, njqj = jdof*ndof, jdof*ndof + ndof
        # assemble global stiffness matrix, quadrant 1 to 4
        aa[niqi:niqj, niqi:niqj] += keg[:6, :6]         # 2nd
        aa[niqi:niqj, njqi:njqj] += keg[:6, 6:12]       # 1st
        aa[njqi:njqj, niqi:niqj] += keg[6:12, :6]       # 3rd
        aa[njqi:njqj, njqi:njqj] += keg[6:12, 6:12]     # 4th
    #
    return aa
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
# ---------------------------------------------