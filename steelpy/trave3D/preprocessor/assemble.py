# 
# Copyright (c) 2009-2022 steelpy
#
# Python stdlib imports
import itertools as it
#from itertools import chain
import math as math
import pickle
import time
from typing import Dict, List, Tuple # , ClassVar, Iterable, Union
#from multiprocessing import Process, Manager

# package imports
from steelpy.process.math.operations import zeros, to_matrix #, transposeM
from steelpy.trave3D.processor.static_solver import UDUt
from steelpy.trave3D.processor.operations import get_bandwidth #, shape_cond

#
# -------------------- 
#     assemble
# --------------------
#
#
def assembleX(idof, jdof, jbc, a, aa):
    """
    """
    # set ipv to the positions in the array of the nodes
    ipv = [int(idof + i) for i in range(6)]
    ipv.extend([int(jdof + i) for i in range(6)])
    # store the values for individual array in global array
    for i in range(12):
        try:
            1.0 / (ieqn1 := jbc[ipv[i]])
            for j in range(i, 12):
                try:
                    1.0 / (ieqn2 := jbc[ipv[j]])
                    if ieqn1 > ieqn2:
                        jband = (ieqn1 - ieqn2)
                        aa[ieqn2-1][jband] += a[i][j]
                    else:
                        jband = (ieqn2 - ieqn1)
                        # TODO: check index = jband - 1
                        aa[ieqn1-1][jband] += a[i][j]
                except ZeroDivisionError:
                    continue
        except ZeroDivisionError:
            continue
#
def assemble(ipv:List, a:List, aa:List):
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
                    iband = min(ieqn1, ieqn2) - 1
                    jband = abs(ieqn1 - ieqn2)
                    aa[iband][jband] += a[i][j]
                except ZeroDivisionError:
                    continue
        except ZeroDivisionError:
            continue
#
#
# ---------------------------------------------
# pure python
#
def form_Kmatrix(elements, jbc:List, neq:int, iband:int):
    """
    elements
    jbc : node equations
    neq : number of equations
    iband : bandwidth

    :return
    a : global banded stiffness matrix (neq X iband)
    """
    #start_time = time.time()
    aa = zeros(neq, iband)
    for key, element in elements.items():
        idof, jdof = element.DoF
        a = element.Kglobal
        ipv = jbc[idof] + jbc[jdof]
        assemble(ipv, a, aa)
    #
    #uptime = time.time() - start_time
    #print("** [Kb] assembly Finish Process Time: {:1.4e} sec".format(uptime))
    return aa
#
# ---------------------
#
def get_Kmatrix(elements, nodes, boundaries):
    """ """
    #
    free_nodes = elements.get_free_nodes
    jbc, neq, iband = get_bandwidth(elements=elements, 
                                     nodes=nodes, 
                                     boundaries=boundaries, 
                                     free_nodes=free_nodes)
    print("** From datack: half band = {:}".format(iband))
    #
    # ---------------------------
    # multiprocessing
    #
    #call(["python", "steelpy//frame3D//preprocessor//assemblyMatrix.py"])    
    #
    #with Manager() as manager:
    #    d = manager.dict(elements)
    #    l = manager.list(jbc)
    #    st = manager.list(aa)
    #    p = Process(target=loop_members, args=(d, l, st))
    #    p.start()
    #    p.join()
    #
    # ---------------------------
    # Normal
    stf = form_Kmatrix(elements=elements, jbc=jbc,
                        neq=neq, iband=iband)

    # set matrix
    aa = UDUt(stf)   
    #print("** Finished Processing Global [K] Matrix")
    return aa, jbc
#
#
def assemble_banded_Kmatrix(elements, nodes, boundaries):
    """
    Asseable the element matrices in upper band form;
    call separatly from formstif, formmass, formgeom
    -------------------------------------------------
    aa : stiffness matrix
    jbc : nodes freedom
    """
    print("** Processing Global [K] Banded Matrix")
    start_time = time.time()    
    aa, jbc = get_Kmatrix(elements, nodes, boundaries)
    #
    with open("stfmx.f2u", "wb") as f:
        pickle.dump(jbc, f)
        pickle.dump(aa, f)
    #
    uptime = time.time() - start_time
    print("** [Kb] assembly Finish Process Time: {:1.4e} sec".format(uptime))
#
#
#