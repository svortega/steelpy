#
# Copyright (c) 2009-2021 fem2ufo
#
# Python stdlib imports
# from math import fsum
from array import array
import pickle
from typing import List #, NamedTuple, Union
from itertools import chain
#
# package imports
from steelpy.trave3D.preprocessor.assemble import Rmatrix # get_element_K, 
#from steelpy.trave3D.processor.operations import zeros
#from steelpy.f2uModel.load.operations.basic_load import get_basic_load
from steelpy.trave3D.processor.operations import zeros, matAbd, trns_3Dv
from steelpy.process.math.vector import Vector
#
#
# --------------------
# solver
# --------------------
#
def UDUt(a: List, neq: int, iband: int) -> List:
    """
    Solution of banded system asing UDU decomposition 
    based on Weaver & Gare pp 469-472.  
    """
    print("** Processing UDU")
    imult = 0
    if a[0][0] == 0.0:
        raise RuntimeError("error:  zero diagonal term")
    #
    if neq == 1:
        return a
    #
    for j in range(1, neq):
        j2 = max(j - iband + 1, 0)
        # off-diagonal terms
        for i in range(j2 + 1, j):
            a[i][j - i] -= sum([a[k][i - k] * a[k][j - k]
                                for k in range(j2, i)])
            imult += abs(i - j2)
        # diagonal  terms
        for k in range(j2, j):
            a[j][0] -= a[k][j - k] * a[k][j - k] / a[k][0]
            a[k][j - k] /= a[k][0]
            imult += 2
        # L30:
        try:
            1.0 / a[j][0]
        except ZeroDivisionError:
            raise RuntimeError("error:  zero diagonal term")
    # L10:
    print("** Finished Processing UDU: mults & divs {:}".format(imult))
    return a
#
#
def BAK(a: List, b: List[float], neq: int, iband: int) -> List:
    """
    back substitution
    """
    #print("** Calculating Joint Displacements")
    wk = zeros(neq)
    # imult = 0
    # forward substitution
    for i in range(neq):
        j = max(i - iband + 1, 0)
        wk[i] = b[i]
        if j < i:
            wk[i] -= sum([a[k][i - k] * wk[k]
                          for k in range(j, i)])
            # imult += abs(i - j)
    # middle terms
    wk = array('d', [wk[i] / a[i][0] for i in range(neq)])
    # imult += neq
    #  backward substitution
    for i1 in range(neq):
        i = neq - i1 - 1
        j = min(i + iband, neq)
        k2 = i + 1
        if k2 < j:
            wk[i] -= sum([a[i][k - i] * wk[k]
                          for k in range(k2, j)])
            # imult += abs(j - k2)
    # L50:
    #print("** Finished Calculating Joint Displacements")
    return wk
#
#
# --------------------
#
#
# --------------------
#
#
#
def solve_basic_load(elements, nodes, materials,
                     sections, basic_load):
    """
    """
    file = open( "stfmx.f2u", "rb" )
    neq = pickle.load(file)
    iband = pickle.load(file)
    jbc = pickle.load(file)
    stf = pickle.load(file)
    file.close()
    #
    memb_force = {}
    basicl_res = {}
    for item in basic_load.get_basic_load(elements, nodes, 
                                          materials, sections):
        nodal_load = item.nodal_load
        memb_force[item.title] = item.member_load
        nloads = [nodal_load[i][j] 
                  for i in range(len(jbc)) for j in range(6) 
                  if jbc[i][j] != 0]
        basicl_res[item.title] = Vector(BAK(stf, nloads, neq, iband))
    return basicl_res, memb_force
#
#
#
def beam_force(elements, disp, mbload):
    """
    Determine element beam loads raferred to local s of m convent.
    """
    # global node displacement
    dispp = array('d', list(chain.from_iterable(disp)))
    gndisp = zeros(12)
    member_load = []
    for element in elements.values():
        in1, in2 = element.DoF
        # set ipv to the positions in the array of the nodes
        gndisp[:6] = dispp[in1: in1 + 6]
        gndisp[6:12] = dispp[in2: in2 + 6]
        # solve K matrix
        R = Rmatrix(*element.unit_vector, element.beta)
        K = element.Kmatrix
        # get nodal force global system
        ngforce = matAbd(K, gndisp)
        # get nodal force local system
        nlforce = trns_3Dv(ngforce, R)
        # get nodal force based on beam local system
        index = element.name
        nodes = element.connectivity
        try:
            m_nload = mbload[index]
            nlforce = [(nlforce[i] - m_nload[i]) if i >= 6 
                       else -1*(nlforce[i] - m_nload[i])
                       for i in range(12)]
        except KeyError:
            nlforce = [nlforce[i] if i >= 6 
                       else -1 * nlforce[i] 
                       for i in range(12)]
        #member_load[element.name] = {"global":ngforce, "local":nlforce}
        member_load.append([index, "global", nodes[0], *ngforce[:6], nodes[1], *ngforce[6:]])
        member_load.append([index, "local", nodes[0], *nlforce[:6], nodes[1], *nlforce[6:]])
    return member_load
#
#