#
# Copyright (c) 2009-2019 fem2ufo
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
from steelpy.trave3D.processor.operations import zeros
from steelpy.trave3D.preprocessor.load import get_basic_load
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
    return a  # imult,


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
    return wk  # imult,
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
    #jbc = pickle.load(open( "jbc.f2u", "rb" ))
    memb_force = {}
    basicl_res = {}
    for item in get_basic_load(elements, nodes, materials,
                               sections, basic_load):
        nodal_load = item.nodal_load
        memb_force[item.name] = item.member_load
        #memb_force.update(item.member_load)
        nloads = [nodal_load[i][j] 
                  for i in range(len(jbc)) for j in range(6) 
                  if jbc[i][j] != 0]
        basicl_res[item.name] = Vector(BAK(stf, nloads, neq, iband))
    #
    #pickle.dump( basicl_res, open( "resd.p", "wb" ))
    #pickle.dump( memb_force, open( "elemout.f2u", "wb" ))
    return basicl_res, memb_force
#
def solve_combinations(basic_res, memb_force, load_combination):
    """
    """
    #basic_res = pickle.load(open( "resd.p", "rb" ))
    #memb_force = pickle.load(open( "elemout.f2u", "rb" ))
    comb_res = {}
    memb_comb = {}
    for cname, comb in load_combination.items():
        #comb_res[cname] = zeros_vector(nbl, nnp)
        beam_load = {}
        for bname, factor in comb.basic_load.items():
            try:
                comb_res[cname] += basic_res[bname] * factor
            except KeyError:
                comb_res[cname] = basic_res[bname] * factor
            #
            for mname, member in memb_force[bname].items():
                try:
                    beam_load[mname] += member * factor
                except KeyError:
                    beam_load[mname] =  member * factor
        #
        memb_comb[cname] = beam_load
    #
    #file = open( "elemout.f2u", "ab" )
    #pickle.dump( memb_comb, file)
    #file.close()
    return comb_res, memb_comb
#
#
#
def memload(elements, nodes, materials, sections, 
            dispp, mbload):
    """
    Determine member loads raferred to local s of m convent.
    """
    # global node displacement
    dispp = array('d', list(chain.from_iterable( dispp )))
    gndisp = zeros(12)
    member_load = {}
    #load_global = []
    #load_local = []
    for element in elements.values():
        in1, in2 = element.DoF(nodes)
        # set ipv to the positions in the array of the nodes
        gndisp[:6] = dispp[in1: in1 + 6]
        gndisp[6:12] = dispp[in2: in2 + 6]
        #
        #material = materials[element.material]
        #section = sections[element.section]     
        ## solve K matrix
        #R = Rmatrix(*element.direction_cosines, element.beta)
        R = Rmatrix(*element.unit_vector(nodes), element.beta)
        #K = get_element_K(element, section, material)
        K = element.Kmatrix(nodes, materials, sections)
        # get nodal force global system
        ngforce = matAbd(K, gndisp)
        # get nodal force local system
        nlforce = trns_3Dv(ngforce, R)
        #
        # get nodal force based on beam local system
        index = element.name
        try:
            m_nload = mbload[index]
            nlforce = [(nlforce[i] - m_nload[i]) if i >= 6 
                      else -1*(nlforce[i] - m_nload[i])
                      for i in range(12)]
        except KeyError:
            nlforce = [nlforce[i] if i >= 6 
                       else -1 * nlforce[i] 
                       for i in range(12)]
        #
        member_load[element.name] = {"global":ngforce, "local":nlforce}
    #
    return member_load
    #return member_load, load_global, load_local
#