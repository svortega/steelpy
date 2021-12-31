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
from steelpy.trave3D.processor.operations import matAbd, trns_3Dv, to_matrix, zeros_vector
from steelpy.process.math.vector import Vector
#
#
# --------------------
# solver
# --------------------
#
def UDUt(a: List[float]) -> List:
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
                       for k in range(i + 1, j)])
    #
    # for i in range(neq-1 , -1, -1):
    #    wkk[i] -= sum([a[i][k - i] * wkk[k]
    #                  for k in range(iband-1, i+1, -1)])
    return wk
#
#
# --------------------
#
#
# --------------------
#
#
def solve_displacement(elements, nodes, boundaries,
                       materials, sections,
                       basic_load, load_combination):
    """
    """
    #
    file = open( "stfmx.f2u", "rb" )
    jbc = pickle.load( file )
    stf = pickle.load( file )
    file.close()
    #
    jbcc = list( chain.from_iterable( jbc ) )
    #
    print("** Calculating Joint Displacements")
    print("** reloaded [k] & {p}")
    #
    fmemb = {}
    memf_basic = {}
    basic_res = {}
    for key, item in basic_load.items():
        load_res = basic_load.get_basic_load(key, elements, nodes, boundaries)
        # member end load in local system
        memf_basic[item.title] = load_res.member_load
        # node load in global system
        nodal_load = load_res.nodal_load
        nloads = [nodal_load[i][j] 
                  for i in range(len(jbc)) for j in range(6) 
                  if jbc[i][j] != 0]
        # get displacement in global system
        basic = Vector(BAK(stf, nloads))
        basic_res[item.title] = basic  # Vector(BAK(stf, nloads))
    #
    comb_res, memf_comb = load_combination.solve_combinations(basic_res, memf_basic)
    #
    with open("elemfout.f2u", "wb") as f:
        pickle.dump(memf_basic, f)
        pickle.dump(memf_comb, f)
    #
    # expand displacement to full vector
    print("** assigning node displacement")
    #
    deflection = []
    # basic load
    # TODO: node name instead number
    node_index = {item.index:key for key, item in nodes.items()}
    for load_tile, basic in basic_res.items():
        disp = [basic[ieqnum - 1] if ieqnum != 0 else ieqnum
                for ieqnum in jbcc]
        disp = to_matrix(disp, 6)
        disp = [[node_index[number], load_tile, "global"] + item
                for number, item in enumerate(disp)]
        deflection.extend(disp)
    # load combination
    for load_tile, comb in comb_res.items():
        disp = [comb[ieqnum - 1] if ieqnum != 0 else ieqnum
                for ieqnum in jbcc]
        disp = to_matrix(disp, 6)
        disp = [[node_index[number], load_tile, "global"] + item
                for number, item in enumerate(disp)]
        deflection.extend(disp)
    print("** Finished Calculating Joint Displacements")
    return deflection
#
#
def beam_end_force(elements, ndisp:List[List[float]], mbload:List[Vector]):
    """
    Determine element beam loads raferred to local s of m convent.
    """
    nnp = len(ndisp)
    nreac = zeros_vector(nnp, 6, code='I')
    # global node displacement
    beam_res = {}
    member_load = []
    for element in elements.values():
        mname = element.name
        in1, in2 = element.DoF
        nodes = element.connectivity
        # set ipv to the positions in the array of the nodes
        gndisp = ndisp[in1] + ndisp[in2]
        #gndisp = ndisp[nodes[0]] + ndisp[nodes[1]]
        # solve K matrix
        R = element.R
        K = element.Kmatrix
        # get nodal force global system
        ngforce = matAbd(K, gndisp)
        lcndisp = trns_3Dv(gndisp, R)
        # get nodal force local system
        nlforce = trns_3Dv(ngforce, R)
        # get nodal force based on beam local system
        try:
            try:
                m_nload = mbload[mname].end_force
            except AttributeError:
                m_nload = mbload[mname]
            #nlforce = [(nlforce[i] - m_nload[i]) if i >= 6
            #           else 1*(nlforce[i] - m_nload[i])
            #           for i in range(12)]
            nlforce = [nlforce[i] - m_nload[i] for i in range(12)]
            # [FV,FM,Fw,Ftheta]
            #nlforce2 = [nlforce[i] if i >= 6 else -1*nlforce[i]
            #           for i in range(12)]
            #bres1 = [ -nlforce[ 1 ], nlforce[ 5 ], lcndisp[ 1 ], -lcndisp[ 5 ] ]
            #bres2 = [ nlforce[ 2 ], nlforce[ 4 ],  lcndisp[ 2 ], lcndisp[ 4 ] ]
            #beam_res[mname] = mbload[mname].response([bres1, bres2])
            #1/0
        except KeyError:
            #nlforce = [nlforce[i] if i >= 6 
            #           else 1 * nlforce[i] 
            #           for i in range(12)]
            mbload[mname] = element.beam()
        # get beam forces along length
        bres1 = [ -nlforce[ 1 ], nlforce[ 5 ], lcndisp[ 1 ], -lcndisp[ 5 ] ]
        bres2 = [ nlforce[ 2 ], nlforce[ 4 ],  lcndisp[ 2 ], lcndisp[ 4 ] ]
        try:
            beam_res[mname] = mbload[mname].response([bres1, bres2])
        except AttributeError:
            beam_res[mname] = mbload[mname]
        #
        nreac[in1] += Vector(ngforce[:6])
        nreac[in2] += Vector(ngforce[6:])
        member_load.append([mname, "local", 
                            nodes[0], *nlforce[:6], 
                            nodes[1], *nlforce[6:]])
    return member_load, nreac
#
def solve_forces(elements, displacement):
    """
    """
    #basicitems = {item.title: key for key, item in basic_load.items()}
    start_time = time.time()
    #print(" ")
    print("** Calculating Member Forces")
    #print("** Reloaded Joint Displacements")
    #nnp = len(nodes)
    membf = []
    nres = []
    filem = open("elemfout.f2u", "rb")
    mbload = pickle.load(filem)
    memf_basic = []
    for key, disp in displacement["basic"].items():
        #basic_name = basicitems[key]
        m_load = mbload[key]
        # select [x, y, z, rx, ry, rz] with list in node index order
        ndisp = [item[1:] for item in disp.items]
        #ndisp = {item[0]:item[ 1: ] for item in disp.items }
        force, nreac = beam_end_force(elements, ndisp, m_load)
        nres.append([key, nreac])
        memf_basic.extend([[item[0], key, item[2], 0, item[1], *item[3:9]]
                           for item in force])
        memf_basic.extend([[item[0], key, item[9], 1, item[1], *item[10:]]
                           for item in force])
    membf.extend(memf_basic)
    #
    if displacement["combination"]:
        memf_comb = []
        mcload = pickle.load(filem)
        for key, disp in displacement["combination"].items():
            m_load = mcload[key]
            ndisp = [item[1:] for item in disp.items]
            #force, nreac = beam_force(elements, ndisp, nreac, m_load)
            force, nreac = beam_end_force(elements, ndisp, m_load)
            nres.append([key, nreac])
            memf_comb.extend([[item[0], key, item[2], 0, item[1], *item[3:9]]
                              for item in force])
            memf_comb.extend([[item[0], key, item[9], 1, item[1], *item[10:]]
                              for item in force])
        membf.extend(memf_comb)
    #
    filem.close()
    end_time = time.time()
    uptime = end_time - start_time
    print("** Calculating Member Forces Process Time: {:1.4e} sec".format(uptime))       
    print("** End Calculating Member Forces")
    return membf, nres
#
#