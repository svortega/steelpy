# 
# Copyright (c) 2009-2022 fem2ufo
#
# Python stdlib imports
import itertools as it
from typing import List, ClassVar, Dict, NamedTuple, Union
#from array import array
from itertools import chain
#from copy import  copy
#import math as math
import pickle
import time
#
# package imports
from steelpy.process.math.vector import Vector
from steelpy.process.math.operations  import zeros, to_matrix, zeros_vector, matAbd, trns_3Dv
#
#

# ---------------------
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
    ibndm4 = [0]
    for key, element in elements.items():
        idof, jdof = element.DoF
        bc1 = jbc[idof]
        bc2 = jbc[jdof]
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
def bd_condition(nodes, boundaries):
    """
    set boundary conditions
    bcs: set default = 0 free (1 fix)
    """
    nnp = len(nodes)
    jbc = zeros(nnp, 6, code='I')
    for node_name, bd in boundaries.node.items():
        ind = nodes[bd.node].index
        jbc[ind] = list(bd[:6])
    return jbc
#
def spclbc(elements, nodes, free_nodes, jbc):
    """
    Impose condition for  special GLOBAL shapes
    bcs: set default = 0 free (1 fix)
    """
    # free_nodes = elements.get_free_nodes()
    #
    #for _element in elements.values():
    #for element in elements.iter_elements:
    for key, element in elements.items():
        conn = element.connectivity
        pos_node = set(conn) - set(free_nodes)
        for node_name in pos_node:
            ind = nodes[node_name].index
            #ind
            # jbc[ind][3:] = [1, 1, 1]
            # jbc[ind][3] = 1
            # jbc[ind][4] = 1
        # if _element.type == "truss":
        #    
        #    # end 1
        #    ind = self._nodes._labels.index(conn[0])
        #    jbc[ind][3:] = [0, 0, 0]
        #    # end 2
        #    ind = self._nodes._labels.index(conn[1])
        #    jbc[ind][3:] = [0, 0, 0]
    return jbc
#
def shape_cond(elements, nodes, boundaries, free_nodes):
    """
    jcs: modify default = 0 free (1 fix)
    """
    jbc = bd_condition(nodes, boundaries)
    # TODO : check this module
    #jbc = spclbc(elements, nodes, free_nodes, jbc)
    #
    # Number the equations  in jbc from 1 up to the order.
    # Start assigning equation numbers for zero dof's
    # from 1 up;  only zero given a number.        
    #
    jbc = list(it.chain.from_iterable(jbc))
    counter = it.count(start=1)
    jbc = [next(counter) if _item == 0  else 0
            for _item in jbc]
    neq = max(jbc)
    jbc = to_matrix(jbc, 6)
    return jbc, neq
#
def get_bandwidth(elements, nodes, boundaries, free_nodes):
    """ """
    jbcc, neq = shape_cond(elements=elements, 
                           nodes=nodes, 
                           boundaries=boundaries, 
                           free_nodes=free_nodes)
    #
    iband = max_bandwidth(elements=elements, jbc=jbcc)
    return jbcc, neq, iband
#
#
#
def get_deflection(nodes, basic_res, comb_res):
    """Expand displacement to full vector"""
    # 
    print("** assigning node displacement")
    deflection = []
    # basic load
    # TODO: node name instead number
    node_index = {item.index:key for key, item in nodes.items()}
    for load_tile, basic in basic_res.items():
        disp = to_matrix(basic, 6)
        disp = [[node_index[number], load_tile, "global"] + item.tolist()
                for number, item in enumerate(disp)]
        deflection.extend(disp)
    # load combination
    for load_tile, comb in comb_res.items():
        disp = to_matrix(comb, 6)
        disp = [[node_index[number], load_tile, "global"] + item.tolist()
                for number, item in enumerate(disp)]
        deflection.extend(disp)
    return deflection
#
# ---------------------------------------------
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
        K = element.Kglobal
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