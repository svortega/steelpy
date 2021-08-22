# 
# Copyright (c) 2009-2021 steelpy
#
# Python stdlib imports
from itertools import chain
import math as math
import pickle
import time
from typing import Dict, List, Tuple # , ClassVar, Iterable, Union
#from multiprocessing import Process, Manager

# package imports
from steelpy.trave3D.processor.operations import zeros, to_matrix, transposeM
from steelpy.trave3D.processor.static_solver import UDUt

#
# -------------------- 
#     assemble
# --------------------
#
def beam_stiffness(length: float, 
                   area:float, J:float, Iy:float, Iz:float,
                   Emod:float, Gmod:float):
    """
    Calculate the beam element stiffness matrix
    length, J, Iy, Iz, area, Emod, Gmod
    """
    # initialize all ek elements to zero  
    ek = zeros( 12, 12 )
    # stiffness matrix in local coordinates 
    emlen = Emod / length
    emlen2 = emlen / length
    emlen3 = emlen2 / length
    # 
    ek[ 0 ][ 0 ] = area * emlen
    ek[ 1 ][ 1 ] = 12.0 * emlen3 * Iz
    ek[ 2 ][ 2 ] = 12.0 * emlen3 * Iy
    ek[ 3 ][ 3 ] = Gmod * J / length
    ek[ 4 ][ 4 ] = 4.0 * emlen * Iy
    ek[ 5 ][ 5 ] = 4.0 * emlen * Iz
    #
    ek[ 1 ][ 5 ] = 6.0 * emlen2 * Iz
    ek[ 2 ][ 4 ] = -6.0 * emlen2 * Iy
    #
    ek[ 6 ][ 6 ] = ek[ 0 ][ 0 ]
    ek[ 7 ][ 7 ] = ek[ 1 ][ 1 ]
    ek[ 8 ][ 8 ] = ek[ 2 ][ 2 ]
    ek[ 9 ][ 9 ] = ek[ 3 ][ 3 ]
    ek[ 10 ][ 10 ] = ek[ 4 ][ 4 ]
    ek[ 11 ][ 11 ] = ek[ 5 ][ 5 ]
    #
    ek[ 0 ][ 6 ] = -ek[ 0 ][ 0 ]
    ek[ 1 ][ 7 ] = -ek[ 1 ][ 1 ]
    ek[ 1 ][ 11 ] = ek[ 1 ][ 5 ]
    ek[ 2 ][ 8 ] = -ek[ 2 ][ 2 ]
    ek[ 2 ][ 10 ] = ek[ 2 ][ 4 ]
    ek[ 3 ][ 9 ] = -ek[ 3 ][ 3 ]
    ek[ 4 ][ 8 ] = -ek[ 2 ][ 4 ]
    ek[ 4 ][ 10 ] = ek[ 4 ][ 4 ] / 2.0
    ek[ 5 ][ 7 ] = -ek[ 1 ][ 5 ]
    ek[ 5 ][ 11 ] = ek[ 5 ][ 5 ] / 2.0
    #
    ek[ 7 ][ 11 ] = -ek[ 1 ][ 5 ]
    ek[ 8 ][ 10 ] = -ek[ 2 ][ 4 ]
    #
    # impose the geometry
    for i in range( 12 ):
        for j in range( i, 12 ):
            ek[ j ][ i ] = ek[ i ][ j ]
    # L10:
    return ek
#
def Rmatrix(l:float, m:float, n:float, beta: float = 0):
    """
    Rotation Matrix 
    """
    #l, m, n = unit_vector
    sb = math.sin( beta * math.pi / 180.0 )
    cb = math.cos( beta * math.pi / 180.0 )
    d = math.sqrt( 1 - n**2 )
    r = zeros( 3, 3 )
    # member orientated along global vertical z axis
    if abs(n) > 0.995:
        r[ 0 ][ 0 ] = 0.0
        r[ 0 ][ 1 ] = 0.0
        r[ 0 ][ 2 ] = n
        r[ 1 ][ 0 ] = -n * sb
        r[ 1 ][ 1 ] = cb
        r[ 1 ][ 2 ] = 0.0
        r[ 2 ][ 0 ] = -n * cb
        r[ 2 ][ 1 ] = -sb
        r[ 2 ][ 2 ] = 0.0
    else:
        r[ 0 ][ 0 ] = l
        r[ 0 ][ 1 ] = m
        r[ 0 ][ 2 ] = n
        if abs(beta) <= 0.01:
            r[ 1 ][ 0 ] = -m / d
            r[ 1 ][ 1 ] = l / d
            r[ 1 ][ 2 ] = 0.0
            r[ 2 ][ 0 ] = -l * n / d
            r[ 2 ][ 1 ] = -m * n / d
            r[ 2 ][ 2 ] = d
        else:
            r[ 1 ][ 0 ] = -(m * cb + l * n * sb) / d
            r[ 1 ][ 1 ] = (l * cb - m * n * sb) / d
            r[ 1 ][ 2 ] = d * sb
            r[ 2 ][ 0 ] = (m * sb - l * n * cb) / d
            r[ 2 ][ 1 ] = -(l * sb + m * n * cb) / d
            r[ 2 ][ 2 ] = d * cb
    #
    #for i in range(3):
    #    for j in range(3):
    #        try:
    #            1/r[i][j]
    #        except ZeroDivisionError:
    #            r[ i ][ j ] = 0.0
    #
    return r
#
def trans_3d_beam(ek: List, r_matrix: List):
    """
    Makes 3-d coordinate transformations. 
    ek : stiffness matrix
    r  : rotation matrix 
    """
    # transpose rotation matrix
    rt = transposeM(r_matrix)
    #
    # take [rtrans][k][r] using the nature of [r] for speed.  
    # k is sectioned off into 3x3s : multiplied [rtrans][k][r]
    ktemp = zeros(12, 12)
    for j1 in range(0, 12, 3):
        for j2 in range(0, 12, 3):
            # [k][r]
            for k in range(3):
                for ii in range(3):
                    ktemp[j1 + k][j2 + ii] = math.fsum([ek[j1 + k][j2 + jj] * r_matrix[jj][ii]
                                                        for jj in range(3)])
            # [rtrans][k][r]
            for k in range(3):
                for ii in range(3):
                    ek[j1 + k][j2 + ii] = math.fsum([rt[k][jj] * ktemp[j1 + jj][j2 + ii]
                                                    for jj in range(3)])
            # 24
    # 22
    return ek
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
def beam_Ks(length: float, 
            area:float, J:float, Iy:float, Iz:float,
            Emod:float, Gmod:float, 
            areasy:float, areasz:float):
    """
    Calculate the beam element stiffness matrix
    length, J, Iy, Iz, area, Emod, Gmod
    """
    # initialize all ek elements to zero  
    ek = zeros( 12, 12 )
    # stiffness matrix in local coordinates 
    emlen = Emod / length
    emlen2 = emlen / length
    emlen3 = emlen2 / length
    #
    areasy = areasz = area
    Phiy = 12*Emod*Iz / (Gmod * areasy * length**2)
    Phiz = 12*Emod*Iy / (Gmod * areasz * length**2)    
    # 
    ek[ 0 ][ 0 ] = area * emlen
    ek[ 1 ][ 1 ] = 12.0 * emlen3 * Iz / (1+Phiy)
    ek[ 2 ][ 2 ] = 12.0 * emlen3 * Iy / (1+Phiz)
    ek[ 3 ][ 3 ] = Gmod * J / length
    ek[ 4 ][ 4 ] = emlen * Iy * (4+Phiz)/(1+Phiz)
    ek[ 5 ][ 5 ] = emlen * Iz * (4+Phiy)/(1+Phiy)
    #
    ek[ 1 ][ 5 ] =  6.0 * emlen2 * Iz / (1+Phiy)
    ek[ 2 ][ 4 ] = -6.0 * emlen2 * Iy / (1+Phiz)
    #
    ek[ 6 ][ 6 ] = ek[ 0 ][ 0 ]
    ek[ 7 ][ 7 ] = ek[ 1 ][ 1 ]
    ek[ 8 ][ 8 ] = ek[ 2 ][ 2 ]
    ek[ 9 ][ 9 ] = ek[ 3 ][ 3 ]
    ek[ 10 ][ 10 ] = ek[ 4 ][ 4 ]
    ek[ 11 ][ 11 ] = ek[ 5 ][ 5 ]
    #
    ek[ 0 ][ 6 ] = -ek[ 0 ][ 0 ]
    ek[ 1 ][ 7 ] = -ek[ 1 ][ 1 ]
    ek[ 1 ][ 11 ] = ek[ 1 ][ 5 ]
    ek[ 2 ][ 8 ] = -ek[ 2 ][ 2 ]
    ek[ 2 ][ 10 ] = ek[ 2 ][ 4 ]
    ek[ 3 ][ 9 ] = -ek[ 3 ][ 3 ]
    ek[ 4 ][ 8 ] = -ek[ 2 ][ 4 ]
    ek[ 4 ][ 10 ] = emlen * Iy * (2-Phiz)/(1+Phiz)
    ek[ 5 ][ 7 ] = -ek[ 1 ][ 5 ]
    ek[ 5 ][ 11 ] = emlen * Iz * (2-Phiy)/(1+Phiy)
    #
    ek[ 7 ][ 11 ] = -ek[ 1 ][ 5 ]
    ek[ 8 ][ 10 ] = -ek[ 2 ][ 4 ]
    #
    # impose the geometry
    for i in range(12):
        for j in range(i, 12):
            ek[j][i] = ek[i][j]
    # L10:
    return ek    
#
#
def get_element_K(element:Tuple, section:Tuple, material:Tuple):
    """ """
    # solve K matrix
    R = Rmatrix(*element.direction_cosines, element.beta)
    K = beam_Ks(element.length,
                section.area, section.J, 
                section.Iy, section.Iz,
                material.E, material.G,
                section.area, section.area)
    return trans_3d_beam(K, R)
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
    start_time = time.time()
    aa = zeros(neq, iband)
    for key, element in elements.items():
        idof, jdof = element.DoF
        a = element.Kmatrix
        ipv = jbc[idof] + jbc[jdof]
        assemble(ipv, a, aa)
    end_time = time.time()
    uptime = end_time - start_time
    print("** [K] assembly Finish Process Time: {:1.4e} sec".format(uptime))
    return aa
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
    jbc = list(chain.from_iterable(jbc))
    #
    #neq2 = 0
    #jbc2 = [neq2 := neq2+1 if _item == 0  else 0
    #       for _item in jbc]
    #jbc2 = to_matrix(jbc2, 6)
    #
    neq = 0
    for i, item in enumerate(jbc):
        try:
            1 / item
            jbc[i] = 0
        except ZeroDivisionError:
            neq += 1
            jbc[i] = neq
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
# ---------------------
#
def get_stiffnes_matrix(elements, nodes, boundaries):
    """ """
    free_nodes = elements.get_free_nodes
    #
    print("** Processing Global [K] Matrix")
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
    stf = form_Kmatrix(elements= elements, jbc=jbc,
                       neq=neq, iband=iband)
    # set matrix
    aa = UDUt(stf)
    #aa = udu( stf )
    #print("** Finished Processing Global [K] Matrix")
    return aa, jbc
#
#
def assemble_banded_matrix(elements, nodes, boundaries):
    """
    Asseable the element matrices in upper band form;
    call separatly from formstif, formmass, formgeom
    -------------------------------------------------
    aa : stiffness matrix
    jbc : nodes freedom
    """
    aa, jbc = get_stiffnes_matrix(elements, nodes, boundaries)
    #
    with open("stfmx.f2u", "wb") as f:
        pickle.dump(jbc, f)
        pickle.dump(aa, f)
#