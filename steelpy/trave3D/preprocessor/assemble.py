# 
# Copyright (c) 2009-2021 steelpy
#
# Python stdlib imports
from itertools import chain
import math as math
#import pickle
import time
from typing import Dict, List, Tuple # , ClassVar, Iterable, Union
#from multiprocessing import Process, Manager

# package imports
from steelpy.trave3D.processor.operations import zeros, to_matrix


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
            r[ 2 ][ 2 ] = (m * sb - l * n * cb) / d
            r[ 2 ][ 2 ] = -(l * sb + m * n * cb) / d
            r[ 2 ][ 2 ] = d * cb
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
    # rt = list( map( list, zip( *r_matrix ) ) )
    rt = list(zip(*r_matrix))
    #
    # take [rtrans][k][r] using the nature of [r] for speed.  
    # k is sectioned off into 3x3s : multiplied [rtrans][k][r]
    ktemp = zeros(12, 12)
    for i in range(4):
        for j in range(4):
            j1 = i * 3
            j2 = j * 3
            for k in range(3):
                for ii in range(3):
                    ktemp[j1 + k][j2 + ii] = sum([ek[j1 + k][j2 + jj] * r_matrix[jj][ii]
                                                  for jj in range(3)])
            # 23
            for k in range(3):
                for ii in range(3):
                    ek[j1 + k][j2 + ii] = sum([rt[k][jj] * ktemp[j1 + jj][j2 + ii]
                                               for jj in range(3)])
            # 24
    # 22
    return ek
#
def assemble(idof, jdof, jbc, a, aa):
    """
    """
    # set ipv to the positions in the array of the nodes
    ipv = [int(idof + i) for i in range(6)]
    ipv.extend([int(jdof + i) for i in range(6)])
    # store the values for individual array in global array
    for i in range(12):
        try:
            1.0 / jbc[ipv[i]]
            ieqn1 = jbc[ipv[i]]
            for j in range(i, 12):
                try:
                    1.0 / jbc[ipv[j]]
                    ieqn2 = jbc[ipv[j]]
                    if ieqn1 > ieqn2:
                        jband = (ieqn1 - ieqn2)
                        aa[ieqn2-1][jband] += a[i][j]
                    else:
                        jband = (ieqn2 - ieqn1)
                        # TODO: check index = jband - 1
                        aa[ieqn1-1][jband] += a[i][j]
                except ZeroDivisionError:
                    continue
            # L30:
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
    for i in range( 12 ):
        for j in range( i, 12 ):
            ek[ j ][ i ] = ek[ i ][ j ]
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
def form_Kmatrix(elements, nodes, materials, 
                 sections, jbc, neq, iband):
    """
    """
    start_time = time.time()
    aa = zeros(neq, iband)
    #for element in elements.values():
    #for element in elements.iter_elements:
    for key, element in elements.items():
        idof, jdof = element.DoF #(nodes)
        a = element.Kmatrix # (nodes, materials, sections)
        assemble(idof, jdof, jbc, a, aa)
    end_time = time.time()
    uptime = end_time - start_time
    print("** [K] assembly Finish Process Time: {:1.4e} sec".format(uptime))
    return aa
#
#
# ---------------------
#
def max_bandwidth(elements, nodes, jbc):
    """
    calculate max bandwidth
    ------------------------  
    npi : connectivity end 1
    npj : connectivity end 2
    jbc : nodes freedom
    nel: number of elements
    
    npi ,npj, jbc, nel
    """
    #ibndm3 = [0]
    ibndm4 = [0]
    for key, element in elements.items():
        conn = element.connectivity
        # end 1
        end_1 = nodes[conn[0]].index
        bc1 = jbc[end_1]
        # end 2
        end_2 = nodes[conn[1]].index
        bc2 = jbc[end_2]
        #
        ieqn = bc1 + bc2
        #
        try:
            ibndm4.append(max([abs(ieqn1 - ieqn2)
                               for x, ieqn1 in enumerate(ieqn) if ieqn1 > 0
                               for ieqn2 in ieqn[x+1:] if ieqn2 > 0]))
            #ibndm3.append(max([abs(ieqn1 - ieqn2)
            #                  for ieqn1 in bc1 if ieqn1 > 0
            #                  for ieqn2 in bc2]))
            #ibndm3.append(max([abs(ieqn1 - ieqn2)
            #                  for ieqn1 in bc1 if ieqn1 > 0
            #                  for ieqn2 in bc2 if ieqn2 > 0]))
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
    jbc = spclbc(elements, nodes, free_nodes, jbc)
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
    iband = max_bandwidth(elements=elements, nodes=nodes,
                          jbc=jbcc)
    return jbcc, neq, iband
#    
#
