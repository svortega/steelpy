#
# Copyright (c) 2009-2023 fem2ufo
#
from __future__ import annotations
#
# Python stdlib imports
import math 
#from typing import Dict, List, ClassVar, Tuple, Iterable, Union

#
# package imports
from steelpy.process.math.operations import zeros, transposeM #, matAbd, trns_3Dv, to_matrix, zeros_vector# FIXME
from steelpy.process.math.vector import Vector
#
#import numpy as np
#import numpy.linalg as la
#from numpy.testing import assert_array_almost_equal
#from numpy import sin, cos, arccos, arctan

#
def cross3(a: list, b: list):
    """Cross produc of 2 vectors"""
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return Vector(c)
#
#
def beam3D_K(length: float, 
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
#
def Rmatrix_new(node1, node2, beta: float = 0):
    """
    Rotation Matrix 
    """
    up = "y"
    sb = math.sin( beta * math.pi / 180.0 )
    cb = math.cos( beta * math.pi / 180.0 )
    #
    # note, local y is being set to a global direction
    if up == "y":
        local_y = Vector([0., 1., 0.])
    elif up == "z":
        local_y = Vector([0., 0., 1.])
    #
    local_x = Vector([x - y for x, y in zip(node2, node1)])
    #
    a = Vector([abs(item) for item in local_x]) / abs(local_x)
    b = Vector([abs(item) for item in local_y]) / abs(local_y)
    #
    is_vertical = all(math.isclose(*values, rel_tol=1e-9, abs_tol=0.0) 
                      for values in zip(a, b))
    #
    if is_vertical:    # set to global x
        local_y = Vector([1., 0., 0.])
    #
    local_z = cross3(local_x, local_y)
    # recalculate local y so that it's orthogonal
    local_y = cross3(local_z, local_x)
    #
    r = [local_x.direction(),
         local_y.direction(),
         local_z.direction()]
    #
    #r[1][0] *= cb
    #r[1][2] += sb
    #
    #r[2][0] *= sb
    #r[2][2] += cb    
    #
    return r
    #return [local_x.direction(),
    #        local_y.direction(),
    #        local_z.direction()]
#
#
def Rmatrix_new2(node1, node2, beta: float = 0):
    """
    Rotation Matrix 
    """
    up = "y"
    sb = math.sin( beta * math.pi / 180.0 )
    cb = math.cos( beta * math.pi / 180.0 )
    #
    # note, local y is being set to a global direction
    if up == "y":
        local_y = Vector([0., 1., 0.])
    elif up == "z":
        local_y = Vector([0., 0., 1.])
    #
    local_x = Vector([x - y for x, y in zip(node2, node1)])
    #
    a = Vector([abs(item) for item in local_x]) / abs(local_x)
    b = Vector([abs(item) for item in local_y]) / abs(local_y)
    #
    is_vertical = all(math.isclose(*values, rel_tol=1e-9, abs_tol=0.0) 
                      for values in zip(a, b))
    #
    if is_vertical:    # set to global x
        local_y = Vector([1., 0., 0.])
    #
    local_z = cross3(local_x, local_y)
    # recalculate local y so that it's orthogonal
    local_y = cross3(local_z, local_x)
    #
    r = [local_x.direction(),
         local_y.direction(),
         local_z.direction()]
    #
    #r[1][0] *= cb
    #r[1][2] += sb
    #
    #r[2][0] *= sb
    #r[2][2] += cb    
    #
    return r
    #return [local_x.direction(),
    #        local_y.direction(),
    #        local_z.direction()]
#
#
def trans_3d_beam(ek: list, r_matrix: list):
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
# ---------------------------------------------
#
def beam3D_Klocal(length: float, 
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
    #areasy = areasz = area
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
def get_element_K(element:tuple, section:tuple, material:tuple):
    """ """
    # solve K matrix
    R = Rmatrix(*element.direction_cosines, element.beta)
    K = beam_Klocal(element.length,
                    section.area, section.J, 
                    section.Iy, section.Iz,
                    material.E, material.G,
                    section.area, section.area)
    return trans_3d_beam(K, R)
#
# ---------------------------------------------
#
def Tr(node1, node2):
    """Local to global transformation matrix"""
    tf = np.zeros((12, 12), dtype=np.float64)
    dc = dircos(node1, node2)
    tf[:3, :3] = tf[3:6, 3:6] = tf[6:9, 6:9] = tf[9:12, 9:12] = dc[:, :]
    return tf
#
def trans3Dbeam(ek, node1, node2):
    """ """
    T = Tr(node1, node2)
    #print('here')
    klocal = np.array(ek)
    return T.transpose() @ klocal @ T 
#    
#
def dircos(node1, node2):
    """
    """
    up = "y"
    # note, local y is being set to a global direction
    if up == "y":
        local_y = np.array([0., 1., 0.], dtype=np.float64)
    elif up == "z":
        local_y = np.array([0., 0., 1.], dtype=np.float64)

    from_vert = np.array(node1)
    to_vert = np.array(node2)
    local_x = to_vert - from_vert
    #
    vert = is_vertical(local_y, local_x)
    if vert:    # set to global x
        local_y = np.array([1., 0., 0.], dtype=np.float64)

    local_z = np.cross(local_x, local_y)

    # recalculate local y so that it's orthogonal
    local_y = np.cross(local_z, local_x)

    dc = np.array([local_x/la.norm(local_x),
                   local_y/la.norm(local_y),
                   local_z/la.norm(local_z)], dtype=np.float64)

    return dc
#
def is_vertical(vertical, local_x):
    try:
        assert_array_almost_equal(np.abs(local_x) / la.norm(local_x),
                                  np.abs(vertical) / la.norm(vertical),
                                  decimal=5)
        return True
    except AssertionError:
    #except AssertionError:
        return False
    except AssertionError as error:
        return False
#
#