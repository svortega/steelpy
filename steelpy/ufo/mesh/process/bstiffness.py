#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
#
# Python stdlib imports
import math
from operator import sub

#
# package imports
#from steelpy.utils.math.operations import zeros  #, transposeM, matAbd, trns_3Dv, to_matrix, zeros_vector# FIXME
from steelpy.utils.math.vector import Vector
#
import numpy as np
#import numpy.linalg as la
#from numpy.testing import assert_array_almost_equal
#from numpy import sin, cos, arccos, arctan

#
# ---------------------------------------------
#
def B3D2_Ke(Le: float,
            Ax: float, Asy: float, Asz: float,
            Jx: float, Iy: float, Iz: float,
            Emod: float, Gmod: float,
            shear: bool = True) -> list:
    """
    3D elastic stiffness matrix for frame elements in local coordinates
    including bending and shear deformation (H.P. Gavin).
    
    Input:  
    Le : Effective length
    Ax, Asy, Asz : cross-section areas, including shear
    Jx, Iy, Iz : sections' inertia
    Emod : elastic modulus
    Gmod : shear moduli
    shear : Flag to include shear deformation [True/False]
    
    Return :
    ek : beam local stiffness matrix
    """
    # stiffness matrix in local coordinates
    emlen = Emod / Le
    emlen2 = emlen / Le
    emlen3 = emlen2 / Le
    #
    Phiy = 0
    Phiz = 0
    if shear:
        Phiy = Emod * Iz / (Gmod * Asy * Le**2)
        Phiz = Emod * Iy / (Gmod * Asz * Le**2)
    #
    niy = 1.0 + 12.0 * Phiy
    niz = 1.0 + 12.0 * Phiz
    #
    lay = 1.0 + 3.0 * Phiy
    laz = 1.0 + 3.0 * Phiz
    #
    gay = 1.0 - 6.0 * Phiy
    gaz = 1.0 - 6.0 * Phiz
    #
    # initialize all ek elements to zero
    ek = np.zeros((12, 12))
    #
    ek[0][0] = ek[6][6] = Ax * emlen
    
    ek[1][1] = ek[7][7] = 12.0 * emlen3 * Iz / niy
    ek[2][2] = ek[8][8] = 12.0 * emlen3 * Iy / niz
    
    ek[3][3] = ek[9][9] = Gmod * Jx / Le
    
    ek[4][4] = ek[10][10] = 4.0 * emlen * Iy * laz / niz
    ek[5][5] = ek[11][11] = 4.0 * emlen * Iz * lay / niy
    #
    ek[1][5] = ek[1][11] = 6.0 * emlen2 * Iz / niy
    ek[2][4] = ek[2][10] = -6.0 * emlen2 * Iy / niz
    #
    ek[0][6] = -ek[0][0]
    ek[1][7] = -ek[1][1]
    #
    ek[2][8] = -ek[2][2]
    #
    ek[3][9] = -ek[3][3]
    ek[4][8] = -ek[2][4]
    ek[4][10] = 2.0 * emlen * Iy * gaz / niz
    ek[5][7] = -ek[1][5]
    ek[5][11] = 2.0 * emlen * Iz * gay / niy
    #
    ek[7][11] = -ek[1][5]
    ek[8][10] = -ek[2][4]
    #
    # impose the geometry
    ek += np.triu(ek, k=1).T
    return ek
#
# ---------------------------------------------
#
def beam3D_B3D2(Le: float,
                Ax: float, Asy: float, Asz: float, 
                Jx: float, Iy: float, Iz: float,
                Emod: float, Gmod: float,
                shear=True):
    """
    Beam's stiffness matrix based on:
    An Efficient 3D Timoshenko Beam Element
    with Consistent Shape Functions (Yunhua Luo)
    
    Ax, Asy, Asz : cross section areas, including shear
    """
    # Shear correction factor (k)
    # k usually lies in the range between 2/3 and 1.0,
    # and predictions of first-order shear deformation theory
    # appear to be relative insensitive to the value
    # chosen within this range.
    # k = 5 / 6 # rectangular section
    #
    # initialize all ek elements to zero
    ek = np.zeros((12, 12))
    # stiffness matrix in local coordinates
    emlen = Emod / Le
    #emlen2 = emlen / length
    #emlen3 = emlen2 / length
    #
    Phiy = ((12 * Emod * Iy + Asy * Gmod * Le**2)
            / (12 * Emod * Iy - Asy * Gmod * Le**2)**2)
    #
    Phiz = ((12 * Emod * Iz + Asz * Gmod * Le**2)
            / (12 * Emod * Iz - Asz * Gmod * Le**2)**2)
    #
    #
    ek[0][0] = Ax * emlen
    #
    ek[1][1] = 12.0 * Asy * Gmod * Emod * Iy * Phiy / Le
    ek[1][5] = 6.00 * Asy * Gmod * Emod * Iy * Phiy
    #
    ek[2][2] = 12.0 * Asz * Gmod * Emod * Iz * Phiz / Le
    ek[2][4] = -6.0 * Asz * Gmod * Emod * Iz * Phiz
    #
    ek[3][3] = Gmod * Jx / Le
    #
    ek[4][4] = (4 * Emod * Iz
                * ((Asz * Gmod)**2 * Le**4
                   + (3 * Asz * Gmod * Le**2 * Emod * Iz)
                   + 36 * (Emod * Iz)**2)
                / (Le * (12 * Emod * Iz - Asz * Gmod *  Le**2)**2))

    #ek[4][8] = -ek[2][4]

    ek[4][10] = (-2 * Emod * Iz
                 * (72 * (Emod * Iz)**2
                    - (Asz * Gmod)**2 * Le**4
                    - 30 * Asz * Gmod * Le**2 * Emod * Iz)
                 / (Le * (12 * Emod * Iz - Asz * Gmod * Le**2)**2))
    #
    ek[5][5] = (4 * Emod * Iy
                * ((Asy * Gmod )**2 * Le**4
                 + 3 * Asy * Gmod * Le**2 * Emod * Iy
                   + 36 * (Emod * Iy)**2 )
                / (Le * (12 * Emod * Iy - Asy * Gmod *  Le**2)**2))

    #ek[5][7] = -ek[1][5]
    ek[5][11] = (-2 * Emod * Iy
                 * (-(Asy * Gmod )**2 * Le**4
                  - 30 * Asy * Gmod * Le**2 * Emod * Iy
                    + 72 * (Emod * Iy)**2)
                 / (Le * (12 * Emod * Iy - Asy * Gmod *  Le**2)**2))
    #
    #
    ek[6][6] = ek[0][0]
    ek[7][7] = ek[1][1]
    ek[8][8] = ek[2][2]
    ek[9][9] = ek[3][3]
    ek[10][10] = ek[4][4]
    ek[11][11] = ek[5][5]
    #
    ek[0][6] = -ek[0][0]
    ek[1][7] = -ek[1][1]
    ek[1][11] = ek[1][5]
    ek[2][8] = -ek[2][2]
    ek[2][10] = ek[2][4]
    ek[3][9] = -ek[3][3]
    ek[4][8] = -ek[2][4]
    ek[5][7] = -ek[1][5]
    #
    ek[7][11] = -ek[1][5]
    ek[8][10] = -ek[2][4]
    #
    # impose the geometry
    ek += np.triu(ek, k=1).T
    return ek
#
# ---------------------------------------------
#
#
def cross3(a: list, b: list):
    """Cross produc of 2 vectors"""
    c = [a[1] * b[2] - a[2] * b[1],
         a[2] * b[0] - a[0] * b[2],
         a[0] * b[1] - a[1] * b[0]]

    return Vector(c)
#
#
#def Rmatrix(l: float, m: float, n: float, beta: float = 0):
def RmatrixXX(nodei:list, nodej:list, beta:float):
    """
    Rotation Matrix 
    """
    #l, m, n = unit_vector
    #
    # direction cosines
    L = math.dist(nodei[:3], nodej[:3])
    uv = list(map(sub, nodej[:3], nodei[:3]))
    l, m, n = [item / L for item in uv]
    #
    beta = np.deg2rad(beta)
    sb = np.sin(beta)
    cb = np.cos(beta)    
    #sb = math.sin(beta * math.pi / 180.0)
    #cb = math.cos(beta * math.pi / 180.0)
    d = np.sqrt(1 - n**2)
    #
    r = np.zeros((3, 3))
    # member orientated along global vertical z axis
    if abs(n) > 0.995:
        r[0][0] = 0.0
        r[0][1] = 0.0
        r[0][2] = n
        r[1][0] = -n * sb
        r[1][1] = cb
        r[1][2] = 0.0
        r[2][0] = -n * cb
        r[2][1] = -sb
        r[2][2] = 0.0
    else:
        r[0][0] = l
        r[0][1] = m
        r[0][2] = n
        if abs(beta) <= 0.01:
            r[1][0] = -m / d
            r[1][1] = l / d
            r[1][2] = 0.0
            r[2][0] = -l * n / d
            r[2][1] = -m * n / d
            r[2][2] = d
        else:
            r[1][0] = -(m * cb + l * n * sb) / d
            r[1][1] = (l * cb - m * n * sb) / d
            r[1][2] = d * sb
            r[2][0] = (m * sb - l * n * cb) / d
            r[2][1] = -(l * sb + m * n * cb) / d
            r[2][2] = d * cb
    #
    #
    # Build full transformation matrix
    #transMatrix = np.zeros((12, 12))
    #transMatrix[0:3, 0:3] = r
    #transMatrix[3:6, 3:6] = r
    #transMatrix[6:9, 6:9] = r
    #transMatrix[9:12, 9:12] = r
    #
    #return Tmatrix(r)
    #return transMatrix
    #
    Rc = np.isclose(r, np.array([[0, 0, 0],
                                 [0, 0, 0],
                                 [0, 0, 0]]))
    r[Rc] = float(0)    
    return r
#
#
def Rmatrix(nodei:list, nodej:list, beta:float):
    """
    Rotation Matrix 
    """
    #l, m, n = unit_vector
    #
    # direction cosines
    L = math.dist(nodei[:3], nodej[:3])
    uv = list(map(sub, nodej[:3], nodei[:3]))
    Cx, Cy, Cz = [item / L for item in uv]
    #
    beta = np.deg2rad(beta)
    Sp = np.sin(beta)
    Cp = np.cos(beta)    
    #sb = math.sin(beta * math.pi / 180.0)
    #cb = math.cos(beta * math.pi / 180.0)
    #d = np.sqrt(1 - n**2)
    #
    r = np.zeros((3, 3))
    # member orientated along global vertical z axis
    if abs(Cz) == 1:
        r[0][2] =  Cz
        r[1][0] = -Cz*Sp
        r[1][1] =  Cp
        r[2][0] = -Cz*Cp
        r[2][1] = -Sp
    else:
        den = np.sqrt ( 1.0 - Cz*Cz );
    
        r[0][0] = Cx
        r[0][1] = Cy
        r[0][2] = Cz

        r[1][0] = (-Cx*Cz*Sp - Cy*Cp)/den
        r[1][1] = (-Cy*Cz*Sp + Cx*Cp)/den
        r[1][2] = Sp*den

        r[2][0] = (-Cx*Cz*Cp + Cy*Sp)/den
        r[2][1] = (-Cy*Cz*Cp - Cx*Sp)/den
        r[2][2] = Cp*den
    #
    if abs(Cy) == 1.0:
        r[0][1] =  Cy
        r[1][0] = -Cy*Cp
        r[1][2] =  Sp
        r[2][0] =  Cy*Sp
        r[2][2] =  Cp
    else:
        den = np.sqrt( 1.0 - Cy*Cy );
    
        r[0][0] = Cx
        r[0][1] = Cy
        r[0][2] = Cz

        r[1][0] = (-Cx*Cy*Cp - Cz*Sp)/den
        r[1][1] = den*Cp
        r[1][2] = (-Cy*Cz*Cp + Cx*Sp)/den
    
        r[2][0] = (Cx*Cy*Sp - Cz*Cp)/den
        r[2][1] = -den*Sp
        r[2][2] = (Cy*Cz*Sp + Cx*Cp)/den
    #
    #
    Rc = np.isclose(r, np.array([[0, 0, 0],
                                 [0, 0, 0],
                                 [0, 0, 0]]))
    r[Rc] = float(0)    
    return r
#
#
# ---------------------------------------------
#
# Transformation matrix
#
#
def unitvec_0(nodei: list, nodej: list,
              beta: float,
              auxNode: list | None = None):
    """
    Input:
    nodei: [x,y,z]
    nodej: [x,y,z]
    beta: angle of roll
    auxNode : [x,y,z]|None

    Returns the transformation matrix for the member.
    r : [rxX, rxY, rxZ,
         ryX, ryY, ryZ,
         rzX, rzY, rzZ]
    """
    #
    x1, y1, z1 = nodei
    x2, y2, z2 = nodej
    L = math.dist(nodei, nodej)
    # Calculate the direction cosines for the local x-axis
    # [rxX, rxY, rxZ] = [cos(theta_xX), cos(theta_xY), cos(theta_xZ)]
    x = [(x2 - x1) / L, (y2 - y1) / L, (z2 - z1) / L]
    #
    # Beam cross-section orientation
    beta = np.deg2rad(beta)
    sb = np.sin(beta)
    cb = np.cos(beta)
    Rb = np.array([[1, 0, 0],
                   [0, cb, sb],
                   [0, -sb, cb]])
    #
    # Calculate the direction cosines for the local z-axis based on the auxiliary node
    if auxNode != None:
        xa, ya, za = auxNode
        # Define a vector in the local xz plane using the auxiliary point
        z = [xa - x1, ya - y1, za - z1]
        # Find the direction cosines for the local y-axis
        y = np.cross(z, x)
        y = np.divide(y, (y[0]**2 + y[1]**2 + y[2]**2)**0.5)
        # Ensure the z-axis is perpendicular to the x-axis.
        # If the vector from the i-node to the auxiliary node is not perpendicular
        # to the member, this will ensure the local coordinate system is orthogonal
        z = np.cross(x, y)
        # Turn the z-vector into a unit vector of direction cosines
        z = np.divide(z, (z[0]**2 + z[1]**2 + z[2]**2)**0.5)
    # If no auxiliary node is specified the program will determine the member's
    # local z-axis automatically
    else:
        # Calculate the remaining direction cosines. The local z-axis will be 
        # kept parallel to the global XZ plane in all cases
        # Vertical members
        if np.isclose(x1, x2) and np.isclose(z1, z2):
            # For vertical members, keep the local y-axis in the XY plane 
            # to make 2D problems easier to solve in the XY plane
            if y2 > y1:
                y = [-1, 0, 0]
                z = [0, 0, 1]
            else:
                y = [1, 0, 0]
                z = [0, 0, 1]
            #
            #Rb = np.array([[cb, 0, sb],
            #               [0, 1, 0],
            #               [-sb, 0, cb]])
        # Horizontal members
        elif np.isclose(y1, y2):
            # Find a vector in the direction of the local z-axis by taking the cross-product
            # of the local x-axis and the local y-axis.
            # This vector will be perpendicular to both the local x-axis and the local y-axis.
            y = [0, 1, 0]
            z = np.cross(x, y)
            # Divide the z-vector by its magnitude to produce a unit vector of direction cosines
            z = np.divide(z, (z[0]**2 + z[1]**2 + z[2]**2)**0.5)
            #
            #Rb = np.array([[cb, -sb, 0],
            #               [sb, cb, 0],
            #               [0, 0, 1]])
        # Members neither vertical or horizontal
        else:
            # Find the projection of x on the global XZ plane
            proj = [x2 - x1, 0, z2 - z1]
            # Find a vector in the direction of the local z-axis by taking the cross-product
            # of the local x-axis and its projection on a plane parallel to the XZ plane.
            # This produces a vector perpendicular to both the local x-axis and its projection.
            # This vector will always be horizontal since it's parallel to the XZ plane.
            # The order in which the vectors are 'crossed' has been selected to ensure the y-axis
            # always has an upward component (i.e. the top of the beam is always on top).
            if y2 > y1:
                z = np.cross(proj, x)
            else:
                z = np.cross(x, proj)
            # Divide the z-vector by its magnitude to produce a unit vector of direction cosines
            z = np.divide(z, (z[0]**2 + z[1]**2 + z[2]**2)**0.5)
            # Find the direction cosines for the local y-axis
            y = np.cross(z, x)
            y = np.divide(y, (y[0]**2 + y[1]**2 + y[2]**2)**0.5)
    #
    # Create the direction cosines matrix
    dirCos = np.array([x, y, z])
    #
    # Beam cross-section orientation
    #R = np.dot(Rb, dirCos)
    R = Rb @ dirCos
    Rc = np.isclose(R, np.array([[0, 0, 0],
                                 [0, 0, 0],
                                 [0, 0, 0]]))
    R[Rc] = float(0)
    return R
#
def Tmatrix(dirCos: list):
    """ Build the transformation matrix """
    # Build the transformation matrix
    transMatrix = np.zeros((12, 12))
    transMatrix[0:3, 0:3] = dirCos
    transMatrix[3:6, 3:6] = dirCos
    transMatrix[6:9, 6:9] = dirCos
    transMatrix[9:12, 9:12] = dirCos
    return transMatrix
#
def Rmatrix2(nodei, nodej, auxNode=None):
    """ Returns the transformation matrix for the member."""
    dirCos = unitvec_0(nodei[:3],
                       nodej[:3],
                       auxNode)
    return Tmatrix(dirCos)
#
#
def Rmatrix_new(node1, node2, beta: float = 0):
    """
    Rotation Matrix 
    """
    up = "y"
    sb = math.sin(beta * math.pi / 180.0)
    cb = math.cos(beta * math.pi / 180.0)
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
    sb = math.sin(beta * math.pi / 180.0)
    cb = math.cos(beta * math.pi / 180.0)
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
# ---------------------------------------------
#
def is_vertical(vertical, local_x):
    try:
        np.testing.assert_allclose(np.abs(local_x) / np.linalg.norm(local_x),
                                   np.abs(vertical) / np.linalg.norm(vertical),
                                   rtol=1e-5, atol=0)
        return True
    except AssertionError:
        return False
#
def dircos(node1, node2, beta:float = 0, 
           up:str = "y"):
    """
    """
    #sb = np.sin(beta * math.pi / 180.0)
    #cb = np.cos(beta * math.pi / 180.0)
    #
    #up = "y"
    # note, local y is being set to a global direction
    #if up == "y":
    local_y = np.array([0., 1., 0.], dtype=np.float64)
    if up == "z":
        local_y = np.array([0., 0., 1.], dtype=np.float64)
    #
    local_x = np.array(node2) - np.array(node1)
    #
    vert = is_vertical(local_y, local_x)
    if vert:    # set to global x
        local_y = np.array([1., 0., 0.], dtype=np.float64)

    local_z = np.cross(local_x, local_y)
    #
    # recalculate local y so that it's orthogonal
    local_y = np.cross(local_z, local_x)
    #
    dc = np.array([local_x / np.linalg.norm(local_x),
                   local_y / np.linalg.norm(local_y),
                   local_z / np.linalg.norm(local_z)], dtype=np.float64)
    #
    return dc
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
    klocal = np.array(ek)
    return T.transpose() @ klocal @ T
#
#
# ---------------------------------------------
#
#
