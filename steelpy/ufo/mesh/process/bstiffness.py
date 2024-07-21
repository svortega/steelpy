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
# ---------------------------------------------
#
#
