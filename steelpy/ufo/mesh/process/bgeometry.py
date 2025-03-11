#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
#
# Python stdlib imports
from steelpy.ufo.mesh.process.bstiffness import B3D2_Ke, beam3D_B3D2

# package imports
import numpy as np
from numpy import sin, cos, sinh, cosh, sqrt, zeros, pi, triu

# ---------------------
#
# Geometry
#
def beam_geom(Le, s):
    """
    Calculates the element geometric matrix  
    """
    # initialize all eg elements to zero 
    eg = np.zeros((12,12), dtype=np.float32) 
    emlenz = s / (30.0 * Le)
    alpha = (s / Le) * 1.0e-6

    if abs( alpha ) < 1.0e-10:
        alpha = 1.0e-10
    #
    eg[ 0 ][ 0 ] = alpha
    eg[ 1 ][ 1 ] = 36 * emlenz
    eg[ 2 ][ 2 ] = 36 * emlenz
    eg[ 3 ][ 3 ] = alpha
    eg[ 4 ][ 4 ] = 4.0 * emlenz * Le * Le
    eg[ 5 ][ 5 ] = 4.0 * emlenz * Le * Le
    #
    eg[ 1 ][ 5 ] = 3.0 * emlenz * Le
    eg[ 2 ][ 4 ] = -3.0 * emlenz * Le
    #
    eg[ 6 ][ 6 ] = eg[ 0 ][ 0 ]
    eg[ 7 ][ 7 ] = eg[ 1 ][ 1 ]
    eg[ 8 ][ 8 ] = eg[ 2 ][ 2 ]
    eg[ 9 ][ 9 ] = eg[ 3 ][ 3 ]
    eg[ 10 ][ 10 ] = eg[ 4 ][ 4 ]
    eg[ 11 ][ 11 ] = eg[ 5 ][ 5 ]
    #
    eg[ 0 ][ 6 ] = -eg[ 0 ][ 0 ]
    eg[ 1 ][ 7 ] = -eg[ 1 ][ 1 ]
    eg[ 1 ][ 11 ] = eg[ 1 ][ 5 ]
    eg[ 2 ][ 8 ] = -eg[ 2 ][ 2 ]
    eg[ 2 ][ 10 ] = eg[ 2 ][ 4 ]
    eg[ 3 ][ 9 ] = -eg[ 3 ][ 3 ]
    eg[ 4 ][ 8 ] = -eg[ 2 ][ 4 ]
    eg[ 4 ][ 10 ] = -eg[ 4 ][ 4 ] / 4.0
    eg[ 5 ][ 7 ] = -eg[ 1 ][ 5 ]
    eg[ 5 ][ 11 ] = -eg[ 5 ][ 5 ] / 4.0
    #
    eg[ 7 ][ 11 ] = -eg[ 1 ][ 5 ]
    eg[ 8 ][ 10 ] = -eg[ 2 ][ 4 ]
    # impose the symmetry 
    for i in range( 12 ):
        for j in range( i, 12 ):
            eg[ j ][ i ] = eg[ i ][ j ]
    #
    return eg
#
#
def beam_KG(Le: float,
            Ax: float, Asy: float, Asz: float,
            Jx: float, Iy: float, Iz: float,
            Emod: float, Gmod: float,
            Fb:list,
            shear: bool = True):
    """
    3D geometric stiffness matrix for frame elements in local coordinates
    including axial, bending, shear and torsional warping effects (H.P. Gavin).
    """
    #if len(Fb) == 6:
    #    # ['x', 'y', 'rz']
    #    Fxa, Fya, Mza, Fxb, Fyb, Mzb = Fb
    #    Fb = np.array([Fxa, Fya, 0.0, 0.0, 0.0, Mza,
    #                   Fxb, Fyb, 0.0, 0.0, 0.0, Mzb])
    #else:
    #(Fxa, Fya, Fza, Mxa, Mya, Mza,
    # Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb.q_loc
    #
    P = Fb.Fx
    #
    #ax = min(Asy, Asz)
    Phiy = 0
    Phiz = 0
    if shear:
        Phiy = 12*Emod*Iz / (Gmod * Asy * Le**2)
        Phiz = 12*Emod*Iy / (Gmod * Asz * Le**2)
    #
    gk = np.zeros((12,12), dtype=np.float32) 
    #
    #gk[ 0 ][ 0 ] = 0
    gk[ 1 ][ 1 ] = (6/5 + 2*Phiz + Phiz**2) / (1 + Phiz)**2
    gk[ 2 ][ 2 ] = (6/5 + 2*Phiy + Phiy**2) / (1 + Phiy)**2
    gk[ 3 ][ 3 ] = Jx/Ax
    gk[ 4 ][ 4 ] = (2*Le**2/15 + Le**2*Phiy/6 + Le**2*Phiy**2/12) / (1+Phiy)**2
    gk[ 5 ][ 5 ] = (2*Le**2/15 + Le**2*Phiz/6 + Le**2*Phiz**2/12) / (1+Phiz)**2
    #
    gk[ 1 ][ 5 ] =  Le / (10 * (1+Phiz)**2)
    gk[ 2 ][ 4 ] = -Le / (10 * (1+Phiy)**2)
    #
    gk[4 , 10] = (-Le**2 / 30 - Le**2 * Phiy / 6 - Le**2 * Phiy**2 / 12) / (1+Phiy)**2
    gk[5 , 11] = (-Le**2 / 30 - Le**2 * Phiz / 6 - Le**2 * Phiz**2 / 12) / (1+Phiz)**2
    #
    gk[ 7 ][ 1 ] = - gk[ 1 ][ 1 ]
    #gk[ 7 ][ 1 ] = gk[ 1 ][ 1 ]
    #
    gk[ 6 ][ 6 ]  = -gk[ 0 ][ 0 ]
    gk[ 7 ][ 7 ]  = -gk[ 1 ][ 1 ]
    gk[ 8 ][ 8 ]  = gk[ 2 ][ 2 ]
    gk[ 9 ][ 9 ]  = gk[ 3 ][ 3 ]
    gk[ 10 ][ 10 ] = gk[ 4 ][ 4 ]
    gk[ 11 ][ 11 ] = gk[ 5 ][ 5 ]
    #
    #gk[1 , 7] = -gk[1 , 1]
    gk[1 , 11] = gk[1 , 5]
    gk[2 , 8] = -gk[2 , 2]
    gk[2 , 10] = gk[2 , 4]
    gk[3 , 9] = -gk[3 , 3]
    #gk[9, 3] = gk[3 , 9]
    #
    gk[4 , 8]  = -gk[2 , 4]
    gk[5 , 7]  = -gk[1 , 5]
    gk[7 , 11] = -gk[1 , 5]
    gk[8 , 10] = -gk[2 , 4]    
    #
    # impose the geometry
    gk += np.triu(gk, k=1).T
    #
    return -1 * P * gk / Le
#
#
# ============================ GeoStfBeamTimoshenko ============================ 

def GeoStfBeamTimoshenko(Le:float, 
                         Ax:float, Asy:float, Asz:float,
                         Jx:float, Iy:float, Iz:float,
                         Emod:float, Gmod:float,
                         Fb:list,
                         shear: bool = True,
                         order:int=4):
    """ 
    This method calculates the local geometric stiffness matrix 
    for an element considering the Timoshenko beam theory. 

    Can consider differents number of terms in the matrix. 
    For more precise results, the case 4 is indicated.

    Input:
    order : Number of terms considered in stiffness matrix 
    Le : Current element length                                          
    Emod : Young's modulus                                                 
    Gmod : Shear modulus                                                   
    Ax : Element cross section area                                     
    Asy : Effective shear area in y direct.;                             
    Asz : Effective shear area in z direct.;                             
    Jx : Mom. of inertia wrt "x" axis                                   
    Iy : Mom. of inertia wrt "y" axis                                   
    Iz : Mom. of inertia wrt "z" axis 
    P :    
    Fya, Fza, Mxa, Mya, Mza - Forces and moments (Node start) 
    Fyb, Fzb, Mxb, Myb, Mzb - Forces and moments (Node end)

    Return :
    et : Tangent stiffness matrix
    """
    #
    (Fxa, Fya, Fza, Mxa, Mya, Mza,
     Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb.q_loc
    #
    P = 1 * Fb.Fx
    L2 = Le*Le
    L3 = L2*Le
    L4 = L3*Le
    L5 = L4*Le
    L6 = L5*Le
    #
    Oy = 0.0
    Oz = 0.0
    if shear:
        # Shear to Flexural rig. ratio y
        Oy = (Emod*Iz)/(Gmod*Asy*L2)
        # Shear to Flexural rig. ratio z
        Oz = (Emod*Iy)/(Gmod*Asz*L2)
    #    
    Fxb2 = P*P
    Fxb3 = Fxb2*P
    niy = (1.0 + 12.0*Oy)
    niy2 = niy*niy
    niy3 = niy2*niy
    niy4 = niy3*niy
    niy5 = niy4*niy
    niz = (1.0 + 12.0*Oz)
    niz2 = niz*niz
    niz3 = niz2*niz
    niz4 = niz3*niz
    niz5 = niz4*niz
    #
    Oy2 = Oy*Oy
    Oy3 = Oy2*Oy
    Oy4 = Oy3*Oy
    Oy5 = Oy4*Oy
    Oy6 = Oy5*Oy
    Oy7 = Oy6*Oy
    Oz2 = Oz*Oz
    Oz3 = Oz2*Oz
    Oz4 = Oz3*Oz
    Oz5 = Oz4*Oz
    Oz6 = Oz5*Oz
    Oz7 = Oz6*Oz
    Iy2 = Iy * Iy
    Iy3 = Iy2 * Iy
    Iz2 = Iz * Iz
    Iz3 = Iz2 * Iz
    E2 = Emod * Emod
    E3 = E2 * Emod
    E4 = E3 * Emod
    E5 = E4 * Emod
    E6 = E5 * Emod
    #
    # coefficient K
    K = P*(Iy + Iz) / Ax
    #
    eg = np.zeros((12,12), dtype=np.float32) 
    #
    match order:
        # ========== 2 terms ========== 
        case 2:
            eg[0 , 0] = eg[6 , 6] = P / Le
            eg[1 , 1] = eg[7 , 7] = 6.0*P*(120.0*Oy2 + 20.0*Oy + 1.0) / (5.0*niy2*Le) + 12.0*P*Iz / (Ax*L3*niy2)
            eg[2 , 2] = eg[8 , 8] = 6.0*P*(120.0*Oz2 + 20.0*Oz + 1.0) / (5.0*niz2*Le) + 12.0*P*Iy / (Ax*L3*niz2)
            eg[3 , 3] = eg[9 , 9] = K / Le
            eg[4 , 4] = eg[10 , 10] = 2.0*P*Le*(90.0*Oz2 + 15.0*Oz + 1.0) / (15.0*niz2) + 4.0*P*Iy*(36.0*Oz2 + 6.0*Oz + 1.0) / (Ax*Le*niz2)
            eg[5 , 5] = eg[11 , 11] = 2.0*P*Le*(90.0*Oy2 + 15.0*Oy + 1.0) / (15.0*niy2) + 4.0*P*Iz*(36.0*Oy2 + 6.0*Oy + 1.0) / (Ax*Le*niy2)

            eg[2 , 3] = eg[3 , 2] = Mza / Le
            eg[3 , 4] = eg[4 , 3] = -Mza / 3.0 + Mzb / 6.0
            eg[5 , 6] = eg[6 , 5] = Mza / Le
            eg[8 , 9] = eg[9 , 8] = -Mzb / Le
            eg[9 , 10] = eg[10 , 9] = Mza / 6.0 - Mzb / 3.0

            eg[1 , 3] = eg[3 , 1] = Mya / Le
            eg[2 , 4] = eg[4 , 2] = -P / (10.0*niz2) - 6.0*P*Iy / (Ax*L2*niz2)
            eg[3 , 5] = eg[5 , 3] = Mya / 3.0 - Myb / 6.0
            eg[4 , 6] = eg[6 , 4] = Mya / Le
            eg[5 , 7] = eg[7 , 5] = -P / (10.0*niy2) - 6.0*P*Iz / (Ax*L2*niy2)
            eg[7 , 9] = eg[9 , 7] = -Myb / Le
            eg[8 , 10] = eg[10 , 8] = P / (10.0*niz2) + 6.0*P*Iy / (Ax*L2*niz2)
            eg[9 , 11] = eg[11 , 9] = -Mya / 6.0 + Myb / 3.0

            eg[1 , 4] = eg[4 , 1] = (Mxb / Le)*(1 / niz)
            eg[2 , 5] = eg[5 , 2] = (Mxb / Le)*(1 / niy)
            eg[4 , 7] = eg[7 , 4] = (-Mxb / Le)*(1 / niz)
            eg[5 , 8] = eg[8 , 5] = (-Mxb / Le)*(1 / niy)
            eg[7 , 10] = eg[10 , 7] = (Mxb / Le)*(1 / niz)
            eg[8 , 11] = eg[11 , 8] = (Mxb / Le)*(1 / niy)

            eg[0 , 4] = eg[4 , 0] = -Mya / Le
            eg[1 , 5] = eg[5 , 1] = P / (10.0*niy2) + 6.0*P*Iz / (Ax*L2*niy2)
            eg[3 , 7] = eg[7 , 3] = -Mya / Le
            eg[4 , 8] = eg[8 , 4] = P / (10.0*niz2) + 6.0*P*Iy / (Ax*L2*niz2)
            eg[5 , 9] = eg[9 , 5] = Mya / 6.0 + Myb / 6.0
            eg[6 , 10] = eg[10 , 6] = Myb / Le
            eg[7 , 11] = eg[11 , 7] = -P / (10.0*niy2) - 6.0*P*Iz / (Ax*L2*niy2)

            eg[0 , 5] = eg[5 , 0] = -Mza / Le
            eg[3 , 8] = eg[8 , 3] = -Mza / Le
            eg[4 , 9] = eg[9 , 4] = -Mza / 6.0 - Mzb / 6.0
            eg[5 , 10] = eg[10 , 5] = (-Mxb / 2.0)*((1 - 144 * Oy*Oz) / (niy*niz))
            eg[6 , 11] = eg[11 , 6] = Mzb / Le

            eg[0 , 6] = eg[6 , 0] = -P / Le
            eg[1 , 7] = eg[7 , 1] = -6.0*P*(120.0*Oy2 + 20.0*Oy + 1.0) / (5.0*niy2*Le) - 12.0*P*Iz / (Ax*L3*niy2)
            eg[2 , 8] = eg[8 , 2] = -6.0*P*(120.0*Oz2 + 20.0*Oz + 1.0) / (5.0*niz2*Le) - 12.0*P*Iy / (Ax*L3*niz2)
            eg[3 , 9] = eg[9 , 3] = -K / Le
            eg[4 , 10] = eg[10 , 4] = -P*Le*(360.0*Oz2 + 60.0*Oz + 1.0) / (30.0*niz2) - 2.0*P*Iy*(72.0*Oz2 + 12.0*Oz - 1.0) / (Ax*Le*niz2)
            eg[5 , 11] = eg[11 , 5] = -P*Le*(360.0*Oy2 + 60.0*Oy + 1.0) / (30.0*niy2) - 2.0*P*Iz*(72.0*Oy2 + 12.0*Oy - 1.0) / (Ax*Le*niy2)

            eg[2 , 9] = eg[9 , 2] = Mzb / Le
            eg[3 , 10] = eg[10 , 3] = -Mza / 6.0 - Mzb / 6.0
            eg[4 , 11] = eg[11 , 4] = (Mxb / 2.0)*((1 - 144 * Oy*Oz) / (niy*niz))

            eg[1 , 9] = eg[9 , 1] = Myb / Le
            eg[2 , 10] = eg[10 , 2] = -P / (10.0*niz2) - 6.0*P*Iy / (Ax*L2*niz2)
            eg[3 , 11] = eg[11 , 3] = Mya / 6.0 + Myb / 6.0

            eg[1 , 10] = eg[10 , 1] = (-Mxb / Le)*(1 / niz)
            eg[2 , 11] = eg[11 , 2] = (-Mxb / Le)*(1 / niy)

            eg[0 , 10] = eg[10 , 0] = -Myb / Le
            eg[1 , 11] = eg[11 , 1] = P / (10.0*niy2) + 6.0*P*Iz / (Ax*L2*niy2)

            eg[0 , 11] = eg[11 , 0] = -Mzb / Le

            eg[4 , 5] = eg[5 , 4] = -6 * Mxb*(Oy - Oz) / (niy*niz)
            eg[10 , 11] = eg[11 , 10] = 6 * Mxb*(Oy - Oz) / (niy*niz)

        # ========== 3 terms ========== 
        case 3:
            #eg[0 , 0] = eg[6 , 6] = 0
            eg[1 , 1] = eg[7 , 7] = ((Le*(907200 * Oy3 + 172800 * Oy2 + 8640 * Oy + 45)) / (31500 * Emod*Iz*niy4) - ((Le*(1680 * Oy2 + 180 * Oy + 1)) / (350 * niy3) + (24 * Iz*Oy) / (5 * Ax*Le*niy3)) / (Emod*Iz))*Fxb2
            eg[2 , 2] = eg[8 , 8] = ((Le*(907200 * Oz3 + 172800 * Oz2 + 8640 * Oz + 45)) / (31500 * Emod*Iy*niz4) - ((Le*(1680 * Oz2 + 180 * Oz + 1)) / (350 * niz3) + (24 * Iy*Oz) / (5 * Ax*Le*niz3)) / (Emod*Iy))*Fxb2
            #eg[3 , 3] = eg[9 , 9] = 0
            eg[4 , 4] = eg[10 , 10] = (-((L3 * (907200 * Oz4 + 241920 * Oz3 + 26460 * Oz2 + 1245 * Oz + 11)) / (3150 * niz3) + (6 * Iy*Le*Oz) / (5 * Ax*niz3)) / (Emod*Iy) + (L3 * (163296000 * Oz5 + 57153600 * Oz4 + 8391600 * Oz3 + 621000 * Oz2 + 20655 * Oz + 165)) / (94500 * Emod*Iy*niz4))*Fxb2
            eg[5 , 5] = eg[11 , 11] = (-((L3 * (907200 * Oy4 + 241920 * Oy3 + 26460 * Oy2 + 1245 * Oy + 11)) / (3150 * niy3) + (6 * Iz*Le*Oy) / (5 * Ax*niy3)) / (Emod*Iz) + (L3 * (163296000 * Oy5 + 57153600 * Oy4 + 8391600 * Oy3 + 621000 * Oy2 + 20655 * Oy + 165)) / (94500 * Emod*Iz*niy4))*Fxb2

            #eg[2 , 3] = eg[3 , 2] = 0
            eg[3 , 4] = eg[4 , 3] = ((Mza + Mzb)*((40 * L4) / (E2 * Iy2) + (100800 * L4 * Oz2) / (E2 * Iy2) + (3360 * L4 * Oz) / (E2 * Iy2))*Fxb2) / 604800 - (((1680 * L2) / (Emod*Iy) + (100800 * L2 * Oz) / (Emod*Iy))*(Mza + Mzb)*P) / 604800
            #eg[5 , 6] = eg[6 , 5] = 0
            #eg[8 , 9] = eg[9 , 8] = 0
            eg[9 , 10] = eg[10 , 9] = ((Mza + Mzb)*((40 * L4) / (E2 * Iy2) + (100800 * L4 * Oz2) / (E2 * Iy2) + (3360 * L4 * Oz) / (E2 * Iy2))*Fxb2) / 604800 - (((1680 * L2) / (Emod*Iy) + (100800 * L2 * Oz) / (Emod*Iy))*(Mza + Mzb)*P) / 604800

            #eg[1 , 3] = eg[3 , 1] = 0
            eg[2 , 4] = eg[4 , 2] = (((L2 * (1680 * Oz2 + 180 * Oz + 1)) / (700 * niz3) + (12 * Iy*Oz) / (5 * Ax*niz3)) / (Emod*Iy) - (L2 * (907200 * Oz3 + 172800 * Oz2 + 8640 * Oz + 45)) / (63000 * Emod*Iy*niz4))*Fxb2
            eg[3 , 5] = eg[5 , 3] = -((Mya + Myb)*(100800 * L4 * Oy2 + 3360 * L4 * Oy + 40 * L4)*Fxb2) / (604800 * E2 * Iz2) + ((100800 * L2 * Oy + 1680 * L2)*(Mya + Myb)*P) / (604800 * Emod*Iz)
            #eg[4 , 6] = eg[6 , 4] = 0
            eg[5 , 7] = eg[7 , 5] = (((L2 * (1680 * Oy2 + 180 * Oy + 1)) / (700 * niy3) + (12 * Iz*Oy) / (5 * Ax*niy3)) / (Emod*Iz) - (L2 * (907200 * Oy3 + 172800 * Oy2 + 8640 * Oy + 45)) / (63000 * Emod*Iz*niy4))*Fxb2
            #eg[7 , 9] = eg[9 , 7] = 0
            eg[8 , 10] = eg[10 , 8] = ((L2 * (907200 * Oz3 + 172800 * Oz2 + 8640 * Oz + 45)) / (63000 * Emod*Iy*niz4) - ((L2 * (1680 * Oz2 + 180 * Oz + 1)) / (700 * niz3) + (12 * Iy*Oz) / (5 * Ax*niz3)) / (Emod*Iy))*Fxb2
            eg[9 , 11] = eg[11 , 9] = -((Mya + Myb)*(100800 * L4 * Oy2 + 3360 * L4 * Oy + 40 * L4)*Fxb2) / (604800 * E2 * Iz2) + ((100800 * L2 * Oy + 1680 * L2)*(Mya + Myb)*P) / (604800 * Emod*Iz)

            eg[1 , 4] = eg[4 , 1] = ((Mxb*(40 * Emod*Iy*L4 + 1680 * Emod*Iy*L4 *Oz)) / (100800 * E3 * Iy3 * Le*(12 * Oy + 1)) + (Le*Mxb*(3891888000 * E4 * Iy2 * Iz2 * L2 + 249080832000 * E4 * Iy2 * Iz2 * L2 * Oy + 8966909952000 * E4 * Iy3 * Iz*L2 * Oy2 + 179338199040000 * E4 * Iy3 * Iz*L2 * Oy3 + 941525544960000 * E4 * Iy3 * Iz*L2 * Oy4 + 4296644352000 * E4 * Iy2 * Iz2 * L2 * Oy2 + 22417274880000 * E4 * Iy2 * Iz2 * L2 * Oy3 + 46702656000 * E4 * Iy3 * Iz*L2 * Oy)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*Fxb2 + (-(Le*Mxb) / (60 * Emod*Iy*(12 * Oy + 1)) - (Le*Mxb*(156920924160000 * E5 * Iy3 * Iz2 * Oy2 + 941525544960000 * E5 * Iy3 * Iz2 * Oy3 + 6538371840000 * E5 * Iy3 * Iz2 * Oy)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*P
            eg[2 , 5] = eg[5 , 2] = ((Mxb*(40 * Emod*Iz*L4 + 1680 * Emod*Iz*L4 *Oy)) / (100800 * E3 * Iz3 * Le*(12 * Oz + 1)) + (Le*Mxb*(3891888000 * E4 * Iy2 * Iz2 * L2 + 249080832000 * E4 * Iy2 * Iz2 * L2 * Oz + 8966909952000 * E4 * Iy*Iz3 * L2 * Oz2 + 179338199040000 * E4 * Iy*Iz3 * L2 * Oz3 + 941525544960000 * E4 * Iy*Iz3 * L2 * Oz4 + 4296644352000 * E4 * Iy2 * Iz2 * L2 * Oz2 + 22417274880000 * E4 * Iy2 * Iz2 * L2 * Oz3 + 46702656000 * E4 * Iy*Iz3 * L2 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niz4))*Fxb2 + (-(Le*Mxb) / (60 * Emod*Iz*(12 * Oz + 1)) - (Le*Mxb*(156920924160000 * E5 * Iy2 * Iz3 * Oz2 + 941525544960000 * E5 * Iy2 * Iz3 * Oz3 + 6538371840000 * E5 * Iy2 * Iz3 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niz4))*P
            eg[4 , 7] = eg[7 , 4] = (-(Mxb*(40 * Emod*Iy*L4 + 1680 * Emod*Iy*L4 * Oz)) / (100800 * E3 * Iy3 * Le*(12 * Oy + 1)) - (Le*Mxb*(3891888000 * E4 * Iy2 * Iz2 * L2 + 249080832000 * E4 * Iy2 * Iz2 * L2 * Oy + 8966909952000 * E4 * Iy3 * Iz*L2 * Oy2 + 179338199040000 * E4 * Iy3 * Iz*L2 * Oy3 + 941525544960000 * E4 * Iy3 * Iz*L2 * Oy4 + 4296644352000 * E4 * Iy2 * Iz2 * L2 * Oy2 + 22417274880000 * E4 * Iy2 * Iz2 * L2 * Oy3 + 46702656000 * E4 * Iy3 * Iz*L2 * Oy)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*Fxb2 + ((Le*Mxb) / (60 * Emod*Iy*(12 * Oy + 1)) + (Le*Mxb*(156920924160000 * E5 * Iy3 * Iz2 * Oy2 + 941525544960000 * E5 * Iy3 * Iz2 * Oy3 + 6538371840000 * E5 * Iy3 * Iz2 * Oy)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*P
            eg[5 , 8] = eg[8 , 5] = (-(Mxb*(40 * Emod*Iz*L4 + 1680 * Emod*Iz*L4 * Oy)) / (100800 * E3 * Iz3 * Le*(12 * Oz + 1)) - (Le*Mxb*(3891888000 * E4 * Iy2 * Iz2 * L2 + 249080832000 * E4 * Iy2 * Iz2 * L2 * Oz + 8966909952000 * E4 * Iy*Iz3 * L2 * Oz2 + 179338199040000 * E4 * Iy*Iz3 * L2 * Oz3 + 941525544960000 * E4 * Iy*Iz3 * L2 * Oz4 + 4296644352000 * E4 * Iy2 * Iz2 * L2 * Oz2 + 22417274880000 * E4 * Iy2 * Iz2 * L2 * Oz3 + 46702656000 * E4 * Iy*Iz3 * L2 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niz4))*Fxb2 + ((Le*Mxb) / (60 * Emod*Iz*(12 * Oz + 1)) + (Le*Mxb*(156920924160000 * E5 * Iy2 * Iz3 * Oz2 + 941525544960000 * E5 * Iy2 * Iz3 * Oz3 + 6538371840000 * E5 * Iy2 * Iz3 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niz4))*P
            eg[7 , 10] = eg[10 , 7] = ((Mxb*(40 * Emod*Iy*L4 + 1680 * Emod*Iy*L4 * Oz)) / (100800 * E3 * Iy3 * Le*(12 * Oy + 1)) + (Le*Mxb*(3891888000 * E4 * Iy2 * Iz2 * L2 + 249080832000 * E4 * Iy2 * Iz2 * L2 * Oy + 8966909952000 * E4 * Iy3 * Iz*L2 * Oy2 + 179338199040000 * E4 * Iy3 * Iz*L2 * Oy3 + 941525544960000 * E4 * Iy3 * Iz*L2 * Oy4 + 4296644352000 * E4 * Iy2 * Iz2 * L2 * Oy2 + 22417274880000 * E4 * Iy2 * Iz2 * L2 * Oy3 + 46702656000 * E4 * Iy3 * Iz*L2 * Oy)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*Fxb2 + (-(Le*Mxb) / (60 * Emod*Iy*(12 * Oy + 1)) - (Le*Mxb*(156920924160000 * E5 * Iy3 * Iz2 * Oy2 + 941525544960000 * E5 * Iy3 * Iz2 * Oy3 + 6538371840000 * E5 * Iy3 * Iz2 * Oy)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*P
            eg[8 , 11] = eg[11 , 8] = ((Mxb*(40 * Emod*Iz*L4 + 1680 * Emod*Iz*L4 * Oy)) / (100800 * E3 * Iz3 * Le*(12 * Oz + 1)) + (Le*Mxb*(3891888000 * E4 * Iy2 * Iz2 * L2 + 249080832000 * E4 * Iy2 * Iz2 * L2 * Oz + 8966909952000 * E4 * Iy*Iz3 * L2 * Oz2 + 179338199040000 * E4 * Iy*Iz3 * L2 * Oz3 + 941525544960000 * E4 * Iy*Iz3 * L2 * Oz4 + 4296644352000 * E4 * Iy2 * Iz2 * L2 * Oz2 + 22417274880000 * E4 * Iy2 * Iz2 * L2 * Oz3 + 46702656000 * E4 * Iy*Iz3 * L2 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niz4))*Fxb2 + (-(Le*Mxb) / (60 * Emod*Iz*(12 * Oz + 1)) - (Le*Mxb*(156920924160000 * E5 * Iy2 * Iz3 * Oz2 + 941525544960000 * E5 * Iy2 * Iz3 * Oz3 + 6538371840000 * E5 * Iy2 * Iz3 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niz4))*P

            #eg[0 , 4] = eg[4 , 0] = 0
            eg[1 , 5] = eg[5 , 1] = ((L2 * (907200 * Oy3 + 172800 * Oy2 + 8640 * Oy + 45)) / (63000 * Emod*Iz*niy4) - ((L2 * (1680 * Oy2 + 180 * Oy + 1)) / (700 * niy3) + (12 * Iz*Oy) / (5 * Ax*niy3)) / (Emod*Iz))*Fxb2
            #eg[3 , 7] = eg[7 , 3] = 0
            eg[4 , 8] = eg[8 , 4] = ((L2 * (907200 * Oz3 + 172800 * Oz2 + 8640 * Oz + 45)) / (63000 * Emod*Iy*niz4) - ((L2 * (1680 * Oz2 + 180 * Oz + 1)) / (700 * niz3) + (12 * Iy*Oz) / (5 * Ax*niz3)) / (Emod*Iy))*Fxb2
            eg[5 , 9] = eg[9 , 5] = ((Mya + Myb)*(100800 * L4 * Oy2 + 3360 * L4 * Oy + 40 * L4)*Fxb2) / (604800 * E2 * Iz2) - ((100800 * L2 * Oy + 1680 * L2)*(Mya + Myb)*P) / (604800 * Emod*Iz)
            eg[6 , 10] = eg[10 , 6] = 0
            eg[7 , 11] = eg[11 , 7] = (((L2 * (1680 * Oy2 + 180 * Oy + 1)) / (700 * niy3) + (12 * Iz*Oy) / (5 * Ax*niy3)) / (Emod*Iz) - (L2 * (907200 * Oy3 + 172800 * Oy2 + 8640 * Oy + 45)) / (63000 * Emod*Iz*niy4))*Fxb2

            #eg[0 , 5] = eg[5 , 0] = 0
            #eg[3 , 8] = eg[8 , 3] = 0
            eg[4 , 9] = eg[9 , 4] = -((Mza + Mzb)*((40 * L4) / (E2 * Iy2) + (100800 * L4 * Oz2) / (E2 * Iy2) + (3360 * L4 * Oz) / (E2 * Iy2))*Fxb2) / 604800 + (((1680 * L2) / (Emod*Iy) + (100800 * L2 * Oz) / (Emod*Iy))*(Mza + Mzb)*P) / 604800
            eg[5 , 10] = eg[10 , 5] = -(L4 * Mxb*(52254720 * Iy2 * Oy4 * Oz2 + 8709120 * Iy2 * Oy4 * Oz + 362880 * Iy2 * Oy4 + 52254720 * Iy2 * Oy3 * Oz3 + 27371520 * Iy2 * Oy3 * Oz2 + 3473280 * Iy2 * Oy3 * Oz + 129600 * Iy2 * Oy3 + 5598720 * Iy2 * Oy2 * Oz3 + 2799360 * Iy2 * Oy2 * Oz2 + 349920 * Iy2 * Oy2 * Oz + 12960 * Iy2 * Oy2 + 31104 * Iy2 * Oy*Oz3 + 63936 * Iy2 * Oy*Oz2 + 10008 * Iy2 * Oy*Oz + 408 * Iy2 * Oy + 720 * Iy2 * Oz2 + 120 * Iy2 * Oz + 5 * Iy2 + 1244160 * Iy*Iz*Oy3 * Oz2 + 134784 * Iy*Iz*Oy3 * Oz + 2592 * Iy*Iz*Oy3 + 1244160 * Iy*Iz*Oy2 * Oz3 + 622080 * Iy*Iz*Oy2 * Oz2 + 59616 * Iy*Iz*Oy2 * Oz + 1368 * Iy*Iz*Oy2 + 134784 * Iy*Iz*Oy*Oz3 + 59616 * Iy*Iz*Oy*Oz2 + 5616 * Iy*Iz*Oy*Oz + 132 * Iy*Iz*Oy + 2592 * Iy*Iz*Oz3 + 1368 * Iy*Iz*Oz2 + 132 * Iy*Iz*Oz + 3 * Iy*Iz + 52254720 * Iz2 * Oy3 * Oz3 + 5598720 * Iz2 * Oy3 * Oz2 + 31104 * Iz2 * Oy3 * Oz + 52254720 * Iz2 * Oy2 * Oz4 + 27371520 * Iz2 * Oy2 * Oz3 + 2799360 * Iz2 * Oy2 * Oz2 + 63936 * Iz2 * Oy2 * Oz + 720 * Iz2 * Oy2 + 8709120 * Iz2 * Oy*Oz4 + 3473280 * Iz2 * Oy*Oz3 + 349920 * Iz2 * Oy*Oz2 + 10008 * Iz2 * Oy*Oz + 120 * Iz2 * Oy + 362880 * Iz2 * Oz4 + 129600 * Iz2 * Oz3 + 12960 * Iz2 * Oz2 + 408 * Iz2 * Oz + 5 * Iz2)*Fxb2) / (25200 * E2 * Iy2 * Iz2 * niy3 * niz3) + (L2 * Mxb*(Iy + Iz + 144 * Iy*Oy2 + 144 * Iz*Oz2 + 36 * Iy*Oy + 12 * Iy*Oz + 12 * Iz*Oy + 36 * Iz*Oz + 576 * Iy*Oy*Oz + 576 * Iz*Oy*Oz + 1728 * Iy*Oy*Oz2 + 1728 * Iy*Oy2 * Oz + 1728 * Iz*Oy*Oz2 + 1728 * Iz*Oy2 * Oz)*P) / (120 * Emod*Iy*Iz*niy2 * niz2)
            #eg[6 , 11] = eg[11 , 6] = 0

            #eg[0 , 6] = eg[6 , 0] = 0
            eg[1 , 7] = eg[7 , 1] = (((Le*(1680 * Oy2 + 180 * Oy + 1)) / (350 * niy3) + (24 * Iz*Oy) / (5 * Ax*Le*niy3)) / (Emod*Iz) - (Le*(907200 * Oy3 + 172800 * Oy2 + 8640 * Oy + 45)) / (31500 * Emod*Iz*niy4))*Fxb2
            eg[2 , 8] = eg[8 , 2] = (((Le*(1680 * Oz2 + 180 * Oz + 1)) / (350 * niz3) + (24 * Iy*Oz) / (5 * Ax*Le*niz3)) / (Emod*Iy) - (Le*(907200 * Oz3 + 172800 * Oz2 + 8640 * Oz + 45)) / (31500 * Emod*Iy*niz4))*Fxb2
            #eg[3 , 9] = eg[9 , 3] = 0
            eg[4 , 10] = eg[10 , 4] = (((L3 * (1814400 * Oz4 + 483840 * Oz3 + 37800 * Oz2 + 870 * Oz + 13)) / (6300 * niz3) - (6 * Iy*Le*Oz) / (5 * Ax*niz3)) / (Emod*Iy) - (L3 * (326592000 * Oz5 + 114307200 * Oz4 + 14061600 * Oz3 + 723600 * Oz2 + 15390 * Oz + 195)) / (189000 * Emod*Iy*niz4))*Fxb2
            eg[5 , 11] = eg[11 , 5] = (((L3 * (1814400 * Oy4 + 483840 * Oy3 + 37800 * Oy2 + 870 * Oy + 13)) / (6300 * niy3) - (6 * Iz*Le*Oy) / (5 * Ax*niy3)) / (Emod*Iz) - (L3 * (326592000 * Oy5 + 114307200 * Oy4 + 14061600 * Oy3 + 723600 * Oy2 + 15390 * Oy + 195)) / (189000 * Emod*Iz*niy4))*Fxb2

            #eg[2 , 9] = eg[9 , 2] = 0
            eg[3 , 10] = eg[10 , 3] = -((Mza + Mzb)*((40 * L4) / (E2 * Iy2) + (100800 * L4 * Oz2) / (E2 * Iy2) + (3360 * L4 * Oz) / (E2 * Iy2))*Fxb2) / 604800 + (((1680 * L2) / (Emod*Iy) + (100800 * L2 * Oz) / (Emod*Iy))*(Mza + Mzb)*P) / 604800
            eg[4 , 11] = eg[11 , 4] = (L4 * Mxb*(52254720 * Iy2 * Oy4 * Oz2 + 8709120 * Iy2 * Oy4 * Oz + 362880 * Iy2 * Oy4 + 52254720 * Iy2 * Oy3 * Oz3 + 27371520 * Iy2 * Oy3 * Oz2 + 3473280 * Iy2 * Oy3 * Oz + 129600 * Iy2 * Oy3 + 5598720 * Iy2 * Oy2 * Oz3 + 2799360 * Iy2 * Oy2 * Oz2 + 349920 * Iy2 * Oy2 * Oz + 12960 * Iy2 * Oy2 + 31104 * Iy2 * Oy*Oz3 + 63936 * Iy2 * Oy*Oz2 + 10008 * Iy2 * Oy*Oz + 408 * Iy2 * Oy + 720 * Iy2 * Oz2 + 120 * Iy2 * Oz + 5 * Iy2 + 1244160 * Iy*Iz*Oy3 * Oz2 + 134784 * Iy*Iz*Oy3 * Oz + 2592 * Iy*Iz*Oy3 + 1244160 * Iy*Iz*Oy2 * Oz3 + 622080 * Iy*Iz*Oy2 * Oz2 + 59616 * Iy*Iz*Oy2 * Oz + 1368 * Iy*Iz*Oy2 + 134784 * Iy*Iz*Oy*Oz3 + 59616 * Iy*Iz*Oy*Oz2 + 5616 * Iy*Iz*Oy*Oz + 132 * Iy*Iz*Oy + 2592 * Iy*Iz*Oz3 + 1368 * Iy*Iz*Oz2 + 132 * Iy*Iz*Oz + 3 * Iy*Iz + 52254720 * Iz2 * Oy3 * Oz3 + 5598720 * Iz2 * Oy3 * Oz2 + 31104 * Iz2 * Oy3 * Oz + 52254720 * Iz2 * Oy2 * Oz4 + 27371520 * Iz2 * Oy2 * Oz3 + 2799360 * Iz2 * Oy2 * Oz2 + 63936 * Iz2 * Oy2 * Oz + 720 * Iz2 * Oy2 + 8709120 * Iz2 * Oy*Oz4 + 3473280 * Iz2 * Oy*Oz3 + 349920 * Iz2 * Oy*Oz2 + 10008 * Iz2 * Oy*Oz + 120 * Iz2 * Oy + 362880 * Iz2 * Oz4 + 129600 * Iz2 * Oz3 + 12960 * Iz2 * Oz2 + 408 * Iz2 * Oz + 5 * Iz2)*Fxb2) / (25200 * E2 * Iy2 * Iz2 * niy3 * niz3) - (L2 * Mxb*(Iy + Iz + 144 * Iy*Oy2 + 144 * Iz*Oz2 + 36 * Iy*Oy + 12 * Iy*Oz + 12 * Iz*Oy + 36 * Iz*Oz + 576 * Iy*Oy*Oz + 576 * Iz*Oy*Oz + 1728 * Iy*Oy*Oz2 + 1728 * Iy*Oy2 * Oz + 1728 * Iz*Oy*Oz2 + 1728 * Iz*Oy2 * Oz)*P) / (120 * Emod*Iy*Iz*niy2 * niz2)

            #eg[1 , 9] = eg[9 , 1] = 0
            eg[2 , 10] = eg[10 , 2] = (((L2 * (1680 * Oz2 + 180 * Oz + 1)) / (700 * niz3) + (12 * Iy*Oz) / (5 * Ax*niz3)) / (Emod*Iy) - (L2 * (907200 * Oz3 + 172800 * Oz2 + 8640 * Oz + 45)) / (63000 * Emod*Iy*niz4))*Fxb2
            eg[3 , 11] = eg[11 , 3] = ((Mya + Myb)*(100800 * L4 * Oy2 + 3360 * L4 * Oy + 40 * L4)*Fxb2) / (604800 * E2 * Iz2) - ((100800 * L2 * Oy + 1680 * L2)*(Mya + Myb)*P) / (604800 * Emod*Iz)

            eg[1 , 10] = eg[10 , 1] = (-(Mxb*(40 * Emod*Iy*L4 + 1680 * Emod*Iy*L4 * Oz)) / (100800 * E3 * Iy3 * Le*(12 * Oy + 1)) - (Le*Mxb*(3891888000 * E4 * Iy2 * Iz2 * L2 + 249080832000 * E4 * Iy2 * Iz2 * L2 * Oy + 8966909952000 * E4 * Iy3 * Iz*L2 * Oy2 + 179338199040000 * E4 * Iy3 * Iz*L2 * Oy3 + 941525544960000 * E4 * Iy3 * Iz*L2 * Oy4 + 4296644352000 * E4 * Iy2 * Iz2 * L2 * Oy2 + 22417274880000 * E4 * Iy2 * Iz2 * L2 * Oy3 + 46702656000 * E4 * Iy3 * Iz*L2 * Oy)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*Fxb2 + ((Le*Mxb) / (60 * Emod*Iy*(12 * Oy + 1)) + (Le*Mxb*(156920924160000 * E5 * Iy3 * Iz2 * Oy2 + 941525544960000 * E5 * Iy3 * Iz2 * Oy3 + 6538371840000 * E5 * Iy3 * Iz2 * Oy)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*P
            eg[2 , 11] = eg[11 , 2] = (-(Mxb*(40 * Emod*Iz*L4 + 1680 * Emod*Iz*L4 * Oy)) / (100800 * E3 * Iz3 * Le*(12 * Oz + 1)) - (Le*Mxb*(3891888000 * E4 * Iy2 * Iz2 * L2 + 249080832000 * E4 * Iy2 * Iz2 * L2 * Oz + 8966909952000 * E4 * Iy*Iz3 * L2 * Oz2 + 179338199040000 * E4 * Iy*Iz3 * L2 * Oz3 + 941525544960000 * E4 * Iy*Iz3 * L2 * Oz4 + 4296644352000 * E4 * Iy2 * Iz2 * L2 * Oz2 + 22417274880000 * E4 * Iy2 * Iz2 * L2 * Oz3 + 46702656000 * E4 * Iy*Iz3 * L2 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niz4))*Fxb2 + ((Le*Mxb) / (60 * Emod*Iz*(12 * Oz + 1)) + (Le*Mxb*(156920924160000 * E5 * Iy2 * Iz3 * Oz2 + 941525544960000 * E5 * Iy2 * Iz3 * Oz3 + 6538371840000 * E5 * Iy2 * Iz3 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niz4))*P

            #eg[0 , 10] = eg[10 , 0] = 0
            eg[1 , 11] = eg[11 , 1] = ((L2 * (907200 * Oy3 + 172800 * Oy2 + 8640 * Oy + 45)) / (63000 * Emod*Iz*niy4) - ((L2 * (1680 * Oy2 + 180 * Oy + 1)) / (700 * niy3) + (12 * Iz*Oy) / (5 * Ax*niy3)) / (Emod*Iz))*Fxb2

            #eg[0 , 11] = eg[11 , 0] = 0

            eg[4 , 5] = eg[5 , 4] = -(L4 * Mxb*(52254720 * Iy2 * Oy4 * Oz2 + 8709120 * Iy2 * Oy4 * Oz + 362880 * Iy2 * Oy4 - 52254720 * Iy2 * Oy3 * Oz3 + 1244160 * Iy2 * Oy3 * Oz2 + 1296000 * Iy2 * Oy3 * Oz + 69120 * Iy2 * Oy3 - 5598720 * Iy2 * Oy2 * Oz3 + 116640 * Iy2 * Oy2 * Oz + 6480 * Iy2 * Oy2 - 31104 * Iy2 * Oy*Oz3 + 48384 * Iy2 * Oy*Oz2 + 8712 * Iy2 * Oy*Oz + 372 * Iy2 * Oy + 720 * Iy2 * Oz2 + 120 * Iy2 * Oz + 5 * Iy2 + 1244160 * Iy*Iz*Oy3 * Oz2 + 134784 * Iy*Iz*Oy3 * Oz + 2592 * Iy*Iz*Oy3 - 1244160 * Iy*Iz*Oy2 * Oz3 + 7776 * Iy*Iz*Oy2 * Oz - 72 * Iy*Iz*Oy2 - 134784 * Iy*Iz*Oy*Oz3 - 7776 * Iy*Iz*Oy*Oz2 - 24 * Iy*Iz*Oy - 2592 * Iy*Iz*Oz3 + 72 * Iy*Iz*Oz2 + 24 * Iy*Iz*Oz + 52254720 * Iz2 * Oy3 * Oz3 + 5598720 * Iz2 * Oy3 * Oz2 + 31104 * Iz2 * Oy3 * Oz - 52254720 * Iz2 * Oy2 * Oz4 - 1244160 * Iz2 * Oy2 * Oz3 - 48384 * Iz2 * Oy2 * Oz - 720 * Iz2 * Oy2 - 8709120 * Iz2 * Oy*Oz4 - 1296000 * Iz2 * Oy*Oz3 - 116640 * Iz2 * Oy*Oz2 - 8712 * Iz2 * Oy*Oz - 120 * Iz2 * Oy - 362880 * Iz2 * Oz4 - 69120 * Iz2 * Oz3 - 6480 * Iz2 * Oz2 - 372 * Iz2 * Oz - 5 * Iz2)*Fxb2) / (25200 * E2 * Iy2 * Iz2 * niy3 * niz3) + (L2 * Mxb*(Iy - Iz + 144 * Iy*Oy2 - 144 * Iz*Oz2 + 12 * Iy*Oy + 12 * Iy*Oz - 12 * Iz*Oy - 12 * Iz*Oz - 1728 * Iy*Oy*Oz2 + 1728 * Iy*Oy2 * Oz - 1728 * Iz*Oy*Oz2 + 1728 * Iz*Oy2 * Oz)*P) / (120 * Emod*Iy*Iz*niy2 * niz2)
            eg[10 , 11] = eg[11 , 10] = (L4 * Mxb*(52254720 * Iy2 * Oy4 * Oz2 + 8709120 * Iy2 * Oy4 * Oz + 362880 * Iy2 * Oy4 - 52254720 * Iy2 * Oy3 * Oz3 + 1244160 * Iy2 * Oy3 * Oz2 + 1296000 * Iy2 * Oy3 * Oz + 69120 * Iy2 * Oy3 - 5598720 * Iy2 * Oy2 * Oz3 + 116640 * Iy2 * Oy2 * Oz + 6480 * Iy2 * Oy2 - 31104 * Iy2 * Oy*Oz3 + 48384 * Iy2 * Oy*Oz2 + 8712 * Iy2 * Oy*Oz + 372 * Iy2 * Oy + 720 * Iy2 * Oz2 + 120 * Iy2 * Oz + 5 * Iy2 + 1244160 * Iy*Iz*Oy3 * Oz2 + 134784 * Iy*Iz*Oy3 * Oz + 2592 * Iy*Iz*Oy3 - 1244160 * Iy*Iz*Oy2 * Oz3 + 7776 * Iy*Iz*Oy2 * Oz - 72 * Iy*Iz*Oy2 - 134784 * Iy*Iz*Oy*Oz3 - 7776 * Iy*Iz*Oy*Oz2 - 24 * Iy*Iz*Oy - 2592 * Iy*Iz*Oz3 + 72 * Iy*Iz*Oz2 + 24 * Iy*Iz*Oz + 52254720 * Iz2 * Oy3 * Oz3 + 5598720 * Iz2 * Oy3 * Oz2 + 31104 * Iz2 * Oy3 * Oz - 52254720 * Iz2 * Oy2 * Oz4 - 1244160 * Iz2 * Oy2 * Oz3 - 48384 * Iz2 * Oy2 * Oz - 720 * Iz2 * Oy2 - 8709120 * Iz2 * Oy*Oz4 - 1296000 * Iz2 * Oy*Oz3 - 116640 * Iz2 * Oy*Oz2 - 8712 * Iz2 * Oy*Oz - 120 * Iz2 * Oy - 362880 * Iz2 * Oz4 - 69120 * Iz2 * Oz3 - 6480 * Iz2 * Oz2 - 372 * Iz2 * Oz - 5 * Iz2)*Fxb2) / (25200 * E2 * Iy2 * Iz2 * niy3 * niz3) - (L2 * Mxb*(Iy - Iz + 144 * Iy*Oy2 - 144 * Iz*Oz2 + 12 * Iy*Oy + 12 * Iy*Oz - 12 * Iz*Oy - 12 * Iz*Oz - 1728 * Iy*Oy*Oz2 + 1728 * Iy*Oy2 * Oz - 1728 * Iz*Oy*Oz2 + 1728 * Iz*Oy2 * Oz)*P) / (120 * Emod*Iy*Iz*niy2 * niz2)

        # ========== 4 terms ========== */
        case 4:
            #eg[0 , 0] = eg[6 , 6] = 0
            eg[1 , 1] = eg[7 , 7] = (((L3 * (21772800 * Oy5 + 6480000 * Oy4 + 665280 * Oy3 + 25920 * Oy2 + 252 * Oy + 1)) / (21000 * niy5) + (Iz*Le*(483840 * Oy4 + 97920 * Oy3 + 5376 * Oy2 + 60 * Oy + 1)) / (700 * Ax*niy5)) / (E2 * Iz2) - (Le*(1814400 * L2 * Oy4 + 388800 * L2 * Oy3 + 23040 * L2 * Oy2 + 240 * L2 * Oy + L2)) / (31500 * E2 * Iz2 * niy4))*Fxb3
            eg[2 , 2] = eg[8 , 8] = (((L3 * (21772800 * Oz5 + 6480000 * Oz4 + 665280 * Oz3 + 25920 * Oz2 + 252 * Oz + 1)) / (21000 * niz5) + (Iy*Le*(483840 * Oz4 + 97920 * Oz3 + 5376 * Oz2 + 60 * Oz + 1)) / (700 * Ax*niz5)) / (E2 * Iy2) - (Le*(1814400 * L2 * Oz4 + 388800 * L2 * Oz3 + 23040 * L2 * Oz2 + 240 * L2 * Oz + L2)) / (31500 * E2 * Iy2 * niz4))*Fxb3
            #eg[3 , 3] = eg[9 , 9] = 0
            eg[4 , 4] = eg[10 , 10] = (((L5 * (3919104000 * Oz7 + 1763596800 * Oz6 + 344476800 * Oz5 + 37260000 * Oz4 + 2307960 * Oz3 + 75690 * Oz2 + 1089 * Oz + 7)) / (63000 * niz5) + (Iy*L3 * (2177280 * Oz5 + 1995840 * Oz4 + 371520 * Oz3 + 24696 * Oz2 + 660 * Oz + 11)) / (6300 * Ax*niz5)) / (E2 * Iy2) - (L3 * (326592000 * L2 * Oz6 + 119750400 * L2 * Oz5 + 18727200 * L2 * Oz4 + 1544400 * L2 * Oz3 + 63630 * L2 * Oz2 + 1005 * L2 * Oz + 7 * L2)) / (94500 * E2 * Iy2 * niz4))*Fxb3
            eg[5 , 5] = eg[11 , 11] = (((L5 * (3919104000 * Oy7 + 1763596800 * Oy6 + 344476800 * Oy5 + 37260000 * Oy4 + 2307960 * Oy3 + 75690 * Oy2 + 1089 * Oy + 7)) / (63000 * niy5) + (Iz*L3 * (2177280 * Oy5 + 1995840 * Oy4 + 371520 * Oy3 + 24696 * Oy2 + 660 * Oy + 11)) / (6300 * Ax*niy5)) / (E2 * Iz2) - (L3 * (326592000 * L2 * Oy6 + 119750400 * L2 * Oy5 + 18727200 * L2 * Oy4 + 1544400 * L2 * Oy3 + 63630 * L2 * Oy2 + 1005 * L2 * Oy + 7 * L2)) / (94500 * E2 * Iz2 * niy4))*Fxb3

            #eg[2 , 3] = eg[3 , 2] = 0
            eg[3 , 4] = eg[4 , 3] = -((Mza + Mzb)*(L6 / (E3 * Iy3) + (5040 * L6 * Oz2) / (E3 * Iy3) + (100800 * L6 * Oz3) / (E3 * Iy3) + (120 * L6 * Oz) / (E3 * Iy3))*Fxb3) / 604800
            #eg[5 , 6] = eg[6 , 5] = 0
            #eg[8 , 9] = eg[9 , 8] = 0
            eg[9 , 10] = eg[10 , 9] = -((Mza + Mzb)*(L6 / (E3 * Iy3) + (5040 * L6 * Oz2) / (E3 * Iy3) + (100800 * L6 * Oz3) / (E3 * Iy3) + (120 * L6 * Oz) / (E3 * Iy3))*Fxb3) / 604800

            #eg[1 , 3] = eg[3 , 1] = 0
            eg[2 , 4] = eg[4 , 2] = (-((L4 * (21772800 * Oz5 + 6480000 * Oz4 + 665280 * Oz3 + 25920 * Oz2 + 252 * Oz + 1)) / (42000 * niz5) + (Iy*L2 * (483840 * Oz4 + 97920 * Oz3 + 5376 * Oz2 + 60 * Oz + 1)) / (1400 * Ax*niz5)) / (E2 * Iy2) + (L2 * (1814400 * L2 * Oz4 + 388800 * L2 * Oz3 + 23040 * L2 * Oz2 + 240 * L2 * Oz + L2)) / (63000 * E2 * Iy2 * niz4))*Fxb3
            eg[3 , 5] = eg[5 , 3] = ((Mya + Myb)*(100800 * L6 * Oy3 + 5040 * L6 * Oy2 + 120 * L6 * Oy + L6)*Fxb3) / (604800 * E3 * Iz3)
            #eg[4 , 6] = eg[6 , 4] = 0
            eg[5 , 7] = eg[7 , 5] = (-((L4 * (21772800 * Oy5 + 6480000 * Oy4 + 665280 * Oy3 + 25920 * Oy2 + 252 * Oy + 1)) / (42000 * niy5) + (Iz*L2 * (483840 * Oy4 + 97920 * Oy3 + 5376 * Oy2 + 60 * Oy + 1)) / (1400 * Ax*niy5)) / (E2 * Iz2) + (L2 * (1814400 * L2 * Oy4 + 388800 * L2 * Oy3 + 23040 * L2 * Oy2 + 240 * L2 * Oy + L2)) / (63000 * E2 * Iz2 * niy4))*Fxb3
            #eg[7 , 9] = eg[9 , 7] = 0
            eg[8 , 10] = eg[10 , 8] = (((L4 * (21772800 * Oz5 + 6480000 * Oz4 + 665280 * Oz3 + 25920 * Oz2 + 252 * Oz + 1)) / (42000 * niz5) + (Iy*L2 * (483840 * Oz4 + 97920 * Oz3 + 5376 * Oz2 + 60 * Oz + 1)) / (1400 * Ax*niz5)) / (E2 * Iy2) - (L2 * (1814400 * L2 * Oz4 + 388800 * L2 * Oz3 + 23040 * L2 * Oz2 + 240 * L2 * Oz + L2)) / (63000 * E2 * Iy2 * niz4))*Fxb3
            eg[9 , 11] = eg[11 , 9] = ((Mya + Myb)*(100800 * L6 * Oy3 + 5040 * L6 * Oy2 + 120 * L6 * Oy + L6)*Fxb3) / (604800 * E3 * Iz3)

            eg[1 , 4] = eg[4 , 1] = (-(Mxb*(1680 * L6 * Oz2 + 80 * L6 * Oz + L6)) / (100800 * E3 * Iy3 * Le*(12 * Oy + 1)) - (Le*Mxb*(124540416000 * E3 * Iy3 * L4 * Oy2 + 11955879936000 * E3 * Iy3 * L4 * Oy3 + 201755473920000 * E3 * Iy3 * L4 * Oy4 + 941525544960000 * E3 * Iy3 * L4 * Oy5 + 108108000 * E3 * Iy*Iz2 * L4 + 43243200 * E3 * Iy2 * Iz*L4 + 518918400 * E3 * Iy3 * L4 * Oy + 108972864000 * E3 * Iy*Iz2 * L4 * Oy2 + 326918592000 * E3 * Iy2 * Iz*L4 * Oy2 + 560431872000 * E3 * Iy*Iz2 * L4 * Oy3 + 4857076224000 * E3 * Iy2 * Iz*L4 * Oy3 + 22417274880000 * E3 * Iy2 * Iz*L4 * Oy4 + 6486480000 * E3 * Iy*Iz2 * L4 * Oy + 7005398400 * E3 * Iy2 * Iz*L4 * Oy + 3891888000 * E3 * Iy*Iz2 * L4 * Oz + 249080832000 * E3 * Iy*Iz2 * L4 * Oy*Oz + 4296644352000 * E3 * Iy*Iz2 * L4 * Oy2 * Oz + 22417274880000 * E3 * Iy*Iz2 * L4 * Oy3 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*Fxb3
            eg[2 , 5] = eg[5 , 2] = (-(Mxb*(1680 * L6 * Oy2 + 80 * L6 * Oy + L6)) / (100800 * E3 * Iz3 * Le*(12 * Oz + 1)) - (Le*Mxb*(124540416000 * E3 * Iz3 * L4 * Oz2 + 11955879936000 * E3 * Iz3 * L4 * Oz3 + 201755473920000 * E3 * Iz3 * L4 * Oz4 + 941525544960000 * E3 * Iz3 * L4 * Oz5 + 43243200 * E3 * Iy*Iz2 * L4 + 108108000 * E3 * Iy2 * Iz*L4 + 518918400 * E3 * Iz3 * L4 * Oz + 326918592000 * E3 * Iy*Iz2 * L4 * Oz2 + 108972864000 * E3 * Iy2 * Iz*L4 * Oz2 + 4857076224000 * E3 * Iy*Iz2 * L4 * Oz3 + 560431872000 * E3 * Iy2 * Iz*L4 * Oz3 + 22417274880000 * E3 * Iy*Iz2 * L4 * Oz4 + 3891888000 * E3 * Iy2 * Iz*L4 * Oy + 7005398400 * E3 * Iy*Iz2 * L4 * Oz + 6486480000 * E3 * Iy2 * Iz*L4 * Oz + 249080832000 * E3 * Iy2 * Iz*L4 * Oy*Oz + 4296644352000 * E3 * Iy2 * Iz*L4 * Oy*Oz2 + 22417274880000 * E3 * Iy2 * Iz*L4 * Oy*Oz3)) / (32691859200000 * E6 * Iy3 * Iz3 * niz4))*Fxb3
            eg[4 , 7] = eg[7 , 4] = ((Mxb*(1680 * L6 * Oz2 + 80 * L6 * Oz + L6)) / (100800 * E3 * Iy3 * Le*(12 * Oy + 1)) + (Le*Mxb*(124540416000 * E3 * Iy3 * L4 * Oy2 + 11955879936000 * E3 * Iy3 * L4 * Oy3 + 201755473920000 * E3 * Iy3 * L4 * Oy4 + 941525544960000 * E3 * Iy3 * L4 * Oy5 + 108108000 * E3 * Iy*Iz2 * L4 + 43243200 * E3 * Iy2 * Iz*L4 + 518918400 * E3 * Iy3 * L4 * Oy + 108972864000 * E3 * Iy*Iz2 * L4 * Oy2 + 326918592000 * E3 * Iy2 * Iz*L4 * Oy2 + 560431872000 * E3 * Iy*Iz2 * L4 * Oy3 + 4857076224000 * E3 * Iy2 * Iz*L4 * Oy3 + 22417274880000 * E3 * Iy2 * Iz*L4 * Oy4 + 6486480000 * E3 * Iy*Iz2 * L4 * Oy + 7005398400 * E3 * Iy2 * Iz*L4 * Oy + 3891888000 * E3 * Iy*Iz2 * L4 * Oz + 249080832000 * E3 * Iy*Iz2 * L4 * Oy*Oz + 4296644352000 * E3 * Iy*Iz2 * L4 * Oy2 * Oz + 22417274880000 * E3 * Iy*Iz2 * L4 * Oy3 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*Fxb3
            eg[5 , 8] = eg[8 , 5] = ((Mxb*(1680 * L6 * Oy2 + 80 * L6 * Oy + L6)) / (100800 * E3 * Iz3 * Le*(12 * Oz + 1)) + (Le*Mxb*(124540416000 * E3 * Iz3 * L4 * Oz2 + 11955879936000 * E3 * Iz3 * L4 * Oz3 + 201755473920000 * E3 * Iz3 * L4 * Oz4 + 941525544960000 * E3 * Iz3 * L4 * Oz5 + 43243200 * E3 * Iy*Iz2 * L4 + 108108000 * E3 * Iy2 * Iz*L4 + 518918400 * E3 * Iz3 * L4 * Oz + 326918592000 * E3 * Iy*Iz2 * L4 * Oz2 + 108972864000 * E3 * Iy2 * Iz*L4 * Oz2 + 4857076224000 * E3 * Iy*Iz2 * L4 * Oz3 + 560431872000 * E3 * Iy2 * Iz*L4 * Oz3 + 22417274880000 * E3 * Iy*Iz2 * L4 * Oz4 + 3891888000 * E3 * Iy2 * Iz*L4 * Oy + 7005398400 * E3 * Iy*Iz2 * L4 * Oz + 6486480000 * E3 * Iy2 * Iz*L4 * Oz + 249080832000 * E3 * Iy2 * Iz*L4 * Oy*Oz + 4296644352000 * E3 * Iy2 * Iz*L4 * Oy*Oz2 + 22417274880000 * E3 * Iy2 * Iz*L4 * Oy*Oz3)) / (32691859200000 * E6 * Iy3 * Iz3 * niz4))*Fxb3
            eg[7 , 10] = eg[10 , 7] = (-(Mxb*(1680 * L6 * Oz2 + 80 * L6 * Oz + L6)) / (100800 * E3 * Iy3 * Le*(12 * Oy + 1)) - (Le*Mxb*(124540416000 * E3 * Iy3 * L4 * Oy2 + 11955879936000 * E3 * Iy3 * L4 * Oy3 + 201755473920000 * E3 * Iy3 * L4 * Oy4 + 941525544960000 * E3 * Iy3 * L4 * Oy5 + 108108000 * E3 * Iy*Iz2 * L4 + 43243200 * E3 * Iy2 * Iz*L4 + 518918400 * E3 * Iy3 * L4 * Oy + 108972864000 * E3 * Iy*Iz2 * L4 * Oy2 + 326918592000 * E3 * Iy2 * Iz*L4 * Oy2 + 560431872000 * E3 * Iy*Iz2 * L4 * Oy3 + 4857076224000 * E3 * Iy2 * Iz*L4 * Oy3 + 22417274880000 * E3 * Iy2 * Iz*L4 * Oy4 + 6486480000 * E3 * Iy*Iz2 * L4 * Oy + 7005398400 * E3 * Iy2 * Iz*L4 * Oy + 3891888000 * E3 * Iy*Iz2 * L4 * Oz + 249080832000 * E3 * Iy*Iz2 * L4 * Oy*Oz + 4296644352000 * E3 * Iy*Iz2 * L4 * Oy2 * Oz + 22417274880000 * E3 * Iy*Iz2 * L4 * Oy3 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*Fxb3
            eg[8 , 11] = eg[11 , 8] = (-(Mxb*(1680 * L6 * Oy2 + 80 * L6 * Oy + L6)) / (100800 * E3 * Iz3 * Le*(12 * Oz + 1)) - (Le*Mxb*(124540416000 * E3 * Iz3 * L4 * Oz2 + 11955879936000 * E3 * Iz3 * L4 * Oz3 + 201755473920000 * E3 * Iz3 * L4 * Oz4 + 941525544960000 * E3 * Iz3 * L4 * Oz5 + 43243200 * E3 * Iy*Iz2 * L4 + 108108000 * E3 * Iy2 * Iz*L4 + 518918400 * E3 * Iz3 * L4 * Oz + 326918592000 * E3 * Iy*Iz2 * L4 * Oz2 + 108972864000 * E3 * Iy2 * Iz*L4 * Oz2 + 4857076224000 * E3 * Iy*Iz2 * L4 * Oz3 + 560431872000 * E3 * Iy2 * Iz*L4 * Oz3 + 22417274880000 * E3 * Iy*Iz2 * L4 * Oz4 + 3891888000 * E3 * Iy2 * Iz*L4 * Oy + 7005398400 * E3 * Iy*Iz2 * L4 * Oz + 6486480000 * E3 * Iy2 * Iz*L4 * Oz + 249080832000 * E3 * Iy2 * Iz*L4 * Oy*Oz + 4296644352000 * E3 * Iy2 * Iz*L4 * Oy*Oz2 + 22417274880000 * E3 * Iy2 * Iz*L4 * Oy*Oz3)) / (32691859200000 * E6 * Iy3 * Iz3 * niz4))*Fxb3

            #eg[0 , 4] = eg[4 , 0] = 0
            eg[1 , 5] = eg[5 , 1] = (((L4 * (21772800 * Oy5 + 6480000 * Oy4 + 665280 * Oy3 + 25920 * Oy2 + 252 * Oy + 1)) / (42000 * niy5) + (Iz*L2 * (483840 * Oy4 + 97920 * Oy3 + 5376 * Oy2 + 60 * Oy + 1)) / (1400 * Ax*niy5)) / (E2 * Iz2) - (L2 * (1814400 * L2 * Oy4 + 388800 * L2 * Oy3 + 23040 * L2 * Oy2 + 240 * L2 * Oy + L2)) / (63000 * E2 * Iz2 * niy4))*Fxb3
            #eg[3 , 7] = eg[7 , 3] = 0
            eg[4 , 8] = eg[8 , 4] = (((L4 * (21772800 * Oz5 + 6480000 * Oz4 + 665280 * Oz3 + 25920 * Oz2 + 252 * Oz + 1)) / (42000 * niz5) + (Iy*L2 * (483840 * Oz4 + 97920 * Oz3 + 5376 * Oz2 + 60 * Oz + 1)) / (1400 * Ax*niy5)) / (E2 * Iy2) - (L2 * (1814400 * L2 * Oz4 + 388800 * L2 * Oz3 + 23040 * L2 * Oz2 + 240 * L2 * Oz + L2)) / (63000 * E2 * Iy2 * niz4))*Fxb3
            eg[5 , 9] = eg[9 , 5] = -((Mya + Myb)*(100800 * L6 * Oy3 + 5040 * L6 * Oy2 + 120 * L6 * Oy + L6)*Fxb3) / (604800 * E3 * Iz3)
            #eg[6 , 10] = eg[10 , 6] = 0
            eg[7 , 11] = eg[11 , 7] = (-((L4 * (21772800 * Oy5 + 6480000 * Oy4 + 665280 * Oy3 + 25920 * Oy2 + 252 * Oy + 1)) / (42000 * niy5) + (Iz*L2 * (483840 * Oy4 + 97920 * Oy3 + 5376 * Oy2 + 60 * Oy + 1)) / (1400 * Ax*niy5)) / (E2 * Iz2) + (L2 * (1814400 * L2 * Oy4 + 388800 * L2 * Oy3 + 23040 * L2 * Oy2 + 240 * L2 * Oy + L2)) / (63000 * E2 * Iz2 * niy4))*Fxb3

            #eg[0 , 5] = eg[5 , 0] = 0
            #eg[3 , 8] = eg[8 , 3] = 0
            eg[4 , 9] = eg[9 , 4] = ((Mza + Mzb)*(L6 / (E3 * Iy3) + (5040 * L6 * Oz2) / (E3 * Iy3) + (100800 * L6 * Oz3) / (E3 * Iy3) + (120 * L6 * Oz) / (E3 * Iy3))*Fxb3) / 604800
            eg[5 , 10] = eg[10 , 5] = ((L6 * Mxb*(648 * Oz + 23040 * Oy*Oz2 + 2211840 * Oy*Oz3 + 37324800 * Oy*Oz4 + 174182400 * Oy*Oz5 + 33840 * Oz2 + 967680 * Oz3 + 13236480 * Oz4 + 80870400 * Oz5 + 174182400 * Oz6 + 96 * Oy*Oz + 5)) / (1008000 * E3 * Iy3 * niy*niz4) + (L6 * Mxb*(648 * Oy + 23040 * Oy2 * Oz + 2211840 * Oy3 * Oz + 37324800 * Oy4 * Oz + 174182400 * Oy5 * Oz + 33840 * Oy2 + 967680 * Oy3 + 13236480 * Oy4 + 80870400 * Oy5 + 174182400 * Oy6 + 96 * Oy*Oz + 5)) / (1008000 * E3 * Iz3 * niy4 * niz) + (L6 * Mxb*(12441600 * Oy4 * Oz + 311040 * Oy4 + 12441600 * Oy3 * Oz2 + 5495040 * Oy3 * Oz + 172800 * Oy3 + 1658880 * Oy2 * Oz2 + 613440 * Oy2 * Oz + 20160 * Oy2 + 43200 * Oy*Oz2 + 20880 * Oy*Oz + 660 * Oy + 288 * Oz2 + 228 * Oz + 7)) / (3024000 * E3 * Iy*Iz2 * niy3 * niz2) + (L6 * Mxb*(12441600 * Oy2 * Oz3 + 1658880 * Oy2 * Oz2 + 43200 * Oy2 * Oz + 288 * Oy2 + 12441600 * Oy*Oz4 + 5495040 * Oy*Oz3 + 613440 * Oy*Oz2 + 20880 * Oy*Oz + 228 * Oy + 311040 * Oz4 + 172800 * Oz3 + 20160 * Oz2 + 660 * Oz + 7)) / (3024000 * E3 * Iy2 * Iz*niy2 * niz3)) * Fxb3
            #eg[6 , 11] = eg[11 , 6] = 0

            #eg[0 , 6] = eg[6 , 0] = 0
            eg[1 , 7] = eg[7 , 1] = (-((L3 * (21772800 * Oy5 + 6480000 * Oy4 + 665280 * Oy3 + 25920 * Oy2 + 252 * Oy + 1)) / (21000 * niy5) + (Iz*Le*(483840 * Oy4 + 97920 * Oy3 + 5376 * Oy2 + 60 * Oy + 1)) / (700 * Ax*niy5)) / (E2 * Iz2) + (Le*(1814400 * L2 * Oy4 + 388800 * L2 * Oy3 + 23040 * L2 * Oy2 + 240 * L2 * Oy + L2)) / (31500 * E2 * Iz2 * niy4))*Fxb3
            eg[2 , 8] = eg[8 , 2] = (-((L3 * (21772800 * Oz5 + 6480000 * Oz4 + 665280 * Oz3 + 25920 * Oz2 + 252 * Oz + 1)) / (21000 * niz5) + (Iy*Le*(483840 * Oz4 + 97920 * Oz3 + 5376 * Oz2 + 60 * Oz + 1)) / (700 * Ax*niz5)) / (E2 * Iy2) + (Le*(1814400 * L2 * Oz4 + 388800 * L2 * Oz3 + 23040 * L2 * Oz2 + 240 * L2 * Oz + L2)) / (31500 * E2 * Iy2 * niz4))*Fxb3
            #eg[3 , 9] = eg[9 , 3] = 0
            eg[4 , 10] = eg[10 , 4] = (-((L5 * (7838208000 * Oz7 + 3527193600 * Oz6 + 623635200 * Oz5 + 55080000 * Oz4 + 2620080 * Oz3 + 73620 * Oz2 + 1422 * Oz + 11)) / (126000 * niz5) + (Iy*L3 * (4354560 * Oz5 - 362880 * Oz4 - 138240 * Oz3 + 1008 * Oz2 + 780 * Oz + 13)) / (12600 * Ax*niz5)) / (E2 * Iy2) + (L3 * (653184000 * L2 * Oz6 + 239500800 * L2 * Oz5 + 32011200 * L2 * Oz4 + 1922400 * L2 * Oz3 + 58140 * L2 * Oz2 + 1290 * L2 * Oz + 11 * L2)) / (189000 * E2 * Iy2 * niz4))*Fxb3
            eg[5 , 11] = eg[11 , 5] = (-((L5 * (7838208000 * Oy7 + 3527193600 * Oy6 + 623635200 * Oy5 + 55080000 * Oy4 + 2620080 * Oy3 + 73620 * Oy2 + 1422 * Oy + 11)) / (126000 * niy5) + (Iz*L3 * (4354560 * Oy5 - 362880 * Oy4 - 138240 * Oy3 + 1008 * Oy2 + 780 * Oy + 13)) / (12600 * Ax*niy5)) / (E2 * Iz2) + (L3 * (653184000 * L2 * Oy6 + 239500800 * L2 * Oy5 + 32011200 * L2 * Oy4 + 1922400 * L2 * Oy3 + 58140 * L2 * Oy2 + 1290 * L2 * Oy + 11 * L2)) / (189000 * E2 * Iz2 * niy4))*Fxb3

            #eg[2 , 9] = eg[9 , 2] = 0
            eg[3 , 10] = eg[10 , 3] = ((Mza + Mzb)*(L6 / (E3 * Iy3) + (5040 * L6 * Oz2) / (E3 * Iy3) + (100800 * L6 * Oz3) / (E3 * Iy3) + (120 * L6 * Oz) / (E3 * Iy3))*Fxb3) / 604800
            eg[4 , 11] = eg[11 , 4] = -((L6 * Mxb*(648 * Oz + 23040 * Oy*Oz2 + 2211840 * Oy*Oz3 + 37324800 * Oy*Oz4 + 174182400 * Oy*Oz5 + 33840 * Oz2 + 967680 * Oz3 + 13236480 * Oz4 + 80870400 * Oz5 + 174182400 * Oz6 + 96 * Oy*Oz + 5)) / (1008000 * E3 * Iy3 * niy*niz4) + (L6 * Mxb*(648 * Oy + 23040 * Oy2 * Oz + 2211840 * Oy3 * Oz + 37324800 * Oy4 * Oz + 174182400 * Oy5 * Oz + 33840 * Oy2 + 967680 * Oy3 + 13236480 * Oy4 + 80870400 * Oy5 + 174182400 * Oy6 + 96 * Oy*Oz + 5)) / (1008000 * E3 * Iz3 * niy4 * niz) + (L6 * Mxb*(12441600 * Oy4 * Oz + 311040 * Oy4 + 12441600 * Oy3 * Oz2 + 5495040 * Oy3 * Oz + 172800 * Oy3 + 1658880 * Oy2 * Oz2 + 613440 * Oy2 * Oz + 20160 * Oy2 + 43200 * Oy*Oz2 + 20880 * Oy*Oz + 660 * Oy + 288 * Oz2 + 228 * Oz + 7)) / (3024000 * E3 * Iy*Iz2 * niy3 * niz2) + (L6 * Mxb*(12441600 * Oy2 * Oz3 + 1658880 * Oy2 * Oz2 + 43200 * Oy2 * Oz + 288 * Oy2 + 12441600 * Oy*Oz4 + 5495040 * Oy*Oz3 + 613440 * Oy*Oz2 + 20880 * Oy*Oz + 228 * Oy + 311040 * Oz4 + 172800 * Oz3 + 20160 * Oz2 + 660 * Oz + 7)) / (3024000 * E3 * Iy2 * Iz*niy2 * niz3))*Fxb3

            #eg[1 , 9] = eg[9 , 1] = 0
            eg[2 , 10] = eg[10 , 2] = (-((L4 * (21772800 * Oz5 + 6480000 * Oz4 + 665280 * Oz3 + 25920 * Oz2 + 252 * Oz + 1)) / (42000 * niz5) + (Iy*L2 * (483840 * Oz4 + 97920 * Oz3 + 5376 * Oz2 + 60 * Oz + 1)) / (1400 * Ax*niz5)) / (E2 * Iy2) + (L2 * (1814400 * L2 * Oz4 + 388800 * L2 * Oz3 + 23040 * L2 * Oz2 + 240 * L2 * Oz + L2)) / (63000 * E2 * Iy2 * niz4))*Fxb3
            eg[3 , 11] = eg[11 , 3] = -((Mya + Myb)*(100800 * L6 * Oy3 + 5040 * L6 * Oy2 + 120 * L6 * Oy + L6)*Fxb3) / (604800 * E3 * Iz3)

            eg[1 , 10] = eg[10 , 1] = ((Mxb*(1680 * L6 * Oz2 + 80 * L6 * Oz + L6)) / (100800 * E3 * Iy3 * Le*(12 * Oy + 1)) + (Le*Mxb*(124540416000 * E3 * Iy3 * L4 * Oy2 + 11955879936000 * E3 * Iy3 * L4 * Oy3 + 201755473920000 * E3 * Iy3 * L4 * Oy4 + 941525544960000 * E3 * Iy3 * L4 * Oy5 + 108108000 * E3 * Iy*Iz2 * L4 + 43243200 * E3 * Iy2 * Iz*L4 + 518918400 * E3 * Iy3 * L4 * Oy + 108972864000 * E3 * Iy*Iz2 * L4 * Oy2 + 326918592000 * E3 * Iy2 * Iz*L4 * Oy2 + 560431872000 * E3 * Iy*Iz2 * L4 * Oy3 + 4857076224000 * E3 * Iy2 * Iz*L4 * Oy3 + 22417274880000 * E3 * Iy2 * Iz*L4 * Oy4 + 6486480000 * E3 * Iy*Iz2 * L4 * Oy + 7005398400 * E3 * Iy2 * Iz*L4 * Oy + 3891888000 * E3 * Iy*Iz2 * L4 * Oz + 249080832000 * E3 * Iy*Iz2 * L4 * Oy*Oz + 4296644352000 * E3 * Iy*Iz2 * L4 * Oy2 * Oz + 22417274880000 * E3 * Iy*Iz2 * L4 * Oy3 * Oz)) / (32691859200000 * E6 * Iy3 * Iz3 * niy4))*Fxb3
            eg[2 , 11] = eg[11 , 2] = ((Mxb*(1680 * L6 * Oy2 + 80 * L6 * Oy + L6)) / (100800 * E3 * Iz3 * Le*(12 * Oz + 1)) + (Le*Mxb*(124540416000 * E3 * Iz3 * L4 * Oz2 + 11955879936000 * E3 * Iz3 * L4 * Oz3 + 201755473920000 * E3 * Iz3 * L4 * Oz4 + 941525544960000 * E3 * Iz3 * L4 * Oz5 + 108108000 * E3 * Iz*Iy2 * L4 + 43243200 * E3 * Iz2 * Iy*L4 + 518918400 * E3 * Iz3 * L4 * Oz + 108972864000 * E3 * Iz*Iy2 * L4 * Oz2 + 326918592000 * E3 * Iz2 * Iy*L4 * Oz2 + 560431872000 * E3 * Iz*Iy2 * L4 * Oz3 + 4857076224000 * E3 * Iz2 * Iy*L4 * Oz3 + 22417274880000 * E3 * Iz2 * Iy*L4 * Oz4 + 6486480000 * E3 * Iz*Iy2 * L4 * Oz + 7005398400 * E3 * Iz2 * Iy*L4 * Oz + 3891888000 * E3 * Iz*Iy2 * L4 * Oz + 249080832000 * E3 * Iz*Iy2 * L4 * Oz*Oy + 4296644352000 * E3 * Iz*Iy2 * L4 * Oz2 * Oy + 22417274880000 * E3 * Iz*Iy2 * L4 * Oz3 * Oy)) / (32691859200000 * E6 * Iz3 * Iy3 * niz4))*Fxb3

            #eg[0 , 10] = eg[10 , 0] = 0
            eg[1 , 11] = eg[11 , 1] = (((L4 * (21772800 * Oy5 + 6480000 * Oy4 + 665280 * Oy3 + 25920 * Oy2 + 252 * Oy + 1)) / (42000 * niy5) + (Iz*L2 * (483840 * Oy4 + 97920 * Oy3 + 5376 * Oy2 + 60 * Oy + 1)) / (1400 * Ax*niy5)) / (E2 * Iz2) - (L2 * (1814400 * L2 * Oy4 + 388800 * L2 * Oy3 + 23040 * L2 * Oy2 + 240 * L2 * Oy + L2)) / (63000 * E2 * Iz2 * niy4))*Fxb3

            #eg[0 , 11] = eg[11 , 0] = 0

            eg[4 , 5] = eg[5 , 4] = ((Iy*L6 * Mxb*Fxb3 * (49766400 * Oy4 * Oz2 + 5391360 * Oy4 * Oz + 103680 * Oy4 - 49766400 * Oy3 * Oz3 + 1244160 * Oy3 * Oz2 + 449280 * Oy3 * Oz - 6635520 * Oy2 * Oz3 - 311040 * Oy2 * Oz2 + 8640 * Oy2 * Oz - 960 * Oy2 - 172800 * Oy*Oz3 + 11520 * Oy*Oz2 + 2400 * Oy*Oz + 20 * Oy - 1152 * Oz3 + 432 * Oz2 + 56 * Oz + 1)) / 1008000 - (Iz*L6 * Mxb*Fxb3 * (-49766400 * Oy3 * Oz3 - 6635520 * Oy3 * Oz2 - 172800 * Oy3 * Oz - 1152 * Oy3 + 49766400 * Oy2 * Oz4 + 1244160 * Oy2 * Oz3 - 311040 * Oy2 * Oz2 + 11520 * Oy2 * Oz + 432 * Oy2 + 5391360 * Oy*Oz4 + 449280 * Oy*Oz3 + 8640 * Oy*Oz2 + 2400 * Oy*Oz + 56 * Oy + 103680 * Oz4 - 960 * Oz2 + 20 * Oz + 1)) / 1008000) / (E3 * Iy2 * Iz2 * niy3 * niz3) + (L6 * Mxb*Fxb3 * (632 * Oy - 23040 * Oy2 * Oz - 2211840 * Oy3 * Oz - 37324800 * Oy4 * Oz - 174182400 * Oy5 * Oz + 30000 * Oy2 + 599040 * Oy3 + 7015680 * Oy4 + 51840000 * Oy5 + 174182400 * Oy6 - 96 * Oy*Oz + 5)) / (1008000 * E3 * Iz3 * niy4 * niz) - (L6 * Mxb*Fxb3 * (-25082265600 * Oy3 * Oz5 - 5374771200 * Oy3 * Oz4 - 318504960 * Oy3 * Oz3 - 3317760 * Oy3 * Oz2 - 13824 * Oy3 * Oz + 25082265600 * Oy2 * Oz6 + 3284582400 * Oy2 * Oz5 + 114462720 * Oy2 * Oz4 + 33177600 * Oy2 * Oz3 + 3767040 * Oy2 * Oz2 + 88704 * Oy2 * Oz + 720 * Oy2 + 4180377600 * Oy*Oz6 + 1069977600 * Oy*Oz5 + 131051520 * Oy*Oz4 + 12165120 * Oy*Oz3 + 696960 * Oy*Oz2 + 15072 * Oy*Oz + 120 * Oy + 174182400 * Oz6 + 51840000 * Oz5 + 7015680 * Oz4 + 599040 * Oz3 + 30000 * Oz2 + 632 * Oz + 5)) / (1008000 * E3 * Iy3 * niy3 * niz4)
            eg[10 , 11] = eg[11 , 10] = (L6 * Mxb*Fxb3 * (-25082265600 * Oy3 * Oz5 - 5374771200 * Oy3 * Oz4 - 318504960 * Oy3 * Oz3 - 3317760 * Oy3 * Oz2 - 13824 * Oy3 * Oz + 25082265600 * Oy2 * Oz6 + 3284582400 * Oy2 * Oz5 + 114462720 * Oy2 * Oz4 + 33177600 * Oy2 * Oz3 + 3767040 * Oy2 * Oz2 + 88704 * Oy2 * Oz + 720 * Oy2 + 4180377600 * Oy*Oz6 + 1069977600 * Oy*Oz5 + 131051520 * Oy*Oz4 + 12165120 * Oy*Oz3 + 696960 * Oy*Oz2 + 15072 * Oy*Oz + 120 * Oy + 174182400 * Oz6 + 51840000 * Oz5 + 7015680 * Oz4 + 599040 * Oz3 + 30000 * Oz2 + 632 * Oz + 5)) / (1008000 * E3 * Iy3 * niy3 * niz4) - (L6 * Mxb*Fxb3 * (632 * Oy - 23040 * Oy2 * Oz - 2211840 * Oy3 * Oz - 37324800 * Oy4 * Oz - 174182400 * Oy5 * Oz + 30000 * Oy2 + 599040 * Oy3 + 7015680 * Oy4 + 51840000 * Oy5 + 174182400 * Oy6 - 96 * Oy*Oz + 5)) / (1008000 * E3 * Iz3 * niy4 * niz) - ((Iy*L6 * Mxb*Fxb3 * (49766400 * Oy4 * Oz2 + 5391360 * Oy4 * Oz + 103680 * Oy4 - 49766400 * Oy3 * Oz3 + 1244160 * Oy3 * Oz2 + 449280 * Oy3 * Oz - 6635520 * Oy2 * Oz3 - 311040 * Oy2 * Oz2 + 8640 * Oy2 * Oz - 960 * Oy2 - 172800 * Oy*Oz3 + 11520 * Oy*Oz2 + 2400 * Oy*Oz + 20 * Oy - 1152 * Oz3 + 432 * Oz2 + 56 * Oz + 1)) / 1008000 - (Iz*L6 * Mxb*Fxb3 * (-49766400 * Oy3 * Oz3 - 6635520 * Oy3 * Oz2 - 172800 * Oy3 * Oz - 1152 * Oy3 + 49766400 * Oy2 * Oz4 + 1244160 * Oy2 * Oz3 - 311040 * Oy2 * Oz2 + 11520 * Oy2 * Oz + 432 * Oy2 + 5391360 * Oy*Oz4 + 449280 * Oy*Oz3 + 8640 * Oy*Oz2 + 2400 * Oy*Oz + 56 * Oy + 103680 * Oz4 - 960 * Oz2 + 20 * Oz + 1)) / 1008000) / (E3 * Iy2 * Iz2 * niy3 * niz3)

    return eg


#
#
#
#
def B3D2_Kt(Le:float, 
            Ax:float, Asy:float, Asz:float,
            Jx:float, Iy:float, Iz:float,
            Emod:float, Gmod:float,
            Fb:list,
            shear: bool = True,
            toler: float = 0.01):
    """
    This method calculates the local tangent stiffness matrix 
    for an element considering the Timoshenko beam theory.
    
    This matrix is obtained directly from the solution of the 
    equilibrium differential equation of a deformed infinitesimal
    element (Rodrigues, 2019).
    
    Input: 
    Le : Current element length
    Ax : Element cross-section area
    Asy : Effective shear area in y direct.
    Asz : Effective shear area in z direct.
    Emod : Young's modulus                                                 
    Gmod : Shear modulus
    Jx : Mom. of inertia wrt "x" axis                                   
    Iy : Mom. of inertia wrt "y" axis                                   
    Iz : Mom. of inertia wrt "z" axis
    Fb : Forces and moments at beam's end nodes [Fya, Fza, Mxa, Mya, Mza, Fyb, Fzb, Mxb, Myb, Mzb]
    toler : tolerance for checking null axial force
    
    Return :
    et : Tangent stiffness matrix
    """
    #if len(Fb) == 6:
    #    # ['x', 'y', 'rz']
    #    Fxa, Fya, Mza, Fxb, Fyb, Mzb = Fb
    #    Fb = np.array([Fxa, Fya, 0.0, 0.0, 0.0, Mza,
    #                   Fxb, Fyb, 0.0, 0.0, 0.0, Mzb])
    #else:
    #    (Fxa, Fya, Fza, Mxa, Mya, Mza,
    #     Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb
    #
    #P = 0.50 * (Fb[6] - Fb[0])
    P = 1 * Fb.Fx
    #
    #1 / 0
    #
    #A = Ax
    L2 = Le*Le
    #L3 = L2*Le
    pi = np.pi
    #
    Oy = 0.0
    Oz = 0.0
    if shear:
        # Shear to Flexural rig. ratio y
        Oy = (Emod*Iz) / (Gmod*Asy*L2)
        # Shear to Flexural rig. ratio z
        Oz = (Emod*Iy) / (Gmod*Asz*L2)
    #
    # minimum cross-section mom. inertia
    Imin = min(Iy, Iz)
    # Euler critical load 
    PE = pi*pi*Emod*Imin / L2
    # factor axial force to Euler load
    # Verifies if it is important consider axial load in the problem 
    # Tries to avoid numerical instability
    fac = P / PE
    #
    mu_y = sqrt(abs(P) / (Emod * Iz))
    mu_z = sqrt(abs(P) / (Emod * Iy))   
    #
    # Tension
    # TODO : Check tension is not required
    if fac > toler:
        print('Positive axial force')
        #Kt1 = Kt_tension(Le=Le, Ax=Ax,
        #                 Jx=Jx, Iy=Iy, Iz=Iz,
        #                 Emod=Emod, Gmod=Gmod,
        #                 Oy=Oy, Oz=Oz, Fb=Fb)
        #
        #Kt = StfBeamTimoshenkoT(L=Le, A=Ax,
        #                        Jx=Jx, Iy=Iy, Iz=Iz,
        #                        E=Emod, G=Gmod,
        #                        Oy=Oy, Oz=Oz, Fb=Fb)
        #
        Ke5 = Ke_new(Le=Le, Ax=Ax,
                     Jx=Jx, Iy=Iy, Iz=Iz,
                     Emod=Emod, Gmod=Gmod,
                     Oy=Oy, Oz=Oz,
                     mu_y=mu_y, mu_z=mu_z,
                     Fb=Fb, tension=True)
        #
        Kg5 = Kg_new(Le=Le, Ax=Ax,
                     Jx=Jx, Iy=Iy, Iz=Iz,
                     Emod=Emod, Gmod=Gmod,
                     Oy=Oy, Oz=Oz,
                     mu_y=mu_y, mu_z=mu_z,
                     Fb=Fb, tension=True)
        #
        Kint5 = Kint_new(Le=Le, Ax=Ax,
                         Jx=Jx, Iy=Iy, Iz=Iz,
                         Emod=Emod, Gmod=Gmod,
                         Oy=Oy, Oz=Oz,
                         mu_y=mu_y, mu_z=mu_z,
                         Fb=Fb, tension=True)
        #
        Kt = Ke5 + Kg5 + Kint5
    # Compression
    elif fac < -toler:
        print('Negative axial force')
        #
        #Kt = Kt_compression(Le=Le, Ax=Ax,
        #                    Jx=Jx, Iy=Iy, Iz=Iz,
        #                    Emod=Emod, Gmod=Gmod,
        #                    Oy=Oy, Oz=Oz, Fb=Fb)
        #
        #Kt = StfBeamTimoshenkoC(L=Le, A=Ax,
        #                        Jx=Jx, Iy=Iy, Iz=Iz,
        #                        E=Emod, G=Gmod,
        #                        Oy=Oy, Oz=Oz, Fb=Fb)
        #
        Ke5 = Ke_new(Le=Le, Ax=Ax,
                     Jx=Jx, Iy=Iy, Iz=Iz,
                     Emod=Emod, Gmod=Gmod,
                     Oy=Oy, Oz=Oz,
                     mu_y=mu_y, mu_z=mu_z,
                     Fb=Fb, tension=False)
        #
        Kg5 = Kg_new(Le=Le,Ax=Ax,
                     Jx=Jx, Iy=Iy, Iz=Iz,
                     Emod=Emod, Gmod=Gmod,
                     Oy=Oy, Oz=Oz,
                     mu_y=mu_y, mu_z=mu_z,
                     Fb=Fb, tension=False)
        #
        Kint5 = Kint_new(Le=Le, Ax=Ax,
                         Jx=Jx, Iy=Iy, Iz=Iz,
                         Emod=Emod, Gmod=Gmod,
                         Oy=Oy, Oz=Oz,
                         mu_y=mu_y, mu_z=mu_z,
                         Fb=Fb, tension=False)
        #
        Kt = Ke5 + Kg5 + Kint5
    else:
        #print('NoAxial')
        #
        Kt = Kt_NoAxial(Le=Le, Ax=Ax,
                        Jx=Jx, Iy=Iy, Iz=Iz,
                        Emod=Emod, Gmod=Gmod,
                        Oy=Oy, Oz=Oz, Fb=Fb)
        #
        #Kg = GeoStfBeamTimoshenko(Le=Le, 
        #                          Ax=Ax, Asy=Asy, Asz=Asz,
        #                          Jx=Jx, Iy=Iy, Iz=Iz,
        #                          Emod=Emod, Gmod=Gmod,
        #                          Fb=Fb,
        #                          shear= True,
        #                          order=2)
        #
        #        
        #
        #
        #
        #Kg = beam_KG(Le=Le, 
        #             Ax=Ax, Asy=Asy, Asz=Asz,
        #             Jx=Jx, Iy=Iy, Iz=Iz,
        #             Emod=Emod, Gmod=Gmod,
        #             Fb=Fb,
        #             shear=True)   
        #
        #Ke = B3D2_Ke(Le=Le, 
        #             Ax=Ax, Asy=Asy, Asz=Asz,
        #             Jx=Jx, Iy=Iy, Iz=Iz,
        #             Emod=Emod, Gmod=Gmod,
        #             shear=True)
        #
        #Ke1 = beam3D_B3D2(Le=Le, 
        #                 Ax=Ax, Asy=Asz, Asz=Asy,
        #                 Jx=Jx, Iy=Iz, Iz=Iy,
        #                 Emod=Emod, Gmod=Gmod,
        #                 shear = False)
        #Kt1 = Ke + Kg
    #
    #Ke = B3D2_Ke(Le=Le, 
    #             Ax=Ax, Asy=Asy, Asz=Asz,
    #             Jx=Jx, Iy=Iy, Iz=Iz,
    #             Emod=Emod, Gmod=Gmod,
    #             shear=True)
    #
    #Ke1 = beam3D_B3D2(Le=Le, 
    #                 Ax=Ax, Asy=Asz, Asz=Asy,
    #                 Jx=Jx, Iy=Iz, Iz=Iy,
    #                 Emod=Emod, Gmod=Gmod,
    #                 shear = False)
    
    #Kg = beam_KG(Le=Le, 
    #             Ax=Ax, Asy=Asy, Asz=Asz,
    #             Jx=Jx, Iy=Iy, Iz=Iz,
    #             Emod=Emod, Gmod=Gmod,
    #             Fb=Fb,
    #             shear=True)
    #
    return Kt
#
#
def Kt_NoAxial(Le:float, Ax:float,
               Jx:float, Iy:float, Iz:float,
               Emod:float, Gmod:float,
               Oy:float, Oz:float,
               Fb:list[float]):
    """
    Null axial force: complete stiffness matrix is the elastic matrix
    + geometric matrix with 1 term"""
    #
    (Fxa, Fya, Fza, Mxa, Mya, Mza,
     Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb.q_loc
    #
    #Fxb = 0
    #
    A = Ax
    L = Le
    E = Emod
    G = Gmod    
    #
    #Oy = (E*Iz)/(G*Ay*L**2)      # Shear to Flexural rig. ratio y                         
    #Oz = (E*Iy)/(G*Az*L**2)      # Shear to Flexural rig. ratio z    
    #
    niy = 1.0 + 12.0 * Oy
    niz = 1.0 + 12.0 * Oz
    #
    lay = 1.0 + 3.0 * Oy
    laz = 1.0 + 3.0 * Oz
    #
    gay = 1.0 - 6.0 * Oy
    gaz = 1.0 - 6.0 * Oz
    #
    L2 = L*L
    L3 = L2*L
    #pi = np.pi
    #P = 0.0
    kvvy = 1.0 / niy
    kvry = 1.0 / niy
    krry = lay / niy
    krty = gay / niy
    kvvz = 1.0 / niz
    kvrz = 1.0 / niz
    krrz = laz / niz
    krtz = gaz / niz
    #
    #
    et = np.zeros((12,12), dtype=np.float32) 
    #
    #
    et[0 , 0] = E*A / L  #+ Mxb / L
    
    et[1 , 1] = (12.0*E*Iz / L3)*kvvy
    et[2 , 2] = (12.0*E*Iy / L3)*kvvz
    
    et[3 , 3] = G*Jx / L
    
    et[4 , 4] = (4.0*E*Iy / L)*krrz
    et[5 , 5] = (4.0*E*Iz / L)*krry
    #
    et[1 , 5] = (6.0*E*Iz / L2)*kvry
    et[2 , 4] = -(6.0*E*Iy / L2)*kvrz
    #
    et[4 , 10] = (2.0*E*Iy / L)*krtz
    et[5 , 11] = (2.0*E*Iz / L)*krty
    #
    #et[3, 9] = -et[3 , 3]
    #et[4 , 2] = et[2 , 4]
    #et[4 , 8] = -et[2 , 4]
    #et[5 , 1] = et[1 , 5]
    #et[5 , 7] = -et[1 , 5]
    #et[6 , 0] = -et[0 , 0]
    #et[7 , 1] = -et[1 , 1]
    #et[7 , 5] = -et[1 , 5]
    #et[8 , 2] = -et[2 , 2]
    #et[8 , 4] = -et[2 , 4]
    #
    et[1, 3] = Mya / L
    et[1, 9] = Myb / L
    
    et[2, 3] = Mza / L
    et[2, 9] = Mzb / L
    #
    et[3, 1] = Mya / L
    et[3, 2] = Mza / L    

    #et[2 , 3] = Mza / L
    et[3 , 4] = -Mza / 3.0 + Mzb / 6.0
    et[5 , 6] = Mza / L
    #et[8 , 9] = -Mzb / L
    et[9 , 10] = Mza / 6.0 - Mzb / 3.0

    #et[1 , 3] = Mya / L
    et[3 , 5] = Mya / 3.0 - Myb / 6.0
    et[4 , 6] = Mya / L
    #et[7 , 9] = -Myb / L
    et[9 , 11] = -Mya / 6.0 + Myb / 3.0

    et[1 , 4] = (Mxb / L)*(1 / niz)
    et[2 , 5] = (Mxb / L)*(1 / niy)
    et[4 , 7] = (-Mxb / L)*(1 / niz)
    #et[5 , 8] = (-Mxb / L)*(1 / niy)
    #et[7 , 10] = (Mxb / L)*(1 / niz)
    #et[8 , 11] = (Mxb / L)*(1 / niy)
    et[5, 2] = (Mxb / L)*(1 / niy)

    et[4 , 0] = -Mya / L
    et[0 , 4] = -Mya / L
    #et[3 , 7] = -Mya / L
    et[5 , 9] = Mya / 6.0 + Myb / 6.0
    et[6 , 10] = Myb / L
    et[10 , 6] = Myb / L

    et[0 , 5] = -Mza / L
    #et[3 , 8] = -Mza / L
    et[4 , 9] = -Mza / 6.0 - Mzb / 6.0
    et[5 , 10] = (-Mxb / 2.0)*((1 - 144 * Oy*Oz) / (niy*niz))
    et[6 , 11] = Mzb / L

    #et[2 , 9] = Mzb / L
    et[3 , 10] = -Mza / 6.0 - Mzb / 6.0
    et[4 , 11] = (Mxb / 2.0)*((1 - 144 * Oy*Oz) / (niy*niz))

    #et[1 , 9] = Myb / L
    et[3 , 11] = Mya / 6.0 + Myb / 6.0

    #et[1 , 10] = (-Mxb / L)*(1 / niz)
    #et[2 , 11] = (-Mxb / L)*(1 / niy)

    et[0 , 10] = -Myb / L
    et[0 , 11] = -Mzb / L
    
    et[6 , 4] = -et[4 , 0]

    #et[0 , 11] = -Mzb / L

    et[4 , 5] = -6 * Mxb*(Oy - Oz) / (niy*niz)
    et[10 , 11] = 6 * Mxb*(Oy - Oz) / (niy*niz)
    #
    #
    et[5, 4] = et[4, 5]
    et[5, 8] = -et[5, 2]
    et[7, 10] = - et[4, 7]
    et[8, 11] = et[2, 5]
    #
    et[4, 1] = et[1, 4]
    #et[4, 11] = -et[5, 10]
    #
    et[1, 10] = -et[1, 4]
    et[2, 11] = -et[2, 5]    
    #
    # ---------------
    # impose the symmetry
    #
    et[6, 6] = et[0, 0]
    et[7, 7] = et[1, 1]
    et[8, 8] = et[2, 2]
    et[9, 9] = et[3, 3]
    et[10, 10] = et[4, 4]
    et[11, 11] = et[5, 5] 
    #    
    #
    et[0, 6] = -et[0 , 0]
    et[1, 7] = -et[1 , 1]
    
    et[2, 8] = -et[2 , 2]
    
    et[3, 9] = -et[3 , 3]
    et[4, 8] = -et[2, 4]
    et[5, 7] = -et[1, 5]
    
    et[7, 11] = -et[1, 5]
    et[8, 10] = -et[2, 4]    
    
    et[1 , 11] = et[1 , 5]
    et[2 , 10] = et[2 , 4]
    #   
    # ---------------------------------------------
    # Filling the rest of the matrix
    # impose the geometry
    et += np.triu(et, k=1).T    
    #
    return et
#
#
#
def Kt_compression(Le:float, Ax:float,
                   Jx:float, Iy:float, Iz:float,
                   Emod:float, Gmod:float,
                   Oy:float, Oz:float,
                   Fb:list[float]):
    """Negative axial force"""
    (Fxa, Fya, Fza, Mxa, Mya, Mza,
     Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb.q_loc
    #
    P = Fb.Fx
    A = Ax
    L2 = Le*Le
    #L3 = L2*Le
    #pi = np.pi
    #
    power = np.power
    sin = np.sin
    cos = np.cos
    sqrt = np.sqrt 
    #
    mu_y = sqrt(np.abs(P) / (Emod*Iz))
    mu_z = sqrt(np.abs(P) / (Emod*Iy))
    #
    Vy = mu_y / (sqrt(1 - Oy*mu_y*mu_y*L2))
    Vz = mu_z / (sqrt(1 - Oz*mu_z*mu_z*L2))

    K = P*(Iy + Iz) / Ax
    
    ay = ((4.0 * (sin(Le * Vy / 2.0) * sin(Le * Vy / 2.0))
           - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sin(Le * Vy / 2.0) * sin(Le * Vy / 2.0)))
    az = ((4.0 * (sin(Le * Vz / 2.0) * sin(Le * Vz / 2.0))
           - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sin(Le * Vz / 2.0) * sin(Le * Vz / 2.0)))
    
    by = Oy * (Le * Le) * (Vy * Vy) + 1.0
    bz = Oz * (Le * Le) * (Vz * Vz) + 1.0
    #
    cy = ((((2.0 * cos(Le * Vy) - 2.0 * (Le * Le) * (Vy * Vy) * Oy) + Le * Vy * sin(Le * Vy))
           + 2.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) - 2.0)
    
    cz = ((((2.0 * cos(Le * Vz) - 2.0 * (Le * Le) * (Vz * Vz) * Oz) + Le * Vz * sin(Le * Vz))
           + 2.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) - 2.0)
    #
    dy = 4.0 * A * (Le * Le) * (Vy * Vy) * Oy * P
    dz = 4.0 * A * (Le * Le) * (Vz * Vz) * Oz * P

    ey = ((((((((((((((((((((2.0 * A * P * sin(Le * Vy)
                             - A * P * sin(2.0 * Le * Vy))
                            - 2.0 * Iz * (Vy * Vy) * P * sin(Le * Vy))
                           + Iz * (Vy * Vy) * P * sin(2.0 * Le * Vy)) - 6.0 * A * Le * Vy * P)
                         - Iz * Le * power(Vy, 3.0) * P) - A * Emod * Iz * Le * power(Vy, 3.0))
                       - 2.0 * A * Emod * Iz * (Vy * Vy) * sin(Le * Vy))
                      + A * Emod * Iz * (Vy * Vy) * sin(2.0 * Le * Vy))
                     + 2.0 * Iz * Le * power(Vy, 3.0) * P * cos(Le * Vy))
                    - Iz * Le * power(Vy, 3.0) * P * cos(2.0 * Le * Vy))
                   - 17.0 * A * power(Le, 3.0) * power(Vy, 3.0) * Oy * P)
                  - 3.0 * Iz * power(Le, 3.0) * power(Vy, 5.0) * Oy * P)
                 - A * power(Le, 3.0) * power(Vy, 3.0) * P * cos(Le * Vy))
                - Iz * power(Le, 3.0) * power(Vy, 5.0) * P * cos(Le * Vy))
               + 3.0 * A * (Le * Le) * (Vy * Vy) * P * sin(Le * Vy))
              + Iz * (Le * Le) * power(Vy, 4.0) * P * sin(Le * Vy))
             + 6.0 * A * Le * Vy * P * cos(Le * Vy))
            - 18.0 * A * power(Le, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P)
           - 9.0 * A * power(Le, 7.0) * power(Vy, 7.0) * power(Oy, 3.0) * P)
          - 2.0 * A * power(Le, 9.0) * power(Vy, 9.0) * power(Oy, 4.0) * P)

    ez = ((((((((((((((((((((2.0 * A * P * sin(Le * Vz)
                             - A * P * sin(2.0 * Le * Vz))
                            - 2.0 * Iy * (Vz * Vz) * P * sin(Le * Vz))
                           + Iy * (Vz * Vz) * P * sin(2.0 * Le * Vz))
                          - 6.0 * A * Le * Vz * P) - Iy * Le * power(Vz, 3.0) * P)
                        - A * Emod * Iy * Le * power(Vz, 3.0))
                       - 2.0 * A * Emod * Iy * (Vz * Vz) * sin(Le * Vz))
                      + A * Emod * Iy * (Vz * Vz) * sin(2.0 * Le * Vz))
                     + 2.0 * Iy * Le * power(Vz, 3.0) * P * cos(Le * Vz))
                    - Iy * Le * power(Vz, 3.0) * P * cos(2.0 * Le * Vz))
                   - 17.0 * A * power(Le, 3.0) * power(Vz, 3.0) * Oz * P)
                  - 3.0 * Iy * power(Le, 3.0) * power(Vz, 5.0) * Oz * P)
                 - A * power(Le, 3.0) * power(Vz, 3.0) * P * cos(Le * Vz))
                - Iy * power(Le, 3.0) * power(Vz, 5.0) * P * cos(Le * Vz))
               + 3.0 * A * (Le * Le) * (Vz * Vz) * P * sin(Le * Vz))
              + Iy * (Le * Le) * power(Vz, 4.0) * P * sin(Le * Vz))
             + 6.0 * A * Le * Vz * P * cos(Le * Vz))
            - 18.0 * A * power(Le, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P)
           - 9.0 * A * power(Le, 7.0) * power(Vz, 7.0) * power(Oz, 3.0) * P)
          - 2.0 * A * power(Le, 9.0) * power(Vz, 9.0) * power(Oz, 4.0) * P)

    #
    fy = 2.0 * Iz * (Le * Le) * power(Vy, 4.0) * Oy
    fz = 2.0 * Iy * (Le * Le) * power(Vz, 4.0) * Oz
    #
    sy = sin(Le * Vy / 2.0)
    sz = sin(Le * Vz / 2.0)
    #
    c_y = cos(Le * Vy / 2.0)
    c_z = cos(Le * Vz / 2.0)
    #
    #
    et = np.zeros((12,12), dtype=np.float32)
    #
    #
    et[0 , 0] = et[6 , 6] = P / Le + A * Emod / Le

    et[1 , 1] = et[7 , 7] = (((((sy * sy * (2.0 * Le * (Vy * Vy)
                                            * P * (((A * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)
                                                     + 2.0 * A * (Le * Le) * (Vy * Vy) * Oy)
                                                    + Iz * (Vy * Vy)) + A)
                                            - Vy * P * sin(Le * Vy) * (((-2.0 * A * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)
                                                                            + 4.0 * A * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Iz * (Vy * Vy)) + 6.0 * A))
                                 + A * Le * (Vy * Vy) * P * (sin(Le * Vy) * sin(Le * Vy))) / (A * (ay * ay))
                                - 2.0 * Emod * Iz * power(Vy, 3.0) * (sy * sy) * (sin(Le * Vy) - Le * Vy) / (ay * ay))
                               + Emod * Iz * Vy * ((((((Le * Vy * (sin(Le * Vy) * sin(Le * Vy))
                                                        - 6.0 * sin(Le * Vy) * (sy * sy)) + 2.0 * Le * Vy * (sy * sy))
                                                      + 2.0 * power(Le, 5.0) * power(Vy, 5.0) * (Oy * Oy) * (sy * sy))
                                                     + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * (sy * sy))
                                                    - 4.0 * (Le * Le) * (Vy * Vy) * Oy * sin(Le * Vy) * (sy * sy))
                                                   + 2.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * sin(Le * Vy) * (sy * sy)) / (Le * Le * Oy * (ay * ay)))
                              + 2.0 * Emod * Iz * Vy * ((((2.0 * Le * Vy - 3.0 * sin(Le * Vy))
                                                          + power(Le, 3.0) * power(Vy, 3.0) * Oy)
                                                         + Le * Vy * cos(Le * Vy))
                                                        - Le * Le * (Vy * Vy) * Oy * sin(Le * Vy))
                              / (Le * Le * Oy * (((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy))
                                                        - 8.0 * (Le * Le) * (Vy * Vy) * Oy)
                                                       - Le * Le * (Vy * Vy) * cos(Le * Vy))
                                                      + 4.0 * Le * Vy * sin(Le * Vy))
                                                     - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy))
                                                    + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy))
                                                   + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy))
                                                  + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) - 4.0)))
                             + 2.0 * Emod * Iz * Vy * (sy * sy) * ((2.0 * Le * Vy - 3.0 * sin(Le * Vy))
                                                                   + Le * Vy * cos(Le * Vy)) / (Le * Le * Oy * (ay * ay)))

    et[2 , 2] = et[8 , 8] = (((((sz * sz * (2.0 * Le * (Vz * Vz) * P * (((A * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)
                                                                          + 2.0 * A * (Le * Le) * (Vz * Vz) * Oz)
                                                                         + Iy * (Vz * Vz)) + A) - Vz * P * sin(Le * Vz)
                                            * (((-2.0 * A * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)
                                                 + 4.0 * A * (Le * Le) * (Vz * Vz) * Oz) + 2.0 * Iy * (Vz * Vz))
                                               + 6.0 * A)) + A * Le * (Vz * Vz) * P * (sin(Le * Vz) * sin(Le * Vz)))
                                / (A * (az * az)) - 2.0 * Emod * Iy * power(Vz, 3.0) * (sz * sz) * (sin(Le * Vz) - Le * Vz) / (az * az))
                               + Emod * Iy * Vz * ((((((Le * Vz * (sin(Le * Vz) * sin(Le * Vz))
                                                        - 6.0 * sin(Le * Vz) * (sz * sz)) + 2.0 * Le * Vz * (sz * sz))
                                                      + 2.0 * power(Le, 5.0) * power(Vz, 5.0) * (Oz * Oz) * (sz * sz))
                                                     + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * (sz * sz))
                                                    - 4.0 * (Le * Le) * (Vz * Vz) * Oz * sin(Le * Vz) * (sz * sz))
                                                   + 2.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * sin(Le * Vz) * (sz * sz)) / (Le * Le * Oz * (az * az)))
                              + 2.0 * Emod * Iy * Vz * ((((2.0 * Le * Vz - 3.0 * sin(Le * Vz))
                                                          + power(Le, 3.0) * power(Vz, 3.0) * Oz)
                                                         + Le * Vz * cos(Le * Vz)) - Le * Le * (Vz * Vz) * Oz * sin(Le * Vz))
                              / (Le * Le * Oz * (((((((((4.0 * cos(Le * Vz) - Le * Le * (Vz * Vz))
                                                        - 8.0 * (Le * Le) * (Vz * Vz) * Oz)
                                                       - Le * Le * (Vz * Vz) * cos(Le * Vz))
                                                      + 4.0 * Le * Vz * sin(Le * Vz))
                                                     - 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz))
                                                    + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz))
                                                   + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz))
                                                  + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) - 4.0)))
                             + 2.0 * Emod * Iy * Vz * (sz * sz) * ((2.0 * Le * Vz - 3.0 * sin(Le * Vz))
                                                                   + Le * Vz * cos(Le * Vz)) / (Le * Le * Oz * (az * az)))

    et[3 , 3] = et[9 , 9] = K / Le + Gmod*Jx / Le

    et[4 , 4] = et[10 , 10] = (((((P * (bz * bz) * (((((((((((Vz * sin(2.0 * Le * Vz) / 2.0 + power(Le, 3.0) * power(Vz, 4.0) / 2.0) + 5.0 * power(Le, 3.0) * power(Vz, 4.0) * Oz / 2.0) - Le * Le * power(Vz, 3.0) * sin(2.0 * Le * Vz) / 4.0) + power(Le, 5.0) * power(Vz, 6.0) * (Oz * Oz)) + Le * (Vz * Vz) * cos(Le * Vz)) - Le * (Vz * Vz) * cos(2.0 * Le * Vz)) - 2.0 * power(Le, 3.0) * power(Vz, 4.0) * Oz * cos(Le * Vz)) - power(Le, 3.0) * power(Vz, 4.0) * Oz * cos(2.0 * Le * Vz) / 2.0) + Le * Le * power(Vz, 3.0) * Oz * sin(2.0 * Le * Vz)) - power(Le, 5.0) * power(Vz, 6.0) * (Oz * Oz) * cos(Le * Vz)) + power(Le, 4.0) * power(Vz, 5.0) * (Oz * Oz) * sin(2.0 * Le * Vz) / 2.0) - P * sin(Le * Vz) * (bz * bz) * ((((power(Le, 4.0) * power(Vz, 5.0) * (Oz * Oz) + power(Le, 4.0) * power(Vz, 5.0) * Oz) + 2.0 * (Le * Le) * power(Vz, 3.0) * Oz) + Le * Le * power(Vz, 3.0)) + Vz)) / (Vz * Vz * (cz * cz)) - (((((Emod * Iy * (sin(Le * Vz) - sin(2.0 * Le * Vz) / 2.0) / Vz - Emod * Iy * Le * (cos(Le * Vz) - cos(2.0 * Le * Vz))) - Emod * Iy * power(Le, 3.0) * (Vz * Vz) * (((5.0 * Oz - 4.0 * Oz * cos(Le * Vz)) - Oz * cos(2.0 * Le * Vz)) + 1.0) / 2.0) + Emod * Iy * (Le * Le) * Vz * (((4.0 * sin(Le * Vz) + sin(2.0 * Le * Vz)) + 8.0 * Oz * sin(Le * Vz)) - 4.0 * Oz * sin(2.0 * Le * Vz)) / 4.0) - Emod * Iy * power(Le, 5.0) * power(Vz, 4.0) * (Oz * Oz) * ((cos(2.0 * Le * Vz) / 2.0 - 3.0 * cos(Le * Vz)) + 2.5)) + Emod * Iy * power(Le, 4.0) * power(Vz, 3.0) * Oz * ((2.0 * sin(Le * Vz) + 2.0 * Oz * sin(Le * Vz)) - Oz * sin(2.0 * Le * Vz)) / 2.0) / (Le * Le * Oz * (cz * cz))) - Emod * Iy * (((((((((Le * (Vz * Vz) / 2.0 - cos(2.0 * Le * Vz) * (Oz * power(Le, 3.0) * power(Vz, 4.0) + Le * (Vz * Vz)) / 2.0) + power(Le, 3.0) * power(Vz, 4.0) * Oz / 2.0) - Le * (Vz * Vz) * (bz * bz) / 2.0) + Vz * sin(2.0 * Le * Vz) * (((power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz)) + 1.0) / 4.0) + Le * sin(Le * Vz) * (Oz * power(Le, 3.0) * power(Vz, 5.0) + Le * power(Vz, 3.0))) - Le * (Vz * Vz) * (((power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) + Le * Le * (Vz * Vz)) + 1.0) / 2.0) - Vz * sin(Le * Vz) * (bz * bz)) + Vz * sin(2.0 * Le * Vz) * (bz * bz) / 4.0) + Le * (Vz * Vz) * cos(Le * Vz) * (bz * bz)) / (cz * cz)) - Iy * P * (((((((((Le * (Vz * Vz) / 2.0 + power(Le, 3.0) * power(Vz, 4.0) * Oz / 2.0) - Le * (Vz * Vz) * (bz * bz) / 2.0) + Vz * sin(2.0 * Le * Vz) * (((power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz)) + 1.0) / 4.0) + Le * sin(Le * Vz) * (Oz * power(Le, 3.0) * power(Vz, 5.0) + Le * power(Vz, 3.0))) - Le * (Vz * Vz) * (((power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) + Le * Le * (Vz * Vz)) + 1.0) / 2.0) - Vz * sin(Le * Vz) * (bz * bz)) + Vz * sin(2.0 * Le * Vz) * (bz * bz) / 4.0) - Le * (Vz * Vz) * cos(2.0 * Le * Vz) * (Oz * (Le * Le) * (Vz * Vz) + 1.0) / 2.0) + Le * (Vz * Vz) * cos(Le * Vz) * (bz * bz)) / (A * (cz * cz))) + (Emod * Iy * (bz * bz) * (((((((((((Vz * sin(2.0 * Le * Vz) / 2.0 + power(Le, 3.0) * power(Vz, 4.0) / 2.0) + 5.0 * power(Le, 3.0) * power(Vz, 4.0) * Oz / 2.0) - Le * Le * power(Vz, 3.0) * sin(2.0 * Le * Vz) / 4.0) + power(Le, 5.0) * power(Vz, 6.0) * (Oz * Oz)) + Le * (Vz * Vz) * cos(Le * Vz)) - Le * (Vz * Vz) * cos(2.0 * Le * Vz)) - 2.0 * power(Le, 3.0) * power(Vz, 4.0) * Oz * cos(Le * Vz)) - power(Le, 3.0) * power(Vz, 4.0) * Oz * cos(2.0 * Le * Vz) / 2.0) + Le * Le * power(Vz, 3.0) * Oz * sin(2.0 * Le * Vz)) - power(Le, 5.0) * power(Vz, 6.0) * (Oz * Oz) * cos(Le * Vz)) + power(Le, 4.0) * power(Vz, 5.0) * (Oz * Oz) * sin(2.0 * Le * Vz) / 2.0) - Emod * Iy * sin(Le * Vz) * (bz * bz) * ((((power(Le, 4.0) * power(Vz, 5.0) * (Oz * Oz) + power(Le, 4.0) * power(Vz, 5.0) * Oz) + 2.0 * (Le * Le) * power(Vz, 3.0) * Oz) + Le * Le * power(Vz, 3.0)) + Vz)) / (Le * Le * (Vz * Vz) * Oz * (cz * cz))) + 2.0 * Emod * Iy * (Oz * (Le * Le) * (Vz * Vz) + 1.0) * ((((((((((((((((2.0 * sin(Le * Vz) - sin(2.0 * Le * Vz)) - power(Le, 3.0) * power(Vz, 3.0)) - 5.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz) + 2.0 * (Le * Le) * (Vz * Vz) * sin(Le * Vz)) + Le * Le * (Vz * Vz) * sin(2.0 * Le * Vz) / 2.0) - 2.0 * Le * Vz * cos(Le * Vz)) + 2.0 * Le * Vz * cos(2.0 * Le * Vz)) - 2.0 * power(Le, 5.0) * power(Vz, 5.0) * (Oz * Oz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * cos(Le * Vz)) + power(Le, 3.0) * power(Vz, 3.0) * Oz * cos(2.0 * Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * sin(Le * Vz)) - 2.0 * (Le * Le) * (Vz * Vz) * Oz * sin(2.0 * Le * Vz)) + 2.0 * power(Le, 4.0) * power(Vz, 4.0) * Oz * sin(Le * Vz)) + 2.0 * power(Le, 5.0) * power(Vz, 5.0) * (Oz * Oz) * cos(Le * Vz)) + 2.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * sin(Le * Vz)) - power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * sin(2.0 * Le * Vz)) / (Le * Le * Vz * Oz * ((((((((((((((4.0 * cos(2.0 * Le * Vz) - 16.0 * cos(Le * Vz)) + Le * Le * (Vz * Vz)) + 24.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * cos(2.0 * Le * Vz)) - 8.0 * Le * Vz * sin(Le * Vz)) + 4.0 * Le * Vz * sin(2.0 * Le * Vz)) + 12.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)) - 32.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(2.0 * Le * Vz)) - 8.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(2.0 * Le * Vz)) - 16.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(2.0 * Le * Vz)) + 12.0))
    et[5 , 5] = et[11 , 11] = (((((P * (by * by) * (((((((((((Vy * sin(2.0 * Le * Vy) / 2.0 + power(Le, 3.0) * power(Vy, 4.0) / 2.0) + 5.0 * power(Le, 3.0) * power(Vy, 4.0) * Oy / 2.0) - Le * Le * power(Vy, 3.0) * sin(2.0 * Le * Vy) / 4.0) + power(Le, 5.0) * power(Vy, 6.0) * (Oy * Oy)) + Le * (Vy * Vy) * cos(Le * Vy)) - Le * (Vy * Vy) * cos(2.0 * Le * Vy)) - 2.0 * power(Le, 3.0) * power(Vy, 4.0) * Oy * cos(Le * Vy)) - power(Le, 3.0) * power(Vy, 4.0) * Oy * cos(2.0 * Le * Vy) / 2.0) + Le * Le * power(Vy, 3.0) * Oy * sin(2.0 * Le * Vy)) - power(Le, 5.0) * power(Vy, 6.0) * (Oy * Oy) * cos(Le * Vy)) + power(Le, 4.0) * power(Vy, 5.0) * (Oy * Oy) * sin(2.0 * Le * Vy) / 2.0) - P * sin(Le * Vy) * (by * by) * ((((power(Le, 4.0) * power(Vy, 5.0) * (Oy * Oy) + power(Le, 4.0) * power(Vy, 5.0) * Oy) + 2.0 * (Le * Le) * power(Vy, 3.0) * Oy) + Le * Le * power(Vy, 3.0)) + Vy)) / (Vy * Vy * (cy * cy)) - (((((Emod * Iz * (sin(Le * Vy) - sin(2.0 * Le * Vy) / 2.0) / Vy - Emod * Iz * Le * (cos(Le * Vy) - cos(2.0 * Le * Vy))) - Emod * Iz * power(Le, 3.0) * (Vy * Vy) * (((5.0 * Oy - 4.0 * Oy * cos(Le * Vy)) - Oy * cos(2.0 * Le * Vy)) + 1.0) / 2.0) + Emod * Iz * (Le * Le) * Vy * (((4.0 * sin(Le * Vy) + sin(2.0 * Le * Vy)) + 8.0 * Oy * sin(Le * Vy)) - 4.0 * Oy * sin(2.0 * Le * Vy)) / 4.0) - Emod * Iz * power(Le, 5.0) * power(Vy, 4.0) * (Oy * Oy) * ((cos(2.0 * Le * Vy) / 2.0 - 3.0 * cos(Le * Vy)) + 2.5)) + Emod * Iz * power(Le, 4.0) * power(Vy, 3.0) * Oy * ((2.0 * sin(Le * Vy) + 2.0 * Oy * sin(Le * Vy)) - Oy * sin(2.0 * Le * Vy)) / 2.0) / (Le * Le * Oy * (cy * cy))) - Emod * Iz * (((((((((Le * (Vy * Vy) / 2.0 - cos(2.0 * Le * Vy) * (Oy * power(Le, 3.0) * power(Vy, 4.0) + Le * (Vy * Vy)) / 2.0) + power(Le, 3.0) * power(Vy, 4.0) * Oy / 2.0) - Le * (Vy * Vy) * (by * by) / 2.0) + Vy * sin(2.0 * Le * Vy) * (((power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy)) + 1.0) / 4.0) + Le * sin(Le * Vy) * (Oy * power(Le, 3.0) * power(Vy, 5.0) + Le * power(Vy, 3.0))) - Le * (Vy * Vy) * (((power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) + Le * Le * (Vy * Vy)) + 1.0) / 2.0) - Vy * sin(Le * Vy) * (by * by)) + Vy * sin(2.0 * Le * Vy) * (by * by) / 4.0) + Le * (Vy * Vy) * cos(Le * Vy) * (by * by)) / (cy * cy)) - Iz * P * (((((((((Le * (Vy * Vy) / 2.0 + power(Le, 3.0) * power(Vy, 4.0) * Oy / 2.0) - Le * (Vy * Vy) * (by * by) / 2.0) + Vy * sin(2.0 * Le * Vy) * (((power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy)) + 1.0) / 4.0) + Le * sin(Le * Vy) * (Oy * power(Le, 3.0) * power(Vy, 5.0) + Le * power(Vy, 3.0))) - Le * (Vy * Vy) * (((power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) + Le * Le * (Vy * Vy)) + 1.0) / 2.0) - Vy * sin(Le * Vy) * (by * by)) + Vy * sin(2.0 * Le * Vy) * (by * by) / 4.0) - Le * (Vy * Vy) * cos(2.0 * Le * Vy) * (Oy * (Le * Le) * (Vy * Vy) + 1.0) / 2.0) + Le * (Vy * Vy) * cos(Le * Vy) * (by * by)) / (A * (cy * cy))) + (Emod * Iz * (by * by) * (((((((((((Vy * sin(2.0 * Le * Vy) / 2.0 + power(Le, 3.0) * power(Vy, 4.0) / 2.0) + 5.0 * power(Le, 3.0) * power(Vy, 4.0) * Oy / 2.0) - Le * Le * power(Vy, 3.0) * sin(2.0 * Le * Vy) / 4.0) + power(Le, 5.0) * power(Vy, 6.0) * (Oy * Oy)) + Le * (Vy * Vy) * cos(Le * Vy)) - Le * (Vy * Vy) * cos(2.0 * Le * Vy)) - 2.0 * power(Le, 3.0) * power(Vy, 4.0) * Oy * cos(Le * Vy)) - power(Le, 3.0) * power(Vy, 4.0) * Oy * cos(2.0 * Le * Vy) / 2.0) + Le * Le * power(Vy, 3.0) * Oy * sin(2.0 * Le * Vy)) - power(Le, 5.0) * power(Vy, 6.0) * (Oy * Oy) * cos(Le * Vy)) + power(Le, 4.0) * power(Vy, 5.0) * (Oy * Oy) * sin(2.0 * Le * Vy) / 2.0) - Emod * Iz * sin(Le * Vy) * (by * by) * ((((power(Le, 4.0) * power(Vy, 5.0) * (Oy * Oy) + power(Le, 4.0) * power(Vy, 5.0) * Oy) + 2.0 * (Le * Le) * power(Vy, 3.0) * Oy) + Le * Le * power(Vy, 3.0)) + Vy)) / (Le * Le * (Vy * Vy) * Oy * (cy * cy))) + 2.0 * Emod * Iz * (Oy * (Le * Le) * (Vy * Vy) + 1.0) * ((((((((((((((((2.0 * sin(Le * Vy) - sin(2.0 * Le * Vy)) - power(Le, 3.0) * power(Vy, 3.0)) - 5.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy) + 2.0 * (Le * Le) * (Vy * Vy) * sin(Le * Vy)) + Le * Le * (Vy * Vy) * sin(2.0 * Le * Vy) / 2.0) - 2.0 * Le * Vy * cos(Le * Vy)) + 2.0 * Le * Vy * cos(2.0 * Le * Vy)) - 2.0 * power(Le, 5.0) * power(Vy, 5.0) * (Oy * Oy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(Le * Vy)) + power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(2.0 * Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * sin(Le * Vy)) - 2.0 * (Le * Le) * (Vy * Vy) * Oy * sin(2.0 * Le * Vy)) + 2.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * sin(Le * Vy)) + 2.0 * power(Le, 5.0) * power(Vy, 5.0) * (Oy * Oy) * cos(Le * Vy)) + 2.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * sin(Le * Vy)) - power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * sin(2.0 * Le * Vy)) / (Le * Le * Vy * Oy * ((((((((((((((4.0 * cos(2.0 * Le * Vy) - 16.0 * cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 24.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * sin(Le * Vy)) + 4.0 * Le * Vy * sin(2.0 * Le * Vy)) + 12.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)) - 32.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(2.0 * Le * Vy)) - 8.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(2.0 * Le * Vy)) - 16.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(2.0 * Le * Vy)) + 12.0))

    et[2 , 3] = et[3 , 2] = Mza / Le
    et[3 , 4] = et[4 , 3] = (Oz * (Le * Le) * (Vz * Vz) + 1.0) * (Mza + Mzb) * ((((((Le * Le * (Vz * Vz) + 8.0 * (sz * sz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * (2.0 * (sz * sz) - 1.0)) - 4.0 * Le * Vz * sin(Le * Vz)) - 2.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (2.0 * (sz * sz) - 1.0)) / (Le * Le * (Vz * Vz) * (((4.0 * (sz * sz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Vz * sin(Le * Vz)) + 2.0 * (Le * Le) * (Vz * Vz) * Oz * (2.0 * (sz * sz) - 1.0))) - Mza / 2.0
    et[5 , 6] = et[6 , 5] = Mza / Le
    et[8 , 9] = et[9 , 8] = -Mzb / Le
    et[9 , 10] = et[10 , 9] = (Oz * (Le * Le) * (Vz * Vz) + 1.0) * (Mza + Mzb) * ((((((Le * Le * (Vz * Vz) + 8.0 * (sz * sz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * (2.0 * (sz * sz) - 1.0)) - 4.0 * Le * Vz * sin(Le * Vz)) - 2.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (2.0 * (sz * sz) - 1.0)) / (Le * Le * (Vz * Vz) * (((4.0 * (sz * sz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Vz * sin(Le * Vz)) + 2.0 * (Le * Le) * (Vz * Vz) * Oz * (2.0 * (sz * sz) - 1.0))) - Mzb / 2.0

    et[1 , 3] = et[3 , 1] = Mya / Le
    et[2 , 4] = et[4 , 2] = ((P * (bz * bz) * (((4.0 * cos(Le * Vz) + Le * Le * (Vz * Vz)) + Le * Vz * sin(Le * Vz)) - 4.0) / (2.0 * (((((((((4.0 * cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * cos(Le * Vz)) + 4.0 * Le * Vz * sin(Le * Vz)) - 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) - 4.0)) + Emod * Iy * Le * power(Vz, 3.0) * (sz * sz) * (sin(Le * Vz) - Le * Vz) / (az * az)) + Emod * Iy * power(Le, 3.0) * power(Vz, 5.0) * Oz * (sin(Le * Vz) + Le * Vz) / (2.0 * (((((((((4.0 * cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * cos(Le * Vz)) + 4.0 * Le * Vz * sin(Le * Vz)) - 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) - 4.0))) + Iy * Le * power(Vz, 3.0) * P * (sz * sz) * (sin(Le * Vz) - Le * Vz) / (A * (az * az))
    et[3 , 5] = et[5 , 3] = Mya / 2.0 - (Oy * (Le * Le) * (Vy * Vy) + 1.0) * (Mya + Myb) * ((((((Le * Le * (Vy * Vy) + 8.0 * (sy * sy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * (2.0 * (sy * sy) - 1.0)) - 4.0 * Le * Vy * sin(Le * Vy)) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (2.0 * (sy * sy) - 1.0)) / (Le * Le * (Vy * Vy) * (((4.0 * (sy * sy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Vy * sin(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * (2.0 * (sy * sy) - 1.0)))
    et[4 , 6] = et[6 , 4] = Mya / Le
    et[5 , 7] = et[7 , 5] = ((P * (by * by) * (((4.0 * cos(Le * Vy) + Le * Le * (Vy * Vy)) + Le * Vy * sin(Le * Vy)) - 4.0) / (2.0 * (((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * cos(Le * Vy)) + 4.0 * Le * Vy * sin(Le * Vy)) - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) - 4.0)) + Emod * Iz * Le * power(Vy, 3.0) * (sy * sy) * (sin(Le * Vy) - Le * Vy) / (ay * ay)) + Emod * Iz * power(Le, 3.0) * power(Vy, 5.0) * Oy * (sin(Le * Vy) + Le * Vy) / (2.0 * (((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * cos(Le * Vy)) + 4.0 * Le * Vy * sin(Le * Vy)) - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) - 4.0))) + Iz * Le * power(Vy, 3.0) * P * (sy * sy) * (sin(Le * Vy) - Le * Vy) / (A * (ay * ay))
    et[7 , 9] = et[9 , 7] = -Myb / Le
    et[8 , 10] = et[10 , 8] = ((-(P * (bz * bz) * (((4.0 * cos(Le * Vz) + Le * Le * (Vz * Vz)) + Le * Vz * sin(Le * Vz)) - 4.0)) / (2.0 * (((((((((4.0 * cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * cos(Le * Vz)) + 4.0 * Le * Vz * sin(Le * Vz)) - 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) - 4.0)) - Emod * Iy * Le * power(Vz, 3.0) * (sz * sz) * (sin(Le * Vz) - Le * Vz) / (az * az)) - Emod * Iy * power(Le, 3.0) * power(Vz, 5.0) * Oz * (sin(Le * Vz) + Le * Vz) / (2.0 * (((((((((4.0 * cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * cos(Le * Vz)) + 4.0 * Le * Vz * sin(Le * Vz)) - 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) - 4.0))) - Iy * Le * power(Vz, 3.0) * P * (sz * sz) * (sin(Le * Vz) - Le * Vz) / (A * (az * az))
    et[9 , 11] = et[11 , 9] = Myb / 2.0 - (Oy * (Le * Le) * (Vy * Vy) + 1.0) * (Mya + Myb) * ((((((Le * Le * (Vy * Vy) + 8.0 * (sy * sy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * (2.0 * (sy * sy) - 1.0)) - 4.0 * Le * Vy * sin(Le * Vy)) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (2.0 * (sy * sy) - 1.0)) / (Le * Le * (Vy * Vy) * (((4.0 * (sy * sy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Vy * sin(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * (2.0 * (sy * sy) - 1.0)))

    et[0 , 4] = et[4 , 0] = -Mya / Le
    et[1 , 5] = et[5 , 1] = ((-(P * (by * by) * (((4.0 * cos(Le * Vy) + Le * Le * (Vy * Vy)) + Le * Vy * sin(Le * Vy)) - 4.0)) / (2.0 * (((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * cos(Le * Vy)) + 4.0 * Le * Vy * sin(Le * Vy)) - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) - 4.0)) - Emod * Iz * Le * power(Vy, 3.0) * (sy * sy) * (sin(Le * Vy) - Le * Vy) / (ay * ay)) - Emod * Iz * power(Le, 3.0) * power(Vy, 5.0) * Oy * (sin(Le * Vy) + Le * Vy) / (2.0 * (((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * cos(Le * Vy)) + 4.0 * Le * Vy * sin(Le * Vy)) - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) - 4.0))) - Iz * Le * power(Vy, 3.0) * P * (sy * sy) * (sin(Le * Vy) - Le * Vy) / (A * (ay * ay))
    et[3 , 7] = et[7 , 3] = -Mya / Le
    et[4 , 8] = et[8 , 4] = ((-(P * (bz * bz) * (((4.0 * cos(Le * Vz) + Le * Le * (Vz * Vz)) + Le * Vz * sin(Le * Vz)) - 4.0)) / (2.0 * (((((((((4.0 * cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * cos(Le * Vz)) + 4.0 * Le * Vz * sin(Le * Vz)) - 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) - 4.0)) - Emod * Iy * Le * power(Vz, 3.0) * (sz * sz) * (sin(Le * Vz) - Le * Vz) / (az * az)) - Emod * Iy * power(Le, 3.0) * power(Vz, 5.0) * Oz * (sin(Le * Vz) + Le * Vz) / (2.0 * (((((((((4.0 * cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * cos(Le * Vz)) + 4.0 * Le * Vz * sin(Le * Vz)) - 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) - 4.0))) - Iy * Le * power(Vz, 3.0) * P * (sz * sz) * (sin(Le * Vz) - Le * Vz) / (A * (az * az))
    et[5 , 9] = et[9 , 5] = 2.0 * (Oy * (Le * Le) * (Vy * Vy) + 1.0) * (Mya + Myb) * (((((Le * Le * (Vy * Vy) + 4.0 * (sy * sy)) - 2.0 * Le * Vy * sin(Le * Vy)) - Le * Le * (Vy * Vy) * (sy * sy)) - power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) / (Le * Le * (Vy * Vy) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)))
    et[6 , 10] = et[10 , 6] = Myb / Le
    et[7 , 11] = et[11 , 7] = ((P * (by * by) * (((4.0 * cos(Le * Vy) + Le * Le * (Vy * Vy)) + Le * Vy * sin(Le * Vy)) - 4.0) / (2.0 * (((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * cos(Le * Vy)) + 4.0 * Le * Vy * sin(Le * Vy)) - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) - 4.0)) + Emod * Iz * Le * power(Vy, 3.0) * (sy * sy) * (sin(Le * Vy) - Le * Vy) / (ay * ay)) + Emod * Iz * power(Le, 3.0) * power(Vy, 5.0) * Oy * (sin(Le * Vy) + Le * Vy) / (2.0 * (((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * cos(Le * Vy)) + 4.0 * Le * Vy * sin(Le * Vy)) - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) - 4.0))) + Iz * Le * power(Vy, 3.0) * P * (sy * sy) * (sin(Le * Vy) - Le * Vy) / (A * (ay * ay))

    et[0 , 5] = et[5 , 0] = -Mza / Le
    et[3 , 8] = et[8 , 3] = -Mza / Le
    et[4 , 9] = et[9 , 4] = -(2.0 * (Oz * (Le * Le) * (Vz * Vz) + 1.0) * (Mza + Mzb) * (((((Le * Le * (Vz * Vz) + 4.0 * (sz * sz)) - 2.0 * Le * Vz * sin(Le * Vz)) - Le * Le * (Vz * Vz) * (sz * sz)) - power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz))) / (Le * Le * (Vz * Vz) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    et[6 , 11] = et[11 , 6] = Mzb / Le

    et[0 , 6] = et[6 , 0] = -P / Le - A * Emod / Le  
    et[1 , 7] = et[7 , 1] = (((2.0 * Emod * Iz * power(Vy, 3.0) * (sy * sy) * (sin(Le * Vy) - Le * Vy) / (ay * ay) - (sy * sy * (2.0 * Le * (Vy * Vy) * P * (((A * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) + 2.0 * A * (Le * Le) * (Vy * Vy) * Oy) + Iz * (Vy * Vy)) + A) - Vy * P * sin(Le * Vy) * (((-2.0 * A * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) + 4.0 * A * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Iz * (Vy * Vy)) + 6.0 * A)) + A * Le * (Vy * Vy) * P * (sin(Le * Vy) * sin(Le * Vy))) / (A * (ay * ay))) - Emod * Iz * Vy * ((((((Le * Vy * (sin(Le * Vy) * sin(Le * Vy)) - 6.0 * sin(Le * Vy) * (sy * sy)) + 2.0 * Le * Vy * (sy * sy)) + 2.0 * power(Le, 5.0) * power(Vy, 5.0) * (Oy * Oy) * (sy * sy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * (sy * sy)) - 4.0 * (Le * Le) * (Vy * Vy) * Oy * sin(Le * Vy) * (sy * sy)) + 2.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * sin(Le * Vy) * (sy * sy)) / (Le * Le * Oy * (ay * ay))) - 2.0 * Emod * Iz * Vy * ((((2.0 * Le * Vy - 3.0 * sin(Le * Vy)) + power(Le, 3.0) * power(Vy, 3.0) * Oy) + Le * Vy * cos(Le * Vy)) - Le * Le * (Vy * Vy) * Oy * sin(Le * Vy)) / (Le * Le * Oy * (((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * cos(Le * Vy)) + 4.0 * Le * Vy * sin(Le * Vy)) - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) - 4.0))) - 2.0 * Emod * Iz * Vy * (sy * sy) * ((2.0 * Le * Vy - 3.0 * sin(Le * Vy)) + Le * Vy * cos(Le * Vy)) / (Le * Le * Oy * (ay * ay))
    et[2 , 8] = et[8 , 2] = (((2.0 * Emod * Iy * power(Vz, 3.0) * (sz * sz) * (sin(Le * Vz) - Le * Vz) / (az * az) - (sz * sz * (2.0 * Le * (Vz * Vz) * P * (((A * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) + 2.0 * A * (Le * Le) * (Vz * Vz) * Oz) + Iy * (Vz * Vz)) + A) - Vz * P * sin(Le * Vz) * (((-2.0 * A * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) + 4.0 * A * (Le * Le) * (Vz * Vz) * Oz) + 2.0 * Iy * (Vz * Vz)) + 6.0 * A)) + A * Le * (Vz * Vz) * P * (sin(Le * Vz) * sin(Le * Vz))) / (A * (az * az))) - Emod * Iy * Vz * ((((((Le * Vz * (sin(Le * Vz) * sin(Le * Vz)) - 6.0 * sin(Le * Vz) * (sz * sz)) + 2.0 * Le * Vz * (sz * sz)) + 2.0 * power(Le, 5.0) * power(Vz, 5.0) * (Oz * Oz) * (sz * sz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * (sz * sz)) - 4.0 * (Le * Le) * (Vz * Vz) * Oz * sin(Le * Vz) * (sz * sz)) + 2.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * sin(Le * Vz) * (sz * sz)) / (Le * Le * Oz * (az * az))) - 2.0 * Emod * Iy * Vz * ((((2.0 * Le * Vz - 3.0 * sin(Le * Vz)) + power(Le, 3.0) * power(Vz, 3.0) * Oz) + Le * Vz * cos(Le * Vz)) - Le * Le * (Vz * Vz) * Oz * sin(Le * Vz)) / (Le * Le * Oz * (((((((((4.0 * cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * cos(Le * Vz)) + 4.0 * Le * Vz * sin(Le * Vz)) - 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) - 4.0))) - 2.0 * Emod * Iy * Vz * (sz * sz) * ((2.0 * Le * Vz - 3.0 * sin(Le * Vz)) + Le * Vz * cos(Le * Vz)) / (Le * Le * Oz * (az * az))
    et[3 , 9] = et[9 , 3] = -(K / Le) - (Gmod*Jx / Le)
    et[4 , 10] = et[10 , 4] = (((((((((((((((((((((((((((((((((((((((((((((((((ez - 2.0 * Iy * power(Le, 5.0) * power(Vz, 7.0) * (Oz * Oz) * P) + 2.0 * dz * sin(Le * Vz)) - dz * sin(2.0 * Le * Vz)) + 8.0 * A * power(Le, 4.0) * power(Vz, 4.0) * Oz * P * sin(Le * Vz)) - 2.0 * fz * P * sin(Le * Vz)) + fz * P * sin(2.0 * Le * Vz)) + 2.0 * Iy * power(Le, 4.0) * power(Vz, 6.0) * Oz * P * sin(Le * Vz)) - 7.0 * A * Emod * Iy * power(Le, 5.0) * power(Vz, 7.0) * (Oz * Oz)) - 2.0 * A * Emod * Iy * power(Le, 7.0) * power(Vz, 9.0) * power(Oz, 3.0)) + 16.0 * A * power(Le, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P * cos(Le * Vz)) + 2.0 * A * power(Le, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P * cos(2.0 * Le * Vz)) - A * power(Le, 7.0) * power(Vz, 7.0) * (Oz * Oz) * P * cos(Le * Vz)) + 8.0 * A * power(Le, 7.0) * power(Vz, 7.0) * power(Oz, 3.0) * P * cos(Le * Vz)) + A * power(Le, 7.0) * power(Vz, 7.0) * power(Oz, 3.0) * P * cos(2.0 * Le * Vz)) + 2.0 * A * power(Le, 9.0) * power(Vz, 9.0) * power(Oz, 4.0) * P * cos(Le * Vz)) + 2.0 * Iy * power(Le, 5.0) * power(Vz, 7.0) * (Oz * Oz) * P * cos(Le * Vz)) + 12.0 * A * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * sin(Le * Vz)) - 6.0 * A * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * sin(2.0 * Le * Vz)) + 7.0 * A * power(Le, 6.0) * power(Vz, 6.0) * (Oz * Oz) * P * sin(Le * Vz)) + 8.0 * A * power(Le, 6.0) * power(Vz, 6.0) * power(Oz, 3.0) * P * sin(Le * Vz)) - 4.0 * A * power(Le, 6.0) * power(Vz, 6.0) * power(Oz, 3.0) * P * sin(2.0 * Le * Vz)) + 2.0 * A * power(Le, 8.0) * power(Vz, 8.0) * power(Oz, 3.0) * P * sin(Le * Vz)) + 2.0 * A * power(Le, 8.0) * power(Vz, 8.0) * power(Oz, 4.0) * P * sin(Le * Vz)) - A * power(Le, 8.0) * power(Vz, 8.0) * power(Oz, 4.0) * P * sin(2.0 * Le * Vz)) - 2.0 * Iy * power(Le, 4.0) * power(Vz, 6.0) * (Oz * Oz) * P * sin(Le * Vz)) + Iy * power(Le, 4.0) * power(Vz, 6.0) * (Oz * Oz) * P * sin(2.0 * Le * Vz)) + 2.0 * A * Emod * Iy * Le * power(Vz, 3.0) * cos(Le * Vz)) - A * Emod * Iy * Le * power(Vz, 3.0) * cos(2.0 * Le * Vz)) - 6.0 * A * Emod * Iy * power(Le, 3.0) * power(Vz, 5.0) * Oz) - A * Emod * Iy * power(Le, 3.0) * power(Vz, 5.0) * cos(Le * Vz)) + A * Emod * Iy * (Le * Le) * power(Vz, 4.0) * sin(Le * Vz)) + 16.0 * A * power(Le, 3.0) * power(Vz, 3.0) * Oz * P * cos(Le * Vz)) + A * power(Le, 3.0) * power(Vz, 3.0) * Oz * P * cos(2.0 * Le * Vz)) - 2.0 * A * power(Le, 5.0) * power(Vz, 5.0) * Oz * P * cos(Le * Vz)) + 4.0 * Iy * power(Le, 3.0) * power(Vz, 5.0) * Oz * P * cos(Le * Vz)) - Iy * power(Le, 3.0) * power(Vz, 5.0) * Oz * P * cos(2.0 * Le * Vz)) + 6.0 * A * Emod * Iy * power(Le, 3.0) * power(Vz, 5.0) * Oz * cos(Le * Vz)) - A * Emod * Iy * power(Le, 5.0) * power(Vz, 7.0) * Oz * cos(Le * Vz)) - 2.0 * A * Emod * Iy * (Le * Le) * power(Vz, 4.0) * Oz * sin(Le * Vz)) + A * Emod * Iy * (Le * Le) * power(Vz, 4.0) * Oz * sin(2.0 * Le * Vz)) + 5.0 * A * Emod * Iy * power(Le, 4.0) * power(Vz, 6.0) * Oz * sin(Le * Vz)) + 6.0 * A * Emod * Iy * power(Le, 5.0) * power(Vz, 7.0) * (Oz * Oz) * cos(Le * Vz)) + A * Emod * Iy * power(Le, 5.0) * power(Vz, 7.0) * (Oz * Oz) * cos(2.0 * Le * Vz)) + 2.0 * A * Emod * Iy * power(Le, 7.0) * power(Vz, 9.0) * power(Oz, 3.0) * cos(Le * Vz)) + 2.0 * A * Emod * Iy * power(Le, 4.0) * power(Vz, 6.0) * (Oz * Oz) * sin(Le * Vz)) - A * Emod * Iy * power(Le, 4.0) * power(Vz, 6.0) * (Oz * Oz) * sin(2.0 * Le * Vz)) + 2.0 * A * Emod * Iy * power(Le, 6.0) * power(Vz, 8.0) * (Oz * Oz) * sin(Le * Vz)) + 2.0 * A * Emod * Iy * power(Le, 6.0) * power(Vz, 8.0) * power(Oz, 3.0) * sin(Le * Vz)) - A * Emod * Iy * power(Le, 6.0) * power(Vz, 8.0) * power(Oz, 3.0) * sin(2.0 * Le * Vz)) / (A * Vz * ((((((((((((((4.0 * cos(2.0 * Le * Vz) - 16.0 * cos(Le * Vz)) + Le * Le * (Vz * Vz)) + 24.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * cos(2.0 * Le * Vz)) - 8.0 * Le * Vz * sin(Le * Vz)) + 4.0 * Le * Vz * sin(2.0 * Le * Vz)) + 12.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)) - 32.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(2.0 * Le * Vz)) - 8.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(2.0 * Le * Vz)) - 16.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(2.0 * Le * Vz)) + 12.0))    
    et[5 , 11] = et[11 , 5] = (((((((((((((((((((((((((((((((((((((((((((((((((ey - 2.0 * Iz * power(Le, 5.0) * power(Vy, 7.0) * (Oy * Oy) * P) + 2.0 * dy * sin(Le * Vy)) - dy * sin(2.0 * Le * Vy)) + 8.0 * A * power(Le, 4.0) * power(Vy, 4.0) * Oy * P * sin(Le * Vy)) - 2.0 * fy * P * sin(Le * Vy)) + fy * P * sin(2.0 * Le * Vy)) + 2.0 * Iz * power(Le, 4.0) * power(Vy, 6.0) * Oy * P * sin(Le * Vy)) - 7.0 * A * Emod * Iz * power(Le, 5.0) * power(Vy, 7.0) * (Oy * Oy)) - 2.0 * A * Emod * Iz * power(Le, 7.0) * power(Vy, 9.0) * power(Oy, 3.0)) + 16.0 * A * power(Le, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P * cos(Le * Vy)) + 2.0 * A * power(Le, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P * cos(2.0 * Le * Vy)) - A * power(Le, 7.0) * power(Vy, 7.0) * (Oy * Oy) * P * cos(Le * Vy)) + 8.0 * A * power(Le, 7.0) * power(Vy, 7.0) * power(Oy, 3.0) * P * cos(Le * Vy)) + A * power(Le, 7.0) * power(Vy, 7.0) * power(Oy, 3.0) * P * cos(2.0 * Le * Vy)) + 2.0 * A * power(Le, 9.0) * power(Vy, 9.0) * power(Oy, 4.0) * P * cos(Le * Vy)) + 2.0 * Iz * power(Le, 5.0) * power(Vy, 7.0) * (Oy * Oy) * P * cos(Le * Vy)) + 12.0 * A * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * sin(Le * Vy)) - 6.0 * A * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * sin(2.0 * Le * Vy)) + 7.0 * A * power(Le, 6.0) * power(Vy, 6.0) * (Oy * Oy) * P * sin(Le * Vy)) + 8.0 * A * power(Le, 6.0) * power(Vy, 6.0) * power(Oy, 3.0) * P * sin(Le * Vy)) - 4.0 * A * power(Le, 6.0) * power(Vy, 6.0) * power(Oy, 3.0) * P * sin(2.0 * Le * Vy)) + 2.0 * A * power(Le, 8.0) * power(Vy, 8.0) * power(Oy, 3.0) * P * sin(Le * Vy)) + 2.0 * A * power(Le, 8.0) * power(Vy, 8.0) * power(Oy, 4.0) * P * sin(Le * Vy)) - A * power(Le, 8.0) * power(Vy, 8.0) * power(Oy, 4.0) * P * sin(2.0 * Le * Vy)) - 2.0 * Iz * power(Le, 4.0) * power(Vy, 6.0) * (Oy * Oy) * P * sin(Le * Vy)) + Iz * power(Le, 4.0) * power(Vy, 6.0) * (Oy * Oy) * P * sin(2.0 * Le * Vy)) + 2.0 * A * Emod * Iz * Le * power(Vy, 3.0) * cos(Le * Vy)) - A * Emod * Iz * Le * power(Vy, 3.0) * cos(2.0 * Le * Vy)) - 6.0 * A * Emod * Iz * power(Le, 3.0) * power(Vy, 5.0) * Oy) - A * Emod * Iz * power(Le, 3.0) * power(Vy, 5.0) * cos(Le * Vy)) + A * Emod * Iz * (Le * Le) * power(Vy, 4.0) * sin(Le * Vy)) + 16.0 * A * power(Le, 3.0) * power(Vy, 3.0) * Oy * P * cos(Le * Vy)) + A * power(Le, 3.0) * power(Vy, 3.0) * Oy * P * cos(2.0 * Le * Vy)) - 2.0 * A * power(Le, 5.0) * power(Vy, 5.0) * Oy * P * cos(Le * Vy)) + 4.0 * Iz * power(Le, 3.0) * power(Vy, 5.0) * Oy * P * cos(Le * Vy)) - Iz * power(Le, 3.0) * power(Vy, 5.0) * Oy * P * cos(2.0 * Le * Vy)) + 6.0 * A * Emod * Iz * power(Le, 3.0) * power(Vy, 5.0) * Oy * cos(Le * Vy)) - A * Emod * Iz * power(Le, 5.0) * power(Vy, 7.0) * Oy * cos(Le * Vy)) - 2.0 * A * Emod * Iz * (Le * Le) * power(Vy, 4.0) * Oy * sin(Le * Vy)) + A * Emod * Iz * (Le * Le) * power(Vy, 4.0) * Oy * sin(2.0 * Le * Vy)) + 5.0 * A * Emod * Iz * power(Le, 4.0) * power(Vy, 6.0) * Oy * sin(Le * Vy)) + 6.0 * A * Emod * Iz * power(Le, 5.0) * power(Vy, 7.0) * (Oy * Oy) * cos(Le * Vy)) + A * Emod * Iz * power(Le, 5.0) * power(Vy, 7.0) * (Oy * Oy) * cos(2.0 * Le * Vy)) + 2.0 * A * Emod * Iz * power(Le, 7.0) * power(Vy, 9.0) * power(Oy, 3.0) * cos(Le * Vy)) + 2.0 * A * Emod * Iz * power(Le, 4.0) * power(Vy, 6.0) * (Oy * Oy) * sin(Le * Vy)) - A * Emod * Iz * power(Le, 4.0) * power(Vy, 6.0) * (Oy * Oy) * sin(2.0 * Le * Vy)) + 2.0 * A * Emod * Iz * power(Le, 6.0) * power(Vy, 8.0) * (Oy * Oy) * sin(Le * Vy)) + 2.0 * A * Emod * Iz * power(Le, 6.0) * power(Vy, 8.0) * power(Oy, 3.0) * sin(Le * Vy)) - A * Emod * Iz * power(Le, 6.0) * power(Vy, 8.0) * power(Oy, 3.0) * sin(2.0 * Le * Vy)) / (A * Vy * ((((((((((((((4.0 * cos(2.0 * Le * Vy) - 16.0 * cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 24.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * sin(Le * Vy)) + 4.0 * Le * Vy * sin(2.0 * Le * Vy)) + 12.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)) - 32.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(2.0 * Le * Vy)) - 8.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(2.0 * Le * Vy)) - 16.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(2.0 * Le * Vy)) + 12.0))

    et[2 , 9] = et[9 , 2] = Mzb / Le
    et[3 , 10] = et[10 , 3] = -(2.0 * (Oz * (Le * Le) * (Vz * Vz) + 1.0) * (Mza + Mzb) * (((((Le * Le * (Vz * Vz) + 4.0 * (sz * sz)) - 2.0 * Le * Vz * sin(Le * Vz)) - Le * Le * (Vz * Vz) * (sz * sz)) - power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz))) / (Le * Le * (Vz * Vz) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    
    et[1 , 9] = et[9 , 1] = Myb / Le
    et[2 , 10] = et[10 , 2] = ((P * (bz * bz) * (((4.0 * cos(Le * Vz) + Le * Le * (Vz * Vz)) + Le * Vz * sin(Le * Vz)) - 4.0) / (2.0 * (((((((((4.0 * cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * cos(Le * Vz)) + 4.0 * Le * Vz * sin(Le * Vz)) - 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) - 4.0)) + Emod * Iy * Le * power(Vz, 3.0) * (sz * sz) * (sin(Le * Vz) - Le * Vz) / (az * az)) + Emod * Iy * power(Le, 3.0) * power(Vz, 5.0) * Oz * (sin(Le * Vz) + Le * Vz) / (2.0 * (((((((((4.0 * cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * cos(Le * Vz)) + 4.0 * Le * Vz * sin(Le * Vz)) - 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * cos(Le * Vz)) + 4.0 * power(Le, 3.0) * power(Vz, 3.0) * Oz * sin(Le * Vz)) + 4.0 * power(Le, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cos(Le * Vz)) - 4.0))) + Iy * Le * power(Vz, 3.0) * P * (sz * sz) * (sin(Le * Vz) - Le * Vz) / (A * (az * az))
    et[3 , 11] = et[11 , 3] = 2.0 * (Oy * (Le * Le) * (Vy * Vy) + 1.0) * (Mya + Myb) * (((((Le * Le * (Vy * Vy) + 4.0 * (sy * sy)) - 2.0 * Le * Vy * sin(Le * Vy)) - Le * Le * (Vy * Vy) * (sy * sy)) - power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) / (Le * Le * (Vy * Vy) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)))

    et[0 , 10] = et[10 , 0] = -Myb / Le
    et[1 , 11] = et[11 , 1] = ((-(P * (by * by) * (((4.0 * cos(Le * Vy) + Le * Le * (Vy * Vy)) + Le * Vy * sin(Le * Vy)) - 4.0)) / (2.0 * (((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * cos(Le * Vy)) + 4.0 * Le * Vy * sin(Le * Vy)) - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) - 4.0)) - Emod * Iz * Le * power(Vy, 3.0) * (sy * sy) * (sin(Le * Vy) - Le * Vy) / (ay * ay)) - Emod * Iz * power(Le, 3.0) * power(Vy, 5.0) * Oy * (sin(Le * Vy) + Le * Vy) / (2.0 * (((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * cos(Le * Vy)) + 4.0 * Le * Vy * sin(Le * Vy)) - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cos(Le * Vy)) - 4.0))) - Iz * Le * power(Vy, 3.0) * P * (sy * sy) * (sin(Le * Vy) - Le * Vy) / (A * (ay * ay))

    et[0 , 11] = et[11 , 0] = -Mzb / Le
    #
    #
    #
    #
    # Torsion for symmetric sections
    if Iy == Iz:
        et[1 , 4] = et[4 , 1] = Vy * Mxb * (sin(Le * Vy) - Le * Vy) / ((((4.0 * cos(Le * Vy) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) - 4.0)
        et[2 , 5] = et[5 , 2] = -(Vy * Mxb * (((((((((2.0 * sin(Le * Vy) - sin(2.0 * Le * Vy)) - 5.0 * Le * Vy / 2.0) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy) + Le * Le * (Vy * Vy) * sin(Le * Vy)) + 2.0 * Le * Vy * cos(Le * Vy)) + Le * Vy * cos(2.0 * Le * Vy) / 2.0) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * sin(Le * Vy)) - Le * Le * (Vy * Vy) * Oy * sin(2.0 * Le * Vy))) / (((((((((((((((((((4.0 * cos(2.0 * Le * Vy) - 16.0 * cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy) + 12.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * sin(Le * Vy)) + 4.0 * Le * Vy * sin(2.0 * Le * Vy)) + 12.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz) - 16.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) - 16.0 * (Le * Le) * (Vy * Vy) * Oz * cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * cos(2.0 * Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * cos(2.0 * Le * Vy)) - 4.0 * power(Le,3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) - 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * sin(Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(2.0 * Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * sin(2.0 * Le * Vy)) - 16.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * cos(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * cos(2.0 * Le * Vy)) + 12.0)
        et[4 , 7] = et[7 , 4] = -(Vy * Mxb * (sin(Le * Vy) - Le * Vy)) / ((((4.0 * cos(Le * Vy) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) - 4.0)
        et[5 , 8] = et[8 , 5] = Vy * Mxb * (((((((((2.0 * sin(Le * Vy) - sin(2.0 * Le * Vy)) - 5.0 * Le * Vy / 2.0) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy) + Le * Le * (Vy * Vy) * sin(Le * Vy)) + 2.0 * Le * Vy * cos(Le * Vy)) + Le * Vy * cos(2.0 * Le * Vy) / 2.0) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * sin(Le * Vy)) - Le * Le * (Vy * Vy) * Oy * sin(2.0 * Le * Vy)) / (((((((((((((((((((4.0 * cos(2.0 * Le * Vy) - 16.0 * cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy) + 12.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * sin(Le * Vy)) + 4.0 * Le * Vy * sin(2.0 * Le * Vy)) + 12.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz) - 16.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) - 16.0 * (Le * Le) * (Vy * Vy) * Oz * cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * cos(2.0 * Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * cos(2.0 * Le * Vy)) - 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) - 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * sin(Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(2.0 * Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * sin(2.0 * Le * Vy)) - 16.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * cos(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * cos(2.0 * Le * Vy)) + 12.0)
        et[7 , 10] = et[10 , 7] = Vy * Mxb * (sin(Le * Vy) - Le * Vy) / ((((4.0 * cos(Le * Vy) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) - 4.0)
        et[8 , 11] = et[11 , 8] = -(Vy * Mxb * (((((((((2.0 * sin(Le * Vy) - sin(2.0 * Le * Vy)) - 5.0 * Le * Vy / 2.0) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy) + Le * Le * (Vy * Vy) * sin(Le * Vy)) + 2.0 * Le * Vy * cos(Le * Vy)) + Le * Vy * cos(2.0 * Le * Vy) / 2.0) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * sin(Le * Vy)) - Le * Le * (Vy * Vy) * Oy * sin(2.0 * Le * Vy))) / (((((((((((((((((((4.0 * cos(2.0 * Le * Vy) - 16.0 * cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy) + 12.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * sin(Le * Vy)) + 4.0 * Le * Vy * sin(2.0 * Le * Vy)) + 12.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz) - 16.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) - 16.0 * (Le * Le) * (Vy * Vy) * Oz * cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * cos(2.0 * Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * cos(2.0 * Le * Vy)) - 4.0 * power(Le,3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) - 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * sin(Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(2.0 * Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * sin(2.0 * Le * Vy)) - 16.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * cos(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * cos(2.0 * Le * Vy)) + 12.0)

        et[5 , 10] = et[10 , 5] = -(Mxb * (((((((((((((((((4.0 * sin(3.0 * Le * Vy / 2.0) - 12.0 * sin(Le * Vy / 2.0)) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * cos(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * sin(Le * Vy / 2.0)) + 2.0 * Le * Vy * cos(Le * Vy / 2.0)) - 2.0 * Le * Vy * cos(3.0 * Le * Vy / 2.0)) + power(Le,3.0) * power(Vy, 3.0) * Oy * cos(Le * Vy / 2.0)) + power(Le, 3.0) * power(Vy, 3.0) * Oz * cos(Le * Vy / 2.0)) - power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(3.0 * Le * Vy / 2.0)) - power(Le, 3.0) * power(Vy, 3.0) * Oz * cos(3.0 * Le * Vy / 2.0)) - 12.0 * (Le * Le) * (Vy * Vy) * Oy * sin(Le * Vy / 2.0)) - 12.0 * (Le * Le) * (Vy * Vy) * Oz * sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * sin(3.0 * Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * sin(3.0 * Le * Vy / 2.0)) + 2.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * sin(Le * Vy / 2.0)) + 2.0 * power(Le, 4.0) * power(Vy, 4.0) * Oz * sin(Le * Vy / 2.0)) - 12.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * sin(Le * Vy / 2.0)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * sin(3.0 * Le * Vy / 2.0))) / (2.0 * (((((((((((((((12.0 * sin(Le * Vy / 2.0) - 4.0 * sin(3.0 * Le * Vy / 2.0)) + Le * Le * (Vy * Vy) * sin(Le * Vy / 2.0)) + Le * Le * (Vy * Vy) * sin(3.0 * Le * Vy / 2.0)) - 4.0 * Le * Vy * cos(Le * Vy / 2.0)) + 4.0 * Le * Vy * cos(3.0 * Le * Vy / 2.0)) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(Le * Vy / 2.0)) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * cos(Le * Vy / 2.0)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(3.0 * Le * Vy / 2.0)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * cos(3.0 * Le * Vy / 2.0)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy * sin(Le * Vy / 2.0)) + 12.0 * (Le * Le) * (Vy * Vy) * Oz * sin(Le * Vy / 2.0)) - 4.0 * (Le * Le) * (Vy * Vy) * Oy * sin(3.0 * Le * Vy / 2.0)) - 4.0 * (Le * Le) * (Vy * Vy) * Oz * sin(3.0 * Le * Vy / 2.0)) + 12.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * sin(Le * Vy / 2.0)) - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * sin(3.0 * Le * Vy / 2.0)))
        et[4 , 11] = et[11 , 4] = Mxb * (((((((((((((((((4.0 * sin(3.0 * Le * Vy / 2.0) - 12.0 * sin(Le * Vy / 2.0)) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * cos(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * sin(Le * Vy / 2.0)) + 2.0 * Le * Vy * cos(Le * Vy / 2.0)) - 2.0 * Le * Vy * cos(3.0 * Le * Vy / 2.0)) + power(Le,3.0) * power(Vy, 3.0) * Oy * cos(Le * Vy / 2.0)) + power(Le, 3.0) * power(Vy, 3.0) * Oz * cos(Le * Vy / 2.0)) - power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(3.0 * Le * Vy / 2.0)) - power(Le, 3.0) * power(Vy, 3.0) * Oz * cos(3.0 * Le * Vy / 2.0)) - 12.0 * (Le * Le) * (Vy * Vy) * Oy * sin(Le * Vy / 2.0)) - 12.0 * (Le * Le) * (Vy * Vy) * Oz * sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * sin(3.0 * Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * sin(3.0 * Le * Vy / 2.0)) + 2.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * sin(Le * Vy / 2.0)) + 2.0 * power(Le, 4.0) * power(Vy, 4.0) * Oz * sin(Le * Vy / 2.0)) - 12.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * sin(Le * Vy / 2.0)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * sin(3.0 * Le * Vy / 2.0)) / (2.0 * (((((((((((((((12.0 * sin(Le * Vy / 2.0) - 4.0 * sin(3.0 * Le * Vy / 2.0)) + Le * Le * (Vy * Vy) * sin(Le * Vy / 2.0)) + Le * Le * (Vy * Vy) * sin(3.0 * Le * Vy / 2.0)) - 4.0 * Le * Vy * cos(Le * Vy / 2.0)) + 4.0 * Le * Vy * cos(3.0 * Le * Vy / 2.0)) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(Le * Vy / 2.0)) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * cos(Le * Vy / 2.0)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(3.0 * Le * Vy / 2.0)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * cos(3.0 * Le * Vy / 2.0)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy * sin(Le * Vy / 2.0)) + 12.0 * (Le * Le) * (Vy * Vy) * Oz * sin(Le * Vy / 2.0)) - 4.0 * (Le * Le) * (Vy * Vy) * Oy * sin(3.0 * Le * Vy / 2.0)) - 4.0 * (Le * Le) * (Vy * Vy) * Oz * sin(3.0 * Le * Vy / 2.0)) + 12.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * sin(Le * Vy / 2.0)) - 4.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * sin(3.0 * Le * Vy / 2.0)))

        et[1 , 10] = et[10 , 1] = -(Vy * Mxb * (sin(Le * Vy) - Le * Vy)) / ((((4.0 * cos(Le * Vy) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) - 4.0)
        et[2 , 11] = et[11 , 2] = Vy * Mxb * (((((((((2.0 * sin(Le * Vy) - sin(2.0 * Le * Vy)) - 5.0 * Le * Vy / 2.0) - 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy) + Le * Le * (Vy * Vy) * sin(Le * Vy)) + 2.0 * Le * Vy * cos(Le * Vy)) + Le * Vy * cos(2.0 * Le * Vy) / 2.0) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * cos(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * sin(Le * Vy)) - Le * Le * (Vy * Vy) * Oy * sin(2.0 * Le * Vy)) / (((((((((((((((((((4.0 * cos(2.0 * Le * Vy) - 16.0 * cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy) + 12.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * sin(Le * Vy)) + 4.0 * Le * Vy * sin(2.0 * Le * Vy)) + 12.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz) - 16.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) - 16.0 * (Le * Le) * (Vy * Vy) * Oz * cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * cos(2.0 * Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * cos(2.0 * Le * Vy)) - 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) - 4.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * sin(Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(2.0 * Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * sin(2.0 * Le * Vy)) - 16.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * cos(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * cos(2.0 * Le * Vy)) + 12.0)

        et[4 , 5] = et[5 , 4] = -(power(Le, 3.0) * power(Vy, 3.0) * Mxb * (Oy - Oz) * (sin(Le * Vy) - Le * Vy)) / (2.0 * ((((((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy)) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) - 4.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * cos(Le * Vy)) + 4.0 * Le * Vy * sin(Le * Vy)) - 4.0 * power(Le,4.0) * power(Vy, 4.0) * Oy * Oz) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * cos(Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * sin(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * cos(Le * Vy)) - 4.0))
        et[10 , 11] = et[11 , 10] = power(Le, 3.0) * power(Vy, 3.0) * Mxb * (Oy - Oz) * (sin(Le * Vy) - Le * Vy) / (2.0 * ((((((((((((4.0 * cos(Le * Vy) - Le * Le * (Vy * Vy)) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) - 4.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * cos(Le * Vy)) + 4.0 * Le * Vy * sin(Le * Vy)) - 4.0 * power(Le,4.0) * power(Vy, 4.0) * Oy * Oz) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * cos(Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oy * sin(Le * Vy)) + 2.0 * power(Le, 3.0) * power(Vy, 3.0) * Oz * sin(Le * Vy)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * Oy * Oz * cos(Le * Vy)) - 4.0))
    else:
        et[1 , 4] = et[4 , 1] = (Vy * Vy * ((4.0 * Vz * Mxb * (sin(Le * Vz) - 2.0 * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0)) - 4.0 * Le * (Vz * Vz) * Mxb * (c_z * c_z - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * power(Vz, 3.0) * Mxb * (Oz * sin(Le * Vz) - 2.0 * Oz * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0))) - power(Vy, 3.0) * ((4.0 * Mxb * (sin(Le * Vy) - 2.0 * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vz * Vz) * Mxb * (Oz * sin(Le * Vy) - 2.0 * Oz * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0))) - Le * Vz * Mxb * sin(Le * Vy) * sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
        et[2 , 5] = et[5 , 2] = -(Vz * Vz * ((4.0 * Vy * Mxb * (sin(Le * Vy) - 2.0 * cos(Le * Vy / 2.0) * (c_y * c_y) * sin(Le * Vy / 2.0)) - 4.0 * Le * (Vy * Vy) * Mxb * (c_y * c_y - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * power(Vy, 3.0) * Mxb * (Oy * sin(Le * Vy) - 2.0 * Oy * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0))) - power(Vz, 3.0) * ((4.0 * Mxb * (sin(Le * Vz) - 2.0 * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Mxb * (Oy * sin(Le * Vz) - 2.0 * Oy * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0))) - Le * Vy * Mxb * sin(Le * Vy) * sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))       
        et[4 , 7] = et[7 , 4] = -(Vy * Vy * ((4.0 * Vz * Mxb * (sin(Le * Vz) - 2.0 * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0)) - 4.0 * Le * (Vz * Vz) * Mxb * (c_z * c_z - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * power(Vz, 3.0) * Mxb * (Oz * sin(Le * Vz) - 2.0 * Oz * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0))) - power(Vy, 3.0) * ((4.0 * Mxb * (sin(Le * Vy) - 2.0 * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vz * Vz) * Mxb * (Oz * sin(Le * Vy) - 2.0 * Oz * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0))) - Le * Vz * Mxb * sin(Le * Vy) * sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
        et[5 , 8] = et[8 , 5] = (Vz * Vz * ((4.0 * Vy * Mxb * (sin(Le * Vy) - 2.0 * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0)) - 4.0 * Le * (Vy * Vy) * Mxb * (c_y * c_y - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * power(Vy, 3.0) * Mxb * (Oy * sin(Le * Vy) - 2.0 * Oy * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0))) - power(Vz, 3.0) * ((4.0 * Mxb * (sin(Le * Vz) - 2.0 * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Mxb * (Oy * sin(Le * Vz) - 2.0 * Oy * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0))) - Le * Vy * Mxb * sin(Le * Vy) * sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
        et[7 , 10] = et[10 , 7] = (Vy * Vy * ((4.0 * Vz * Mxb * (sin(Le * Vz) - 2.0 * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0)) - 4.0 * Le * (Vz * Vz) * Mxb * (c_z * c_z - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * power(Vz, 3.0) * Mxb * (Oz * sin(Le * Vz) - 2.0 * Oz * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0))) - power(Vy, 3.0) * ((4.0 * Mxb * (sin(Le * Vy) - 2.0 * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vz * Vz) * Mxb * (Oz * sin(Le * Vy) - 2.0 * Oz * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0))) - Le * Vz * Mxb * sin(Le * Vy) * sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
        et[8 , 11] = et[11 , 8] = -(Vz * Vz * ((4.0 * Vy * Mxb * (sin(Le * Vy) - 2.0 * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0)) - 4.0 * Le * (Vy * Vy) * Mxb * (c_y * c_y - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * power(Vy, 3.0) * Mxb * (Oy * sin (Le * Vy) - 2.0 * Oy * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0))) - power(Vz, 3.0) * ((4.0 * Mxb * (sin(Le * Vz) - 2.0 * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Mxb * (Oy * sin(Le * Vz) - 2.0 * Oy * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0))) - Le * Vy * Mxb * sin(Le * Vy) * sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
        #
        et[5 , 10] = et[10 , 5] = 2.0 * Mxb * (((((((((((((((4.0 * (Vy * Vy) * (sy * sy) * (sz * sz) - 4.0 * (Vz * Vz) * (sy * sy) * (sz * sz)) + Le * Le * (Vy * Vy) * (Vz * Vz) * (sy * sy)) - Le * Le * (Vy * Vy) * (Vz * Vz) * (sz * sz)) + 4.0 * (Le * Le) * power(Vy, 4.0) * Oy * (sy * sy) * (sz * sz)) - 4.0 * (Le * Le) * power(Vz, 4.0) * Oz * (sy * sy) * (sz * sz)) + 2.0 * Le * Vy * (Vz * Vz) * sin(Le * Vy) * (sz * sz)) - 2.0 * Le * (Vy * Vy) * Vz * sin(Le * Vz) * (sy * sy)) - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy * (sy * sy) * (sz * sz)) + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz * (sy * sy) * (sz * sz)) - power(Le, 3.0) * power(Vy, 4.0) * Vz * Oy * sin(Le * Vz) * (sy * sy)) + power(Le, 3.0) * Vy * power(Vz, 4.0) * Oz * sin(Le * Vy) * (sz * sz)) + power(Le, 3.0) * power(Vy, 3.0) * (Vz * Vz) * Oy * sin(Le * Vy) * (sz * sz)) - power(Le, 3.0) * (Vy * Vy) * power(Vz, 3.0) * Oz * sin(Le * Vz) * (sy * sy)) - 4.0 * power(Le, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz * (sy * sy) * (sz * sz)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (sy * sy) * (sz * sz)) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
        et[4 , 11] = et[11 , 4] = -(2.0 * Mxb * (((((((((((((((4.0 * (Vy * Vy) * (sy * sy) * (sz * sz) - 4.0 * (Vz * Vz) * (sy * sy) * (sz * sz)) + Le * Le * (Vy * Vy) * (Vz * Vz) * (sy * sy)) - Le * Le * (Vy * Vy) * (Vz * Vz) * (sz * sz)) + 4.0 * (Le * Le) * power(Vy, 4.0) * Oy * (sy * sy) * (sz * sz)) - 4.0 * (Le * Le) * power(Vz, 4.0) * Oz * (sy * sy) * (sz * sz)) + 2.0 * Le * Vy * (Vz * Vz) * sin(Le * Vy) * (sz * sz)) - 2.0 * Le * (Vy * Vy) * Vz * sin(Le * Vz) * (sy * sy)) - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy * (sy * sy) * (sz * sz)) + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz * (sy * sy) * (sz * sz)) - power(Le, 3.0) * power(Vy, 4.0) * Vz * Oy * sin(Le * Vz) * (sy * sy)) + power(Le, 3.0) * Vy * power(Vz, 4.0) * Oz * sin(Le * Vy) * (sz * sz)) + power(Le, 3.0) * power(Vy, 3.0) * (Vz * Vz) * Oy * sin(Le * Vy) * (sz * sz)) - power(Le, 3.0) * (Vy * Vy) * power(Vz, 3.0) * Oz * sin(Le * Vz) * (sy * sy)) - 4.0 * power(Le, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz * (sy * sy) * (sz * sz)) + 4.0 * power(Le, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (sy * sy) * (sz * sz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
        #
        et[1 , 10] = et[10 , 1] = -(Vy * Vy * ((4.0 * Vz * Mxb * (sin(Le * Vz) - 2.0 * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0)) - 4.0 * Le * (Vz * Vz) * Mxb * (c_z * c_z - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * power(Vz, 3.0) * Mxb * (Oz * sin(Le * Vz) - 2.0 * Oz * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0))) - power(Vy, 3.0) * ((4.0 * Mxb * (sin(Le * Vy) - 2.0 * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vz * Vz) * Mxb * (Oz * sin(Le * Vy) - 2.0 * Oz * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0))) - Le * Vz * Mxb * sin(Le * Vy) * sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
        et[2 , 11] = et[11 , 2] = (Vz * Vz * ((4.0 * Vy * Mxb * (sin(Le * Vy) - 2.0 * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0)) - 4.0 * Le * (Vy * Vy) * Mxb * (c_y * c_y - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * power(Vy, 3.0) * Mxb * (Oy * sin(Le * Vy) - 2.0 * Oy * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0))) - power(Vz, 3.0) * ((4.0 * Mxb * (sin(Le * Vz) - 2.0 * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Mxb * (Oy * sin(Le * Vz) - 2.0 * Oy * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0))) - Le * Vy * Mxb * sin(Le * Vy) * sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
        #
        et[4 , 5] = et[5 , 4] = ((power(Le, 3.0) * (power(Vy, 3.0) * (Vz * Vz) * (Mxb * (4.0 * Oy * sin(Le * Vy) - 4.0 * Oz * sin(Le * Vy)) / 2.0 - Mxb * sin(Le * Vy / 2.0) * (8.0 * Oy * cos(Le * Vy / 2.0) * (c_z * c_z) - 8.0 * Oz * cos(Le * Vy / 2.0) * (c_z * c_z)) / 2.0) - Vy * Vy * power(Vz, 3.0) * (Mxb * (4.0 * Oy * sin(Le * Vz) - 4.0 * Oz * sin(Le * Vz)) / 2.0 - Mxb * sin(Le * Vz / 2.0) * (8.0 * Oy * (c_y * c_y) * cos(Le * Vz / 2.0) - 8.0 * Oz * (c_y * c_y) * cos(Le * Vz / 2.0)) / 2.0)) - Le * (((power(Vy, 3.0) * (2.0 * Mxb * sin(Le * Vy) - 4.0 * Mxb * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0)) + power(Vz, 3.0) * (2.0 * Mxb * sin(Le * Vz) - 4.0 * Mxb * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0))) - Vy * (Vz * Vz) * (2.0 * Mxb * sin(Le * Vy) - 4.0 * Mxb * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0))) - Vy * Vy * Vz * (2.0 * Mxb * sin(Le * Vz) - 4.0 * Mxb * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0)))) + Le * Le * ((Vy * power(Vz, 3.0) * Mxb * sin(Le * Vy) * sin(Le * Vz) / 2.0 - Vy * Vy * (Vz * Vz) * Mxb * ((4.0 * (c_y * c_y) + 4.0 * (c_z * c_z)) - 8.0 * (c_y * c_y) * (c_z * c_z)) / 2.0) + power(Vy, 3.0) * Vz * Mxb * sin(Le * Vy) * sin(Le * Vz) / 2.0)) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
        et[10 , 11] = et[11 , 10] = -((power(Le, 3.0) * (power(Vy, 3.0) * (Vz * Vz) * (Mxb * (4.0 * Oy * sin(Le * Vy) - 4.0 * Oz * sin(Le * Vy)) / 2.0 - Mxb * sin(Le * Vy / 2.0) * (8.0 * Oy * cos(Le * Vy / 2.0) * (c_z * c_z) - 8.0 * Oz * cos(Le * Vy / 2.0) * (c_z * c_z)) / 2.0) - Vy * Vy * power(Vz, 3.0) * (Mxb * (4.0 * Oy * sin(Le * Vz) - 4.0 * Oz * sin(Le * Vz)) / 2.0 - Mxb * sin(Le * Vz / 2.0) * (8.0 * Oy * (c_y * c_y) * cos(Le * Vz / 2.0) - 8.0 * Oz * (c_y * c_y) * cos(Le * Vz / 2.0)) / 2.0)) - Le * (((power(Vy, 3.0) * (2.0 * Mxb * sin(Le * Vy) - 4.0 * Mxb * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0)) + power(Vz, 3.0) * (2.0 * Mxb * sin(Le * Vz) - 4.0 * Mxb * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0))) - Vy * (Vz * Vz) * (2.0 * Mxb * sin(Le * Vy) - 4.0 * Mxb * cos(Le * Vy / 2.0) * (c_z * c_z) * sin(Le * Vy / 2.0))) - Vy * Vy * Vz * (2.0 * Mxb * sin(Le * Vz) - 4.0 * Mxb * (c_y * c_y) * cos(Le * Vz / 2.0) * sin(Le * Vz / 2.0)))) + Le * Le * ((Vy * power(Vz, 3.0) * Mxb * sin(Le * Vy) * sin(Le * Vz) / 2.0 - Vy * Vy * (Vz * Vz) * Mxb * ((4.0 * (c_y * c_y) + 4.0 * (c_z * c_z)) - 8.0 * (c_y * c_y) * (c_z * c_z)) / 2.0) + power(Vy, 3.0) * Vz * Mxb * sin(Le * Vy) * sin(Le * Vz) / 2.0)) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    #
    #
    return et
#
def Kt_tension(Le:float, Ax:float,
               Jx:float, Iy:float, Iz:float,
               Emod:float, Gmod:float,
               Oy:float, Oz:float,
               Fb:list[float]):
    """Negative axial force"""
    (Fxa, Fya, Fza, Mxa, Mya, Mza,
     Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb.q_loc
    #
    P = Fb.Fx
    A = Ax
    E = Emod
    G = Gmod
    L = Le
    L2 = Le*Le
    #L3 = L2*Le
    #pi = np.pi
    Ix = Jx
    #
    power = np.power
    sinh = np.sinh
    cosh = np.cosh
    sqrt = np.sqrt 
    #
    #/* --------------------------------------------- Positive axial force --------------------------------------------- */
    mu_y = sqrt(P / (E*Iz))
    mu_z = sqrt(P / (E*Iy))
    Vy = mu_y / (sqrt(1 + Oy*mu_y*mu_y*L2))
    Vz = mu_z / (sqrt(1 + Oz*mu_z*mu_z*L2))

    K = P*(Iy + Iz) / Ax
    ay = (L * Vy * cosh(L * Vy / 2.0) - 2.0 * sinh(L * Vy / 2.0)) + 2.0 * (L * L) * (Vy * Vy) * Oy * sinh(L * Vy / 2.0)
    az = (L * Vz * cosh(L * Vz / 2.0) - 2.0 * sinh(L * Vz / 2.0)) + 2.0 * (L * L) * (Vz * Vz) * Oz * sinh(L * Vz / 2.0)
    by = ((L * Vy - 2.0 * sinh(L * Vy / 2.0)) + 2.0 * L * Vy * (sinh(L * Vy / 4.0) * sinh(L * Vy / 4.0))) + 2.0 * (L * L) * (Vy * Vy) * Oy * sinh(L * Vy / 2.0)
    bz = ((L * Vz - 2.0 * sinh(L * Vz / 2.0)) + 2.0 * L * Vz * (sinh(L * Vz / 4.0) * sinh(L * Vz / 4.0))) + 2.0 * (L * L) * (Vz * Vz) * Oz * sinh(L * Vz / 2.0)
    cy = 2.0 * Iz * power(L, 5.0) * power(Vy, 7.0) * (Oy * Oy) * P
    cz = 2.0 * Iy * power(L, 5.0) * power(Vz, 7.0) * (Oz * Oz) * P
    dy = 4.0 * A * (L * L) * (Vy * Vy) * Oy * P
    dz = 4.0 * A * (L * L) * (Vz * Vz) * Oz * P
    ey = A * power(L, 4.0) * power(Vy, 4.0) * Oy * P
    ez = A * power(L, 4.0) * power(Vz, 4.0) * Oz * P
    fy = 2.0 * Iz * (L * L) * power(Vy, 4.0) * Oy * P
    fz = 2.0 * Iy * (L * L) * power(Vz, 4.0) * Oz * P
    sy = sinh(L * Vy / 2.0)
    sz = sinh(L * Vz / 2.0)
    c_y = cosh(L * Vy / 2.0)
    c_z = cosh(L * Vz / 2.0)

    Ket = np.zeros((12,12), dtype=np.float32)

    Ket[0 , 0] = Ket[6 , 6] = P / L + Ax * E / L
    Ket[1 , 1] = Ket[7 , 7] = Vy * ((((((((((((Iz * (Vy * Vy) * P * sinh(L * Vy) - 3.0 * A * P * sinh(L * Vy)) + 2.0 * A * L * Vy * P) - Iz * L * power(Vy, 3.0) * P) - A * E * Iz * L * power(Vy, 3.0)) + A * E * Iz * (Vy * Vy) * sinh(L * Vy)) - 2.0 * A * power(L, 3.0) * power(Vy, 3.0) * Oy * P) + A * L * Vy * P * cosh(L * Vy)) + A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P) + 2.0 * A * (L * L) * (Vy * Vy) * Oy * P * sinh(L * Vy)) + A * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * sinh(L * Vy)) + A * E * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy) + A * E * Iz * (L * L) * power(Vy, 4.0) * Oy * sinh(L * Vy)) / (2.0 * A * (ay * ay))
    Ket[2 , 2] = Ket[8 , 8] = Vz * ((((((((((((Iy * (Vz * Vz) * P * sinh(L * Vz) - 3.0 * A * P * sinh(L * Vz)) + 2.0 * A * L * Vz * P) - Iy * L * power(Vz, 3.0) * P) - A * E * Iy * L * power(Vz, 3.0)) + A * E * Iy * (Vz * Vz) * sinh(L * Vz)) - 2.0 * A * power(L, 3.0) * power(Vz, 3.0) * Oz * P) + A * L * Vz * P * cosh(L * Vz)) + A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P) + 2.0 * A * (L * L) * (Vz * Vz) * Oz * P * sinh(L * Vz)) + A * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * sinh(L * Vz)) + A * E * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz) + A * E * Iy * (L * L) * power(Vz, 4.0) * Oz * sinh(L * Vz)) / (2.0 * A * (az * az))
    Ket[3 , 3] = Ket[9 , 9] = K / L + G*Ix / L
    x = ((((((((((((((((((((((((A * P * sinh(2.0 * L * Vz) - 2.0 * A * P * sinh(L *	Vz)) - 2.0 * Iy * (Vz * Vz) * P * sinh(L * Vz)) + Iy * (Vz * Vz) * P * sinh(2.0 * L * Vz)) - A * power(L, 3.0) * power(Vz, 3.0) * P) + Iy * power(L, 3.0) * power(Vz, 5.0) * P) - Iy * L * power(Vz, 3.0) * P) - A * E * Iy * L * power(Vz, 3.0)) - 2.0 * A * E * Iy * (Vz * Vz) * sinh(L * Vz)) + A * E * Iy * (Vz * Vz) * sinh(2.0 * L * Vz)) + 2.0 * Iy * L * power(Vz, 3.0) * P * cosh(L * Vz)) - Iy * L * power(Vz, 3.0) * P * cosh(2.0 * L * Vz)) + A * E * Iy * power(L, 3.0) * power(Vz, 5.0)) - 5.0 * A * power(L, 3.0) * power(Vz, 3.0) * Oz * P) + 2.0 * A * power(L, 5.0) * power(Vz, 5.0) * Oz * P) + 3.0 * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * P) + 2.0 * A * (L * L) * (Vz * Vz) * P * sinh(L * Vz)) + A * (L * L) * (Vz * Vz) * P * sinh(2.0 * L * Vz) / 2.0) - 2.0 * Iy * (L * L) * power(Vz, 4.0) * P * sinh(L * Vz)) + Iy * (L * L) * power(Vz, 4.0) * P * sinh(2.0 * L * Vz) / 2.0) + 2.0 * A * L * Vz * P * cosh(L * Vz)) - 2.0 * A * L * Vz * P * cosh(2.0 * L * Vz)) + 12.0 * A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P) - A * power(L, 7.0) * power(Vz, 7.0) * (Oz * Oz) * P) - 9.0 * A * power(L, 7.0) * power(Vz, 7.0) * power(Oz, 3.0) * P) + 2.0 * A * power(L, 9.0) * power(Vz, 9.0) * power(Oz, 4.0) * P
    Ket[4 , 4] = Ket[10 , 10] = (((((((((((((((((((((((((((((((((((((((((((((((((x - cz) + 2.0 * dz * sinh(L * Vz)) - dz * sinh(2.0 * L * Vz)) - 6.0 * ez * sinh(L * Vz)) - ez * sinh(2.0 * L * Vz)) + 2.0 * fz * sinh(L * Vz)) - fz * sinh(2.0 * L * Vz)) + 2.0 * Iy * power(L, 4.0) * power(Vz, 6.0) * Oz * P * sinh(L * Vz)) - 7.0 * A * E * Iy * power(L, 5.0) * power(Vz, 7.0) * (Oz * Oz)) + 2.0 * A * E * Iy * power(L, 7.0) * power(Vz, 9.0) * power(Oz, 3.0)) - 8.0 * A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P * cosh(L * Vz)) - 4.0 * A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P * cosh(2.0 * L * Vz)) + 8.0 * A * power(L, 7.0) * power(Vz, 7.0) * power(Oz, 3.0) * P * cosh(L * Vz)) + A * power(L, 7.0) * power(Vz, 7.0) * power(Oz, 3.0) * P * cosh(2.0 * L * Vz)) - 2.0 * A * power(L, 9.0) * power(Vz, 9.0) * power(Oz, 4.0) * P * cosh(L * Vz)) + 2.0 * Iy * power(L, 5.0) * power(Vz, 7.0) * (Oz * Oz) * P * cosh(L * Vz)) - 12.0 * A * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * sinh(L * Vz)) + 6.0 * A * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * sinh(2.0 * L * Vz)) + 6.0 * A * power(L, 6.0) * power(Vz, 6.0) * (Oz * Oz) * P * sinh(L * Vz)) + 8.0 * A * power(L, 6.0) * power(Vz, 6.0) * power(Oz, 3.0) * P * sinh(L * Vz)) + A * power(L, 6.0) * power(Vz, 6.0) * (Oz * Oz) * P * sinh(2.0 * L * Vz) / 2.0) - 4.0 * A * power(L, 6.0) * power(Vz, 6.0) * power(Oz, 3.0) * P * sinh(2.0 * L * Vz)) - 2.0 * A * power(L, 8.0) * power(Vz, 8.0) * power(Oz, 3.0) * P * sinh(L * Vz)) - 2.0 * A * power(L, 8.0) * power(Vz, 8.0) * power(Oz, 4.0) * P * sinh(L * Vz)) + A * power(L, 8.0) * power(Vz, 8.0) * power(Oz, 4.0) * P * sinh(2.0 * L * Vz)) - 2.0 * Iy * power(L, 4.0) * power(Vz, 6.0) * (Oz * Oz) * P * sinh(L * Vz)) + Iy * power(L, 4.0) * power(Vz, 6.0) * (Oz * Oz) * P * sinh(2.0 * L * Vz)) + 2.0 * A * E * Iy * L * power(Vz, 3.0) * cosh(L * Vz)) - A * E * Iy * L * power(Vz, 3.0) * cosh(2.0 * L * Vz)) + 6.0 * A * E * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz) - A * E * Iy * power(L, 5.0) * power(Vz, 7.0) * Oz) - 2.0 * A * E * Iy * (L * L) * power(Vz, 4.0) * sinh(L * Vz)) + A * E * Iy * (L * L) * power(Vz, 4.0) * sinh(2.0 * L * Vz) / 2.0) + 5.0 * A * power(L, 3.0) * power(Vz, 3.0) * Oz * P * cosh(2.0 * L * Vz)) - 4.0 * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * P * cosh(L * Vz)) + Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * P * cosh(2.0 * L * Vz)) - 6.0 * A * E * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * cosh(L * Vz)) + 2.0 * A * E * Iy * (L * L) * power(Vz, 4.0) * Oz * sinh(L * Vz)) - A * E * Iy * (L * L) * power(Vz, 4.0) * Oz * sinh(2.0 * L * Vz)) + 4.0 * A * E * Iy * power(L, 4.0) * power(Vz, 6.0) * Oz * sinh(L * Vz)) + A * E * Iy * power(L, 4.0) * power(Vz, 6.0) * Oz * sinh(2.0 * L * Vz) / 2.0) + 6.0 * A * E * Iy * power(L, 5.0) * power(Vz, 7.0) * (Oz * Oz) * cosh(L * Vz)) + A * E * Iy * power(L, 5.0) * power(Vz, 7.0) * (Oz * Oz) * cosh(2.0 * L * Vz)) - 2.0 * A * E * Iy * power(L, 7.0) * power(Vz, 9.0) * power(Oz, 3.0) * cosh(L * Vz)) + 2.0 * A * E * Iy * power(L, 4.0) * power(Vz, 6.0) * (Oz * Oz) * sinh(L * Vz)) - A * E * Iy * power(L, 4.0) * power(Vz, 6.0) * (Oz * Oz) * sinh(2.0 * L * Vz)) - 2.0 * A * E * Iy * power(L, 6.0) * power(Vz, 8.0) * (Oz * Oz) * sinh(L * Vz)) - 2.0 * A * E * Iy * power(L, 6.0) * power(Vz, 8.0) * power(Oz, 3.0) * sinh(L * Vz)) + A * E * Iy * power(L, 6.0) * power(Vz, 8.0) * power(Oz, 3.0) * sinh(2.0 * L * Vz)) / (2.0 * A * Vz * (cosh(L * Vz) - 1.0) * (((((((((4.0 * cosh(L * Vz) + L * L * (Vz * Vz)) + 8.0 * (L * L) * (Vz * Vz) * Oz) + L * L * (Vz * Vz) * cosh(L * Vz)) - 4.0 * L * Vz * sinh(L * Vz)) - 4.0 * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz)) - 8.0 * (L * L) * (Vz * Vz) * Oz * cosh(L * Vz)) + 4.0 * power(L, 3.0)	* power(Vz, 3.0) * Oz * sinh(L * Vz)) + 4.0 * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cosh(L * Vz)) - 4.0))
    x = ((((((((((((((((((((((((A * P * sinh(2.0 * L * Vy) - 2.0 * A * P * sinh(L * Vy)) - 2.0 * Iz * (Vy * Vy) * P * sinh(L * Vy)) + Iz * (Vy * Vy) * P * sinh(2.0 * L * Vy)) - A * power(L, 3.0) * power(Vy, 3.0) * P) + Iz * power(L, 3.0) * power(Vy, 5.0) * P) - Iz * L * power(Vy, 3.0) * P) - A * E * Iz * L * power(Vy, 3.0)) - 2.0 * A * E * Iz * (Vy * Vy) * sinh(L * Vy)) + A * E * Iz * (Vy * Vy) * sinh(2.0 * L * Vy)) + 2.0 * Iz * L * power(Vy, 3.0) * P * cosh(L * Vy)) - Iz * L * power(Vy, 3.0) * P * cosh(2.0 * L * Vy)) + A * E * Iz * power(L, 3.0) * power(Vy, 5.0)) - 5.0 * A * power(L, 3.0) * power(Vy, 3.0) * Oy * P) + 2.0 * A * power(L, 5.0) * power(Vy, 5.0) * Oy * P) + 3.0 * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * P) + 2.0 * A * (L * L) * (Vy * Vy) * P * sinh(L * Vy)) + A * (L * L) * (Vy * Vy) * P * sinh(2.0 * L * Vy) / 2.0) - 2.0 * Iz * (L * L) * power(Vy, 4.0) * P * sinh(L * Vy)) + Iz * (L * L) * power(Vy, 4.0) * P * sinh(2.0 * L * Vy) / 2.0) + 2.0 * A * L * Vy * P * cosh(L * Vy)) - 2.0 * A * L * Vy * P * cosh(2.0 * L * Vy)) + 12.0 * A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P) - A * power(L, 7.0) * power(Vy, 7.0) * (Oy * Oy) * P) - 9.0 * A * power(L, 7.0) * power(Vy, 7.0) * power(Oy, 3.0) * P) + 2.0 * A * power(L, 9.0) * power(Vy, 9.0) * power(Oy, 4.0) * P
    Ket[5 , 5] = Ket[11 , 11] = (((((((((((((((((((((((((((((((((((((((((((((((((x - cy) + 2.0 * dy * sinh(L * Vy)) - dy * sinh(2.0 * L * Vy)) - 6.0 * ey * sinh(L * Vy)) - ey * sinh(2.0 * L * Vy)) + 2.0 * fy * sinh(L * Vy)) - fy * sinh(2.0 * L * Vy)) + 2.0 * Iz * power(L, 4.0) * power(Vy, 6.0) * Oy * P * sinh(L * Vy)) - 7.0 * A * E * Iz * power(L, 5.0) * power(Vy, 7.0) * (Oy * Oy)) + 2.0 * A * E * Iz * power(L, 7.0) * power(Vy, 9.0) * power(Oy, 3.0)) - 8.0 * A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P * cosh(L * Vy)) - 4.0 * A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P * cosh(2.0 * L * Vy)) + 8.0 * A * power(L, 7.0) * power(Vy, 7.0) * power(Oy, 3.0) * P * cosh(L * Vy)) + A * power(L, 7.0) * power(Vy, 7.0) * power(Oy, 3.0) * P * cosh(2.0 * L * Vy)) - 2.0 * A * power(L, 9.0) * power(Vy, 9.0) * power(Oy, 4.0) * P * cosh(L * Vy)) + 2.0 * Iz * power(L, 5.0) * power(Vy, 7.0) * (Oy * Oy) * P * cosh(L * Vy)) - 12.0 * A * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * sinh(L * Vy)) + 6.0 * A * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * sinh(2.0 * L * Vy)) + 6.0 * A * power(L, 6.0) * power(Vy, 6.0) * (Oy * Oy) * P * sinh(L * Vy)) + 8.0 * A * power(L, 6.0) * power(Vy, 6.0) * power(Oy, 3.0) * P * sinh(L * Vy)) + A * power(L, 6.0) * power(Vy, 6.0) * (Oy * Oy) * P * sinh(2.0 * L * Vy) / 2.0) - 4.0 * A * power(L, 6.0) * power(Vy, 6.0) * power(Oy, 3.0) * P * sinh(2.0 * L * Vy)) - 2.0 * A * power(L, 8.0) * power(Vy, 8.0) * power(Oy, 3.0) * P * sinh(L * Vy)) - 2.0 * A * power(L, 8.0) * power(Vy, 8.0) * power(Oy, 4.0) * P * sinh(L * Vy)) + A * power(L, 8.0) * power(Vy, 8.0) * power(Oy, 4.0) * P * sinh(2.0 * L * Vy)) - 2.0 * Iz * power(L, 4.0) * power(Vy, 6.0) * (Oy * Oy) * P * sinh(L * Vy)) + Iz * power(L, 4.0) * power(Vy, 6.0) * (Oy * Oy) * P * sinh(2.0 * L * Vy)) + 2.0 * A * E * Iz * L * power(Vy, 3.0) * cosh(L * Vy)) - A * E * Iz * L * power(Vy, 3.0) * cosh(2.0 * L * Vy)) + 6.0 * A * E * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy) - A * E * Iz * power(L, 5.0) * power(Vy, 7.0) * Oy) - 2.0 * A * E * Iz * (L * L) * power(Vy, 4.0) * sinh(L * Vy)) + A * E * Iz * (L * L) * power(Vy, 4.0) * sinh(2.0 * L * Vy) / 2.0) + 5.0 * A * power(L, 3.0) * power(Vy, 3.0) * Oy * P * cosh(2.0 * L * Vy)) - 4.0 * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * P * cosh(L * Vy)) + Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * P * cosh(2.0 * L * Vy)) - 6.0 * A * E * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * cosh(L * Vy)) + 2.0 * A * E * Iz * (L * L) * power(Vy, 4.0) * Oy * sinh(L * Vy)) - A * E * Iz * (L * L) * power(Vy, 4.0) * Oy * sinh(2.0 * L * Vy)) + 4.0 * A * E * Iz * power(L, 4.0) * power(Vy, 6.0) * Oy * sinh(L * Vy)) + A * E * Iz * power(L, 4.0) * power(Vy, 6.0) * Oy * sinh(2.0 * L * Vy) / 2.0) + 6.0 * A * E * Iz * power(L, 5.0) * power(Vy, 7.0) * (Oy * Oy) * cosh(L * Vy)) + A * E * Iz * power(L, 5.0) * power(Vy, 7.0) * (Oy * Oy) * cosh(2.0 * L * Vy)) - 2.0 * A * E * Iz * power(L, 7.0) * power(Vy, 9.0) * power(Oy, 3.0) * cosh(L * Vy)) + 2.0 * A * E * Iz * power(L, 4.0) * power(Vy, 6.0) * (Oy * Oy) * sinh(L * Vy)) - A * E * Iz * power(L, 4.0) * power(Vy, 6.0) * (Oy * Oy) * sinh(2.0 * L * Vy)) - 2.0 * A * E * Iz * power(L, 6.0) * power(Vy, 8.0) * (Oy * Oy) * sinh(L * Vy)) - 2.0 * A * E * Iz * power(L, 6.0) * power(Vy, 8.0) * power(Oy, 3.0) * sinh(L * Vy)) + A * E * Iz * power(L, 6.0) * power(Vy, 8.0) * power(Oy, 3.0) * sinh(2.0 * L * Vy)) / (2.0 * A * Vy * (cosh(L * Vy) - 1.0) * (((((((((4.0 * cosh(L * Vy) + L * L * (Vy * Vy)) + 8.0 * (L * L) * (Vy * Vy) * Oy) + L * L * (Vy * Vy) * cosh(L * Vy)) - 4.0 * L * Vy * sinh(L * Vy)) - 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy)) - 8.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) + 4.0 * power(L, 3.0) * power(Vy, 3.0) * Oy * sinh(L * Vy)) + 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cosh(L * Vy)) - 4.0))

    Ket[2 , 3] = Ket[3 , 2] = Mza / L
    Ket[3 , 4] = Ket[4 , 3] = -Mza / 2.0 - (Oz * (L * L) * (Vz * Vz) - 1.0) * (Mza + Mzb) * ((((L * Vz - 2.0 * sinh(L * Vz)) - 2.0 * cosh(L * Vz)) + L * Vz * (cosh(L * Vz) + sinh(L * Vz))) + 2.0) / (L * L * (Vz * Vz) * ((cosh(L * Vz) + sinh(L * Vz)) - 1.0))
    Ket[5 , 6] = Ket[6 , 5] = Mza / L
    Ket[8 , 9] = Ket[9 , 8] = -Mzb / L
    Ket[9 , 10] = Ket[10 , 9] = -Mzb / 2.0 - (Oz * (L * L) * (Vz * Vz) - 1.0) * (Mza + Mzb) * ((((L * Vz - 2.0 * sinh(L * Vz)) - 2.0 * cosh(L * Vz)) + L * Vz * (cosh(L * Vz) + sinh(L * Vz))) + 2.0) / (L * L * (Vz * Vz) * ((cosh(L * Vz) + sinh(L * Vz)) - 1.0))

    Ket[1 , 3] = Ket[3 , 1] = Mya / L
    Ket[2 , 4] = Ket[4 , 2] = -((((((((((((((A * (L * L) * (Vz * Vz) * P - 8.0 * A * P * (sz * sz)) - Iy * (L * L) * power(Vz, 4.0) * P) + Iy * L * power(Vz, 3.0) * P * sinh(L * Vz)) - A * E * Iy * (L * L) * power(Vz, 4.0)) - 2.0 * A * power(L, 4.0) * power(Vz, 4.0) * Oz * P) + A * L * Vz * P * sinh(L * Vz)) + A * power(L, 6.0) * power(Vz, 6.0) * (Oz * Oz) * P) - 2.0 * A * power(L, 3.0) * power(Vz, 3.0) * Oz * P * sinh(L * Vz)) + 16.0 * A * (L * L) * (Vz * Vz) * Oz * P * (sz * sz)) + A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P * sinh(L * Vz)) + A * E * Iy * L * power(Vz, 3.0) * sinh(L * Vz)) - 8.0 * A * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * (sz * sz)) + A * E * Iy * power(L, 4.0) * power(Vz, 6.0) * Oz) + A * E * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * sinh(L * Vz)) / (4.0 * A * (bz * bz))
    Ket[3 , 5] = Ket[5 , 3] = Mya / 2.0 + (Oy * (L * L) * (Vy * Vy) - 1.0) * (Mya + Myb) * ((((L * Vy - 2.0 * sinh(L * Vy)) - 2.0 * cosh(L * Vy)) + L * Vy * (cosh(L * Vy) + sinh(L * Vy))) + 2.0) / (L * L * (Vy * Vy) * ((cosh(L * Vy) + sinh(L * Vy)) - 1.0))
    Ket[4 , 6] = Ket[6 , 4] = Mya / L
    Ket[5 , 7] = Ket[7 , 5] = -((((((((((((((A * (L * L) * (Vy * Vy) * P - 8.0 * A * P * (sy * sy)) - Iz * (L * L) * power(Vy, 4.0) * P) + Iz * L * power(Vy, 3.0) * P * sinh(L * Vy)) - A * E * Iz * (L * L) * power(Vy, 4.0)) - 2.0 * A * power(L, 4.0) * power(Vy, 4.0) * Oy * P) + A * L * Vy * P * sinh(L * Vy)) + A * power(L, 6.0) * power(Vy, 6.0) * (Oy * Oy) * P) - 2.0 * A * power(L, 3.0) * power(Vy, 3.0) * Oy * P * sinh(L * Vy)) + 16.0 * A * (L * L) * (Vy * Vy) * Oy * P * (sy * sy)) + A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P * sinh(L * Vy)) + A * E * Iz * L * power(Vy, 3.0) * sinh(L * Vy)) - 8.0 * A * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * (sy * sy)) + A * E * Iz * power(L, 4.0) * power(Vy, 6.0) * Oy) + A * E * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * sinh(L * Vy)) / (4.0 * A * (by * by))
    Ket[7 , 9] = Ket[9 , 7] = -Myb / L
    Ket[8 , 10] = Ket[10 , 8] = ((((((((((((((A * (L * L) * (Vz * Vz) * P - 8.0 * A * P * (sz * sz)) - Iy * (L * L) * power(Vz, 4.0) * P) + Iy * L * power(Vz, 3.0) * P * sinh(L * Vz)) - A * E * Iy * (L * L) * power(Vz, 4.0)) - 2.0 * A * power(L, 4.0) * power(Vz, 4.0) * Oz * P) + A * L * Vz * P * sinh(L * Vz)) + A * power(L, 6.0) * power(Vz, 6.0) * (Oz * Oz) * P) - 2.0 * A * power(L, 3.0) * power(Vz, 3.0) * Oz * P * sinh(L * Vz)) + 16.0 * A * (L * L) * (Vz * Vz) * Oz * P * (sz * sz)) + A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P * sinh(L * Vz)) + A * E * Iy * L * power(Vz, 3.0) * sinh(L * Vz)) - 8.0 * A * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * (sz * sz)) + A * E * Iy * power(L, 4.0) * power(Vz, 6.0) * Oz) + A * E * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * sinh(L * Vz)) / (4.0 * A * (bz * bz))
    Ket[9 , 11] = Ket[11 , 9] = Myb / 2.0 + (Oy * (L * L) * (Vy * Vy) - 1.0) * (Mya + Myb) * ((((L * Vy - 2.0 * sinh(L * Vy)) - 2.0 * cosh(L * Vy)) + L * Vy * (cosh(L * Vy) + sinh(L * Vy))) + 2.0) / (L * L * (Vy * Vy) * ((cosh(L * Vy) + sinh(L * Vy)) - 1.0))

    #Ket[1 , 4] = Ket[4 , 1] = Vy * Vy * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0)) / (sinh(L * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vy * cosh(L * Vy / 2.0) - 2.0 * sinh(L * Vy / 2.0)) + 2.0 * (L * L) * (Vy * Vy) * Oy * sinh(L * Vy / 2.0)))
    #Ket[2 , 5] = Ket[5 , 2] = Vz * Vz * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0)) / (sinh(L * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vz * cosh(L * Vz / 2.0) - 2.0 * sinh(L * Vz / 2.0)) + 2.0 * (L * L) * (Vz * Vz) * Oz * sinh(L * Vz / 2.0)))
    #Ket[4 , 7] = Ket[7 , 4] = -(Vy * Vy * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0))) / (sinh(L * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vy * cosh(L * Vy / 2.0) - 2.0 * sinh(L * Vy / 2.0)) + 2.0 * (L * L) * (Vy * Vy) * Oy * sinh(L * Vy / 2.0)))
    #Ket[5 , 8] = Ket[8 , 5] = -(Vz * Vz * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0))) / (sinh(L * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vz * cosh(L * Vz / 2.0) - 2.0 * sinh(L * Vz / 2.0)) + 2.0 * (L * L) * (Vz * Vz) * Oz * sinh(L * Vz / 2.0)))
    #Ket[7 , 10] = Ket[10 , 7] = Vy * Vy * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0)) / (sinh(L * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vy * cosh(L * Vy / 2.0) - 2.0 * sinh(L * Vy / 2.0)) + 2.0 * (L * L) * (Vy * Vy) * Oy * sinh(L * Vy / 2.0)))
    #Ket[8 , 11] = Ket[11 , 8] = Vz * Vz * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0)) / (sinh(L * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vz * cosh(L * Vz / 2.0) - 2.0 * sinh(L * Vz / 2.0)) + 2.0 * (L * L) * (Vz * Vz) * Oz * sinh(L * Vz / 2.0)))

    Ket[0 , 4] = Ket[4 , 0] = -Mya / L
    Ket[1 , 5] = Ket[5 , 1] = ((((((((((((((A * (L * L) * (Vy * Vy) * P - 8.0 * A * P * (sy * sy)) - Iz * (L * L) * power(Vy, 4.0) * P) + Iz * L * power(Vy, 3.0) * P * sinh(L * Vy)) - A * E * Iz * (L * L) * power(Vy, 4.0)) - 2.0 * A * power(L, 4.0) * power(Vy, 4.0) * Oy * P) + A * L * Vy * P * sinh(L * Vy)) + A * power(L, 6.0) * power(Vy, 6.0) * (Oy * Oy) * P) - 2.0 * A * power(L, 3.0) * power(Vy, 3.0) * Oy * P * sinh(L * Vy)) + 16.0 * A * (L * L) * (Vy * Vy) * Oy * P * (sy * sy)) + A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P * sinh(L * Vy)) + A * E * Iz * L * power(Vy, 3.0) * sinh(L * Vy)) - 8.0 * A * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * (sy * sy)) + A * E * Iz * power(L, 4.0) * power(Vy, 6.0) * Oy) + A * E * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * sinh(L * Vy)) / (4.0 * A * (by * by))
    Ket[3 , 7] = Ket[7 , 3] = -Mya / L
    Ket[4 , 8] = Ket[8 , 4] = ((((((((((((((A * (L * L) * (Vz * Vz) * P - 8.0 * A * P * (sz * sz)) - Iy * (L * L) * power(Vz, 4.0) * P) + Iy * L * power(Vz, 3.0) * P * sinh(L * Vz)) - A * E * Iy * (L * L) * power(Vz, 4.0)) - 2.0 * A * power(L, 4.0) * power(Vz, 4.0) * Oz * P) + A * L * Vz * P * sinh(L * Vz)) + A * power(L, 6.0) * power(Vz, 6.0) * (Oz * Oz) * P) - 2.0 * A * power(L, 3.0) * power(Vz, 3.0) * Oz * P * sinh(L * Vz)) + 16.0 * A * (L * L) * (Vz * Vz) * Oz * P * (sz * sz)) + A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P * sinh(L * Vz)) + A * E * Iy * L * power(Vz, 3.0) * sinh(L * Vz)) - 8.0 * A * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * (sz * sz)) + A * E * Iy * power(L, 4.0) * power(Vz, 6.0) * Oz) + A * E * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * sinh(L * Vz)) / (4.0 * A * (bz * bz))
    Ket[5 , 9] = Ket[9 , 5] = -((Oy * (L * L) * (Vy * Vy) - 1.0) * (L * Vy * cosh(L * Vy / 2.0) / sinh(L * Vy / 2.0) - 2.0) * (Mya + Myb)) / (L * L * (Vy * Vy))
    Ket[6 , 10] = Ket[10 , 6] = Myb / L
    Ket[7 , 11] = Ket[11 , 7] = -((((((((((((((A * (L * L) * (Vy * Vy) * P - 8.0 * A * P * (sy * sy)) - Iz * (L * L) * power(Vy, 4.0) * P) + Iz * L * power(Vy, 3.0) * P * sinh(L * Vy)) - A * E * Iz * (L * L) * power(Vy, 4.0)) - 2.0 * A * power(L, 4.0) * power(Vy, 4.0) * Oy * P) + A * L * Vy * P * sinh(L * Vy)) + A * power(L, 6.0) * power(Vy, 6.0) * (Oy * Oy) * P) - 2.0 * A * power(L, 3.0) * power(Vy, 3.0) * Oy * P * sinh(L * Vy)) + 16.0 * A * (L * L) * (Vy * Vy) * Oy * P * (sy * sy)) + A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P * sinh(L * Vy)) + A * E * Iz * L * power(Vy, 3.0) * sinh(L * Vy)) - 8.0 * A * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * (sy * sy)) + A * E * Iz * power(L, 4.0) * power(Vy, 6.0) * Oy) + A * E * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * sinh(L * Vy)) / (4.0 * A * (by * by))

    Ket[0 , 5] = Ket[5 , 0] = -Mza / L
    Ket[3 , 8] = Ket[8 , 3] = -Mza / L
    Ket[4 , 9] = Ket[9 , 4] = (Oz * (L * L) * (Vz * Vz) - 1.0) * (L * Vz * cosh(L * Vz / 2.0) / sinh(L * Vz / 2.0) - 2.0) * (Mza + Mzb) / (L * L * (Vz * Vz))
    #Ket[5 , 10] = Ket[10 , 5] = 2.0 * Mxb * (((((((((((((((((((((((((((((((((((((((((((((4.0 * Vy * Vy - 4.0 * Vz * Vz) - 4.0 * (Vy * Vy) * c_y * c_y) - 4.0 * Vy * Vy * (c_z * c_z)) + 4.0 * Vz * Vz * (c_y * c_y)) + 4.0 * (Vz * Vz) * (c_z * c_z)) - 4.0 * L * L * power(Vy, 4.0) * Oy) + 4.0 * (L * L) * power(Vz, 4.0) * Oz) + 4.0 * (Vy * Vy) * (c_y * c_y) * (c_z * c_z)) - 4.0 * (Vz * Vz) * (c_y * c_y) * (c_z * c_z)) - 2.0 * L * Vy * (Vz * Vz) * sinh(L * Vy)) + 2.0 * L * (Vy * Vy) * Vz * sinh(L * Vz)) + L * L * (Vy * Vy) * (Vz * Vz) * (c_y * c_y)) - L * L * (Vy * Vy) * (Vz * Vz) * (c_z * c_z)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz) + 4.0 * (L * L) * power(Vy, 4.0) * Oy * (c_y * c_y)) + 4.0 * (L * L) * power(Vy, 4.0) * Oy * (c_z * c_z)) - 4.0 * (L * L) * power(Vz, 4.0) * Oz * (c_y * c_y)) - 4.0 * (L * L) * power(Vz, 4.0) * Oz * (c_z * c_z)) - 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz) + 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz) + power(L, 3.0) * power(Vy, 3.0) * (Vz * Vz) * Oy * sinh(L * Vy)) - power(L, 3.0) * (Vy * Vy) * power(Vz, 3.0) * Oz * sinh(L * Vz)) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy * (c_y * c_y)) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy * (c_z * c_z)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz * (c_y * c_y)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz * (c_z * c_z)) - 4.0 * (L * L) * power(Vy, 4.0) * Oy * (c_y * c_y) * (c_z * c_z)) + 4.0 * (L * L) * power(Vz, 4.0) * Oz * (c_y * c_y) * (c_z * c_z)) + power(L, 3.0) * Vy * power(Vz, 4.0) * Oz * sinh(L * Vy)) - power(L, 3.0) * power(Vy, 4.0) * Vz * Oy * sinh(L * Vz)) + 4.0 * L * Vy * (Vz * Vz) * cosh(L * Vy / 2.0) * (c_z * c_z) * sinh(L * Vy / 2.0)) - 4.0 * L * (Vy * Vy) * Vz * (c_y * c_y) * cosh(L * Vz / 2.0) * sinh(L * Vz / 2.0)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy * (c_y * c_y) * (c_z * c_z)) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz * (c_y * c_y) * (c_z * c_z)) + 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz * (c_y * c_y)) - 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_y * c_y)) + 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz * (c_z * c_z)) - 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_z * c_z)) - 2.0 * power(L, 3.0) * power(Vy, 3.0) * (Vz * Vz) * Oy * cosh(L * Vy / 2.0) * (c_z * c_z) * sinh(L * Vy / 2.0)) + 2.0 * power(L, 3.0) * (Vy * Vy) * power(Vz, 3.0) * Oz * (c_y * c_y) * cosh(L * Vz / 2.0) * sinh(L * Vz / 2.0)) - 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz * (c_y * c_y) * (c_z * c_z)) + 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_y * c_y) * (c_z * c_z)) - 2.0 * power(L, 3.0) * Vy * power(Vz, 4.0) * Oz * cosh(L * Vy / 2.0) * (c_z * c_z) * sinh(L * Vy / 2.0)) + 2.0 * power(L, 3.0) * power(Vy, 4.0) * Vz * Oy * (c_y * c_y) * cosh(L * Vz / 2.0) * sinh(L * Vz / 2.0)) / ((Vy * Vy - Vz * Vz) * ((((L * Vy * sinh(L * Vy) - 2.0 * (L * L) * (Vy * Vy) * Oy) - 2.0 * cosh(L * Vy)) + 2.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) + 2.0) * ((((L * Vz * sinh(L * Vz) - 2.0 * (L * L) * (Vz * Vz) * Oz) - 2.0 * cosh(L * Vz)) + 2.0 * (L * L) * (Vz * Vz) * Oz * cosh(L * Vz)) + 2.0))
    Ket[6 , 11] = Ket[11 , 6] = -(2.0 * Mxb * (((((((((((((((((((((((((((((((((((((((((((((4.0 * (Vy * Vy) - 4.0 * (Vz * Vz)) - 4.0 * (Vy * Vy) * c_y * c_y) - 4.0 * (Vy * Vy) * c_z * c_z) + 4.0 * Vz * Vz * (c_y * c_y)) + 4.0 * Vz * Vz * (c_z * c_z)) - 4.0 * L * L * power(Vy, 4.0) * Oy) + 4.0 * L * L * power(Vz, 4.0) * Oz) + 4.0 * (Vy * Vy) * (c_y * c_y) * (c_z * c_z)) - 4.0 * (Vz * Vz) * (c_y * c_y) * (c_z * c_z)) - 2.0 * L * Vy * (Vz * Vz) * sinh(L * Vy)) + 2.0 * L * (Vy * Vy) * Vz * sinh(L * Vz)) + L * L * (Vy * Vy) * (Vz * Vz) * (c_y * c_y)) - L * L * (Vy * Vy) * (Vz * Vz) * (c_z * c_z)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz) + 4.0 * (L * L) * power(Vy, 4.0) * Oy * (c_y * c_y)) + 4.0 * (L * L) * power(Vy, 4.0) * Oy * (c_z * c_z)) - 4.0 * (L * L) * power(Vz, 4.0) * Oz * (c_y * c_y)) - 4.0 * (L * L) * power(Vz, 4.0) * Oz * (c_z * c_z)) - 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz) + 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz) + power(L, 3.0) * power(Vy, 3.0) * (Vz * Vz) * Oy * sinh(L * Vy)) - power(L, 3.0) * (Vy * Vy) * power(Vz, 3.0) * Oz * sinh(L * Vz)) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy * (c_y * c_y)) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy * (c_z * c_z)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz * (c_y * c_y)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz * (c_z * c_z)) - 4.0 * (L * L) * power(Vy, 4.0) * Oy * (c_y * c_y) * (c_z * c_z)) + 4.0 * (L * L) * power(Vz, 4.0) * Oz * (c_y * c_y) * (c_z * c_z)) + power(L, 3.0) * Vy * power(Vz, 4.0) * Oz * sinh(L * Vy)) - power(L, 3.0) * power(Vy, 4.0) * Vz * Oy * sinh(L * Vz)) + 4.0 * L * Vy * (Vz * Vz) * cosh(L * Vy / 2.0) * (c_z * c_z) * sinh(L * Vy / 2.0)) - 4.0 * L * (Vy * Vy) * Vz * (c_y * c_y) * cosh(L * Vz / 2.0) * sinh(L * Vz / 2.0)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy * (c_y * c_y) * (c_z * c_z)) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz * (c_y * c_y) * (c_z * c_z)) + 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz * (c_y * c_y)) - 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_y * c_y)) + 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz * (c_z * c_z)) - 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_z * c_z)) - 2.0 * power(L, 3.0) * power(Vy, 3.0) * (Vz * Vz) * Oy * cosh(L * Vy / 2.0) * (c_z * c_z) * sinh(L * Vy / 2.0)) + 2.0 * power(L, 3.0) * (Vy * Vy) * power(Vz, 3.0) * Oz * (c_y * c_y) * cosh(L * Vz / 2.0) * sinh(L * Vz / 2.0)) - 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz * (c_y * c_y) * (c_z * c_z)) + 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_y * c_y) * (c_z * c_z)) - 2.0 * power(L, 3.0) * Vy * power(Vz, 4.0) * Oz * cosh(L * Vy / 2.0) * (c_z * c_z) * sinh(L * Vy / 2.0)) + 2.0 * power(L, 3.0) * power(Vy, 4.0) * Vz * Oy * (c_y * c_y) * cosh(L * Vz / 2.0) * sinh(L * Vz / 2.0))) / ((Vy * Vy - Vz * Vz) * ((((L * Vy * sinh(L * Vy) - 2.0 * (L * L) * (Vy * Vy) * Oy) - 2.0 * cosh(L * Vy)) + 2.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) + 2.0) * ((((L * Vz * sinh(L * Vz) - 2.0 * (L * L) * (Vz * Vz) * Oz) - 2.0 * cosh(L * Vz)) + 2.0 * (L * L) * (Vz * Vz) * Oz * cosh(L * Vz)) + 2.0))

    Ket[0 , 6] = Ket[6 , 0] = -P / L - Ax * E / L
    Ket[1 , 7] = Ket[7 , 1] = -(Vy * ((((((((((((Iz * (Vy * Vy) * P * sinh(L * Vy) - 3.0 * A * P * sinh(L * Vy)) + 2.0 * A * L * Vy * P) - Iz * L * power(Vy, 3.0) * P) - A * E * Iz * L * power(Vy, 3.0)) + A * E * Iz * (Vy * Vy) * sinh(L * Vy)) - 2.0 * A * power(L, 3.0) * power(Vy, 3.0) * Oy * P) + A * L * Vy * P * cosh(L * Vy)) + A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P) + 2.0 * A * (L * L) * (Vy * Vy) * Oy * P * sinh(L * Vy)) + A * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * sinh(L * Vy)) + A * E * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy) + A * E * Iz * (L * L) * power(Vy, 4.0) * Oy * sinh(L * Vy))) / (2.0 * A * (ay * ay))
    Ket[2 , 8] = Ket[8 , 2] = -(Vz * ((((((((((((Iy * (Vz * Vz) * P * sinh(L * Vz) - 3.0 * A * P * sinh(L * Vz)) + 2.0 * A * L * Vz * P) - Iy * L * power(Vz, 3.0) * P) - A * E * Iy * L * power(Vz, 3.0)) + A * E * Iy * (Vz * Vz) * sinh(L * Vz)) - 2.0 * A * power(L, 3.0) * power(Vz, 3.0) * Oz * P) + A * L * Vz * P * cosh(L * Vz)) + A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P) + 2.0 * A * (L * L) * (Vz * Vz) * Oz * P * sinh(L * Vz)) + A * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * sinh(L * Vz)) + A * E * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz) + A * E * Iy * (L * L) * power(Vz, 4.0) * Oz * sinh(L * Vz))) / (2.0 * A * (az * az))
    Ket[3 , 9] = Ket[9 , 3] = -(K / L) - (G*Ix / L)
    x = (((((((((((((((((((2.0 * A * P * sinh(L * Vz) - A * P * sinh(2.0 * L * Vz)) + 2.0 * Iy * (Vz * Vz) * P * sinh(L * Vz)) - Iy * (Vz * Vz) * P * sinh(2.0 * L * Vz)) - 6.0 * A * L * Vz * P) + Iy * L * power(Vz, 3.0) * P) + A * E * Iy * L * power(Vz, 3.0)) + 2.0 * A * E * Iy * (Vz * Vz) * sinh(L * Vz)) - A * E * Iy * (Vz * Vz) * sinh(2.0 * L * Vz)) - 2.0 * Iy * L * power(Vz, 3.0) * P * cosh(L * Vz)) + Iy * L * power(Vz, 3.0) * P * cosh(2.0 * L * Vz)) + 17.0 * A * power(L, 3.0) * power(Vz, 3.0) * Oz * P) - 3.0 * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * P) + A * power(L, 3.0) * power(Vz, 3.0) * P * cosh(L * Vz)) - Iy * power(L, 3.0) * power(Vz, 5.0) * P * cosh(L * Vz)) - 3.0 * A * (L * L) * (Vz * Vz) * P * sinh(L * Vz)) + Iy * (L * L) * power(Vz, 4.0) * P * sinh(L * Vz)) + 6.0 * A * L * Vz * P * cosh(L * Vz)) - 18.0 * A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P) + 9.0 * A * power(L, 7.0) * power(Vz, 7.0) * power(Oz, 3.0) * P) - 2.0 * A * power(L, 9.0) * power(Vz, 9.0) * power(Oz, 4.0) * P
    Ket[4 , 10] = Ket[10 , 4] = (((((((((((((((((((((((((((((((((((((((((((((((((x + cz) - 2.0 * dz * sinh(L * Vz)) + dz * sinh(2.0 * L * Vz)) + 8.0 * ez * sinh(L * Vz)) - 2.0 * fz * sinh(L * Vz)) + fz * sinh(2.0 * L * Vz)) - 2.0 * Iy * power(L, 4.0) * power(Vz, 6.0) * Oz * P * sinh(L * Vz)) + 7.0 * A * E * Iy * power(L, 5.0) * power(Vz, 7.0) * (Oz * Oz)) - 2.0 * A * E * Iy * power(L, 7.0) * power(Vz, 9.0) * power(Oz, 3.0)) + 16.0 * A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P * cosh(L * Vz)) + 2.0 * A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P * cosh(2.0 * L * Vz)) + A * power(L, 7.0) * power(Vz, 7.0) * (Oz * Oz) * P * cosh(L * Vz)) - 8.0 * A * power(L, 7.0) * power(Vz, 7.0) * power(Oz, 3.0) * P * cosh(L * Vz)) - A * power(L, 7.0) * power(Vz, 7.0) * power(Oz, 3.0) * P * cosh(2.0 * L * Vz)) + 2.0 * A * power(L, 9.0) * power(Vz, 9.0) * power(Oz, 4.0) * P * cosh(L * Vz)) - 2.0 * Iy * power(L, 5.0) * power(Vz, 7.0) * (Oz * Oz) * P * cosh(L * Vz)) + 12.0 * A * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * sinh(L * Vz)) - 6.0 * A * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * sinh(2.0 * L * Vz)) - 7.0 * A * power(L, 6.0) * power(Vz, 6.0) * (Oz * Oz) * P * sinh(L * Vz)) - 8.0 * A * power(L, 6.0) * power(Vz, 6.0) * power(Oz, 3.0) * P * sinh(L * Vz)) + 4.0 * A * power(L, 6.0) * power(Vz, 6.0) * power(Oz, 3.0) * P * sinh(2.0 * L * Vz)) + 2.0 * A * power(L, 8.0) * power(Vz, 8.0) * power(Oz, 3.0) * P * sinh(L * Vz)) + 2.0 * A * power(L, 8.0) * power(Vz, 8.0) * power(Oz, 4.0) * P * sinh(L * Vz)) - A * power(L, 8.0) * power(Vz, 8.0) * power(Oz, 4.0) * P * sinh(2.0 * L * Vz)) + 2.0 * Iy * power(L, 4.0) * power(Vz, 6.0) * (Oz * Oz) * P * sinh(L * Vz)) - Iy * power(L, 4.0) * power(Vz, 6.0) * (Oz * Oz) * P * sinh(2.0 * L * Vz)) - 2.0 * A * E * Iy * L * power(Vz, 3.0) * cosh(L * Vz)) + A * E * Iy * L * power(Vz, 3.0) * cosh(2.0 * L * Vz)) - 6.0 * A * E * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz) - A * E * Iy * power(L, 3.0) * power(Vz, 5.0) * cosh(L * Vz)) + A * E * Iy * (L * L) * power(Vz, 4.0) * sinh(L * Vz)) - 16.0 * A * power(L, 3.0) * power(Vz, 3.0) * Oz * P * cosh(L * Vz)) - A * power(L, 3.0) * power(Vz, 3.0) * Oz * P * cosh(2.0 * L * Vz)) - 2.0 * A * power(L, 5.0) * power(Vz, 5.0) * Oz * P * cosh(L * Vz)) + 4.0 * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * P * cosh(L * Vz)) - Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * P * cosh(2.0 * L * Vz)) + 6.0 * A * E * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * cosh(L * Vz)) + A * E * Iy * power(L, 5.0) * power(Vz, 7.0) * Oz * cosh(L * Vz)) - 2.0 * A * E * Iy * (L * L) * power(Vz, 4.0) * Oz * sinh(L * Vz)) + A * E * Iy * (L * L) * power(Vz, 4.0) * Oz * sinh(2.0 * L * Vz)) - 5.0 * A * E * Iy * power(L, 4.0) * power(Vz, 6.0) * Oz * sinh(L * Vz)) - 6.0 * A * E * Iy * power(L, 5.0) * power(Vz, 7.0) * (Oz * Oz) * cosh(L * Vz)) - A * E * Iy * power(L, 5.0) * power(Vz, 7.0) * (Oz * Oz) * cosh(2.0 * L * Vz)) + 2.0 * A * E * Iy * power(L, 7.0) * power(Vz, 9.0) * power(Oz, 3.0) * cosh(L * Vz)) - 2.0 * A * E * Iy * power(L, 4.0) * power(Vz, 6.0) * (Oz * Oz) * sinh(L * Vz)) + A * E * Iy * power(L, 4.0) * power(Vz, 6.0) * (Oz * Oz) * sinh(2.0 * L * Vz)) + 2.0 * A * E * Iy * power(L, 6.0) * power(Vz, 8.0) * (Oz * Oz) * sinh(L * Vz)) + 2.0 * A * E * Iy * power(L, 6.0) * power(Vz, 8.0) * power(Oz, 3.0) * sinh(L * Vz)) - A * E * Iy * power(L, 6.0) * power(Vz, 8.0) * power(Oz, 3.0) * sinh(2.0 * L * Vz)) / (A * Vz * ((((((((((((((4.0 * cosh(2.0 * L * Vz) - 16.0 * cosh(L * Vz)) - L * L * (Vz * Vz)) - 24.0 * (L * L) * (Vz * Vz) * Oz) + L * L * (Vz * Vz) * cosh(2.0 * L * Vz)) + 8.0 * L * Vz * sinh(L * Vz)) - 4.0 * L * Vz * sinh(2.0 * L * Vz)) + 12.0 * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz)) + 32.0 * (L * L) * (Vz * Vz) * Oz * cosh(L * Vz)) - 8.0 * (L * L) * (Vz * Vz) * Oz * cosh(2.0 * L * Vz)) - 8.0 * power(L, 3.0) * power(Vz, 3.0) * Oz * sinh(L * Vz)) + 4.0 * power(L, 3.0) * power(Vz, 3.0) * Oz * sinh(2.0 * L * Vz)) - 16.0 * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cosh(L * Vz)) + 4.0 * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * cosh(2.0 * L * Vz)) + 12.0))
    x = (((((((((((((((((((2.0 * A * P * sinh(L * Vy) - A * P * sinh(2.0 * L * Vy)) + 2.0 * Iz * (Vy * Vy) * P * sinh(L * Vy)) - Iz * (Vy * Vy) * P * sinh(2.0 * L * Vy)) - 6.0 * A * L * Vy * P) + Iz * L * power(Vy, 3.0) * P) + A * E * Iz * L * power(Vy, 3.0)) + 2.0 * A * E * Iz * (Vy * Vy) * sinh(L * Vy)) - A * E * Iz * (Vy * Vy) * sinh(2.0 * L * Vy)) - 2.0 * Iz * L * power(Vy, 3.0) * P * cosh(L * Vy)) + Iz * L * power(Vy, 3.0) * P * cosh(2.0 * L * Vy)) + 17.0 * A * power(L, 3.0) * power(Vy, 3.0) * Oy * P) - 3.0 * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * P) + A * power(L, 3.0) * power(Vy, 3.0) * P * cosh(L * Vy)) - Iz * power(L, 3.0) * power(Vy, 5.0) * P * cosh(L * Vy)) - 3.0 * A * (L * L) * (Vy * Vy) * P * sinh(L * Vy)) + Iz * (L * L) * power(Vy, 4.0) * P * sinh(L * Vy)) + 6.0 * A * L * Vy * P * cosh(L * Vy)) - 18.0 * A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P) + 9.0 * A * power(L, 7.0) * power(Vy, 7.0) * power(Oy, 3.0) * P) - 2.0 * A * power(L, 9.0) * power(Vy, 9.0) * power(Oy, 4.0) * P
    Ket[5 , 11] = Ket[11 , 5] = (((((((((((((((((((((((((((((((((((((((((((((((((x + cy) - 2.0 * dy * sinh(L * Vy)) + dy * sinh(2.0 * L * Vy)) + 8.0 * ey * sinh(L * Vy)) - 2.0 * fy * sinh(L * Vy)) + fy * sinh(2.0 * L * Vy)) - 2.0 * Iz * power(L, 4.0) * power(Vy, 6.0) * Oy * P * sinh(L * Vy)) + 7.0 * A * E * Iz * power(L, 5.0) * power(Vy, 7.0) * (Oy * Oy)) - 2.0 * A * E * Iz * power(L, 7.0) * power(Vy, 9.0) * power(Oy, 3.0)) + 16.0 * A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P * cosh(L * Vy)) + 2.0 * A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P * cosh(2.0 * L * Vy)) + A * power(L, 7.0) * power(Vy, 7.0) * (Oy * Oy) * P * cosh(L * Vy)) - 8.0 * A * power(L, 7.0) * power(Vy, 7.0) * power(Oy, 3.0) * P * cosh(L * Vy)) - A * power(L, 7.0) * power(Vy, 7.0) * power(Oy, 3.0) * P * cosh(2.0 * L * Vy)) + 2.0 * A * power(L, 9.0) * power(Vy, 9.0) * power(Oy, 4.0) * P * cosh(L * Vy)) - 2.0 * Iz * power(L, 5.0) * power(Vy, 7.0) * (Oy * Oy) * P * cosh(L * Vy)) + 12.0 * A * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * sinh(L * Vy)) - 6.0 * A * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * sinh(2.0 * L * Vy)) - 7.0 * A * power(L, 6.0) * power(Vy, 6.0) * (Oy * Oy) * P * sinh(L * Vy)) - 8.0 * A * power(L, 6.0) * power(Vy, 6.0) * power(Oy, 3.0) * P * sinh(L * Vy)) + 4.0 * A * power(L, 6.0) * power(Vy, 6.0) * power(Oy, 3.0) * P * sinh(2.0 * L * Vy)) + 2.0 * A * power(L, 8.0) * power(Vy, 8.0) * power(Oy, 3.0) * P * sinh(L * Vy)) + 2.0 * A * power(L, 8.0) * power(Vy, 8.0) * power(Oy, 4.0) * P * sinh(L * Vy)) - A * power(L, 8.0) * power(Vy, 8.0) * power(Oy, 4.0) * P * sinh(2.0 * L * Vy)) + 2.0 * Iz * power(L, 4.0) * power(Vy, 6.0) * (Oy * Oy) * P * sinh(L * Vy)) - Iz * power(L, 4.0) * power(Vy, 6.0) * (Oy * Oy) * P * sinh(2.0 * L * Vy)) - 2.0 * A * E * Iz * L * power(Vy, 3.0) * cosh(L * Vy)) + A * E * Iz * L * power(Vy, 3.0) * cosh(2.0 * L * Vy)) - 6.0 * A * E * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy) - A * E * Iz * power(L, 3.0) * power(Vy, 5.0) * cosh(L * Vy)) + A * E * Iz * (L * L) * power(Vy, 4.0) * sinh(L * Vy)) - 16.0 * A * power(L, 3.0) * power(Vy, 3.0) * Oy * P * cosh(L * Vy)) - A * power(L, 3.0) * power(Vy, 3.0) * Oy * P * cosh(2.0 * L * Vy)) - 2.0 * A * power(L, 5.0) * power(Vy, 5.0) * Oy * P * cosh(L * Vy)) + 4.0 * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * P * cosh(L * Vy)) - Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * P * cosh(2.0 * L * Vy)) + 6.0 * A * E * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * cosh(L * Vy)) + A * E * Iz * power(L, 5.0) * power(Vy, 7.0) * Oy * cosh(L * Vy)) - 2.0 * A * E * Iz * (L * L) * power(Vy, 4.0) * Oy * sinh(L * Vy)) + A * E * Iz * (L * L) * power(Vy, 4.0) * Oy * sinh(2.0 * L * Vy)) - 5.0 * A * E * Iz * power(L, 4.0) * power(Vy, 6.0) * Oy * sinh(L * Vy)) - 6.0 * A * E * Iz * power(L, 5.0) * power(Vy, 7.0) * (Oy * Oy) * cosh(L * Vy)) - A * E * Iz * power(L, 5.0) * power(Vy, 7.0) * (Oy * Oy) * cosh(2.0 * L * Vy)) + 2.0 * A * E * Iz * power(L, 7.0) * power(Vy, 9.0) * power(Oy, 3.0) * cosh(L * Vy)) - 2.0 * A * E * Iz * power(L, 4.0) * power(Vy, 6.0) * (Oy * Oy) * sinh(L * Vy)) + A * E * Iz * power(L, 4.0) * power(Vy, 6.0) * (Oy * Oy) * sinh(2.0 * L * Vy)) + 2.0 * A * E * Iz * power(L, 6.0) * power(Vy, 8.0) * (Oy * Oy) * sinh(L * Vy)) + 2.0 * A * E * Iz * power(L, 6.0) * power(Vy, 8.0) * power(Oy, 3.0) * sinh(L * Vy)) - A * E * Iz * power(L, 6.0) * power(Vy, 8.0) * power(Oy, 3.0) * sinh(2.0 * L * Vy)) / (A * Vy * ((((((((((((((4.0 * cosh(2.0 * L * Vy) - 16.0 * cosh(L * Vy)) - L * L * (Vy * Vy)) - 24.0 * (L * L) * (Vy * Vy) * Oy) + L * L * (Vy * Vy) * cosh(2.0 * L * Vy)) + 8.0 * L * Vy * sinh(L * Vy)) - 4.0 * L * Vy * sinh(2.0 * L * Vy)) + 12.0 * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy)) + 32.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) - 8.0 * (L * L) * (Vy * Vy) * Oy * cosh(2.0 * L * Vy)) - 8.0 * power(L, 3.0) * power(Vy, 3.0) * Oy * sinh(L * Vy)) + 4.0 * power(L, 3.0) * power(Vy, 3.0) * Oy * sinh(2.0 * L * Vy)) - 16.0 * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cosh(L * Vy)) + 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * cosh(2.0 * L * Vy)) + 12.0))

    Ket[2 , 9] = Ket[9 , 2] = Mzb / L
    Ket[3 , 10] = Ket[10 , 3] = (Oz * (L * L) * (Vz * Vz) - 1.0) * (L * Vz * cosh(L * Vz / 2.0) / sinh(L * Vz / 2.0) - 2.0) * (Mza + Mzb) / (L * L * (Vz * Vz))
    #Ket[4 , 11] = Ket[11 , 4] = 0

    Ket[1 , 9] = Ket[9 , 1] = Myb / L
    Ket[2 , 10] = Ket[10 , 2] = -((((((((((((((A * (L * L) * (Vz * Vz) * P - 8.0 * A * P * (sz * sz)) - Iy * (L * L) * power(Vz, 4.0) * P) + Iy * L * power(Vz, 3.0) * P * sinh(L * Vz)) - A * E * Iy * (L * L) * power(Vz, 4.0)) - 2.0 * A * power(L, 4.0) * power(Vz, 4.0) * Oz * P) + A * L * Vz * P * sinh(L * Vz)) + A * power(L, 6.0) * power(Vz, 6.0) * (Oz * Oz) * P) - 2.0 * A * power(L, 3.0) * power(Vz, 3.0) * Oz * P * sinh(L * Vz)) + 16.0 * A * (L * L) * (Vz * Vz) * Oz * P * (sz * sz)) + A * power(L, 5.0) * power(Vz, 5.0) * (Oz * Oz) * P * sinh(L * Vz)) + A * E * Iy * L * power(Vz, 3.0) * sinh(L * Vz)) - 8.0 * A * power(L, 4.0) * power(Vz, 4.0) * (Oz * Oz) * P * (sz * sz)) + A * E * Iy * power(L, 4.0) * power(Vz, 6.0) * Oz) + A * E * Iy * power(L, 3.0) * power(Vz, 5.0) * Oz * sinh(L * Vz)) / (4.0 * A * (bz * bz))
    Ket[3 , 11] = Ket[11 , 3] = -((Oy * (L * L) * (Vy * Vy) - 1.0) * (L * Vy * cosh(L * Vy / 2.0) / sinh(L * Vy / 2.0) - 2.0) * (Mya + Myb)) / (L * L * (Vy * Vy))

    #Ket[1 , 10] = Ket[10 , 1] = -(Vy * Vy * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0))) / (sinh(L * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vy * cosh(L * Vy / 2.0) - 2.0 * sinh(L * Vy / 2.0)) + 2.0 * (L * L) * (Vy * Vy) * Oy * sinh(L * Vy / 2.0)))
    #Ket[2 , 11] = Ket[11 , 2] = -(Vz * Vz * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0))) / (sinh(L * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vz * cosh(L * Vz / 2.0) - 2.0 * sinh(L * Vz / 2.0)) + 2.0 * (L * L) * (Vz * Vz) * Oz * sinh(L * Vz / 2.0)))

    Ket[0 , 10] = Ket[10 , 0] = -Myb / L
    Ket[1 , 11] = Ket[11 , 1] = ((((((((((((((A * (L * L) * (Vy * Vy) * P - 8.0 * A * P * (sy * sy)) - Iz * (L * L) * power(Vy, 4.0) * P) + Iz * L * power(Vy, 3.0) * P * sinh(L * Vy)) - A * E * Iz * (L * L) * power(Vy, 4.0)) - 2.0 * A * power(L, 4.0) * power(Vy, 4.0) * Oy * P) + A * L * Vy * P * sinh(L * Vy)) + A * power(L, 6.0) * power(Vy, 6.0) * (Oy * Oy) * P) - 2.0 * A * power(L, 3.0) * power(Vy, 3.0) * Oy * P * sinh(L * Vy)) + 16.0 * A * (L * L) * (Vy * Vy) * Oy * P * (sy * sy)) + A * power(L, 5.0) * power(Vy, 5.0) * (Oy * Oy) * P * sinh(L * Vy)) + A * E * Iz * L * power(Vy, 3.0) * sinh(L * Vy)) - 8.0 * A * power(L, 4.0) * power(Vy, 4.0) * (Oy * Oy) * P * (sy * sy)) + A * E * Iz * power(L, 4.0) * power(Vy, 6.0) * Oy) + A * E * Iz * power(L, 3.0) * power(Vy, 5.0) * Oy * sinh(L * Vy)) / (4.0 * A * (by * by))

    Ket[0 , 11] = Ket[11 , 0] = -Mzb / L

    s_y = sinh(L * Vy / 4.0)
    s_z = sinh(L * Vz / 4.0)
    #Ket[4 , 5] =   Ket[5 , 4] =    ((L * L * ((2.0 * (Vy * Vy) * (Vz * Vz) * Mxb * (((sy * sy - 2.0 * (sy * sy + 1.0) * (sz * sz + 1.0)) + sz * sz) + 2.0) + Vy * power(Vz, 3.0) * Mxb * sinh(L * Vy) * sinh(L * Vz) / 2.0) + power(Vy, 3.0) * Vz * Mxb * sinh(L * Vy) * sinh(L * Vz) / 2.0) + L * (((2.0 * power(Vy, 3.0) * Mxb * (sinh(L * Vy) - 2.0 * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)) + 2.0 * power(Vz, 3.0) * Mxb * (sinh(L * Vz) - 2.0 * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0))) - 2.0 * Vy * (Vz * Vz) * Mxb * (sinh(L * Vy) - 2.0 * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0))) - 2.0 * (Vy * Vy) * Vz * Mxb * (sinh(L * Vz) - 2.0 * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)))) + power(L, 3.0) * (2.0 * power(Vy, 3.0) * (Vz * Vz) * Mxb * (((Oy * sinh(L * Vy) - Oz * sinh(L * Vy)) - 2.0 * Oy * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)) + 2.0 * Oz * sinh(L * Vy / 2.0) * (sy * sy + 1.0) * (2.0 * (s_y * s_y) + 1.0)) - 2.0 * (Vy * Vy) * power(Vz, 3.0) * Mxb * (((Oy * sinh(L * Vz) - Oz * sinh(L * Vz)) - 2.0 * Oy * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)) + 2.0 * Oz * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)))) / ((Vy * Vy - Vz * Vz) * ((L * Vy * sinh(L * Vy) - 4.0 * (sy * sy)) + 4.0 * (L * L) * (Vy * Vy) * Oy * (sy * sy)) * ((L * Vz * sinh(L * Vz) - 4.0 * (sz * sz)) + 4.0 * (L * L) * (Vz * Vz) * Oz * (sz * sz)))
    #Ket[10 , 11] = Ket[11 , 10] = -((L * L * ((2.0 * (Vy * Vy) * (Vz * Vz) * Mxb * (((sy * sy - 2.0 * (sy * sy + 1.0) * (sz * sz + 1.0)) + sz * sz) + 2.0) + Vy * power(Vz, 3.0) * Mxb * sinh(L * Vy) * sinh(L * Vz) / 2.0) + power(Vy, 3.0) * Vz * Mxb * sinh(L * Vy) * sinh(L * Vz) / 2.0) + L * (((2.0 * power(Vy, 3.0) * Mxb * (sinh(L * Vy) - 2.0 * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)) + 2.0 * power(Vz, 3.0) * Mxb * (sinh(L * Vz) - 2.0 * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0))) - 2.0 * Vy * (Vz * Vz) * Mxb * (sinh(L * Vy) - 2.0 * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0))) - 2.0 * (Vy * Vy) * Vz * Mxb * (sinh(L * Vz) - 2.0 * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)))) + power(L, 3.0) * (2.0 * power(Vy, 3.0) * (Vz * Vz) * Mxb * (((Oy * sinh(L * Vy) - Oz * sinh(L * Vy)) - 2.0 * Oy * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)) + 2.0 * Oz * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)) - 2.0 * (Vy * Vy) * power(Vz, 3.0) * Mxb * (((Oy * sinh(L * Vz) - Oz * sinh(L * Vz)) - 2.0 * Oy * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)) + 2.0 * Oz * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)))) / ((Vy * Vy - Vz * Vz) * ((L * Vy * sinh(L * Vy) - 4.0 * (sy * sy)) + 4.0 * (L * L) * (Vy * Vy) * Oy * (sy * sy)) * ((L * Vz * sinh(L * Vz) - 4.0 * (sz * sz)) + 4.0 * (L * L) * (Vz * Vz) * Oz * (sz * sz)))

    #/* --------------------------------------------- Torsion for symetric sections --------------------------------------------- */
    if Iy == Iz:
        #
        Ket[1 , 4] = Ket[4 , 1] = Vy * Mxb * (sinh(L * Vy) - L * Vy) / (2.0 * ((((L * Vy * sinh(L * Vy) - 2.0 * (L * L) * (Vy * Vy) * Oy) - 2.0 * cosh(L * Vy)) + 2.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) + 2.0))
        Ket[2 , 5] = Ket[5 , 2] = Vy * Mxb * (sinh(L * Vy) - L * Vy) / (2.0 * ((((L * Vy * sinh(L * Vy) - 2.0 * (L * L) * (Vy * Vy) * Oz) - 2.0 * cosh(L * Vy)) + 2.0 * (L * L) * (Vy * Vy) * Oz * cosh(L * Vy)) + 2.0))
        Ket[4 , 7] = Ket[7 , 4] = -(Vy * Mxb * (sinh(L * Vy) - L * Vy)) / (2.0 * ((((L * Vy * sinh(L * Vy) - 2.0 * (L * L) * (Vy * Vy) * Oy) - 2.0 * cosh(L * Vy)) + 2.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) + 2.0))
        Ket[5 , 8] = Ket[8 , 5] = -(Vy * Mxb * (sinh(L * Vy) - L * Vy)) / (2.0 * ((((L * Vy * sinh(L * Vy) - 2.0 * (L * L) * (Vy * Vy) * Oz) - 2.0 * cosh(L * Vy)) + 2.0 * (L * L) * (Vy * Vy) * Oz * cosh(L * Vy)) + 2.0))
        Ket[7 , 10] = Ket[10 , 7] = Vy * Mxb * (sinh(L * Vy) - L * Vy) / (2.0 * ((((L * Vy * sinh(L * Vy) - 2.0 * (L * L) * (Vy * Vy) * Oy) - 2.0 * cosh(L * Vy)) + 2.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) + 2.0))
        Ket[8 , 11] = Ket[11 , 8] = Vy * Mxb * (sinh(L * Vy) - L * Vy) / (2.0 * ((((L * Vy * sinh(L * Vy) - 2.0 * (L * L) * (Vy * Vy) * Oz) - 2.0 * cosh(L * Vy)) + 2.0 * (L * L) * (Vy * Vy) * Oz * cosh(L * Vy)) + 2.0))

        Ket[5 , 10] = Ket[10 , 5] = -(Mxb * ((((((((((((((((((((((((6.0 * Vy - 8.0 * Vy * cosh(L * Vy)) + 2.0 * Vy * cosh(2.0 * L * Vy)) + 2.0 * (L * L) * power(Vy, 3.0)) - 6.0 * (L * L) * power(Vy, 3.0) * Oy) - 6.0 * (L * L) * power(Vy, 3.0) * Oz) - power(L, 4.0) * power(Vy, 5.0) * Oy) - power(L,4.0) * power(Vy, 5.0) * Oz) - 2.0 * (L * L) * power(Vy, 3.0) * cosh(L * Vy)) + power(L, 3.0) * power(Vy, 4.0) * sinh(L * Vy)) + 2.0 * L * (Vy * Vy) * sinh(L * Vy)) - L * (Vy * Vy) * sinh(2.0 * L * Vy)) + 6.0 * power(L, 4.0) * power(Vy, 5.0) * Oy * Oz) + 8.0 * (L * L) * power(Vy, 3.0) * Oy * cosh(L * Vy)) + 8.0 * (L * L) * power(Vy,3.0) * Oz * cosh(L * Vy)) - 2.0 * (L * L) * power(Vy, 3.0) * Oy * cosh(2.0 * L * Vy)) - 2.0 * (L * L) * power(Vy, 3.0) * Oz * cosh(2.0 * L * Vy)) + power(L, 4.0) * power(Vy, 5.0) * Oy * cosh(L * Vy)) + power(L, 4.0) * power(Vy, 5.0) * Oz * cosh(L * Vy)) - power(L, 3.0) * power(Vy, 4.0) * Oy * sinh(L * Vy)) - power(L, 3.0) * power(Vy, 4.0) * Oz * sinh(L * Vy)) + power(L, 3.0) * power(Vy, 4.0) * Oy * sinh(2.0 * L * Vy) / 2.0) + power(L, 3.0) * power(Vy, 4.0) * Oz * sinh(2.0 * L * Vy) / 2.0) - 8.0 * power(L, 4.0) * power(Vy, 5.0) * Oy * Oz * cosh(L * Vy)) + 2.0 * power(L, 4.0) * power(Vy, 5.0) * Oy * Oz * cosh(2.0 * L * Vy))) / (Vy * (((((((((((((((((((16.0 * cosh(L * Vy) - 4.0 * cosh(2.0 * L * Vy)) + L * L * (Vy * Vy)) + 12.0 * (L * L) * (Vy * Vy) * Oy) + 12.0 * (L * L) * (Vy * Vy) * Oz) - L * L * (Vy * Vy) * cosh(2.0 * L * Vy)) - 8.0 * L * Vy * sinh(L * Vy)) + 4.0 * L * Vy * sinh(2.0 * L * Vy)) - 12.0 * power(L, 4.0) * power(Vy, 4.0) * Oy * Oz) - 16.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) - 16.0 * (L * L) * (Vy * Vy) * Oz * cosh(L * Vy)) + 4.0 * (L * L) * (Vy * Vy) * Oy * cosh(2.0 * L * Vy)) + 4.0 * (L * L) * (Vy * Vy) * Oz * cosh(2.0 * L * Vy)) + 4.0 * power(L, 3.0) * power(Vy, 3.0) * Oy * sinh(L * Vy)) + 4.0 * power(L, 3.0) * power(Vy, 3.0) * Oz * sinh(L * Vy)) - 2.0 * power(L, 3.0) * power(Vy, 3.0) * Oy * sinh(2.0 * L * Vy)) - 2.0 * power(L, 3.0) * power(Vy, 3.0) * Oz * sinh(2.0 * L * Vy)) + 16.0 * power(L, 4.0) * power(Vy, 4.0) * Oy * Oz * cosh(L * Vy)) - 4.0 * power(L, 4.0) * power(Vy, 4.0) * Oy * Oz * cosh(2.0 * L * Vy)) - 12.0))
        Ket[4 , 11] = Ket[11 , 4] = Mxb * ((((((((((((((((((((((((6.0 * Vy - 8.0 * Vy * cosh(L * Vy)) + 2.0 * Vy * cosh(2.0 * L * Vy)) + 2.0 * (L * L) * power(Vy, 3.0)) - 6.0 * (L * L) * power(Vy, 3.0) * Oy) - 6.0 * (L * L) * power(Vy, 3.0) * Oz) - power(L, 4.0) * power(Vy, 5.0) * Oy) - power(L,4.0) * power(Vy, 5.0) * Oz) - 2.0 * (L * L) * power(Vy, 3.0) * cosh(L * Vy)) + power(L, 3.0) * power(Vy, 4.0) * sinh(L * Vy)) + 2.0 * L * (Vy * Vy) * sinh(L * Vy)) - L * (Vy * Vy) * sinh(2.0 * L * Vy)) + 6.0 * power(L, 4.0) * power(Vy, 5.0) * Oy * Oz) + 8.0 * (L * L) * power(Vy, 3.0) * Oy * cosh(L * Vy)) + 8.0 * (L * L) * power(Vy,3.0) * Oz * cosh(L * Vy)) - 2.0 * (L * L) * power(Vy, 3.0) * Oy * cosh(2.0 * L * Vy)) - 2.0 * (L * L) * power(Vy, 3.0) * Oz * cosh(2.0 * L * Vy)) + power(L, 4.0) * power(Vy, 5.0) * Oy * cosh(L * Vy)) + power(L, 4.0) * power(Vy, 5.0) * Oz * cosh(L * Vy)) - power(L, 3.0) * power(Vy, 4.0) * Oy * sinh(L * Vy)) - power(L, 3.0) * power(Vy, 4.0) * Oz * sinh(L * Vy)) + power(L, 3.0) * power(Vy, 4.0) * Oy * sinh(2.0 * L * Vy) / 2.0) + power(L, 3.0) * power(Vy, 4.0) * Oz * sinh(2.0 * L * Vy) / 2.0) - 8.0 * power(L, 4.0) * power(Vy, 5.0) * Oy * Oz * cosh(L * Vy)) + 2.0 * power(L, 4.0) * power(Vy, 5.0) * Oy * Oz * cosh(2.0 * L * Vy)) / (Vy * (((((((((((((((((((16.0 * cosh(L * Vy) - 4.0 * cosh(2.0 * L * Vy)) + L * L * (Vy * Vy)) + 12.0 * (L * L) * (Vy * Vy) * Oy) + 12.0 * (L * L) * (Vy * Vy) * Oz) - L * L * (Vy * Vy) * cosh(2.0 * L * Vy)) - 8.0 * L * Vy * sinh(L * Vy)) + 4.0 * L * Vy * sinh(2.0 * L * Vy)) - 12.0 * power(L, 4.0) * power(Vy, 4.0) * Oy * Oz) - 16.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) - 16.0 * (L * L) * (Vy * Vy) * Oz * cosh(L * Vy)) + 4.0 * (L * L) * (Vy * Vy) * Oy * cosh(2.0 * L * Vy)) + 4.0 * (L * L) * (Vy * Vy) * Oz * cosh(2.0 * L * Vy)) + 4.0 * power(L, 3.0) * power(Vy, 3.0) * Oy * sinh(L * Vy)) + 4.0 * power(L, 3.0) * power(Vy, 3.0) * Oz * sinh(L * Vy)) - 2.0 * power(L, 3.0) * power(Vy, 3.0) * Oy * sinh(2.0 * L * Vy)) - 2.0 * power(L, 3.0) * power(Vy, 3.0) * Oz * sinh(2.0 * L * Vy)) + 16.0 * power(L, 4.0) * power(Vy, 4.0) * Oy * Oz * cosh(L * Vy)) - 4.0 * power(L, 4.0) * power(Vy, 4.0) * Oy * Oz * cosh(2.0 * L * Vy)) - 12.0))

        Ket[1 , 10] = Ket[10 , 1] = -(Vy * Mxb * (sinh(L * Vy) - L * Vy)) / (2.0 * ((((L * Vy * sinh(L * Vy) - 2.0 * (L * L) * (Vy * Vy) * Oy) - 2.0 * cosh(L * Vy)) + 2.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) + 2.0))
        Ket[2 , 11] = Ket[11 , 2] = -(Vy * Mxb * (sinh(L * Vy) - L * Vy)) / (2.0 * ((((L * Vy * sinh(L * Vy) - 2.0 * (L * L) * (Vy * Vy) * Oz) - 2.0 * cosh(L * Vy)) + 2.0 * (L * L) * (Vy * Vy) * Oz * cosh(L * Vy)) + 2.0))

        Ket[4 , 5] = Ket[5 , 4] = -(power(L, 3.0) * power(Vy, 3.0) * Mxb * (Oy - Oz) * (sinh(L * Vy) - L * Vy)) / (2.0 * ((((((((((((4.0 * cosh(L * Vy) + L * L * (Vy * Vy)) + 4.0 * (L * L) * (Vy * Vy) * Oy) + 4.0 * (L * L) * (Vy * Vy) * Oz) + L * L * (Vy * Vy) * cosh(L * Vy)) - 4.0 * L * Vy * sinh(L * Vy)) - 4.0 * power(L, 4.0) * power(Vy, 4.0) * Oy * Oz) - 4.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) - 4.0 * (L * L) * (Vy * Vy) * Oz * cosh(L * Vy)) + 2.0 * power(L, 3.0) * power(Vy, 3.0) * Oy * sinh(L * Vy)) + 2.0 * power(L, 3.0) * power(Vy, 3.0) * Oz * sinh(L * Vy)) + 4.0 * power(L, 4.0) * power(Vy, 4.0) * Oy * Oz * cosh(L * Vy)) - 4.0))
        Ket[10 , 11] = Ket[11 , 10] = power(L, 3.0) * power(Vy, 3.0) * Mxb * (Oy - Oz) * (sinh(L * Vy) - L * Vy) / (2.0 * ((((((((((((4.0 * cosh(L * Vy) + L * L * (Vy * Vy)) + 4.0 * (L * L) * (Vy * Vy) * Oy) + 4.0 * (L * L) * (Vy * Vy) * Oz) + L * L * (Vy * Vy) * cosh(L * Vy)) - 4.0 * L * Vy * sinh(L * Vy)) - 4.0 * power(L, 4.0) * power(Vy, 4.0) * Oy * Oz) - 4.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) - 4.0 * (L * L) * (Vy * Vy) * Oz * cosh(L * Vy)) + 2.0 * power(L, 3.0) * power(Vy, 3.0) * Oy * sinh(L * Vy)) + 2.0 * power(L, 3.0) * power(Vy, 3.0) * Oz * sinh(L * Vy)) + 4.0 * power(L, 4.0) * power(Vy, 4.0) * Oy * Oz * cosh(L * Vy)) - 4.0))
    else:
        Ket[1 , 4] = Ket[4 , 1] = Vy * Vy * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0)) / (sinh(L * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vy * cosh(L * Vy / 2.0) - 2.0 * sinh(L * Vy / 2.0)) + 2.0 * (L * L) * (Vy * Vy) * Oy * sinh(L * Vy / 2.0)))
        Ket[2 , 5] = Ket[5 , 2] = Vz * Vz * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0)) / (sinh(L * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vz * cosh(L * Vz / 2.0) - 2.0 * sinh(L * Vz / 2.0)) + 2.0 * (L * L) * (Vz * Vz) * Oz * sinh(L * Vz / 2.0)))
        Ket[4 , 7] = Ket[7 , 4] = -(Vy * Vy * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0))) / (sinh(L * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vy * cosh(L * Vy / 2.0) - 2.0 * sinh(L * Vy / 2.0)) + 2.0 * (L * L) * (Vy * Vy) * Oy * sinh(L * Vy / 2.0)))
        Ket[5 , 8] = Ket[8 , 5] = -(Vz * Vz * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0))) / (sinh(L * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vz * cosh(L * Vz / 2.0) - 2.0 * sinh(L * Vz / 2.0)) + 2.0 * (L * L) * (Vz * Vz) * Oz * sinh(L * Vz / 2.0)))
        Ket[7 , 10] = Ket[10 , 7] = Vy * Vy * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0)) / (sinh(L * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vy * cosh(L * Vy / 2.0) - 2.0 * sinh(L * Vy / 2.0)) + 2.0 * (L * L) * (Vy * Vy) * Oy * sinh(L * Vy / 2.0)))
        Ket[8 , 11] = Ket[11 , 8] = Vz * Vz * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0)) / (sinh(L * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vz * cosh(L * Vz / 2.0) - 2.0 * sinh(L * Vz / 2.0)) + 2.0 * (L * L) * (Vz * Vz) * Oz * sinh(L * Vz / 2.0)))

        Ket[5 , 10] = Ket[10 , 5] = 2.0 * Mxb * (((((((((((((((((((((((((((((((((((((((((((((4.0 * Vy * Vy - 4.0 * Vz * Vz) - 4.0 * (Vy * Vy) * c_y * c_y) - 4.0 * Vy * Vy * (c_z * c_z)) + 4.0 * Vz * Vz * (c_y * c_y)) + 4.0 * (Vz * Vz) * (c_z * c_z)) - 4.0 * L * L * power(Vy, 4.0) * Oy) + 4.0 * (L * L) * power(Vz, 4.0) * Oz) + 4.0 * (Vy * Vy) * (c_y * c_y) * (c_z * c_z)) - 4.0 * (Vz * Vz) * (c_y * c_y) * (c_z * c_z)) - 2.0 * L * Vy * (Vz * Vz) * sinh(L * Vy)) + 2.0 * L * (Vy * Vy) * Vz * sinh(L * Vz)) + L * L * (Vy * Vy) * (Vz * Vz) * (c_y * c_y)) - L * L * (Vy * Vy) * (Vz * Vz) * (c_z * c_z)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz) + 4.0 * (L * L) * power(Vy, 4.0) * Oy * (c_y * c_y)) + 4.0 * (L * L) * power(Vy, 4.0) * Oy * (c_z * c_z)) - 4.0 * (L * L) * power(Vz, 4.0) * Oz * (c_y * c_y)) - 4.0 * (L * L) * power(Vz, 4.0) * Oz * (c_z * c_z)) - 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz) + 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz) + power(L, 3.0) * power(Vy, 3.0) * (Vz * Vz) * Oy * sinh(L * Vy)) - power(L, 3.0) * (Vy * Vy) * power(Vz, 3.0) * Oz * sinh(L * Vz)) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy * (c_y * c_y)) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy * (c_z * c_z)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz * (c_y * c_y)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz * (c_z * c_z)) - 4.0 * (L * L) * power(Vy, 4.0) * Oy * (c_y * c_y) * (c_z * c_z)) + 4.0 * (L * L) * power(Vz, 4.0) * Oz * (c_y * c_y) * (c_z * c_z)) + power(L, 3.0) * Vy * power(Vz, 4.0) * Oz * sinh(L * Vy)) - power(L, 3.0) * power(Vy, 4.0) * Vz * Oy * sinh(L * Vz)) + 4.0 * L * Vy * (Vz * Vz) * cosh(L * Vy / 2.0) * (c_z * c_z) * sinh(L * Vy / 2.0)) - 4.0 * L * (Vy * Vy) * Vz * (c_y * c_y) * cosh(L * Vz / 2.0) * sinh(L * Vz / 2.0)) + 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oy * (c_y * c_y) * (c_z * c_z)) - 4.0 * (L * L) * (Vy * Vy) * (Vz * Vz) * Oz * (c_y * c_y) * (c_z * c_z)) + 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz * (c_y * c_y)) - 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_y * c_y)) + 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz * (c_z * c_z)) - 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_z * c_z)) - 2.0 * power(L, 3.0) * power(Vy, 3.0) * (Vz * Vz) * Oy * cosh(L * Vy / 2.0) * (c_z * c_z) * sinh(L * Vy / 2.0)) + 2.0 * power(L, 3.0) * (Vy * Vy) * power(Vz, 3.0) * Oz * (c_y * c_y) * cosh(L * Vz / 2.0) * sinh(L * Vz / 2.0)) - 4.0 * power(L, 4.0) * (Vy * Vy) * power(Vz, 4.0) * Oy * Oz * (c_y * c_y) * (c_z * c_z)) + 4.0 * power(L, 4.0) * power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_y * c_y) * (c_z * c_z)) - 2.0 * power(L, 3.0) * Vy * power(Vz, 4.0) * Oz * cosh(L * Vy / 2.0) * (c_z * c_z) * sinh(L * Vy / 2.0)) + 2.0 * power(L, 3.0) * power(Vy, 4.0) * Vz * Oy * (c_y * c_y) * cosh(L * Vz / 2.0) * sinh(L * Vz / 2.0)) / ((Vy * Vy - Vz * Vz) * ((((L * Vy * sinh(L * Vy) - 2.0 * (L * L) * (Vy * Vy) * Oy) - 2.0 * cosh(L * Vy)) + 2.0 * (L * L) * (Vy * Vy) * Oy * cosh(L * Vy)) + 2.0) * ((((L * Vz * sinh(L * Vz) - 2.0 * (L * L) * (Vz * Vz) * Oz) - 2.0 * cosh(L * Vz)) + 2.0 * (L * L) * (Vz * Vz) * Oz * cosh(L * Vz)) + 2.0))
        Ket[4 , 11] = Ket[11 , 4] = 0

        Ket[1 , 10] = Ket[10 , 1] = -(Vy * Vy * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0))) / (sinh(L * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vy * cosh(L * Vy / 2.0) - 2.0 * sinh(L * Vy / 2.0)) + 2.0 * (L * L) * (Vy * Vy) * Oy * sinh(L * Vy / 2.0)))
        Ket[2 , 11] = Ket[11 , 2] = -(Vz * Vz * Mxb * (Vy * cosh(L * Vy / 2.0) * sinh(L * Vz / 2.0) - Vz * cosh(L * Vz / 2.0) * sinh(L * Vy / 2.0))) / (sinh(L * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((L * Vz * cosh(L * Vz / 2.0) - 2.0 * sinh(L * Vz / 2.0)) + 2.0 * (L * L) * (Vz * Vz) * Oz * sinh(L * Vz / 2.0)))

        Ket[4 , 5] =   Ket[5 , 4] =    ((L * L * ((2.0 * (Vy * Vy) * (Vz * Vz) * Mxb * (((sy * sy - 2.0 * (sy * sy + 1.0) * (sz * sz + 1.0)) + sz * sz) + 2.0) + Vy * power(Vz, 3.0) * Mxb * sinh(L * Vy) * sinh(L * Vz) / 2.0) + power(Vy, 3.0) * Vz * Mxb * sinh(L * Vy) * sinh(L * Vz) / 2.0) + L * (((2.0 * power(Vy, 3.0) * Mxb * (sinh(L * Vy) - 2.0 * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)) + 2.0 * power(Vz, 3.0) * Mxb * (sinh(L * Vz) - 2.0 * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0))) - 2.0 * Vy * (Vz * Vz) * Mxb * (sinh(L * Vy) - 2.0 * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0))) - 2.0 * (Vy * Vy) * Vz * Mxb * (sinh(L * Vz) - 2.0 * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)))) + power(L, 3.0) * (2.0 * power(Vy, 3.0) * (Vz * Vz) * Mxb * (((Oy * sinh(L * Vy) - Oz * sinh(L * Vy)) - 2.0 * Oy * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)) + 2.0 * Oz * sinh(L * Vy / 2.0) * (sy * sy + 1.0) * (2.0 * (s_y * s_y) + 1.0)) - 2.0 * (Vy * Vy) * power(Vz, 3.0) * Mxb * (((Oy * sinh(L * Vz) - Oz * sinh(L * Vz)) - 2.0 * Oy * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)) + 2.0 * Oz * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)))) / ((Vy * Vy - Vz * Vz) * ((L * Vy * sinh(L * Vy) - 4.0 * (sy * sy)) + 4.0 * (L * L) * (Vy * Vy) * Oy * (sy * sy)) * ((L * Vz * sinh(L * Vz) - 4.0 * (sz * sz)) + 4.0 * (L * L) * (Vz * Vz) * Oz * (sz * sz)))
        Ket[10 , 11] = Ket[11 , 10] = -((L * L * ((2.0 * (Vy * Vy) * (Vz * Vz) * Mxb * (((sy * sy - 2.0 * (sy * sy + 1.0) * (sz * sz + 1.0)) + sz * sz) + 2.0) + Vy * power(Vz, 3.0) * Mxb * sinh(L * Vy) * sinh(L * Vz) / 2.0) + power(Vy, 3.0) * Vz * Mxb * sinh(L * Vy) * sinh(L * Vz) / 2.0) + L * (((2.0 * power(Vy, 3.0) * Mxb * (sinh(L * Vy) - 2.0 * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)) + 2.0 * power(Vz, 3.0) * Mxb * (sinh(L * Vz) - 2.0 * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0))) - 2.0 * Vy * (Vz * Vz) * Mxb * (sinh(L * Vy) - 2.0 * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0))) - 2.0 * (Vy * Vy) * Vz * Mxb * (sinh(L * Vz) - 2.0 * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)))) + power(L, 3.0) * (2.0 * power(Vy, 3.0) * (Vz * Vz) * Mxb * (((Oy * sinh(L * Vy) - Oz * sinh(L * Vy)) - 2.0 * Oy * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)) + 2.0 * Oz * sinh(L * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)) - 2.0 * (Vy * Vy) * power(Vz, 3.0) * Mxb * (((Oy * sinh(L * Vz) - Oz * sinh(L * Vz)) - 2.0 * Oy * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)) + 2.0 * Oz * sinh(L * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)))) / ((Vy * Vy - Vz * Vz) * ((L * Vy * sinh(L * Vy) - 4.0 * (sy * sy)) + 4.0 * (L * L) * (Vy * Vy) * Oy * (sy * sy)) * ((L * Vz * sinh(L * Vz) - 4.0 * (sz * sz)) + 4.0 * (L * L) * (Vz * Vz) * Oz * (sz * sz)))
    #
    return Ket
#
def StfBeamTimoshenkoT(L:float, A:float,
                      Jx:float, Iy:float, Iz:float,
                       E:float, G:float,
                       Oy:float, Oz:float,
                       Fb:list[float]):
    """ """
    (Fxa, Fya, Fza, Mxa, Mya, Mza,
     Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb.q_loc
    #
    P = Fb.Fx    
    #
    sinh = np.sinh
    cosh = np.cosh
    sqrt =np.sqrt
    #
    #Oy = (E*Iz)/(G*Ay*L**2);     # Shear to Flexural rig. ratio y                         
    #Oz = (E*Iy)/(G*Az*L**2);     # Shear to Flexural rig. ratio z
    mu_y = sqrt(P / (E*Iz))
    mu_z = sqrt(P / (E*Iy))
    Vy = mu_y / (sqrt(1 + Oy*mu_y*mu_y*L**2))
    Vz = mu_z / (sqrt(1 + Oz*mu_z*mu_z*L**2))
    Jp = P *(Iy + Iz) / A

    elm_tangstf = np.zeros((12,12,), dtype=np.float32) 

    elm_tangstf[0,0] = P/L + (A*E)/L
    #elm_tangstf[0,6]  = - P/L - (A*E)/L
    elm_tangstf[0,6] = -elm_tangstf[0][0]
    elm_tangstf[6,6] = elm_tangstf[0][0]    
    
    #eq1_num = 0
    eq1_den = (2*A*(L*Vy*cosh((L*Vy)/2) - 2*sinh((L*Vy)/2) + 2*L**2*Vy**2*Oy*sinh((L*Vy)/2))**2)
    #eq12_num = 0
    eq12_den = (2*A*(L*Vz*cosh((L*Vz)/2) - 2*sinh((L*Vz)/2) + 2*L**2*Vz**2*Oz*sinh((L*Vz)/2))**2)

    elm_tangstf[1,1] = ((Vy*(Iz*Vy**2*P*sinh(L*Vy) - 3*A*P*sinh(L*Vy) + 2*A*L*Vy*P - Iz*L*Vy**3*P
                             - A*E*Iz*L*Vy**3 + A*E*Iz*Vy**2*sinh(L*Vy) - 2*A*L**3*Vy**3*Oy*P + A*L*Vy*P*cosh(L*Vy)
                             + A*L**5*Vy**5*Oy**2*P + 2*A*L**2*Vy**2*Oy*P*sinh(L*Vy) + A*L**4*Vy**4*Oy**2*P*sinh(L*Vy)
                             + A*E*Iz*L**3*Vy**5*Oy + A*E*Iz*L**2*Vy**4*Oy*sinh(L*Vy))) / eq1_den)

    #elm_tangstf[1,7]  = -((Vy*(Iz*Vy**2*P*sinh(L*Vy) - 3*A*P*sinh(L*Vy) + 2*A*L*Vy*P - Iz*L*Vy**3*P - A*E*Iz*L*Vy**3
    #                           + A*E*Iz*Vy**2*sinh(L*Vy) - 2*A*L**3*Vy**3*Oy*P + A*L*Vy*P*cosh(L*Vy) + A*L**5*Vy**5*Oy**2*P
    #                           + 2*A*L**2*Vy**2*Oy*P*sinh(L*Vy) + A*L**4*Vy**4*Oy**2*P*sinh(L*Vy) + A*E*Iz*L**3*Vy**5*Oy
    #                           + A*E*Iz*L**2*Vy**4*Oy*sinh(L*Vy))) / eq1_den)
    elm_tangstf[1,7] = -elm_tangstf[1,1]
    elm_tangstf[7,7] = elm_tangstf[1,1]
    

    elm_tangstf[2,2] = ((Vz*(Iy*Vz**2*P*sinh(L*Vz) - 3*A*P*sinh(L*Vz) + 2*A*L*Vz*P - Iy*L*Vz**3*P
                             - A*E*Iy*L*Vz**3 + A*E*Iy*Vz**2*sinh(L*Vz) - 2*A*L**3*Vz**3*Oz*P + A*L*Vz*P*cosh(L*Vz)
                             + A*L**5*Vz**5*Oz**2*P + 2*A*L**2*Vz**2*Oz*P*sinh(L*Vz) + A*L**4*Vz**4*Oz**2*P*sinh(L*Vz)
                             + A*E*Iy*L**3*Vz**5*Oz + A*E*Iy*L**2*Vz**4*Oz*sinh(L*Vz))) /eq12_den)

    #elm_tangstf[2,8]  = (-(Vz*(Iy*Vz**2*P*sinh(L*Vz) - 3*A*P*sinh(L*Vz) + 2*A*L*Vz*P - Iy*L*Vz**3*P - A*E*Iy*L*Vz**3 +
    #                           A*E*Iy*Vz**2*sinh(L*Vz) - 2*A*L**3*Vz**3*Oz*P + A*L*Vz*P*cosh(L*Vz) + A*L**5*Vz**5*Oz**2*P
    #                           + 2*A*L**2*Vz**2*Oz*P*sinh(L*Vz) + A*L**4*Vz**4*Oz**2*P*sinh(L*Vz) + A*E*Iy*L**3*Vz**5*Oz
    #                           + A*E*Iy*L**2*Vz**4*Oz*sinh(L*Vz))) / eq12_den)
    #
    elm_tangstf[2, 8] = -elm_tangstf[2, 2]
    elm_tangstf[8, 8] = elm_tangstf[2, 2]


    elm_tangstf[3,3] = (Jp*P)/(A*L)
    #elm_tangstf[3,9]  = -(Jp*P)/(A*L)
    elm_tangstf[3, 9] = -elm_tangstf[3, 3]
    elm_tangstf[9, 9] = elm_tangstf[3, 3]    

    elm_tangstf[4,4] = ((A*P*sinh(2*L*Vz) - 2*A*P*sinh(L*Vz) - 2*Iy*Vz**2*P*sinh(L*Vz) + Iy*Vz**2*P*sinh(2*L*Vz)
                         - A*L**3*Vz**3*P + Iy*L**3*Vz**5*P - Iy*L*Vz**3*P - A*E*Iy*L*Vz**3 - 2*A*E*Iy*Vz**2*sinh(L*Vz)
                         + A*E*Iy*Vz**2*sinh(2*L*Vz) + 2*Iy*L*Vz**3*P*cosh(L*Vz) - Iy*L*Vz**3*P*cosh(2*L*Vz)
                         + A*E*Iy*L**3*Vz**5 - 5*A*L**3*Vz**3*Oz*P + 2*A*L**5*Vz**5*Oz*P + 3*Iy*L**3*Vz**5*Oz*P
                         + 2*A*L**2*Vz**2*P*sinh(L*Vz) + (A*L**2*Vz**2*P*sinh(2*L*Vz))/2 - 2*Iy*L**2*Vz**4*P*sinh(L*Vz)
                         + (Iy*L**2*Vz**4*P*sinh(2*L*Vz))/2 + 2*A*L*Vz*P*cosh(L*Vz) - 2*A*L*Vz*P*cosh(2*L*Vz)
                         + 12*A*L**5*Vz**5*Oz**2*P - A*L**7*Vz**7*Oz**2*P - 9*A*L**7*Vz**7*Oz**3*P + 2*A*L**9*Vz**9*Oz**4*P
                         - 2*Iy*L**5*Vz**7*Oz**2*P + 8*A*L**2*Vz**2*Oz*P*sinh(L*Vz) - 4*A*L**2*Vz**2*Oz*P*sinh(2*L*Vz)
                         - 6*A*L**4*Vz**4*Oz*P*sinh(L*Vz) - A*L**4*Vz**4*Oz*P*sinh(2*L*Vz) + 4*Iy*L**2*Vz**4*Oz*P*sinh(L*Vz)
                         - 2*Iy*L**2*Vz**4*Oz*P*sinh(2*L*Vz) + 2*Iy*L**4*Vz**6*Oz*P*sinh(L*Vz) - 7*A*E*Iy*L**5*Vz**7*Oz**2
                         + 2*A*E*Iy*L**7*Vz**9*Oz**3 - 8*A*L**5*Vz**5*Oz**2*P*cosh(L*Vz) - 4*A*L**5*Vz**5*Oz**2*P*cosh(2*L*Vz)
                         + 8*A*L**7*Vz**7*Oz**3*P*cosh(L*Vz) + A*L**7*Vz**7*Oz**3*P*cosh(2*L*Vz) - 2*A*L**9*Vz**9*Oz**4*P*cosh(L*Vz)
                         + 2*Iy*L**5*Vz**7*Oz**2*P*cosh(L*Vz) - 12*A*L**4*Vz**4*Oz**2*P*sinh(L*Vz) + 6*A*L**4*Vz**4*Oz**2*P*sinh(2*L*Vz)
                         + 6*A*L**6*Vz**6*Oz**2*P*sinh(L*Vz) + 8*A*L**6*Vz**6*Oz**3*P*sinh(L*Vz) + (A*L**6*Vz**6*Oz**2*P*sinh(2*L*Vz))/2
                         - 4*A*L**6*Vz**6*Oz**3*P*sinh(2*L*Vz) - 2*A*L**8*Vz**8*Oz**3*P*sinh(L*Vz) - 2*A*L**8*Vz**8*Oz**4*P*sinh(L*Vz)
                         + A*L**8*Vz**8*Oz**4*P*sinh(2*L*Vz) - 2*Iy*L**4*Vz**6*Oz**2*P*sinh(L*Vz) + Iy*L**4*Vz**6*Oz**2*P*sinh(2*L*Vz)
                         + 2*A*E*Iy*L*Vz**3*cosh(L*Vz) - A*E*Iy*L*Vz**3*cosh(2*L*Vz) + 6*A*E*Iy*L**3*Vz**5*Oz - A*E*Iy*L**5*Vz**7*Oz
                         - 2*A*E*Iy*L**2*Vz**4*sinh(L*Vz) + (A*E*Iy*L**2*Vz**4*sinh(2*L*Vz))/2 + 5*A*L**3*Vz**3*Oz*P*cosh(2*L*Vz)
                         - 4*Iy*L**3*Vz**5*Oz*P*cosh(L*Vz) + Iy*L**3*Vz**5*Oz*P*cosh(2*L*Vz) - 6*A*E*Iy*L**3*Vz**5*Oz*cosh(L*Vz)
                         + 2*A*E*Iy*L**2*Vz**4*Oz*sinh(L*Vz) - A*E*Iy*L**2*Vz**4*Oz*sinh(2*L*Vz) + 4*A*E*Iy*L**4*Vz**6*Oz*sinh(L*Vz)
                         + (A*E*Iy*L**4*Vz**6*Oz*sinh(2*L*Vz))/2 + 6*A*E*Iy*L**5*Vz**7*Oz**2*cosh(L*Vz) + A*E*Iy*L**5*Vz**7*Oz**2*cosh(2*L*Vz)
                         - 2*A*E*Iy*L**7*Vz**9*Oz**3*cosh(L*Vz) + 2*A*E*Iy*L**4*Vz**6*Oz**2*sinh(L*Vz) - A*E*Iy*L**4*Vz**6*Oz**2*sinh(2*L*Vz)
                         - 2*A*E*Iy*L**6*Vz**8*Oz**2*sinh(L*Vz) - 2*A*E*Iy*L**6*Vz**8*Oz**3*sinh(L*Vz) + A*E*Iy*L**6*Vz**8*Oz**3*sinh(2*L*Vz))
                        /(2*A*Vz*(cosh(L*Vz) - 1)*(4*cosh(L*Vz) + L**2*Vz**2 + 8*L**2*Vz**2*Oz + L**2*Vz**2*cosh(L*Vz)
                                                   - 4*L*Vz*sinh(L*Vz) - 4*L**4*Vz**4*Oz**2 - 8*L**2*Vz**2*Oz*cosh(L*Vz)
                                                   + 4*L**3*Vz**3*Oz*sinh(L*Vz) + 4*L**4*Vz**4*Oz**2*cosh(L*Vz) - 4)))
    elm_tangstf[10, 10] = elm_tangstf[4, 4]
    #
    elm_tangstf[5,5] = ((A*P*sinh(2*L*Vy) - 2*A*P*sinh(L*Vy) - 2*Iz*Vy**2*P*sinh(L*Vy) + Iz*Vy**2*P*sinh(2*L*Vy) - A*L**3*Vy**3*P
                         + Iz*L**3*Vy**5*P - Iz*L*Vy**3*P - A*E*Iz*L*Vy**3 - 2*A*E*Iz*Vy**2*sinh(L*Vy) + A*E*Iz*Vy**2*sinh(2*L*Vy)
                         + 2*Iz*L*Vy**3*P*cosh(L*Vy) - Iz*L*Vy**3*P*cosh(2*L*Vy) + A*E*Iz*L**3*Vy**5 - 5*A*L**3*Vy**3*Oy*P
                         + 2*A*L**5*Vy**5*Oy*P + 3*Iz*L**3*Vy**5*Oy*P + 2*A*L**2*Vy**2*P*sinh(L*Vy) + (A*L**2*Vy**2*P*sinh(2*L*Vy))/2
                         - 2*Iz*L**2*Vy**4*P*sinh(L*Vy) + (Iz*L**2*Vy**4*P*sinh(2*L*Vy))/2 + 2*A*L*Vy*P*cosh(L*Vy) - 2*A*L*Vy*P*cosh(2*L*Vy)
                         + 12*A*L**5*Vy**5*Oy**2*P - A*L**7*Vy**7*Oy**2*P - 9*A*L**7*Vy**7*Oy**3*P + 2*A*L**9*Vy**9*Oy**4*P - 2*Iz*L**5*Vy**7*Oy**2*P
                         + 8*A*L**2*Vy**2*Oy*P*sinh(L*Vy) - 4*A*L**2*Vy**2*Oy*P*sinh(2*L*Vy) - 6*A*L**4*Vy**4*Oy*P*sinh(L*Vy) - A*L**4*Vy**4*Oy*P*sinh(2*L*Vy)
                         + 4*Iz*L**2*Vy**4*Oy*P*sinh(L*Vy) - 2*Iz*L**2*Vy**4*Oy*P*sinh(2*L*Vy) + 2*Iz*L**4*Vy**6*Oy*P*sinh(L*Vy) - 7*A*E*Iz*L**5*Vy**7*Oy**2
                         + 2*A*E*Iz*L**7*Vy**9*Oy**3 - 8*A*L**5*Vy**5*Oy**2*P*cosh(L*Vy) - 4*A*L**5*Vy**5*Oy**2*P*cosh(2*L*Vy)
                         + 8*A*L**7*Vy**7*Oy**3*P*cosh(L*Vy) + A*L**7*Vy**7*Oy**3*P*cosh(2*L*Vy) - 2*A*L**9*Vy**9*Oy**4*P*cosh(L*Vy)
                         + 2*Iz*L**5*Vy**7*Oy**2*P*cosh(L*Vy) - 12*A*L**4*Vy**4*Oy**2*P*sinh(L*Vy) + 6*A*L**4*Vy**4*Oy**2*P*sinh(2*L*Vy)
                         + 6*A*L**6*Vy**6*Oy**2*P*sinh(L*Vy) + 8*A*L**6*Vy**6*Oy**3*P*sinh(L*Vy) + (A*L**6*Vy**6*Oy**2*P*sinh(2*L*Vy))/2
                         - 4*A*L**6*Vy**6*Oy**3*P*sinh(2*L*Vy) - 2*A*L**8*Vy**8*Oy**3*P*sinh(L*Vy) - 2*A*L**8*Vy**8*Oy**4*P*sinh(L*Vy)
                         + A*L**8*Vy**8*Oy**4*P*sinh(2*L*Vy) - 2*Iz*L**4*Vy**6*Oy**2*P*sinh(L*Vy) + Iz*L**4*Vy**6*Oy**2*P*sinh(2*L*Vy)
                         + 2*A*E*Iz*L*Vy**3*cosh(L*Vy) - A*E*Iz*L*Vy**3*cosh(2*L*Vy) + 6*A*E*Iz*L**3*Vy**5*Oy - A*E*Iz*L**5*Vy**7*Oy
                         - 2*A*E*Iz*L**2*Vy**4*sinh(L*Vy) + (A*E*Iz*L**2*Vy**4*sinh(2*L*Vy))/2 + 5*A*L**3*Vy**3*Oy*P*cosh(2*L*Vy)
                         - 4*Iz*L**3*Vy**5*Oy*P*cosh(L*Vy) + Iz*L**3*Vy**5*Oy*P*cosh(2*L*Vy) - 6*A*E*Iz*L**3*Vy**5*Oy*cosh(L*Vy)
                         + 2*A*E*Iz*L**2*Vy**4*Oy*sinh(L*Vy) - A*E*Iz*L**2*Vy**4*Oy*sinh(2*L*Vy) + 4*A*E*Iz*L**4*Vy**6*Oy*sinh(L*Vy)
                         + (A*E*Iz*L**4*Vy**6*Oy*sinh(2*L*Vy))/2 + 6*A*E*Iz*L**5*Vy**7*Oy**2*cosh(L*Vy) + A*E*Iz*L**5*Vy**7*Oy**2*cosh(2*L*Vy)
                         - 2*A*E*Iz*L**7*Vy**9*Oy**3*cosh(L*Vy) + 2*A*E*Iz*L**4*Vy**6*Oy**2*sinh(L*Vy) - A*E*Iz*L**4*Vy**6*Oy**2*sinh(2*L*Vy)
                         - 2*A*E*Iz*L**6*Vy**8*Oy**2*sinh(L*Vy) - 2*A*E*Iz*L**6*Vy**8*Oy**3*sinh(L*Vy) + A*E*Iz*L**6*Vy**8*Oy**3*sinh(2*L*Vy))
                        /(2*A*Vy*(cosh(L*Vy) - 1)*(4*cosh(L*Vy) + L**2*Vy**2 + 8*L**2*Vy**2*Oy + L**2*Vy**2*cosh(L*Vy)
                                                   - 4*L*Vy*sinh(L*Vy) - 4*L**4*Vy**4*Oy**2 - 8*L**2*Vy**2*Oy*cosh(L*Vy)
                                                   + 4*L**3*Vy**3*Oy*sinh(L*Vy) + 4*L**4*Vy**4*Oy**2*cosh(L*Vy) - 4)))
    #
    elm_tangstf[11, 11] = elm_tangstf[5, 5]
    #
    
    #
    elm_tangstf[2,3]  = Mza/L
    elm_tangstf[3,2] = elm_tangstf[2,3]
    #
    #
    eq2_num = ((Oz*L**2*Vz**2 - 1)*(Mza + Mzb)*(L*Vz - 2*sinh(L*Vz) - 2*cosh(L*Vz) + L*Vz*(cosh(L*Vz) + sinh(L*Vz)) + 2))
    eq2_den = (L**2*Vz**2*(cosh(L*Vz) + sinh(L*Vz) - 1))    
    elm_tangstf[3,4]  = -Mza/2 - eq2_num/eq2_den
    elm_tangstf[4,3] = elm_tangstf[3,4]
    #
    #
    elm_tangstf[9,10] = -Mzb/2 - eq2_num/eq2_den
    elm_tangstf[10,9] = elm_tangstf[9,10]
    #
    #
    #
    eq11_num = ((Oy*L**2*Vy**2 - 1)*(Mya + Myb)*(L*Vy - 2*sinh(L*Vy) - 2*cosh(L*Vy) + L*Vy*(cosh(L*Vy) + sinh(L*Vy)) + 2))
    eq11_den = (L**2*Vy**2*(cosh(L*Vy) + sinh(L*Vy) - 1))
    #
    elm_tangstf[3,5]  = Mya/2 + eq11_num/eq11_den
    #elm_tangstf[5,3] =  elm_tangstf[3,5]    
    #
    #
    elm_tangstf[4,6]  = Mya/L
    #elm_tangstf[6,4] = elm_tangstf[4,6]
    #
    #
    eq10_num = (A*L**2*Vy**2*P - 8*A*P*sinh((L*Vy)/2)**2 - Iz*L**2*Vy**4*P + Iz*L*Vy**3*P*sinh(L*Vy) - A*E*Iz*L**2*Vy**4 - 2*A*L**4*Vy**4*Oy*P
                + A*L*Vy*P*sinh(L*Vy) + A*L**6*Vy**6*Oy**2*P - 2*A*L**3*Vy**3*Oy*P*sinh(L*Vy) + 16*A*L**2*Vy**2*Oy*P*sinh((L*Vy)/2)**2
                + A*L**5*Vy**5*Oy**2*P*sinh(L*Vy) + A*E*Iz*L*Vy**3*sinh(L*Vy) - 8*A*L**4*Vy**4*Oy**2*P*sinh((L*Vy)/2)**2
                + A*E*Iz*L**4*Vy**6*Oy + A*E*Iz*L**3*Vy**5*Oy*sinh(L*Vy))
    eq10_den = (4*A*(L*Vy - 2*sinh((L*Vy)/2) + 2*L*Vy*sinh((L*Vy)/4)**2 + 2*L**2*Vy**2*Oy*sinh((L*Vy)/2))**2)
    #
    elm_tangstf[1,5]  = eq10_num/eq10_den
    elm_tangstf[1,11] = elm_tangstf[1, 5]
    #elm_tangstf[5,7]  = -eq10_num/eq10_den
    elm_tangstf[5,7] = -elm_tangstf[1,5]
    #elm_tangstf[7,5] = elm_tangstf[5,7]
    #
    #
    elm_tangstf[7,9]  = -Myb/L
    #elm_tangstf[9,7] = elm_tangstf[7,9]
    #
    #
    elm_tangstf[8,10] = -elm_tangstf[2, 4]
    #elm_tangstf[8,10] = -elm_tangstf[2,4] #eq3_num/eq3_den
    #elm_tangstf[10,8] = elm_tangstf[8,10]    
    #
    #
    elm_tangstf[9,11] = Myb/2 + eq11_num/eq11_den
    #elm_tangstf[11,9] = elm_tangstf[9,11]
    #    
    #
    eq3_num = (A*L**2*Vz**2*P - 8*A*P*sinh((L*Vz)/2)**2 - Iy*L**2*Vz**4*P + Iy*L*Vz**3*P*sinh(L*Vz) - A*E*Iy*L**2*Vz**4 - 2*A*L**4*Vz**4*Oz*P
               + A*L*Vz*P*sinh(L*Vz) + A*L**6*Vz**6*Oz**2*P - 2*A*L**3*Vz**3*Oz*P*sinh(L*Vz) + 16*A*L**2*Vz**2*Oz*P*sinh((L*Vz)/2)**2
               + A*L**5*Vz**5*Oz**2*P*sinh(L*Vz) + A*E*Iy*L*Vz**3*sinh(L*Vz) - 8*A*L**4*Vz**4*Oz**2*P*sinh((L*Vz)/2)**2
               + A*E*Iy*L**4*Vz**6*Oz + A*E*Iy*L**3*Vz**5*Oz*sinh(L*Vz))

    eq3_den = (4*A*(L*Vz - 2*sinh((L*Vz)/2) + 2*L*Vz*sinh((L*Vz)/4)**2 + 2*L**2*Vz**2*Oz*sinh((L*Vz)/2))**2)
    #
    elm_tangstf[2,4]  = -eq3_num/eq3_den
    #elm_tangstf[2,10] = -((A*L**2*Vz**2*P - 8*A*P*sinh((L*Vz)/2)**2 - Iy*L**2*Vz**4*P + Iy*L*Vz**3*P*sinh(L*Vz)
    #                       - A*E*Iy*L**2*Vz**4 - 2*A*L**4*Vz**4*Oz*P + A*L*Vz*P*sinh(L*Vz) + A*L**6*Vz**6*Oz**2*P
    #                       - 2*A*L**3*Vz**3*Oz*P*sinh(L*Vz) + 16*A*L**2*Vz**2*Oz*P*sinh((L*Vz)/2)**2 + A*L**5*Vz**5*Oz**2*P*sinh(L*Vz)
    #                       + A*E*Iy*L*Vz**3*sinh(L*Vz) - 8*A*L**4*Vz**4*Oz**2*P*sinh((L*Vz)/2)**2 + A*E*Iy*L**4*Vz**6*Oz
    #                       + A*E*Iy*L**3*Vz**5*Oz*sinh(L*Vz))
    #                      /(4*A*(L*Vz - 2*sinh((L*Vz)/2) + 2*L*Vz*sinh((L*Vz)/4)**2 + 2*L**2*Vz**2*Oz*sinh((L*Vz)/2))**2))
    #
    elm_tangstf[4,2] = elm_tangstf[2,4]
    elm_tangstf[2,10] = elm_tangstf[2, 4]  
    #
    #        
    #
    #
    elm_tangstf[0,4]  = -Mya/L
    elm_tangstf[3,7]  = elm_tangstf[0,4] #-Mya/L
    #
    #
    elm_tangstf[4,8]  = -elm_tangstf[2,4] #eq3_num/eq3_den  
    #
    # ??   
    elm_tangstf[5,9]  = (-((Oy*L**2*Vy**2 - 1)*((L*Vy*cosh((L*Vy)/2))/sinh((L*Vy)/2) - 2)*(Mya + Myb))
                         /(L**2*Vy**2))
    #elm_tangstf[3,11] = (-((Oy*L**2*Vy**2 - 1)*((L*Vy*cosh((L*Vy)/2))/sinh((L*Vy)/2) - 2)*(Mya + Myb))
    #                     /(L**2*Vy**2))
    elm_tangstf[3,11] = elm_tangstf[5,9]
    #
    #
    elm_tangstf[6,10] = Myb/L
    elm_tangstf[1,9]  = elm_tangstf[6,10]
    elm_tangstf[0,10] = -elm_tangstf[6,10] #-Myb/L 
    #
    #
    elm_tangstf[7,11] = -elm_tangstf[1, 5]
    #elm_tangstf[7,11] = -eq10_num/eq10_den     
    #
    #
    elm_tangstf[0,5]  = -Mza/L
    elm_tangstf[3,8]  = elm_tangstf[0,5] #-Mza/L 
    elm_tangstf[5,6]  = -elm_tangstf[0,5] # Mza/L
    #
    # ?
    elm_tangstf[4,9]  = (((Oz*L**2*Vz**2 - 1)*((L*Vz*cosh((L*Vz)/2))/sinh((L*Vz)/2) - 2)*(Mza + Mzb))
                         /(L**2*Vz**2))
    #elm_tangstf[3,10] = ((Oz*L**2*Vz**2 - 1)*((L*Vz*cosh((L*Vz)/2))/sinh((L*Vz)/2) - 2)*(Mza + Mzb))/(L**2*Vz**2) 
    elm_tangstf[3,10] =  elm_tangstf[4,9]
    #
    #
    elm_tangstf[2,9]  = Mzb/L
    elm_tangstf[6,11] = Mzb/L
    elm_tangstf[8,9]  = -Mzb/L
    elm_tangstf[0,11] = -Mzb/L
    #
    #
    elm_tangstf[4,10] = ((2*A*P*sinh(L*Vz) - A*P*sinh(2*L*Vz) + 2*Iy*Vz**2*P*sinh(L*Vz) - Iy*Vz**2*P*sinh(2*L*Vz) - 6*A*L*Vz*P
                          + Iy*L*Vz**3*P + A*E*Iy*L*Vz**3 + 2*A*E*Iy*Vz**2*sinh(L*Vz) - A*E*Iy*Vz**2*sinh(2*L*Vz) - 2*Iy*L*Vz**3*P*cosh(L*Vz)
                          + Iy*L*Vz**3*P*cosh(2*L*Vz) + 17*A*L**3*Vz**3*Oz*P - 3*Iy*L**3*Vz**5*Oz*P + A*L**3*Vz**3*P*cosh(L*Vz)
                          - Iy*L**3*Vz**5*P*cosh(L*Vz) - 3*A*L**2*Vz**2*P*sinh(L*Vz) + Iy*L**2*Vz**4*P*sinh(L*Vz) + 6*A*L*Vz*P*cosh(L*Vz)
                          - 18*A*L**5*Vz**5*Oz**2*P + 9*A*L**7*Vz**7*Oz**3*P - 2*A*L**9*Vz**9*Oz**4*P + 2*Iy*L**5*Vz**7*Oz**2*P
                          - 8*A*L**2*Vz**2*Oz*P*sinh(L*Vz) + 4*A*L**2*Vz**2*Oz*P*sinh(2*L*Vz) + 8*A*L**4*Vz**4*Oz*P*sinh(L*Vz)
                          - 4*Iy*L**2*Vz**4*Oz*P*sinh(L*Vz) + 2*Iy*L**2*Vz**4*Oz*P*sinh(2*L*Vz) - 2*Iy*L**4*Vz**6*Oz*P*sinh(L*Vz)
                          + 7*A*E*Iy*L**5*Vz**7*Oz**2 - 2*A*E*Iy*L**7*Vz**9*Oz**3 + 16*A*L**5*Vz**5*Oz**2*P*cosh(L*Vz)
                          + 2*A*L**5*Vz**5*Oz**2*P*cosh(2*L*Vz) + A*L**7*Vz**7*Oz**2*P*cosh(L*Vz) - 8*A*L**7*Vz**7*Oz**3*P*cosh(L*Vz)
                          - A*L**7*Vz**7*Oz**3*P*cosh(2*L*Vz) + 2*A*L**9*Vz**9*Oz**4*P*cosh(L*Vz) - 2*Iy*L**5*Vz**7*Oz**2*P*cosh(L*Vz)
                          + 12*A*L**4*Vz**4*Oz**2*P*sinh(L*Vz) - 6*A*L**4*Vz**4*Oz**2*P*sinh(2*L*Vz) - 7*A*L**6*Vz**6*Oz**2*P*sinh(L*Vz)
                          - 8*A*L**6*Vz**6*Oz**3*P*sinh(L*Vz) + 4*A*L**6*Vz**6*Oz**3*P*sinh(2*L*Vz) + 2*A*L**8*Vz**8*Oz**3*P*sinh(L*Vz)
                          + 2*A*L**8*Vz**8*Oz**4*P*sinh(L*Vz) - A*L**8*Vz**8*Oz**4*P*sinh(2*L*Vz) + 2*Iy*L**4*Vz**6*Oz**2*P*sinh(L*Vz)
                          - Iy*L**4*Vz**6*Oz**2*P*sinh(2*L*Vz) - 2*A*E*Iy*L*Vz**3*cosh(L*Vz) + A*E*Iy*L*Vz**3*cosh(2*L*Vz) - 6*A*E*Iy*L**3*Vz**5*Oz
                          - A*E*Iy*L**3*Vz**5*cosh(L*Vz) + A*E*Iy*L**2*Vz**4*sinh(L*Vz) - 16*A*L**3*Vz**3*Oz*P*cosh(L*Vz) - A*L**3*Vz**3*Oz*P*cosh(2*L*Vz)
                          - 2*A*L**5*Vz**5*Oz*P*cosh(L*Vz) + 4*Iy*L**3*Vz**5*Oz*P*cosh(L*Vz) - Iy*L**3*Vz**5*Oz*P*cosh(2*L*Vz)
                          + 6*A*E*Iy*L**3*Vz**5*Oz*cosh(L*Vz) + A*E*Iy*L**5*Vz**7*Oz*cosh(L*Vz) - 2*A*E*Iy*L**2*Vz**4*Oz*sinh(L*Vz)
                          + A*E*Iy*L**2*Vz**4*Oz*sinh(2*L*Vz) - 5*A*E*Iy*L**4*Vz**6*Oz*sinh(L*Vz) - 6*A*E*Iy*L**5*Vz**7*Oz**2*cosh(L*Vz)
                          - A*E*Iy*L**5*Vz**7*Oz**2*cosh(2*L*Vz) + 2*A*E*Iy*L**7*Vz**9*Oz**3*cosh(L*Vz) - 2*A*E*Iy*L**4*Vz**6*Oz**2*sinh(L*Vz)
                          + A*E*Iy*L**4*Vz**6*Oz**2*sinh(2*L*Vz) + 2*A*E*Iy*L**6*Vz**8*Oz**2*sinh(L*Vz) + 2*A*E*Iy*L**6*Vz**8*Oz**3*sinh(L*Vz)
                          - A*E*Iy*L**6*Vz**8*Oz**3*sinh(2*L*Vz))
                         /(A*Vz*(4*cosh(2*L*Vz) - 16*cosh(L*Vz) - L**2*Vz**2 - 24*L**2*Vz**2*Oz + L**2*Vz**2*cosh(2*L*Vz) + 8*L*Vz*sinh(L*Vz)
                                 - 4*L*Vz*sinh(2*L*Vz) + 12*L**4*Vz**4*Oz**2 + 32*L**2*Vz**2*Oz*cosh(L*Vz) - 8*L**2*Vz**2*Oz*cosh(2*L*Vz)
                                 - 8*L**3*Vz**3*Oz*sinh(L*Vz) + 4*L**3*Vz**3*Oz*sinh(2*L*Vz) - 16*L**4*Vz**4*Oz**2*cosh(L*Vz)
                                 + 4*L**4*Vz**4*Oz**2*cosh(2*L*Vz) + 12))) 
    #
    #
    elm_tangstf[5,11] = ((2*A*P*sinh(L*Vy) - A*P*sinh(2*L*Vy) + 2*Iz*Vy**2*P*sinh(L*Vy) - Iz*Vy**2*P*sinh(2*L*Vy) - 6*A*L*Vy*P
                          + Iz*L*Vy**3*P + A*E*Iz*L*Vy**3 + 2*A*E*Iz*Vy**2*sinh(L*Vy) - A*E*Iz*Vy**2*sinh(2*L*Vy) - 2*Iz*L*Vy**3*P*cosh(L*Vy)
                          + Iz*L*Vy**3*P*cosh(2*L*Vy) + 17*A*L**3*Vy**3*Oy*P - 3*Iz*L**3*Vy**5*Oy*P + A*L**3*Vy**3*P*cosh(L*Vy)
                          - Iz*L**3*Vy**5*P*cosh(L*Vy) - 3*A*L**2*Vy**2*P*sinh(L*Vy) + Iz*L**2*Vy**4*P*sinh(L*Vy) + 6*A*L*Vy*P*cosh(L*Vy)
                          - 18*A*L**5*Vy**5*Oy**2*P + 9*A*L**7*Vy**7*Oy**3*P - 2*A*L**9*Vy**9*Oy**4*P + 2*Iz*L**5*Vy**7*Oy**2*P
                          - 8*A*L**2*Vy**2*Oy*P*sinh(L*Vy) + 4*A*L**2*Vy**2*Oy*P*sinh(2*L*Vy) + 8*A*L**4*Vy**4*Oy*P*sinh(L*Vy)
                          - 4*Iz*L**2*Vy**4*Oy*P*sinh(L*Vy) + 2*Iz*L**2*Vy**4*Oy*P*sinh(2*L*Vy) - 2*Iz*L**4*Vy**6*Oy*P*sinh(L*Vy)
                          + 7*A*E*Iz*L**5*Vy**7*Oy**2 - 2*A*E*Iz*L**7*Vy**9*Oy**3 + 16*A*L**5*Vy**5*Oy**2*P*cosh(L*Vy)
                          + 2*A*L**5*Vy**5*Oy**2*P*cosh(2*L*Vy) + A*L**7*Vy**7*Oy**2*P*cosh(L*Vy) - 8*A*L**7*Vy**7*Oy**3*P*cosh(L*Vy)
                          - A*L**7*Vy**7*Oy**3*P*cosh(2*L*Vy) + 2*A*L**9*Vy**9*Oy**4*P*cosh(L*Vy) - 2*Iz*L**5*Vy**7*Oy**2*P*cosh(L*Vy)
                          + 12*A*L**4*Vy**4*Oy**2*P*sinh(L*Vy) - 6*A*L**4*Vy**4*Oy**2*P*sinh(2*L*Vy) - 7*A*L**6*Vy**6*Oy**2*P*sinh(L*Vy)
                          - 8*A*L**6*Vy**6*Oy**3*P*sinh(L*Vy) + 4*A*L**6*Vy**6*Oy**3*P*sinh(2*L*Vy) + 2*A*L**8*Vy**8*Oy**3*P*sinh(L*Vy)
                          + 2*A*L**8*Vy**8*Oy**4*P*sinh(L*Vy) - A*L**8*Vy**8*Oy**4*P*sinh(2*L*Vy) + 2*Iz*L**4*Vy**6*Oy**2*P*sinh(L*Vy)
                          - Iz*L**4*Vy**6*Oy**2*P*sinh(2*L*Vy) - 2*A*E*Iz*L*Vy**3*cosh(L*Vy) + A*E*Iz*L*Vy**3*cosh(2*L*Vy) - 6*A*E*Iz*L**3*Vy**5*Oy
                          - A*E*Iz*L**3*Vy**5*cosh(L*Vy) + A*E*Iz*L**2*Vy**4*sinh(L*Vy) - 16*A*L**3*Vy**3*Oy*P*cosh(L*Vy) - A*L**3*Vy**3*Oy*P*cosh(2*L*Vy)
                          - 2*A*L**5*Vy**5*Oy*P*cosh(L*Vy) + 4*Iz*L**3*Vy**5*Oy*P*cosh(L*Vy) - Iz*L**3*Vy**5*Oy*P*cosh(2*L*Vy)
                          + 6*A*E*Iz*L**3*Vy**5*Oy*cosh(L*Vy) + A*E*Iz*L**5*Vy**7*Oy*cosh(L*Vy) - 2*A*E*Iz*L**2*Vy**4*Oy*sinh(L*Vy)
                          + A*E*Iz*L**2*Vy**4*Oy*sinh(2*L*Vy) - 5*A*E*Iz*L**4*Vy**6*Oy*sinh(L*Vy) - 6*A*E*Iz*L**5*Vy**7*Oy**2*cosh(L*Vy)
                          - A*E*Iz*L**5*Vy**7*Oy**2*cosh(2*L*Vy) + 2*A*E*Iz*L**7*Vy**9*Oy**3*cosh(L*Vy) - 2*A*E*Iz*L**4*Vy**6*Oy**2*sinh(L*Vy)
                          + A*E*Iz*L**4*Vy**6*Oy**2*sinh(2*L*Vy) + 2*A*E*Iz*L**6*Vy**8*Oy**2*sinh(L*Vy) + 2*A*E*Iz*L**6*Vy**8*Oy**3*sinh(L*Vy)
                          - A*E*Iz*L**6*Vy**8*Oy**3*sinh(2*L*Vy))
                         /(A*Vy*(4*cosh(2*L*Vy) - 16*cosh(L*Vy) - L**2*Vy**2 - 24*L**2*Vy**2*Oy + L**2*Vy**2*cosh(2*L*Vy) + 8*L*Vy*sinh(L*Vy)
                                 - 4*L*Vy*sinh(2*L*Vy) + 12*L**4*Vy**4*Oy**2 + 32*L**2*Vy**2*Oy*cosh(L*Vy) - 8*L**2*Vy**2*Oy*cosh(2*L*Vy)
                                 - 8*L**3*Vy**3*Oy*sinh(L*Vy) + 4*L**3*Vy**3*Oy*sinh(2*L*Vy) - 16*L**4*Vy**4*Oy**2*cosh(L*Vy)
                                 + 4*L**4*Vy**4*Oy**2*cosh(2*L*Vy) + 12)))
    #
    #    
    #
    elm_tangstf[1,3]  = Mya/L
    #
    #    
    #
    #
    #elm_tangstf[1,11] = ((A*L**2*Vy**2*P - 8*A*P*sinh((L*Vy)/2)**2 - Iz*L**2*Vy**4*P + Iz*L*Vy**3*P*sinh(L*Vy)
    #                      - A*E*Iz*L**2*Vy**4 - 2*A*L**4*Vy**4*Oy*P + A*L*Vy*P*sinh(L*Vy) + A*L**6*Vy**6*Oy**2*P
    #                      - 2*A*L**3*Vy**3*Oy*P*sinh(L*Vy) + 16*A*L**2*Vy**2*Oy*P*sinh((L*Vy)/2)**2 + A*L**5*Vy**5*Oy**2*P*sinh(L*Vy)
    #                      + A*E*Iz*L*Vy**3*sinh(L*Vy) - 8*A*L**4*Vy**4*Oy**2*P*sinh((L*Vy)/2)**2 + A*E*Iz*L**4*Vy**6*Oy
    #                      + A*E*Iz*L**3*Vy**5*Oy*sinh(L*Vy))
    #                     /(4*A*(L*Vy - 2*sinh((L*Vy)/2) + 2*L*Vy*sinh((L*Vy)/4)**2 + 2*L**2*Vy**2*Oy*sinh((L*Vy)/2))**2))
    #
    #
    # ---------------- Torsion for symetric sections --------------------
    if Iy == Iz:
        eq4_num = (Vy*Mxb*(sinh(L*Vy) - L*Vy))
        eq4_den = (2*(L*Vy*sinh(L*Vy) - 2*L**2*Vy**2*Oy - 2*cosh(L*Vy) + 2*L**2*Vy**2*Oy*cosh(L*Vy) + 2))
        eq5_den = (2*(L*Vy*sinh(L*Vy) - 2*L**2*Vy**2*Oz - 2*cosh(L*Vy) + 2*L**2*Vy**2*Oz*cosh(L*Vy) + 2))
        #
        elm_tangstf[1,4]  = eq4_num/eq4_den
        elm_tangstf[2,5]  = eq4_num/eq5_den
        elm_tangstf[4,7]  = -eq4_num/eq4_den
        elm_tangstf[5,8]  = -eq4_num/eq5_den
        elm_tangstf[7,10] = eq4_num/eq4_den
        elm_tangstf[8,11] = eq4_num/eq5_den
        #
        elm_tangstf[1,10] = -eq4_num/eq4_den
        elm_tangstf[2,11] = -eq4_num/eq5_den        
        #
        # ---------------------------------------
        #
        eq6_num = (Mxb*(6*Vy - 8*Vy*cosh(L*Vy) + 2*Vy*cosh(2*L*Vy) + 2*L**2*Vy**3 - 6*L**2*Vy**3*Oy - 6*L**2*Vy**3*Oz - L**4*Vy**5*Oy
                        - L**4*Vy**5*Oz - 2*L**2*Vy**3*cosh(L*Vy) + L**3*Vy**4*sinh(L*Vy) + 2*L*Vy**2*sinh(L*Vy)
                        - L*Vy**2*sinh(2*L*Vy) + 6*L**4*Vy**5*Oy*Oz + 8*L**2*Vy**3*Oy*cosh(L*Vy) + 8*L**2*Vy**3*Oz*cosh(L*Vy)
                        - 2*L**2*Vy**3*Oy*cosh(2*L*Vy) - 2*L**2*Vy**3*Oz*cosh(2*L*Vy) + L**4*Vy**5*Oy*cosh(L*Vy) + L**4*Vy**5*Oz*cosh(L*Vy)
                        - L**3*Vy**4*Oy*sinh(L*Vy) - L**3*Vy**4*Oz*sinh(L*Vy) + (L**3*Vy**4*Oy*sinh(2*L*Vy))/2 + (L**3*Vy**4*Oz*sinh(2*L*Vy))/2
                        - 8*L**4*Vy**5*Oy*Oz*cosh(L*Vy) + 2*L**4*Vy**5*Oy*Oz*cosh(2*L*Vy)))

        eq6_den = (Vy*(16*cosh(L*Vy) - 4*cosh(2*L*Vy) + L**2*Vy**2 + 12*L**2*Vy**2*Oy + 12*L**2*Vy**2*Oz - L**2*Vy**2*cosh(2*L*Vy)
                       - 8*L*Vy*sinh(L*Vy) + 4*L*Vy*sinh(2*L*Vy) - 12*L**4*Vy**4*Oy*Oz - 16*L**2*Vy**2*Oy*cosh(L*Vy)
                       - 16*L**2*Vy**2*Oz*cosh(L*Vy) + 4*L**2*Vy**2*Oy*cosh(2*L*Vy) + 4*L**2*Vy**2*Oz*cosh(2*L*Vy)
                       + 4*L**3*Vy**3*Oy*sinh(L*Vy) + 4*L**3*Vy**3*Oz*sinh(L*Vy) - 2*L**3*Vy**3*Oy*sinh(2*L*Vy)
                       - 2*L**3*Vy**3*Oz*sinh(2*L*Vy) + 16*L**4*Vy**4*Oy*Oz*cosh(L*Vy) - 4*L**4*Vy**4*Oy*Oz*cosh(2*L*Vy) - 12))

        elm_tangstf[5,10] = -eq6_num/eq6_den
        elm_tangstf[4,11] = eq6_num/eq6_den

        #
        # -----------------------------------------
        #
        eq7_num = (L**3*Vy**3*Mxb*(Oy - Oz)*(sinh(L*Vy) - L*Vy))
        eq7_den =  (2*(4*cosh(L*Vy) + L**2*Vy**2 + 4*L**2*Vy**2*Oy + 4*L**2*Vy**2*Oz + L**2*Vy**2*cosh(L*Vy)
                       - 4*L*Vy*sinh(L*Vy) - 4*L**4*Vy**4*Oy*Oz - 4*L**2*Vy**2*Oy*cosh(L*Vy) - 4*L**2*Vy**2*Oz*cosh(L*Vy)
                       + 2*L**3*Vy**3*Oy*sinh(L*Vy) + 2*L**3*Vy**3*Oz*sinh(L*Vy) + 4*L**4*Vy**4*Oy*Oz*cosh(L*Vy) - 4))

        elm_tangstf[4,5]   = -eq7_num/eq7_den
        elm_tangstf[10,11] = eq7_num/eq7_den
    else:
        eq4_num = (Vy**2*Mxb*(Vy*cosh((L*Vy)/2)*sinh((L*Vz)/2) - Vz*cosh((L*Vz)/2)*sinh((L*Vy)/2)))
        eq4_den = (sinh((L*Vz)/2)*(Vy**2 - Vz**2)*(L*Vy*cosh((L*Vy)/2) - 2*sinh((L*Vy)/2) + 2*L**2*Vy**2*Oy*sinh((L*Vy)/2)))

        eq5_num = (Vz**2*Mxb*(Vy*cosh((L*Vy)/2)*sinh((L*Vz)/2) - Vz*cosh((L*Vz)/2)*sinh((L*Vy)/2)))
        eq5_den = (sinh((L*Vy)/2)*(Vy**2 - Vz**2)*(L*Vz*cosh((L*Vz)/2) - 2*sinh((L*Vz)/2) + 2*L**2*Vz**2*Oz*sinh((L*Vz)/2)))
        #
        #
        elm_tangstf[1,4]  = eq4_num/eq4_den
        elm_tangstf[2,5]  = eq5_num/eq5_den
        elm_tangstf[4,7]  = -eq4_num/eq4_den
        elm_tangstf[5,8]  = -eq5_num/eq5_den
        elm_tangstf[7,10] = eq4_num/eq4_den
        elm_tangstf[8,11] = eq5_num/eq5_den
        #
        #
        elm_tangstf[1,10] = -eq4_num/eq4_den
        elm_tangstf[2,11] = -eq5_num/eq5_den         
        #---------------------------------
        #
        eq6_num = (2*Mxb*(4*Vy**2 - 4*Vz**2 - 4*Vy**2*cosh((L*Vy)/2)**2 - 4*Vy**2*cosh((L*Vz)/2)**2
                          + 4*Vz**2*cosh((L*Vy)/2)**2 + 4*Vz**2*cosh((L*Vz)/2)**2 - 4*L**2*Vy**4*Oy
                          + 4*L**2*Vz**4*Oz + 4*Vy**2*cosh((L*Vy)/2)**2*cosh((L*Vz)/2)**2 - 4*Vz**2*cosh((L*Vy)/2)**2*cosh((L*Vz)/2)**2
                          - 2*L*Vy*Vz**2*sinh(L*Vy) + 2*L*Vy**2*Vz*sinh(L*Vz) + L**2*Vy**2*Vz**2*cosh((L*Vy)/2)**2
                          - L**2*Vy**2*Vz**2*cosh((L*Vz)/2)**2 + 4*L**2*Vy**2*Vz**2*Oy - 4*L**2*Vy**2*Vz**2*Oz
                          + 4*L**2*Vy**4*Oy*cosh((L*Vy)/2)**2 + 4*L**2*Vy**4*Oy*cosh((L*Vz)/2)**2 - 4*L**2*Vz**4*Oz*cosh((L*Vy)/2)**2
                          - 4*L**2*Vz**4*Oz*cosh((L*Vz)/2)**2 - 4*L**4*Vy**2*Vz**4*Oy*Oz + 4*L**4*Vy**4*Vz**2*Oy*Oz
                          + L**3*Vy**3*Vz**2*Oy*sinh(L*Vy) - L**3*Vy**2*Vz**3*Oz*sinh(L*Vz) - 4*L**2*Vy**2*Vz**2*Oy*cosh((L*Vy)/2)**2
                          - 4*L**2*Vy**2*Vz**2*Oy*cosh((L*Vz)/2)**2 + 4*L**2*Vy**2*Vz**2*Oz*cosh((L*Vy)/2)**2
                          + 4*L**2*Vy**2*Vz**2*Oz*cosh((L*Vz)/2)**2 - 4*L**2*Vy**4*Oy*cosh((L*Vy)/2)**2*cosh((L*Vz)/2)**2
                          + 4*L**2*Vz**4*Oz*cosh((L*Vy)/2)**2*cosh((L*Vz)/2)**2 + L**3*Vy*Vz**4*Oz*sinh(L*Vy)
                          - L**3*Vy**4*Vz*Oy*sinh(L*Vz) + 4*L*Vy*Vz**2*cosh((L*Vy)/2)*cosh((L*Vz)/2)**2*sinh((L*Vy)/2)
                          - 4*L*Vy**2*Vz*cosh((L*Vy)/2)**2*cosh((L*Vz)/2)*sinh((L*Vz)/2) + 4*L**2*Vy**2*Vz**2*Oy*cosh((L*Vy)/2)**2*cosh((L*Vz)/2)**2
                          - 4*L**2*Vy**2*Vz**2*Oz*cosh((L*Vy)/2)**2*cosh((L*Vz)/2)**2 + 4*L**4*Vy**2*Vz**4*Oy*Oz*cosh((L*Vy)/2)**2
                          - 4*L**4*Vy**4*Vz**2*Oy*Oz*cosh((L*Vy)/2)**2 + 4*L**4*Vy**2*Vz**4*Oy*Oz*cosh((L*Vz)/2)**2
                          - 4*L**4*Vy**4*Vz**2*Oy*Oz*cosh((L*Vz)/2)**2 - 2*L**3*Vy**3*Vz**2*Oy*cosh((L*Vy)/2)*cosh((L*Vz)/2)**2*sinh((L*Vy)/2)
                          + 2*L**3*Vy**2*Vz**3*Oz*cosh((L*Vy)/2)**2*cosh((L*Vz)/2)*sinh((L*Vz)/2) - 4*L**4*Vy**2*Vz**4*Oy*Oz*cosh((L*Vy)/2)**2*cosh((L*Vz)/2)**2
                          + 4*L**4*Vy**4*Vz**2*Oy*Oz*cosh((L*Vy)/2)**2*cosh((L*Vz)/2)**2 - 2*L**3*Vy*Vz**4*Oz*cosh((L*Vy)/2)*cosh((L*Vz)/2)**2*sinh((L*Vy)/2)
                          + 2*L**3*Vy**4*Vz*Oy*cosh((L*Vy)/2)**2*cosh((L*Vz)/2)*sinh((L*Vz)/2)))
        eq6_den = ((Vy**2 - Vz**2)*(L*Vy*sinh(L*Vy) - 2*L**2*Vy**2*Oy - 2*cosh(L*Vy)
                                    + 2*L**2*Vy**2*Oy*cosh(L*Vy) + 2)*(L*Vz*sinh(L*Vz) - 2*L**2*Vz**2*Oz - 2*cosh(L*Vz) + 2*L**2*Vz**2*Oz*cosh(L*Vz) + 2))

        elm_tangstf[5,10] = eq6_num/eq6_den
        elm_tangstf[4,11] = -eq6_num/eq6_den
        #
        # --------------------------------
        #
        eq7_num = (L**2*(2*Vy**2*Vz**2*Mxb*(sinh((L*Vy)/2)**2 - 2*(sinh((L*Vy)/2)**2 + 1)*(sinh((L*Vz)/2)**2 + 1)
                                            + sinh((L*Vz)/2)**2 + 2) + (Vy*Vz**3*Mxb*sinh(L*Vy)*sinh(L*Vz))/2
                         + (Vy**3*Vz*Mxb*sinh(L*Vy)*sinh(L*Vz))/2) + L*(2*Vy**3*Mxb*(sinh(L*Vy) - 2*sinh((L*Vy)/2)*(sinh((L*Vz)/2)**2 + 1)*(2*sinh((L*Vy)/4)**2 + 1))
                                                                        + 2*Vz**3*Mxb*(sinh(L*Vz) - 2*sinh((L*Vz)/2)*(sinh((L*Vy)/2)**2 + 1)*(2*sinh((L*Vz)/4)**2 + 1))
                                                                        - 2*Vy*Vz**2*Mxb*(sinh(L*Vy) - 2*sinh((L*Vy)/2)*(sinh((L*Vz)/2)**2 + 1)*(2*sinh((L*Vy)/4)**2 + 1))
                                                                        - 2*Vy**2*Vz*Mxb*(sinh(L*Vz) - 2*sinh((L*Vz)/2)*(sinh((L*Vy)/2)**2 + 1)*(2*sinh((L*Vz)/4)**2 + 1)))
                   + L**3*(2*Vy**3*Vz**2*Mxb*(Oy*sinh(L*Vy) - Oz*sinh(L*Vy) - 2*Oy*sinh((L*Vy)/2)*(sinh((L*Vz)/2)**2 + 1)*(2*sinh((L*Vy)/4)**2 + 1)
                                              + 2*Oz*sinh((L*Vy)/2)*(sinh((L*Vz)/2)**2 + 1)*(2*sinh((L*Vy)/4)**2 + 1))
                           - 2*Vy**2*Vz**3*Mxb*(Oy*sinh(L*Vz) - Oz*sinh(L*Vz) - 2*Oy*sinh((L*Vz)/2)*(sinh((L*Vy)/2)**2 + 1)*(2*sinh((L*Vz)/4)**2 + 1)
                                                + 2*Oz*sinh((L*Vz)/2)*(sinh((L*Vy)/2)**2 + 1)*(2*sinh((L*Vz)/4)**2 + 1))))
        eq7_den = ((Vy**2 - Vz**2)*(L*Vy*sinh(L*Vy) - 4*sinh((L*Vy)/2)**2 + 4*L**2*Vy**2*Oy*sinh((L*Vy)/2)**2)
                   *(L*Vz*sinh(L*Vz) - 4*sinh((L*Vz)/2)**2 + 4*L**2*Vz**2*Oz*sinh((L*Vz)/2)**2))

        elm_tangstf[4,5]   = eq7_num/eq7_den
        elm_tangstf[10,11] = -eq7_num/eq7_den

    #   
    #
    #end
    # ------------------------------------------------------------------------------------------------------------------------- */ 
    # Filling the rest of the matrix
    # impose the geometry
    elm_tangstf += np.triu(elm_tangstf, k=1).T
    #
    return elm_tangstf
    # =========================================================== End of File */
#
#
#
def Ke_new(Le:float, Ax:float,
           Jx:float, Iy:float, Iz:float,
           Emod:float, Gmod:float,
           Oy:float, Oz:float,
           mu_y: float, mu_z: float, 
           Fb:list[float],
           tension:bool):
    """ """
    #
    # --------------------------------------
    #
    A = Ax
    L = Le
    E = Emod
    G = Gmod
    #
    L2 = Le**2
    #
    # --------------------------------------
    #
    # Tension
    if tension:
        # Positive axial force
        #
        Ay = mu_y / (sqrt(1 + Oy * mu_y**2 * L2))
        Az = mu_z / (sqrt(1 + Oz * mu_z**2 * L2))        
        #
        def Cxy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3             

            return 2 * (sinh(LA) - LA) + Om * (LA3 + LA2 * sinh(LA))
        #
        #Cy = Cxy(A=Ay, Om=Oy, L=Le)
        #Cz = Cxy(A=Az, Om=Oz, L=Le)        
        #
        #
        def Dxy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2

            return (2 * LA * cosh(LA/2) - 4 * sinh(LA/2) + 4 * LA2 * Om * sinh(LA/2))**2
        #
        #Dy = Dxy(A=Ay, Om=Oy, L=Le)
        #Dz = Dxy(A=Az, Om=Oz, L=Le)        
        #
        #
        def PHI_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            #
            Om2 = Om**2

            return (4 * sinh(LA)**2 - 16 * sinh(LA / 2.0)**2 + 4 * LA * sinh(LA) + LA2 * sinh(LA)**2
                    + 4 * LA4 * Om2 * sinh(LA)**2 - 16 * LA4 * Om2 * sinh(LA / 2)**2
                    - 4 * LA3 * Om * sinh(LA) - 4 * LA * cosh(LA) * sinh(LA) - 8 * LA2 * Om * sinh(LA)**2
                    + 32 * LA2 * Om * sinh(LA / 2)**2 + 4 * LA3 * Om * cosh(LA) * sinh(LA))
        #
        #PHI_y = PHI_xy(A=Ay, Om=Oy, L=Le)
        #PHI_z = PHI_xy(A=Az, Om=Oz, L=Le)
        #
        #
        def Bxy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            LA6 = LA**6
            LA7 = LA**7
            #
            Om2 = Om**2
            Om3 = Om**3
            #

            return (LA3 - 2 * sinh(LA) + 2 * cosh(LA) * sinh(LA) - LA5 * Om
                    - 2 * LA * sinh(LA)**2 - 2 * LA2 * sinh(LA) + 4 * LA * sinh(LA / 2)**2 + 2 * LA2 * Om * sinh(LA)
                    + 4 * LA4 * Om * sinh(LA) + LA2 * cosh(LA) * sinh(LA) + 2 * LA4 * Om2 * sinh(LA)
                    - 12 * LA3 * Om * sinh(LA / 2)**2 - 2 * LA6 * Om2 * sinh(LA) - 2 * LA6 * Om3 * sinh(LA)
                    + 2 * LA5 * Om2 * sinh(LA)**2 + 12 * LA5 * Om2 * sinh(LA / 2)**2 - 4 * LA7 * Om3 * sinh(LA / 2)**2
                    - 2 * LA2 * Om * cosh(LA) * sinh(LA) + LA4 * Om * cosh(LA) * sinh(LA)
                    - 2 * LA4 * Om2 * cosh(LA) * sinh(LA) + 2 * LA6 * Om3 * cosh(LA) * sinh(LA))
        #
        #By = Bxy(A=Ay, Om=Oy, L=Le)
        #Bz = Bxy(A=Az, Om=Oz, L=Le)        
        #
        #
        def Fxy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            LA6 = LA**6
            LA7 = LA**7            
            #
            Om2 = Om**2
            Om3 = Om**3

            return (8 * LA7 * Om3 * sinh(LA / 2)**2 - 8 * LA6 * Om3 * sinh(LA) * sinh(LA / 2)**2 + 4 * LA6 * Om2 * sinh(LA)
                    - 4 * LA5 * Om2 * sinh(LA)**2 - 24 * LA5 * Om2 * sinh(LA / 2)**2
                    + 4 * LA5 * Om * sinh(LA / 2)**2 + 2 * LA5 * Om + 8 * LA4 * Om2 * sinh(LA) * sinh(LA / 2)**2
                    - 10 * LA4 * Om * sinh(LA) + 24 * LA3 * Om * sinh(LA / 2)**2 - 4 * LA3 * sinh(LA / 2)**2
                    - 2 * LA3 + 8 * LA2 * Om * sinh(LA) * sinh(LA / 2)**2 + 2 * LA2 * sinh(LA) + 4 * LA * sinh(LA)**2
                    - 8 * LA * sinh(LA / 2)**2 - 8 * sinh(LA) * sinh(LA/2)**2)
        #
        #Fy = Fxy(A=Ay, Om=Oy, L=Le)
        #Fz = Fxy(A=Az, Om=Oz, L=Le)        
        #
    else:
        # Negative axial force
        #
        Ay = mu_y / (sqrt(1 - Oy * mu_y**2 * L2))
        Az = mu_z / (sqrt(1 - Oz * mu_z**2 * L2))        
        #
        def Cxy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3            

            return (sin(LA) - LA - LA3 * Om - LA2 * Om * sin(LA))
        #
        #Cy = Cxy(A=Ay, Om=Oy, L=Le)
        #Cz = Cxy(A=Az, Om=Oz, L=Le)
        #
        #
        def Dxy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            #
            Om2 = Om**2

            return (4 * cos(LA) - LA2 - 8 * LA2 * Om - LA2 * cos(LA) + 4 * LA * sin(LA)
                    - 4 * LA4 * Om2 + 8 * LA2 * Om * cos(LA)
                    + 4 * LA3 * Om * sin(LA) + 4 * LA4 * Om2 * cos(LA) - 4)            
        #
        #Dy = Dxy(A=Ay, Om=Oy, L=Le)
        #Dz = Dxy(A=Az, Om=Oz, L=Le)
        #
        #
        def PHI_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4            
            #
            Om2 = Om**2

            return (4 * cos(2 * LA) - 16 * cos(LA) + LA2 + 24 * LA2 * Om - LA2 * cos(2 * LA)
                    - 8 * LA * sin(LA) + 4 * LA * sin(2 * LA) + 12 * LA4 * Om2 - 32 * LA2 * Om * cos(LA)
                    + 8 * LA2 * Om * cos(2 * LA) - 8 * LA3 * Om * sin(LA) + 4 * LA3 * Om * sin(2 * LA)
                    - 16 * LA4 * Om2 * cos(LA) + 4 * LA4 * Om2 * cos(2 * LA) + 12)
        #
        #PHI_y = PHI_xy(A=Ay, Om=Oy, L=Le)
        #PHI_z = PHI_xy(A=Az, Om=Oz, L=Le)
        #
        #
        def Bxy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            LA6 = LA**6
            LA7 = LA**7
            #
            Om2 = Om**2
            Om3 = Om**3
            #

            return (2 * (-sin(2 * LA) + 2 * sin(LA) + LA3 + LA + 6 * LA3 * Om + LA5 * Om
                         - 2 * LA2 * sin(LA) + (LA2 * sin(2 * LA)) / 2 - 2 * LA * cos(LA) + LA * cos(2 * LA)
                         + 7 * LA5 * Om2 + 2 * LA7 * Om3 - 6 * LA3 * Om * cos(LA) + 2 * LA2 * Om * sin(LA)
                         - LA2 * Om * sin(2 * LA) - 4 * LA4 * Om * sin(LA) - (LA4 * Om * sin(2 * LA)) / 2
                         - 6 * LA5 * Om2 * cos(LA) - LA5 * Om2 * cos(2 * LA) - 2 * LA7 * Om3 * cos(LA)
                         - 2 * LA4 * Om2 * sin(LA) + LA4 * Om2 * sin(2 * LA)
                         - 2 * LA6 * Om2 * sin(LA)
                         - 2 * LA6 * Om3 * sin(LA) + LA6 * Om3 * sin(2 * LA)))
        #
        #By = Bxy(A=Ay, Om=Oy, L=Le)
        #Bz = Bxy(A=Az, Om=Oz, L=Le)
        #
        #
        def Fxy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            LA6 = LA**6
            LA7 = LA**7
            #
            Om2 = Om**2
            Om3 = Om**3
            #

            return 4 * (sin(2 * LA) - 2 * sin(LA) - LA - 6 * LA3 * Om - LA3 * cos(LA)
                        + LA2 * sin(LA) + 2 * LA * cos(LA) - LA * cos(2 * LA) - 7 * LA5 * Om2 - 2 * LA7 * Om3
                        + 6 * LA3 * Om * cos(LA) - LA5 * Om * cos(LA) - 2 * LA2 * Om * sin(LA)
                        + LA2 * Om * sin(2 * LA) + 5 * LA4 * Om * sin(LA) + 6 * LA5 * Om2 * cos(LA)
                        + LA5 * Om2 * cos(2 * LA) + 2 * LA7 * Om3 * cos(LA) + 2 * LA4 * Om2 * sin(LA)
                        - LA4 * Om2 * sin(2 * LA) + 2 * LA6 * Om2 * sin(LA) + 2 * LA6 * Om3 * sin(LA)
                        - LA6 * Om3 * sin(2 * LA))
        #
        #Fy = Fxy(A=Ay, Om=Oy, L=Le)
        #Fz = Fxy(A=Az, Om=Oz, L=Le)
        #
    #
    #
    #
    # --------------------------------------
    #
    Cy = Cxy(A=Ay, Om=Oy, L=Le)
    Cz = Cxy(A=Az, Om=Oz, L=Le)
    #
    Dy = Dxy(A=Ay, Om=Oy, L=Le)
    Dz = Dxy(A=Az, Om=Oz, L=Le)    
    #
    PHI_y = PHI_xy(A=Ay, Om=Oy, L=Le)
    PHI_z = PHI_xy(A=Az, Om=Oz, L=Le)
    #
    By = Bxy(A=Ay, Om=Oy, L=Le)
    Bz = Bxy(A=Az, Om=Oz, L=Le)
    #
    Fy = Fxy(A=Ay, Om=Oy, L=Le)
    Fz = Fxy(A=Az, Om=Oz, L=Le)    
    #
    # --------------------------------------
    #
    ke = zeros((12, 12), dtype=np.float32)

    ke[0, 0] = (A * E / L)

    ke[1, 1] = (12 * E * Iz / L**3) * (L * Ay)**3 * Cy / (12 * Dy)

    ke[2, 2] = (12 * E * Iy / L**3) * (L * Az)**3 * Cz / (12 * Dz)

    ke[3, 3] = (G * Jx / L)
    #
    ke[4, 4] = (4 * E * Iy / L) * (L * Az * Bz / (8 * PHI_z))
    ke[5, 5] = (4 * E * Iz / L) * (L * Ay * By / (8 * PHI_y))
    #
    #
    ke[1, 5] = (6 * E * Iz / L2) * (L * Ay)**3 * Cy / (12 * Dy)
    ke[2, 4] = -(6 * E * Iy / L2) * (L * Az)**3 * Cz / (12 * Dz)
    #
    # ---------------
    # impose the symmetry
    #
    ke[6, 6] = ke[0, 0]
    ke[7, 7] = ke[1, 1]
    ke[8, 8] = ke[2, 2]
    ke[9, 9] = ke[3, 3]
    ke[10, 10] = ke[4, 4]
    ke[11, 11] = ke[5, 5] 
    #    
    #
    ke[0 , 6] = -ke[0 , 0]
    ke[1 , 7] = -ke[1 , 1]
    ke[1 , 11] = ke[1 , 5]
    ke[2 , 8] = -ke[2 , 2]
    ke[2 , 10] = ke[2 , 4]
    ke[3 , 9] = -ke[3 , 3]
    #
    # discrepancy
    ke[4, 10] = (2 * E * Iy / L) * (L * Az * Fz) / (8 * PHI_z)
    ke[5, 11] = (2 * E * Iz / L) * (L * Ay * Fy) / (8 * PHI_y)
    #
    #
    ke[4 , 8]  = -ke[2 , 4]
    ke[5 , 7]  = -ke[1 , 5]
    ke[7 , 11] = -ke[1 , 5]
    ke[8 , 10] = -ke[2 , 4]  
    #
    # impose the geometry
    ke += triu(ke, k=1).T
    #
    return ke   
#
def Kg_new(Le:float, Ax:float,
           Jx:float, Iy:float, Iz:float,
           Emod:float, Gmod:float,
           Oy:float, Oz:float,
           mu_y: float, mu_z: float, 
           Fb:list[float],
           tension:bool):
    """ """
    (Fxa, Fya, Fza, Mxa, Mya, Mza,
     Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb.q_loc
    #
    # --------------------------------------
    #    
    P = abs(Fb.Fx) #+ 0.1
    A = Ax
    L = Le
    L2 = Le**2
    #E = Emod
    #G = Gmod
    #
    #
    # --------------------------------------
    #
    # Tension
    if tension:
        # Positive axial force
        #
        Ay = mu_y / (sqrt(1 + Oy * mu_y**2 * L2))
        Az = mu_z / (sqrt(1 + Oz * mu_z**2 * L2))        
        #
        def Dxy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2            

            return (LA * sinh(LA) - 2 * LA2 * Om - 2 * cosh(LA) + 2 * LA2 * Om * cosh(LA) + 2)
        #
        #
        def Cxy(A: float, L: float):
            """ """
            LA = L * A

            return (cosh(LA) - 1.0) * (sinh(LA) - LA)
        #
        #
        def alpha_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            #
            Om2 = Om**2

            return (A * (3 * LA - 3 * sinh(LA) - 2 * LA3 * Om + LA5 * Om2 + 2 * LA * sinh(LA / 2)**2
                         + 2 * LA2 * Om * sinh(LA) + LA4 * Om2 * sinh(LA)))
        #
        #
        def beta_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            #
            Om2 = Om**2

            return (2 * (LA2 + 4 * sinh(LA / 2)**2 - 2 * LA * sinh(LA) + LA2 * sinh(LA / 2)**2
                         + 4 * LA4 * Om2 * sinh(LA / 2)**2 + 2 * LA3 * Om * sinh(LA)
                         - 8 * LA2 * Om * sinh(LA / 2)**2))
        #
        #
        def gamma_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2

            return ((LA2 * Om - 1)**2 * (LA2 - 8 * sinh(LA / 2)**2 + LA * sinh(LA)))
        #        
        #
        def PHI_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            #
            Om2 = Om**2

            return (4 * LA4 * Om2 * sinh(LA)**2 - 16 * LA4 * Om2 * sinh(LA / 2)**2
                    + 8 * LA3 * Om * sinh(LA) * sinh(LA / 2)**2 - 8 * LA2 * Om * sinh(LA)**2
                    + 32 * LA2 * Om * sinh(LA / 2)**2 + LA2 * sinh(LA)**2 - 8 * LA * sinh(LA) * sinh(LA / 2)**2
                    + 4 * sinh(LA)**2 - 16 * sinh(LA / 2)**2)

        #
        #
        def phi_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            #
            Om2 = Om**2

            return ((LA2 * Om - 1)**2
                    * (-4 * LA5 * Oy**2 * sinh(LA / 2)**2
                       + 4 * LA4 * Om2 * sinh(LA) * sinh(LA / 2)**2 - 2 * LA4 * Om * sinh(LA) + 2 * LA3 * Om * sinh(LA)**2
                       + 8 * LA3 * Om * sinh(LA / 2)**2 - LA3 - 8 * LA2 * Om * sinh(LA) * sinh(LA / 2)**2
                       + 2 * LA2 * sinh(LA) * sinh(LA / 2)**2 + 3 * LA2 * sinh(LA)
                       - 4 * LA * sinh(LA)**2 + 4 * LA * sinh(LA / 2)**2 + 4 * sinh(LA) * sinh(LA / 2)**2))
        #
        #
        def omega_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            #
            Om2 = Om**2
            # TODO : why 2
            return (2 * (4 * LA5 * Om2 * sinh(LA / 2)**2 + 4 * LA4 * Om2 * sinh(LA) * sinh(LA / 2)**2
                         + 2 * LA4 * Om * sinh(LA) + 2 * LA3 * Om * sinh(LA)**2
                         - 8 * LA3 * Om * sinh(LA / 2)**2 + LA3 - 8 * LA2 * Om * sinh(LA) * sinh(LA / 2)**2
                         + 2 * LA2 * sinh(LA) * sinh(LA / 2)**2 - LA2 * sinh(LA)
                         - 2 * LA * sinh(LA)**2 + 4 * LA * sinh(LA / 2)**2 + 4 * sinh(LA) * sinh(LA / 2)**2))
        #
        #
        def psi_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            #
            Om2 = Om**2

            return (- 4 * LA5 * Om2 * sinh(LA / 2)**2 - 4 * LA4 * Om2 * sinh(LA) * sinh(LA / 2)**2
                    - 2 * LA4 * Om * sinh(LA) - 2 * LA3 * Om * sinh(LA)**2
                    + 8 * LA3 * Om * sinh(LA / 2)**2 - 2 * LA3 * sinh(LA / 2)**2 - LA3
                    + 8 * LA2 * Om * sinh(LA) * sinh(LA / 2)**2 + LA2 * sinh(LA)
                    + 2 * LA * sinh(LA)**2 - 4 * LA * sinh(LA / 2)**2 - 4 * sinh(LA) * sinh(LA / 2)**2)
        #
        #
        def eta_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            #
            Om2 = Om**2

            return ((LA2 * Om - 1)**2
                    * (4 * LA5 * Om2 * sinh(LA / 2)**2 - 4 * LA4 * Om2 * sinh(LA) * sinh(LA / 2)**2
                       + 2 * LA4 * Om * sinh(LA) - 2 * LA3 * Om * sinh(LA)**2
                       - 8 * LA3 * Om * sinh(LA / 2)**2 + 2 * LA3 * sinh(LA / 2)**2
                       + LA3 + 8 * LA2 * Om * sinh(LA) * sinh(LA / 2)**2
                       - 3 * LA2 * sinh(LA) + 12 * LA * sinh(LA / 2)**2
                       - 4 * sinh(LA) * sinh(LA / 2)**2))
        #
    else:
        # Negative axial force
        #
        Ay = mu_y / (sqrt(1 - Oy * mu_y**2 * L2))
        Az = mu_z / (sqrt(1 - Oz * mu_z**2 * L2))        
        #
        def Dxy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2

            return (4 * sin(LA / 2)**2 - LA * sin(LA) + 4 * LA2 * Om * sin(LA / 2)**2)
        #
        #
        def Cxy(A: float, L: float):
            """ """
            LA = L * A
            return 2 * sin(LA / 2)**2 * (LA - sin(LA))
        #
        #
        def alpha_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            #
            Om2 = Om**2

            return (A * (3 * LA - 3 * sin(LA) + 2 * LA3 * Om + LA5 * Om2
                         - 2 * LA * sin(LA / 2)**2
                         - 2 * LA2 * Om * sin(LA) + LA4 * Om2 * sin(LA)))
        #
        #
        def beta_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            #LA5 = LA**5
            #
            Om2 = Om**2

            return (2 * (LA2 + 4 * sin(LA / 2)**2 - 2 * LA * sin(LA) - LA2 * sin(LA / 2)**2
                         + 4 * LA4 * Om2 * sin(LA / 2)**2 - 2 * LA3 * Om * sin(LA)
                         + 8 * LA2 * Om * sin(LA / 2)**2))
        #
        #
        def gamma_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2

            return ((Om * LA2 + 1)**2 * (LA2 - 8 * sin(LA / 2)**2 + LA * sin(LA)))
        #        
        #
        def PHI_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            #LA5 = LA**5
            #
            Om2 = Om**2

            return (4 * LA4 * Om2 * sin(LA)**2 - 16 * LA4 * Om2 * sin(LA / 2)**2
                    + 8 * LA3 * Om * sin(LA) * sin(LA / 2)**2
                    + 8 * LA2 * Om * sin(LA)**2 - 32 * LA2 * Om * sin(LA / 2)**2
                    - LA2 * sin(LA)**2 + 8 * LA * sin(LA) * sin(LA / 2)**2
                    + 4 * sin(LA)**2 - 16 * sin(LA / 2)**2)
        #
        #
        def phi_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            #
            Om2 = Om**2            
            #
            return ((Om * LA2 + 1)**2
                    * (-(sin(2 * LA) + LA3 + 5 * LA3 * Om - (LA2 * sin(2 * LA)) / 2
                         + 2 * LA5 * Om2 + 2 * LA * (2 * sin(LA)**2 - 1) + 2 * LA2 * Om * sin(2 * LA)
                         + LA3 * Om * (2 * sin(LA)**2 - 1) + LA4 * Om2 * sin(2 * LA))
                       - (2 * sin(LA / 2)**2 - 1) * (2 * LA5 * Om2 + 4 * LA3 * Om - 2 * LA)
                       + sin(LA) * (2 * LA4 * Om2 + 2 * LA4 * Om + 4 * LA2 * Om + 2 * LA2 + 2)))
        #
        #
        def omega_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            #
            Om2 = Om**2             
            #
            return (sin(LA) * (4 * cos(LA) + 4 * LA2 - 8 * LA2 * Om + 4 * LA4 * Om - 2 * LA2 * cos(LA)
                               - 4 * LA4 * Om2 + 8 * LA2 * Om * cos(LA) + 4 * LA4 * Om2 * cos(LA) - 4)
                    + 4 * LA5 * Om2 * cos(LA) - 4 * LA5 * Om2 - 4 * LA3 * Om * cos(LA)**2 + 8 * LA3 * Om * cos(LA)
                    - 4 * LA3 * Om - 2 * LA3 - 4 * LA * cos(LA)**2 + 4 * LA * cos(LA))
        #
        #
        def psi_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            #
            Om2 = Om**2
            #
            return (-sin(LA) * (2 * cos(LA) + LA2 - 4 * LA2 * Om
                                + 2 * LA4 * Om - 2 * LA4 * Om2
                                + 4 * LA2 * Om * cos(LA) + 2 * LA4 * Om2 * cos(LA) - 2)
                    - 2 * LA5 * Om2 * cos(LA)
                    + 2 * LA5 * Om2 + 2 * LA3 * Om * cos(LA)**2 - 4 * LA3 * Om * cos(LA)
                    + 2 * LA3 * Om + LA3 * cos(LA) + 2 * LA * cos(LA)**2 - 2 * LA * cos(LA))
        #
        #
        def eta_xy(A: float, Om: float,  L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
            LA4 = LA**4
            LA5 = LA**5
            #
            Om2 = Om**2             
            #
            return ((Om * LA2 + 1)**2 *
                    (sin(2 * LA) + 6 * LA + 5 * LA3 * Om
                     + 2 * LA5 * Om2 + 2 * LA2 * Om * sin(2 * LA)
                     + LA3 * Om * (2 * sin(LA)**2 - 1) + LA4 * Om2 * sin(2 * LA)
                     + (2 * sin(LA / 2)**2 - 1) * (2 * LA5 * Om2 + 4 * LA3 * Om - LA3 + 6 * L * A)
                     - sin(LA) * (2 * LA4 * Om2 + 2 * LA4 * Om + 4 * LA2 * Om + 3 * LA2 + 2)))

        #
    #
    # --------------------------------------
    #
    Dy = Dxy(A=Ay, Om=Oy, L=Le)
    Dz = Dxy(A=Az, Om=Oz, L=Le)
    #
    Cy = Cxy(A=Ay, L=Le)
    Cz = Cxy(A=Az, L=Le)
    #
    alpha_y = alpha_xy(A=Ay, Om=Oy, L=Le)
    alpha_z = alpha_xy(A=Az, Om=Oz, L=Le)
    #
    beta_y = beta_xy(A=Ay, Om=Oy, L=Le)
    beta_z = beta_xy(A=Az, Om=Oz, L=Le)
    #
    PHI_y = PHI_xy(A=Ay, Om=Oy, L=Le)
    PHI_z = PHI_xy(A=Az, Om=Oz, L=Le)
    #
    phi_y = phi_xy(A=Ay, Om=Oy, L=Le)
    phi_z = phi_xy(A=Az, Om=Oz, L=Le)
    #
    omega_y = omega_xy(A=Ay, Om=Oy, L=Le)
    omega_z = omega_xy(A=Az, Om=Oz, L=Le)
    #
    psi_y = psi_xy(A=Ay, Om=Oy, L=Le)
    psi_z = psi_xy(A=Az, Om=Oz, L=Le)
    #
    gamma_y = gamma_xy(A=Ay, Om=Oy, L=Le)
    gamma_z = gamma_xy(A=Az, Om=Oz, L=Le)
    #
    eta_y = eta_xy(A=Ay, Om=Oy, L=Le)
    eta_z = eta_xy(A=Az, Om=Oz, L=Le)    
    #
    # --------------------------------------
    #    
    kg = zeros((12, 12), dtype=np.float32)
    #
    # --------------------------------------
    # 
    kg[0 , 0] = P / L
    kg[1 , 1] =  P*(alpha_y/beta_y + Iz*Ay**3*Cy / (A*Dy**2))
    kg[2 , 2] =  P*(alpha_z/beta_z + Iy*Az**3*Cz / (A*Dz**2))
    #kg[3 , 3] = 0
    # TODO : discrepancy
    kg[4 , 4] = P*(phi_z/(2*Az*PHI_z) + Iy*Az*omega_z/(4*A*PHI_z))
    kg[5 , 5] = P*(phi_y/(2*Ay*PHI_y) + Iz*Ay*omega_y/(4*A*PHI_y))
    #
    #
    kg[4 , 0] = -Mya / L
    #kg[2 , 3] = 0 
    #kg[3 , 4] = 0 
    kg[5 , 6] = Mza / L

    #kg[8 , 9] = 0 
    #kg[9 , 10] = 0 

    #kg[1 , 3] = 0 
    kg[2 , 4] = -P*(gamma_z/(2*beta_z) + Iy*Az**3*Cz/(2*A*Dz**2))
    #kg[3 , 5] = 0 
    kg[4 , 6] = -kg[4 , 0]
    #
    #
    kg[0 , 10] = -Myb / L
    #kg[7 , 9] = 0
    #kg[8 , 10] = P*(gamma_z/(2*beta_z) + Iy*Az**3*Cz/(2*A*Dz**2))
    #kg[9 , 11] = 0
    kg[6 , 4] = -kg[4 , 0]
    #
    #
    kg[0 , 4]  = -kg[4 , 6]
    kg[1 , 5]  = P*(gamma_y/(2*beta_y) + Iz*Ay**3*Cy/(2*A*Dy**2))
    #kg[3 , 7] =  0
    #kg[5 , 9] =  0
    kg[10 , 6] =  -kg[0 , 10]
    #
    kg[6 , 10] = Myb / L
    kg[0 , 5]  = - kg[5 , 6]
    #kg[3 , 8] = 0
    #kg[4 , 9] = 0 
    kg[6 , 11]  = Mzb / L
    kg[0 , 11] = -kg[6 , 11]
    #
    # ---------------
    # impose the symmetry
    #
    kg[6, 6] = kg[0, 0]
    kg[7, 7] = kg[1, 1]
    kg[8, 8] = kg[2, 2]
    kg[9, 9] = kg[3, 3]
    kg[10, 10] = kg[4, 4]
    kg[11, 11] = kg[5, 5]    
    #
    kg[0 , 6] = -kg[0 , 0]
    kg[1 , 7] = -kg[1 , 1]
    kg[1 , 11] = kg[1 , 5]
    kg[2 , 8] = -kg[2 , 2]
    kg[2 , 10] = kg[2 , 4]
    #kg[3 , 9] = -kg[3 , 3]
    #kg[9, 3] = kg[3 , 9]
    #
    kg[4 , 10] = P*(eta_z/(2*Az*PHI_z) + Iy*Az*psi_z/(2*A*PHI_z))
    kg[5 , 11] = P*(eta_y/(2*Ay*PHI_y) + Iz*Ay*psi_y/(2*A*PHI_y))
    #
    kg[4 , 8]  = -kg[2 , 4]
    kg[5 , 7]  = -kg[1 , 5]
    kg[7 , 11] = -kg[1 , 5]
    kg[8 , 10] = -kg[2 , 4] 
    #
    #
    # impose the geometry
    kg += triu(kg, k=1).T
    #
    return kg    
#
def Kint_new(Le:float, Ax:float,
             Jx:float, Iy:float, Iz:float,
             Emod:float, Gmod:float,
             Oy:float, Oz:float,
             mu_y: float, mu_z: float, 
             Fb:list[float],
             tension:bool):
    """ """
    (Fxa, Fya, Fza, Mxa, Mya, Mza,
     Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb.q_loc
    #
    # --------------------------------------
    #    
    P = abs(Fb.Fx) #+ 0.1
    A = Ax
    L = Le
    L2 = Le**2
    #E = Emod
    #G = Gmod
    #
    Jp = (Iy + Iz)
    #
    #
    #mu_y = max(1.0, sqrt(abs(P) / (Emod * Iz)))
    #mu_z = max(1.0, sqrt(abs(P) / (Emod * Iy)))
    #
    L2 = Le**2
    #
    #Ay = mu_y / (sqrt(1 + Oy * mu_y**2 * L2))
    #Az = mu_z / (sqrt(1 + Oz * mu_z**2 * L2))
    # 
    # --------------------------------------
    #
    # Tension
    if tension:
        # Positive axial force
        #
        Ay = mu_y / (sqrt(1 + Oy * mu_y**2 * L2))
        Az = mu_z / (sqrt(1 + Oz * mu_z**2 * L2))        
        #
        def Dxy(Ay: float, Az: float, Om: float, L: float):
            """ """
            LAy = L * Ay
            LAy2 = LAy**2
            LAz = L * Az
            
            return (sinh(LAz / 2)
                    * (2 * sinh(LAy / 2) - LAy * cosh(LAy / 2)
                       - 2 * Om * LAy2 * sinh(LAy / 2)))
        #
        Dy = Dxy(Ay=Ay, Az=Az, Om=Oy, L=Le)
        Dz = Dxy(Ay=Az, Az=Ay, Om=Oz, L=Le)
        #
        #
        def D_1(Ay: float, Az: float, Oy: float, Oz: float, L: float):
            """ """
            LAy = L * Ay
            LAy2 = LAy**2
            LAz = L * Az
            LAz2 = LAz**2
            #
            
            return ((LAy * sinh(LAy) - 4 * sinh(LAy / 2)**2 + 4 * LAy2 * Oy * sinh(LAy / 2)**2)
                    * (LAz * sinh(LAz) - 4 * sinh(LAz / 2)**2 + 4 * LAz2 * Oz * sinh(LAz / 2)**2))
        #
        D1 = D_1(Ay=Ay, Az=Az, Oy=Oy, Oz=Oz, L=Le)
        #
        def D_2(Ay: float, Az: float, Oy: float, Oz: float, L: float):
            """ """
            LAy = L * Ay
            LAy2 = LAy**2
            LAz = L * Az
            LAz2 = LAz**2
            #
            Ay2 = Ay**2
            Az2 = Az**2
            
            return ((Ay2 - Az2)
                    * ((LAy * sinh(LAy) - 2 * LAy2 * Oy - 2 * cosh(LAy) + 2 * LAy2 * Oy * cosh(LAy) + 2)
                       * (LAz * sinh(LAz) - 2 * LAz2 * Oz - 2 * cosh(LAz) + 2 * LAz2 * Oz * cosh(LAz) + 2)))
        #
        D2 = D_2(Ay=Ay, Az=Az, Oy=Oy, Oz=Oz, L=Le)
        #
        def psi_xy(Ay: float, Az: float, Om: float, L: float):
            """ """
            LAy = L * Ay
            LAz = L * Az
            
            return (Ay * cosh(LAy / 2) * sinh(LAz / 2) - Az * cosh(LAz / 2) * sinh(LAy / 2))
        #
        #
        def phi_xy(A: float, Om: float, L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            
            return ((1 - Om * LA2)
                    * (LA - 2 * sinh(LA) - 2 * cosh(LA) + LA * (cosh(LA) + sinh(LA)) + 2))
        #
        #
        def PHI_xy(A: float, Om: float, L: float):
            """ """
            LA = L * A
            
            return (cosh(LA) + sinh(LA) - 1.0)
        #
        #
        def phi1(Ay: float, Az: float, Oy: float, Oz: float, L: float):
            """ """
            LAy = L * Ay
            LAz = L * Az
            #
            L2 = L**2
            #
            Ay2 = Ay**2
            Ay3 = Ay**3
            #
            Az2 = Az**2
            Az3 = Az**3
            #
            Ayz = Ay * Az

            return (L**3 * (Ay - Az) * (2 * Ay2 * Az3 * sinh(LAz) * sinh(LAy / 2)**2
                                        - 2 * Ay3 * Az2 * sinh(LAy) * sinh(LAz / 2)**2)
                    + L * (Ay2 - Az2) * (Ay * sinh(LAy) - Az * sinh(LAz) - Ay * sinh(LAy) * (2 * sinh(LAz / 2)**2 + 1.0)
                                         + Az * sinh(LAz) * (2 * sinh(LAy / 2)**2 + 1.0))
                    - (L2 * Ayz * (4 * Ayz * sinh(LAy / 2)**2 + 4 * Ayz * sinh(LAz / 2)**2 - Ay2 * sinh(LAy) * sinh(LAz)
                                   - Az2 * sinh(LAy) * sinh(LAz) + 8 * Ayz * sinh(LAy / 2)**2 * sinh(LAz / 2)**2)) / 2.0)
        #
        #
        def phi2(Ay: float, Az: float, Oy: float, Oz: float, L: float):
            """ """
            LAy = L * Ay
            LAz = L * Az
            #
            L2 = L**2
            L3 = L**3
            L4 = L**4
            #L5 = L**5
            #
            Ay2 = Ay**2
            Ay3 = Ay**3
            Ay4 = Ay**4
            #Ay5 = Ay**5
            #
            Az2 = Az**2
            Az3 = Az**3
            Az4 = Az**4
            
            return (4 * Ay2 - 4 * Az2 - 4 * Ay2 * cosh(LAy / 2)**2 - 4 * Ay2 * cosh(LAz / 2)**2 + 4 * Az2 * cosh(LAy / 2)**2
                    + 4 * Az2 * cosh(LAz / 2)**2 - 4 * L2 * Ay4 * Oy + 4 * L2 * Az4 * Oz + 4 * Ay2 * cosh(LAy / 2)**2 * cosh(LAz / 2)**2
                    - 4 * Az2 * cosh(LAy / 2)**2 * cosh(LAz / 2)**2 - 2 * LAy * Az2 * sinh(LAy) + 2 * L * Ay**2 * Az * sinh(LAz)
                    + L2 * Ay2 * Az2 * cosh(LAy / 2)**2 - L2 * Ay2 * Az2 * cosh(LAz / 2)**2 + 4 * L2 * Ay2 * Az2 * Oy
                    - 4 * L2 * Ay2 * Az2 * Oz + 4 * L2 * Ay4 * Oy * cosh(LAy / 2)**2 + 4 * L2 * Ay4 * Oy * cosh(LAz / 2)**2
                    - 4 * L2 * Az4 * Oz * cosh(LAy / 2)**2 - 4 * L2 * Az4 * Oz * cosh(LAz / 2)**2 - 4 * L4 * Ay2 * Az4 * Oy * Oz + 4 * L4 * Ay4 * Az2 * Oy * Oz
                    + L3 * Ay3 * Az2 * Oy * sinh(LAy) - L3 * Ay2 * Az3 * Oz * sinh(LAz) - 4 * L2 * Ay2 * Az2 * Oy * cosh(LAy / 2)**2
                    - 4 * L2 * Ay2 * Az2 * Oz * cosh(LAy / 2)**2 - 4 * L2 * Ay2 * Az2 * Oy * cosh(LAz / 2)**2
                    + 4 * L2 * Ay2 * Az2 * Oz * cosh(LAz / 2)**2 - 4 * L2 * Ay4 * Oy * cosh(LAy / 2)**2 * cosh(LAz / 2)**2
                    + 4 * L2 * Az4 * Oz * cosh(LAy / 2)**2 * cosh(LAz / 2)**2 + L3 * Ay * Az4 * Oz * sinh(LAy) - L3 * Ay4 * Az * Oy * sinh(LAz)
                    + 4 * L * Ay * Az2 * cosh(LAy / 2) * cosh(LAz / 2)**2 * sinh(LAy / 2) - 4 * L2 * Ay2 * Az2 * Oy * cosh(LAy / 2)**2 * cosh(LAz / 2)**2
                    - 4 * L * Ay2 * Az * cosh(LAy / 2)**2 * cosh(LAz / 2) * sinh(LAz / 2)
                    - 4 * L2 * Ay2 * Az2 * Oz * cosh(LAy / 2)**2 * cosh(LAz / 2)**2 + 4 * L4 * Ay2 * Az4 * Oy * Oz * cosh(LAy / 2)**2
                    - 4 * L4 * Ay4 * Az2 * Oy * Oz * cosh(LAy / 2)**2 + 4 * L4 * Ay2 * Az4 * Oy * Oz * cosh(LAz / 2)**2
                    + 2 * L3 * Ay2 * Az3 * Oz * cosh(LAy / 2)**2 * cosh(LAz / 2) * sinh(LAz / 2)
                    - 2 * L3 * Ay3 * Az2 * Oy * cosh(LAy / 2)* cosh(LAz / 2)**2 * sinh(LAy / 2)
                    - 4 * L4 * Ay4 * Az2 * Oy * Oz * cosh(LAz / 2)**2
                    - 4 * L4 * Ay2 * Az4 * Oy * Oz * cosh(LAy / 2)**2 * cosh(LAz / 2)**2
                    + 4 * L4 * Ay4 * Az2 * Oy * Oz * cosh(LAy / 2)**2 * cosh(LAz / 2)**2
                    - 2 * L3 * Ay * Az4 * Oz * cosh(LAy / 2) * cosh(LAz / 2)**2 * sinh(LAy / 2)
                    + 2 * L3 * Ay4 * Az * Oy * cosh(LAy / 2)**2 * cosh(LAz / 2) * sinh(LAz / 2))
        #
        #
        #         
    else:
        # Negative axial force
        #
        Ay = mu_y / (sqrt(1 - Oy * mu_y**2 * L2))
        Az = mu_z / (sqrt(1 - Oz * mu_z**2 * L2))        
        #
        def Dxy(Ay:float, Az:float, Oy: float, Oz: float, L:float): 
            """ """
            LAy = L * Ay
            LAz = L * Az
            L2 = L**2
            Ay2 = Ay**2
            Az2 = Az**2
            
            return ((4 * sin(LAy / 2)**2 - LAy * sin(LAy) + 4 * L2 * Ay2 * Oy * sin(LAy / 2)**2
                     * (4 * sin(LAz / 2)**2 - LAz * sin(LAz) + 4 * L2 * Az2 * Oz * sin(LAz / 2)**2)))
        #
        Dy = Dxy(Ay=Ay, Az=Az, Oy=Oy, Oz=Oz, L=Le)
        Dz = D1 = D2 = Dy
        #
        def psi_xy(Ay: float, Az: float, Om: float, L: float):
            """ """
            LAy = L * Ay
            LAz = L * Az
            #
            L2 = L**2
            #
            Ay2 = Ay**2
            Ay3 = Ay**3
            #
            Az2 = Az**2
            Az3 = Az**3
        
            return (Az2 * (4 * Ay * (sin(LAy) - 2 * cos(LAy / 2) * cos(LAz / 2)**2 * sin(LAy / 2))
                           - 4 * L * Ay2 * (cos(LAy / 2)**2 - cos(LAy / 2)**2 * cos(LAz / 2)**2)
                           + 4 * L2 * Ay3 * (Om * sin(LAy) - 2 * Om * cos(LAy / 2) * cos(LAz / 2)**2 * sin(LAy / 2)))
                    - Az3 * (4 * (sin(LAz) - 2 * cos(LAy / 2)**2 * cos(LAz / 2) * sin(LAz / 2))
                             + 4 * L2 * Ay2 * (Om * sin(LAz) - 2 * Om * cos(LAy / 2)**2 * cos(LAz / 2) * sin(LAz / 2))
                             - LAy * sin(LAy) * sin(LAz)))
        #  
        #
        #def psi_zz(Ay: float, Az: float, Oz: float, L: float):
        #    LAy = L * Ay
        #    LAz = L * Az
        #    #
        #    L2 = L**2
        #    #
        #    Ay2 = Ay**2
        #    Ay3 = Ay**3
        #    #
        #    Az2 = Az**2
        #    Az3 = Az**3
        #
        #    return (Ay2*(4*Az * (sin(LAz) - 2*cos(LAy/2)**2 * cos(LAz/2) * sin(LAz/2))
        #                   - 4*LAz**2 * (cos(LAz/2)**2 - cos(LAy/2)**2 * cos(LAz/2)**2 )
        #                   + 4*L2 * Az3 * (Oz*sin(LAz) - 2 * Oz * cos(LAy/2)**2 * cos(LAz/2) * sin(LAz/2)))
        #            - Ay3 * (4 *(sin(LAy) - 2 * cos(LAy/2) * cos(LAz/2)**2 * sin(LAy/2))
        #                            + 4*L2 * Az2 * (Oz * sin(LAy) -2*Oz * cos(LAy/2) * cos(LAz/2)**2 * sin(LAy/2))
        #                   - L*Az * sin(LAy) * sin(LAz)))
        #
        #
        def phi_xy(A: float, Om: float, L: float):
            """ """
            LA = L * A
            LA2 = LA**2
            LA3 = LA**3
        
            return (2 * (Om * LA2 + 1)
                    * (LA2 + 4 * sin(LA / 2)**2 - 2 * LA * sin(LA) - LA2 * sin(LA / 2)**2
                       - LA3 * Om * sin(LA) + 4 * LA2 * Om * sin(LA / 2)**2))
        #
        #
        def PHI_xy(A: float, Om: float, L: float):
            """ """
            LA = L * A
            LA2 = LA**2
        
            return ((4 * sin(LA / 2)**2 + 2 * LA2 * Om - LA * sin(LA)
                     + 2 * LA2 * Om * (2 * sin(LA / 2)**2 - 1)))
        #
        #
        def phi1(Ay: float, Az: float, Oy: float, Oz: float, L: float):
            """ """
            LAy = L * Ay
            LAz = L * Az
            #
            L2 = L**2
            L3 = L**3
            #
            Ay2 = Ay**2
            #
            Az2 = Az**2
        
            return (sin(LAz / 2)
                    * (2 * L * (Ay2 - Az2) * (2 * Az * cos(LAz / 2) - 2 * Az * cos(LAy / 2)**2 * cos(LAz / 2))
                       - 2 * L3 * Ay2 * Az2 * (2 * Az * cos(LAz / 2) - 2 * Az * cos(LAy / 2)**2 * cos(LAz / 2)) * (Oy - Oz))
                    - sin(LAy / 2)
                    * (2 * L * (Ay2 - Az2) * (2 * Ay * cos(LAy / 2) - 2 * Ay * cos(LAy / 2) * cos(LAz / 2)**2)
                       - 2 * L3 * Ay2 * Az2 * (2 * Ay * cos(LAy / 2) - 2 * Ay * cos(LAy / 2) * cos(LAz / 2)**2) * (Oy - Oz))
                    + (L2 * Ay * Az
                       * (Ay2 * sin(LAy) * sin(LAz) - 4 * Ay * Az * cos(LAz / 2)**2
                          - 4 * Ay * Az * cos(LAy / 2)**2 + Az2 * sin(LAy) * sin(LAz)
                          + 8 * Ay * Az * cos(LAy / 2)**2 * cos(LAz / 2)**2)) / 2.0)
        #
        #
        def phi2(Ay: float, Az: float, Oy: float, Oz: float, L: float):
            LAy = L * Ay
            LAz = L * Az
            #
            #
            L2 = L**2
            L3 = L**3
            L4 = L**4
            #L5 = L**5
            #
            Ay2 = Ay**2
            Ay3 = Ay**3
            Ay4 = Ay**4
            #Ay5 = Ay**5
            #
            Az2 = Az**2
            Az3 = Az**3
            Az4 = Az**4
            #Az5 = Az**5
        
            return (4 * Ay2 * sin(LAy / 2)**2 * sin(LAz / 2)**2 - 4 * Az2 * sin(LAy / 2)**2 * sin(LAz / 2)**2 + L2 * Ay2 * Az2 * sin(LAy / 2)**2
                    - L2 * Ay2 * Az2 * sin(LAz / 2)**2 + 4 * L2 * Ay4 * Oy * sin(LAy / 2)**2 * sin(LAz / 2)**2
                    - 4 * L2 * Az4 * Oz * sin(LAy / 2)**2 * sin(LAz / 2)**2 + 2 * LAy * Az2 * sin(LAy) * sin(LAz / 2)**2
                    - 2 * L * Ay2 * Az * sin(LAz) * sin(LAy / 2)**2 - 4 * L2 * Ay2 * Az2 * Oy * sin(LAy / 2)**2 * sin(LAz / 2)**2
                    + 4 * L2 * Ay2 * Az2 * Oz * sin(LAy / 2)**2 * sin(LAz / 2)**2 - L3 * Ay4 * Az * Oy * sin(LAz) * sin(LAy / 2)**2
                    + L3 * Ay * Az4 * Oz * sin(LAy) * sin(LAz / 2)**2 + L3 * Ay3 * Az2 * Oy * sin(LAy) * sin(LAz / 2)**2
                    - L3 * Ay2 * Az3 * Oz * sin(LAz) * sin(LAy / 2)**2 - 4 * L4 * Ay2 * Az4 * Oy * Oz * sin(LAy / 2)**2 * sin(LAz / 2)**2
                    + 4 * L4 * Ay4 * Az2 * Oy * Oz * sin(LAy / 2)**2 * sin(LAz / 2)**2)
        #
    #
    # --------------------------------------
    #
    psi_y = psi_xy(Ay=Ay, Az=Az, Om=Oy, L=Le)
    psi_z = psi_xy(Ay=Az, Az=Ay, Om=Oz, L=Le)
    #
    phi_y = phi_xy(A=Ay, Om=Oy, L=Le)
    phi_z = phi_xy(A=Az, Om=Oz, L=Le)
    #
    PHI_y = PHI_xy(A=Ay, Om=Oy, L=Le)
    PHI_z = PHI_xy(A=Az, Om=Oz, L=Le)
    #
    phi_1 = phi1(Ay, Az, Oy, Oz, L=Le)
    phi_2 = phi2(Ay, Az, Oy, Oz, L=Le)
    #
    #
    # --------------------------------------
    #     
    kint = zeros((12, 12), dtype=np.float32)
    # --------------------------------------
    #
    kint[1, 3] = Mya / L
    kint[1, 9] = Myb / L
    
    kint[2, 3] = Mza / L
    kint[2, 9] = Mzb / L
    #
    kint[3, 1] = Mya / L
    kint[3, 2] = Mza / L
    kint[3, 3] = Jp * P / (A * L)
    kint[3, 4] = (Mza + Mzb) * phi_z / (L2 * Az**2 * PHI_z) - Mza / 2.0
    kint[3, 5] = - (Mya + Myb) * phi_y / (L2 * Ay**2 * PHI_y) + Mya / 2.0
    kint[3, 7] = - kint[3, 1]
    kint[3, 8] = - kint[3, 2]
    kint[3, 9] = -kint[3 , 3]
    kint[3, 10] = -(Mza + Mzb) * phi_z / (L2 * Az**2 * PHI_z)
    kint[3, 11] = (Mya + Myb) * phi_y / (L2 * Ay**2 * PHI_y)
    
    kint[4, 3] = kint[3, 4]
    kint[4, 9] = - (Mza + Mzb) * phi_z / (L2 * Az**2 * PHI_z)
    
    kint[5, 3] = -(Mya + Myb) * phi_y / (L2 * Ay**2 * PHI_y) + Mya / 2.0
    
    #kint[5, 9] = (Mya + Myb) * phi_y / (L2 * Ay**2 * PHI_y)
    kint[5, 9] = kint[3, 11]

    kint[8, 9] = - Mzb / L

    kint[9, 10] = ((Mza + Mzb) * phi_z / (L2 * Az**2 * PHI_z) - Mzb / 2.0)


    kint[7, 9] = - Myb / L

    kint[9, 11] = - (Mya + Myb) * phi_y / (L2 * Ay**2 * PHI_y) + Myb / 2.0
    
    #
    # Torsion for symmetric sections
    #
    if Iy == Iz:
        # FIXME : not formulas
        #psi_y = Ay*(sin(L*Ay) - L*Ay)
        #D = (4*cos(L*Ay) - 4*L**2*Ay**2*Oy + 2*L*Ay*sin(L*Ay) + 4*L**2*Ay**2*Oy*cos(L*Ay) - 4) 
        #
        kint[1, 4] = -Mxb * Ay**2 * psi_y / Dy
        kint[2, 5] = -Mxb * Az**2 * psi_z / Dz
        #
        kint[4, 7] = Mxb * Ay**2 * psi_y / Dy
        kint[4, 5] = - Mxb * phi_1 / D1

        kint[5, 2] = Mxb * Az**2 * psi_z /  Dz
        kint[5, 10] = Mxb * phi_2 / D1
        
        #
        #1 / 0
    else:
        Adif = (Ay**2 - Az**2) 
        #
        kint[1, 4] = - Mxb * Ay**2 * psi_y / (Adif * Dy)
        
        kint[2, 5] = - Mxb * Az**2 * psi_z / (Adif * Dz)
        
        kint[4, 7] = Mxb * Ay**2 * psi_y / (Adif * Dy)
        kint[4, 5] = - Mxb * phi_1 / (Adif * D1)
        #
        #kint[5 , 8] = - kint[2 , 5]   # 0 + 0 + Mxb*mu_y**2*psi2_y / ((mu_y**2-mu_z**2)* D)
        #kint[7 , 10] = - kint[4 , 7]  # 0 + 0 - Mxb*mu_y**2*psi2_y / ((mu_y**2-mu_z**2)* D)
        #kint[8 , 11] = kint[2 , 5]    # 0 + 0 - Mxb*mu_z**2*psi2_z / ((mu_y**2-mu_z**2)* D)
        #
        kint[5, 2] = Mxb * Az**2 * psi_z / (Adif * Dz)
        kint[5, 10] = Mxb * phi_2 / (Adif * D2)
        #kint[4 , 11] = -kint[5 , 10] # 0 + 0 - Mxb*phi_2/((mu_y**2 - mu_z**2)*D)
        #
        #kint[1 , 10] = -kint[1 , 4]
        #kint[2 , 11] = -kint[2 , 5]
        #
        
        #kint[10 , 11] = - kint[4 , 5]
    #
    kint[5, 4] = kint[4, 5]
    kint[5, 8] = -kint[5, 2]
    kint[7, 10] = - kint[4, 7]
    kint[8, 11] = kint[2, 5]
    #
    kint[4, 1] = kint[1, 4]
    kint[4, 11] = -kint[5, 10]
    #
    kint[1, 10] = -kint[1, 4]
    kint[2, 11] = -kint[2, 5]
    #
    #kint[10 , 11] = - kint[4 , 5]
    kint[2 , 9] = Mzb / L
    #
    # ---------------
    # impose the symmetry
    #
    kint[6, 6] = kint[0, 0]
    kint[7, 7] = kint[1, 1]
    kint[8, 8] = kint[2, 2]
    kint[9, 9] = kint[3, 3]
    kint[10, 10] = kint[4, 4]
    kint[11, 11] = kint[5, 5]     
    #
    #kint[0 , 6] = -kint[0 , 0]
    #kint[1 , 7] = -kint[1 , 1]
    #kint[1 , 11] = kint[1 , 5]
    #kint[2 , 8] = -kint[2 , 2]
    #kint[2 , 10] = kint[2 , 4]
    #kint[3 , 9] = -kint[3 , 3]
    #
    #kint[4 , 10] = (2*E*Iy/L)*(Az*L*Bz)/(8*PHI_z) + (P*eta_z/(2*Az*PHI2_z) + P*Iy*Az*psi2_z/(2*A*PHI2_z)) + 0
    #kint[5 , 11] = (2*E*Iz/L)*(Ay*L*Fy)/(8*PHI_y) + (P*eta_y/(2*Ay*PHI2_y) + P*Iz*Ay*psi_y/(4*A*PHI2_y)) + 0
    #
    #kint[4 , 8]  = -kint[2 , 4]
    #kint[5 , 7]  = -kint[1 , 5]
    #kint[7 , 11] = -kint[1 , 5]
    #kint[8 , 10] = -kint[2 , 4]
    #
    #
    # impose the geometry
    kint += triu(kint, k=1).T
    #
    return kint         
#