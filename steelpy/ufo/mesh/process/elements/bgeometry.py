#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
#
# Python stdlib imports

# package imports
import numpy as np

# ---------------------
#
# Geometry
#
def beam_geom(Le, s):
    """
    Calculates the element geometric matrix  
    """
    # initialize all eg elements to zero 
    eg = np.zeros( 12, 12 )
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
    if len(Fb) == 6:
        # ['x', 'y', 'rz']
        Fxa, Fya, Mza, Fxb, Fyb, Mzb = Fb
        Fb = np.array([Fxa, Fya, 0.0, 0.0, 0.0, Mza,
                       Fxb, Fyb, 0.0, 0.0, 0.0, Mzb])
    else:
        (Fxa, Fya, Fza, Mxa, Mya, Mza,
         Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb
    #
    P = Fxb
    #
    #ax = min(Asy, Asz)
    Phiy = 0
    Phiz = 0
    if shear:
        Phiy = 12*Emod*Iz / (Gmod * Asy * Le**2)
        Phiz = 12*Emod*Iy / (Gmod * Asz * Le**2)
    #
    gk = np.zeros(( 12, 12 ))
    #
    gk[ 0 ][ 0 ] = 0
    gk[ 1 ][ 1 ] = (6/5 + 2*Phiz + Phiz**2) / (1 + Phiz)**2
    gk[ 2 ][ 2 ] = (6/5 + 2*Phiy + Phiy**2) / (1 + Phiy)**2
    gk[ 3 ][ 3 ] = Jx/Ax
    gk[ 4 ][ 4 ] = (2*Le**2/15 + Le**2*Phiy/6 + Le**2*Phiy**2/12) / (1+Phiy)**2
    gk[ 5 ][ 5 ] = (2*Le**2/15 + Le**2*Phiz/6 + Le**2*Phiz**2/12) / (1+Phiz)**2
    #
    gk[ 1 ][ 5 ] =  Le/10 / (1+Phiz)**2
    gk[ 2 ][ 4 ] = -Le/10 / (1+Phiy)**2
    #
    gk[ 6 ][ 6 ]  = -gk[ 0 ][ 0 ]
    gk[ 7 ][ 7 ]  = -gk[ 1 ][ 1 ]
    gk[ 8 ][ 8 ]  = gk[ 2 ][ 2 ]
    gk[ 9 ][ 9 ]  = gk[ 3 ][ 3 ]
    gk[ 10 ][ 10 ] = gk[ 4 ][ 4 ]
    gk[ 11 ][ 11 ] = gk[ 5 ][ 5 ]
    #
    # impose the geometry
    gk += np.triu(gk, k=1).T
    #
    return -1 * P * gk / Le
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
    if len(Fb) == 6:
        # ['x', 'y', 'rz']
        Fxa, Fya, Mza, Fxb, Fyb, Mzb = Fb
        Fb = np.array([Fxa, Fya, 0.0, 0.0, 0.0, Mza,
                       Fxb, Fyb, 0.0, 0.0, 0.0, Mzb])
    else:
        (Fxa, Fya, Fza, Mxa, Mya, Mza,
         Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb
    #
    P = 1 * Fxb
    #
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
    # Tension
    if fac > toler:
        # Positive axial force
        Kt = Kt_tension(Le=Le, Ax=Ax,
                        Jx=Jx, Iy=Iy, Iz=Iz,
                        Emod=Emod, Gmod=Gmod,
                        Oy=Oy, Oz=Oz, Fb=Fb)
    # Compression
    elif fac < -toler:
        # Negative axial force
        Kt = Kt_compression(Le=Le, Ax=Ax,
                            Jx=Jx, Iy=Iy, Iz=Iz,
                            Emod=Emod, Gmod=Gmod,
                            Oy=Oy, Oz=Oz, Fb=Fb)
    
    else:
        Kt = Kt_NoAxial(Le=Le, Ax=Ax,
                        Jx=Jx, Iy=Iy, Iz=Iz,
                        Emod=Emod, Gmod=Gmod,
                        Oy=Oy, Oz=Oz, Fb=Fb)
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
     Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb
    #
    Fxb = 0
    #
    niy = 1.0 + 12.0 * Oy
    niz = 1.0 + 12.0 * Oz
    lay = 1.0 + 3.0 * Oy
    laz = 1.0 + 3.0 * Oz
    gay = 1.0 - 6.0 * Oy
    gaz = 1.0 - 6.0 * Oz
    #
    L2 = Le*Le
    L3 = L2*Le
    pi = np.pi
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
    et = np.zeros((12, 12))
    #
    #
    et[0 , 0] = et[6 , 6] = Emod*Ax / Le  + Mxb / Le
    et[1 , 1] = et[7 , 7] = (12.0*Emod*Iz / L3)*kvvy
    et[2 , 2] = et[8 , 8] = (12.0*Emod*Iy / L3)*kvvz
    et[3 , 3] = et[9 , 9] = Gmod*Jx / Le
    et[4 , 4] = et[10 , 10] = (4.0*Emod*Iy / Le)*krrz
    et[5 , 5] = et[11 , 11] = (4.0*Emod*Iz / Le)*krry
    et[1 , 5] = et[1 , 11] = (6.0*Emod*Iz / L2)*kvry
    et[2 , 4] = et[2 , 10] = -(6.0*Emod*Iy / L2)*kvrz
    et[3 , 9] = et[9 , 3] = -Gmod*Jx / Le
    et[4 , 10] = et[10 , 4] = (2.0*Emod*Iy / Le)*krtz
    et[5 , 11] = et[11 , 5] = (2.0*Emod*Iz / Le)*krty
    #
    et[4 , 2] = et[10 , 2] = et[2 , 4]
    et[4 , 8] = et[8 , 10] = -et[2 , 4]
    et[5 , 1] = et[11 , 1] = et[1 , 5]
    et[5 , 7] = et[7 , 11] = -et[1 , 5]
    et[6 , 0] = et[0 , 6] = -et[0 , 0]
    et[7 , 1] = et[1 , 7] = -et[1 , 1]
    et[7 , 5] = et[11 , 7] = -et[1 , 5]
    et[8 , 2] = et[2 , 8] = -et[2 , 2]
    et[8 , 4] = et[10 , 8] = -et[2 , 4]

    et[2 , 3] = et[3 , 2] = Mza / Le
    et[3 , 4] = et[4 , 3] = -Mza / 3.0 + Mzb / 6.0
    et[5 , 6] = et[6 , 5] = Mza / Le
    et[8 , 9] = et[9 , 8] = -Mzb / Le
    et[9 , 10] = et[10 , 9] = Mza / 6.0 - Mzb / 3.0

    et[1 , 3] = et[3 , 1] = Mya / Le
    et[3 , 5] = et[5 , 3] = Mya / 3.0 - Myb / 6.0
    et[4 , 6] = et[6 , 4] = Mya / Le
    et[7 , 9] = et[9 , 7] = -Myb / Le
    et[9 , 11] = et[11 , 9] = -Mya / 6.0 + Myb / 3.0

    et[1 , 4] = et[4 , 1] = (Mxb / Le)*(1 / niz)
    et[2 , 5] = et[5 , 2] = (Mxb / Le)*(1 / niy)
    et[4 , 7] = et[7 , 4] = (-Mxb / Le)*(1 / niz)
    et[5 , 8] = et[8 , 5] = (-Mxb / Le)*(1 / niy)
    et[7 , 10] = et[10 , 7] = (Mxb / Le)*(1 / niz)
    et[8 , 11] = et[11 , 8] = (Mxb / Le)*(1 / niy)

    et[0 , 4] = et[4 , 0] = -Mya / Le
    et[3 , 7] = et[7 , 3] = -Mya / Le
    et[5 , 9] = et[9 , 5] = Mya / 6.0 + Myb / 6.0
    et[6 , 10] = et[10 , 6] = Myb / Le

    et[0 , 5] = et[5 , 0] = -Mza / Le
    et[3 , 8] = et[8 , 3] = -Mza / Le
    et[4 , 9] = et[9 , 4] = -Mza / 6.0 - Mzb / 6.0
    et[5 , 10] = et[10 , 5] = (-Mxb / 2.0)*((1 - 144 * Oy*Oz) / (niy*niz))
    et[6 , 11] = et[11 , 6] = Mzb / Le

    et[2 , 9] = et[9 , 2] = Mzb / Le
    et[3 , 10] = et[10 , 3] = -Mza / 6.0 - Mzb / 6.0
    et[4 , 11] = et[11 , 4] = (Mxb / 2.0)*((1 - 144 * Oy*Oz) / (niy*niz))

    et[1 , 9] = et[9 , 1] = Myb / Le
    et[3 , 11] = et[11 , 3] = Mya / 6.0 + Myb / 6.0

    et[1 , 10] = et[10 , 1] = (-Mxb / Le)*(1 / niz)
    et[2 , 11] = et[11 , 2] = (-Mxb / Le)*(1 / niy)

    et[0 , 10] = et[10 , 0] = -Myb / Le

    et[0 , 11] = et[11 , 0] = -Mzb / Le

    et[4 , 5] = et[5 , 4] = -6 * Mxb*(Oy - Oz) / (niy*niz)
    et[10 , 11] = et[11 , 10] = 6 * Mxb*(Oy - Oz) / (niy*niz)
    #
    return et
#
#
def Kt_compression(Le:float, Ax:float,
                   Jx:float, Iy:float, Iz:float,
                   Emod:float, Gmod:float,
                   Oy:float, Oz:float,
                   Fb:list[float]):
    """Negative axial force"""
    (Fxa, Fya, Fza, Mxa, Mya, Mza,
     Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb
    #
    P = Fxb
    A = Ax
    L2 = Le*Le
    #L3 = L2*Le
    #pi = np.pi
    #
    #
    mu_y = np.sqrt(np.abs(P) / (Emod*Iz))
    mu_z = np.sqrt(np.abs(P) / (Emod*Iy))
    Vy = mu_y / (np.sqrt(1 - Oy*mu_y*mu_y*L2))
    Vz = mu_z / (np.sqrt(1 - Oz*mu_z*mu_z*L2))

    K = P*(Iy + Iz) / Ax
    ay = (4.0 * (np.sin(Le * Vy / 2.0) * np.sin(Le * Vy / 2.0)) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (np.sin(Le * Vy / 2.0) * np.sin(Le * Vy / 2.0))
    az = (4.0 * (np.sin(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0)) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (np.sin(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0))
    by = Oy * (Le * Le) * (Vy * Vy) + 1.0
    bz = Oz * (Le * Le) * (Vz * Vz) + 1.0
    cy = (((2.0 * np.cos(Le * Vy) - 2.0 * (Le * Le) * (Vy * Vy) * Oy) + Le * Vy * np.sin(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) - 2.0
    cz = (((2.0 * np.cos(Le * Vz) - 2.0 * (Le * Le) * (Vz * Vz) * Oz) + Le * Vz * np.sin(Le * Vz)) + 2.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) - 2.0
    dy = 4.0 * A * (Le * Le) * (Vy * Vy) * Oy * P
    dz = 4.0 * A * (Le * Le) * (Vz * Vz) * Oz * P
    
    ey = ((((((((((((((((((((2.0 * A * P * np.sin(Le * Vy)
                             - A * P * np.sin(2.0 * Le * Vy))
                            - 2.0 * Iz * (Vy * Vy) * P * np.sin(Le * Vy))
                           + Iz * (Vy * Vy) * P * np.sin(2.0 * Le * Vy)) - 6.0 * A * Le * Vy * P)
                         - Iz * Le * np.power(Vy, 3.0) * P) - A * Emod * Iz * Le * np.power(Vy, 3.0))
                       - 2.0 * A * Emod * Iz * (Vy * Vy) * np.sin(Le * Vy))
                      + A * Emod * Iz * (Vy * Vy) * np.sin(2.0 * Le * Vy))
                     + 2.0 * Iz * Le * np.power(Vy, 3.0) * P * np.cos(Le * Vy))
                    - Iz * Le * np.power(Vy, 3.0) * P * np.cos(2.0 * Le * Vy))
                   - 17.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P)
                  - 3.0 * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * P)
                 - A * np.power(Le, 3.0) * np.power(Vy, 3.0) * P * np.cos(Le * Vy))
                - Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * P * np.cos(Le * Vy))
               + 3.0 * A * (Le * Le) * (Vy * Vy) * P * np.sin(Le * Vy))
              + Iz * (Le * Le) * np.power(Vy, 4.0) * P * np.sin(Le * Vy))
             + 6.0 * A * Le * Vy * P * np.cos(Le * Vy))
            - 18.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P)
           - 9.0 * A * np.power(Le, 7.0) * np.power(Vy, 7.0) * np.power(Oy, 3.0) * P)
          - 2.0 * A * np.power(Le, 9.0) * np.power(Vy, 9.0) * np.power(Oy, 4.0) * P)
    
    ez = ((((((((((((((((((((2.0 * A * P * np.sin(Le * Vz)
                             - A * P * np.sin(2.0 * Le * Vz))
                            - 2.0 * Iy * (Vz * Vz) * P * np.sin(Le * Vz))
                           + Iy * (Vz * Vz) * P * np.sin(2.0 * Le * Vz))
                          - 6.0 * A * Le * Vz * P) - Iy * Le * np.power(Vz, 3.0) * P)
                        - A * Emod * Iy * Le * np.power(Vz, 3.0))
                       - 2.0 * A * Emod * Iy * (Vz * Vz) * np.sin(Le * Vz))
                      + A * Emod * Iy * (Vz * Vz) * np.sin(2.0 * Le * Vz))
                     + 2.0 * Iy * Le * np.power(Vz, 3.0) * P * np.cos(Le * Vz))
                    - Iy * Le * np.power(Vz, 3.0) * P * np.cos(2.0 * Le * Vz))
                   - 17.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P)
                  - 3.0 * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * P)
                 - A * np.power(Le, 3.0) * np.power(Vz, 3.0) * P * np.cos(Le * Vz))
                - Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * P * np.cos(Le * Vz))
               + 3.0 * A * (Le * Le) * (Vz * Vz) * P * np.sin(Le * Vz))
              + Iy * (Le * Le) * np.power(Vz, 4.0) * P * np.sin(Le * Vz))
             + 6.0 * A * Le * Vz * P * np.cos(Le * Vz))
            - 18.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P)
           - 9.0 * A * np.power(Le, 7.0) * np.power(Vz, 7.0) * np.power(Oz, 3.0) * P)
          - 2.0 * A * np.power(Le, 9.0) * np.power(Vz, 9.0) * np.power(Oz, 4.0) * P)
    
    fy = 2.0 * Iz * (Le * Le) * np.power(Vy, 4.0) * Oy
    fz = 2.0 * Iy * (Le * Le) * np.power(Vz, 4.0) * Oz
    sy = np.sin(Le * Vy / 2.0)
    sz = np.sin(Le * Vz / 2.0)
    c_y = np.cos(Le * Vy / 2.0)
    c_z = np.cos(Le * Vz / 2.0)
    #
    #
    et = np.zeros((12, 12))
    #
    #
    et[0 , 0] = et[6 , 6] = P / Le + A * Emod / Le
    
    et[1 , 1] = et[7 , 7] = (((((sy * sy * (2.0 * Le * (Vy * Vy)
                                            * P * (((A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)
                                                     + 2.0 * A * (Le * Le) * (Vy * Vy) * Oy)
                                                    + Iz * (Vy * Vy)) + A)
                                            - Vy * P * np.sin(Le * Vy) * (((-2.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)
                                                                            + 4.0 * A * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Iz * (Vy * Vy)) + 6.0 * A))
                                 + A * Le * (Vy * Vy) * P * (np.sin(Le * Vy) * np.sin(Le * Vy))) / (A * (ay * ay))
                                - 2.0 * Emod * Iz * np.power(Vy, 3.0) * (sy * sy) * (np.sin(Le * Vy) - Le * Vy) / (ay * ay))
                               + Emod * Iz * Vy * ((((((Le * Vy * (np.sin(Le * Vy) * np.sin(Le * Vy))
                                                     - 6.0 * np.sin(Le * Vy) * (sy * sy)) + 2.0 * Le * Vy * (sy * sy))
                                                   + 2.0 * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * (sy * sy))
                                                  + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * (sy * sy))
                                                 - 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(Le * Vy) * (sy * sy))
                                                + 2.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.sin(Le * Vy) * (sy * sy)) / (Le * Le * Oy * (ay * ay)))
                              + 2.0 * Emod * Iz * Vy * ((((2.0 * Le * Vy - 3.0 * np.sin(Le * Vy))
                                                       + np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy)
                                                      + Le * Vy * np.cos(Le * Vy))
                                                     - Le * Le * (Vy * Vy) * Oy * np.sin(Le * Vy))
                              / (Le * Le * Oy * (((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy))
                                                        - 8.0 * (Le * Le) * (Vy * Vy) * Oy)
                                                       - Le * Le * (Vy * Vy) * np.cos(Le * Vy))
                                                      + 4.0 * Le * Vy * np.sin(Le * Vy))
                                                     - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy))
                                                    + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy))
                                                   + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy))
                                                  + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) - 4.0)))
                             + 2.0 * Emod * Iz * Vy * (sy * sy) * ((2.0 * Le * Vy - 3.0 * np.sin(Le * Vy))
                                                                + Le * Vy * np.cos(Le * Vy)) / (Le * Le * Oy * (ay * ay)))
    
    et[2 , 2] = et[8 , 8] = (((((sz * sz * (2.0 * Le * (Vz * Vz) * P * (((A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)
                                                                          + 2.0 * A * (Le * Le) * (Vz * Vz) * Oz)
                                                                         + Iy * (Vz * Vz)) + A) - Vz * P * np.sin(Le * Vz)
                                            * (((-2.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)
                                                 + 4.0 * A * (Le * Le) * (Vz * Vz) * Oz) + 2.0 * Iy * (Vz * Vz))
                                               + 6.0 * A)) + A * Le * (Vz * Vz) * P * (np.sin(Le * Vz) * np.sin(Le * Vz)))
                                / (A * (az * az)) - 2.0 * Emod * Iy * np.power(Vz, 3.0) * (sz * sz) * (np.sin(Le * Vz) - Le * Vz) / (az * az))
                               + Emod * Iy * Vz * ((((((Le * Vz * (np.sin(Le * Vz) * np.sin(Le * Vz))
                                                     - 6.0 * np.sin(Le * Vz) * (sz * sz)) + 2.0 * Le * Vz * (sz * sz))
                                                   + 2.0 * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * (sz * sz))
                                                  + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * (sz * sz))
                                                 - 4.0 * (Le * Le) * (Vz * Vz) * Oz * np.sin(Le * Vz) * (sz * sz))
                                                + 2.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.sin(Le * Vz) * (sz * sz)) / (Le * Le * Oz * (az * az)))
                              + 2.0 * Emod * Iy * Vz * ((((2.0 * Le * Vz - 3.0 * np.sin(Le * Vz))
                                                       + np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz)
                                                      + Le * Vz * np.cos(Le * Vz)) - Le * Le * (Vz * Vz) * Oz * np.sin(Le * Vz))
                              / (Le * Le * Oz * (((((((((4.0 * np.cos(Le * Vz) - Le * Le * (Vz * Vz))
                                                        - 8.0 * (Le * Le) * (Vz * Vz) * Oz)
                                                       - Le * Le * (Vz * Vz) * np.cos(Le * Vz))
                                                      + 4.0 * Le * Vz * np.sin(Le * Vz))
                                                     - 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz))
                                                    + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz))
                                                   + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz))
                                                  + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) - 4.0)))
                             + 2.0 * Emod * Iy * Vz * (sz * sz) * ((2.0 * Le * Vz - 3.0 * np.sin(Le * Vz))
                                                                + Le * Vz * np.cos(Le * Vz)) / (Le * Le * Oz * (az * az)))
    
    et[3 , 3] = et[9 , 9] = K / Le + Gmod*Jx / Le
    
    et[4 , 4] = et[10 , 10] = (((((P * (bz * bz) * (((((((((((Vz * np.sin(2.0 * Le * Vz) / 2.0 + np.power(Le, 3.0) * np.power(Vz, 4.0) / 2.0) + 5.0 * np.power(Le, 3.0) * np.power(Vz, 4.0) * Oz / 2.0) - Le * Le * np.power(Vz, 3.0) * np.sin(2.0 * Le * Vz) / 4.0) + np.power(Le, 5.0) * np.power(Vz, 6.0) * (Oz * Oz)) + Le * (Vz * Vz) * np.cos(Le * Vz)) - Le * (Vz * Vz) * np.cos(2.0 * Le * Vz)) - 2.0 * np.power(Le, 3.0) * np.power(Vz, 4.0) * Oz * np.cos(Le * Vz)) - np.power(Le, 3.0) * np.power(Vz, 4.0) * Oz * np.cos(2.0 * Le * Vz) / 2.0) + Le * Le * np.power(Vz, 3.0) * Oz * np.sin(2.0 * Le * Vz)) - np.power(Le, 5.0) * np.power(Vz, 6.0) * (Oz * Oz) * np.cos(Le * Vz)) + np.power(Le, 4.0) * np.power(Vz, 5.0) * (Oz * Oz) * np.sin(2.0 * Le * Vz) / 2.0) - P * np.sin(Le * Vz) * (bz * bz) * ((((np.power(Le, 4.0) * np.power(Vz, 5.0) * (Oz * Oz) + np.power(Le, 4.0) * np.power(Vz, 5.0) * Oz) + 2.0 * (Le * Le) * np.power(Vz, 3.0) * Oz) + Le * Le * np.power(Vz, 3.0)) + Vz)) / (Vz * Vz * (cz * cz)) - (((((Emod * Iy * (np.sin(Le * Vz) - np.sin(2.0 * Le * Vz) / 2.0) / Vz - Emod * Iy * Le * (np.cos(Le * Vz) - np.cos(2.0 * Le * Vz))) - Emod * Iy * np.power(Le, 3.0) * (Vz * Vz) * (((5.0 * Oz - 4.0 * Oz * np.cos(Le * Vz)) - Oz * np.cos(2.0 * Le * Vz)) + 1.0) / 2.0) + Emod * Iy * (Le * Le) * Vz * (((4.0 * np.sin(Le * Vz) + np.sin(2.0 * Le * Vz)) + 8.0 * Oz * np.sin(Le * Vz)) - 4.0 * Oz * np.sin(2.0 * Le * Vz)) / 4.0) - Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 4.0) * (Oz * Oz) * ((np.cos(2.0 * Le * Vz) / 2.0 - 3.0 * np.cos(Le * Vz)) + 2.5)) + Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 3.0) * Oz * ((2.0 * np.sin(Le * Vz) + 2.0 * Oz * np.sin(Le * Vz)) - Oz * np.sin(2.0 * Le * Vz)) / 2.0) / (Le * Le * Oz * (cz * cz))) - Emod * Iy * (((((((((Le * (Vz * Vz) / 2.0 - np.cos(2.0 * Le * Vz) * (Oz * np.power(Le, 3.0) * np.power(Vz, 4.0) + Le * (Vz * Vz)) / 2.0) + np.power(Le, 3.0) * np.power(Vz, 4.0) * Oz / 2.0) - Le * (Vz * Vz) * (bz * bz) / 2.0) + Vz * np.sin(2.0 * Le * Vz) * (((np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz)) + 1.0) / 4.0) + Le * np.sin(Le * Vz) * (Oz * np.power(Le, 3.0) * np.power(Vz, 5.0) + Le * np.power(Vz, 3.0))) - Le * (Vz * Vz) * (((np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) + Le * Le * (Vz * Vz)) + 1.0) / 2.0) - Vz * np.sin(Le * Vz) * (bz * bz)) + Vz * np.sin(2.0 * Le * Vz) * (bz * bz) / 4.0) + Le * (Vz * Vz) * np.cos(Le * Vz) * (bz * bz)) / (cz * cz)) - Iy * P * (((((((((Le * (Vz * Vz) / 2.0 + np.power(Le, 3.0) * np.power(Vz, 4.0) * Oz / 2.0) - Le * (Vz * Vz) * (bz * bz) / 2.0) + Vz * np.sin(2.0 * Le * Vz) * (((np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz)) + 1.0) / 4.0) + Le * np.sin(Le * Vz) * (Oz * np.power(Le, 3.0) * np.power(Vz, 5.0) + Le * np.power(Vz, 3.0))) - Le * (Vz * Vz) * (((np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) + Le * Le * (Vz * Vz)) + 1.0) / 2.0) - Vz * np.sin(Le * Vz) * (bz * bz)) + Vz * np.sin(2.0 * Le * Vz) * (bz * bz) / 4.0) - Le * (Vz * Vz) * np.cos(2.0 * Le * Vz) * (Oz * (Le * Le) * (Vz * Vz) + 1.0) / 2.0) + Le * (Vz * Vz) * np.cos(Le * Vz) * (bz * bz)) / (A * (cz * cz))) + (Emod * Iy * (bz * bz) * (((((((((((Vz * np.sin(2.0 * Le * Vz) / 2.0 + np.power(Le, 3.0) * np.power(Vz, 4.0) / 2.0) + 5.0 * np.power(Le, 3.0) * np.power(Vz, 4.0) * Oz / 2.0) - Le * Le * np.power(Vz, 3.0) * np.sin(2.0 * Le * Vz) / 4.0) + np.power(Le, 5.0) * np.power(Vz, 6.0) * (Oz * Oz)) + Le * (Vz * Vz) * np.cos(Le * Vz)) - Le * (Vz * Vz) * np.cos(2.0 * Le * Vz)) - 2.0 * np.power(Le, 3.0) * np.power(Vz, 4.0) * Oz * np.cos(Le * Vz)) - np.power(Le, 3.0) * np.power(Vz, 4.0) * Oz * np.cos(2.0 * Le * Vz) / 2.0) + Le * Le * np.power(Vz, 3.0) * Oz * np.sin(2.0 * Le * Vz)) - np.power(Le, 5.0) * np.power(Vz, 6.0) * (Oz * Oz) * np.cos(Le * Vz)) + np.power(Le, 4.0) * np.power(Vz, 5.0) * (Oz * Oz) * np.sin(2.0 * Le * Vz) / 2.0) - Emod * Iy * np.sin(Le * Vz) * (bz * bz) * ((((np.power(Le, 4.0) * np.power(Vz, 5.0) * (Oz * Oz) + np.power(Le, 4.0) * np.power(Vz, 5.0) * Oz) + 2.0 * (Le * Le) * np.power(Vz, 3.0) * Oz) + Le * Le * np.power(Vz, 3.0)) + Vz)) / (Le * Le * (Vz * Vz) * Oz * (cz * cz))) + 2.0 * Emod * Iy * (Oz * (Le * Le) * (Vz * Vz) + 1.0) * ((((((((((((((((2.0 * np.sin(Le * Vz) - np.sin(2.0 * Le * Vz)) - np.power(Le, 3.0) * np.power(Vz, 3.0)) - 5.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz) + 2.0 * (Le * Le) * (Vz * Vz) * np.sin(Le * Vz)) + Le * Le * (Vz * Vz) * np.sin(2.0 * Le * Vz) / 2.0) - 2.0 * Le * Vz * np.cos(Le * Vz)) + 2.0 * Le * Vz * np.cos(2.0 * Le * Vz)) - 2.0 * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.cos(Le * Vz)) + np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.cos(2.0 * Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * np.sin(Le * Vz)) - 2.0 * (Le * Le) * (Vz * Vz) * Oz * np.sin(2.0 * Le * Vz)) + 2.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * Oz * np.sin(Le * Vz)) + 2.0 * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * np.cos(Le * Vz)) + 2.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.sin(Le * Vz)) - np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.sin(2.0 * Le * Vz)) / (Le * Le * Vz * Oz * ((((((((((((((4.0 * np.cos(2.0 * Le * Vz) - 16.0 * np.cos(Le * Vz)) + Le * Le * (Vz * Vz)) + 24.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * np.cos(2.0 * Le * Vz)) - 8.0 * Le * Vz * np.sin(Le * Vz)) + 4.0 * Le * Vz * np.sin(2.0 * Le * Vz)) + 12.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)) - 32.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(2.0 * Le * Vz)) - 8.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(2.0 * Le * Vz)) - 16.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(2.0 * Le * Vz)) + 12.0))
    et[5 , 5] = et[11 , 11] = (((((P * (by * by) * (((((((((((Vy * np.sin(2.0 * Le * Vy) / 2.0 + np.power(Le, 3.0) * np.power(Vy, 4.0) / 2.0) + 5.0 * np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy / 2.0) - Le * Le * np.power(Vy, 3.0) * np.sin(2.0 * Le * Vy) / 4.0) + np.power(Le, 5.0) * np.power(Vy, 6.0) * (Oy * Oy)) + Le * (Vy * Vy) * np.cos(Le * Vy)) - Le * (Vy * Vy) * np.cos(2.0 * Le * Vy)) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy * np.cos(Le * Vy)) - np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy * np.cos(2.0 * Le * Vy) / 2.0) + Le * Le * np.power(Vy, 3.0) * Oy * np.sin(2.0 * Le * Vy)) - np.power(Le, 5.0) * np.power(Vy, 6.0) * (Oy * Oy) * np.cos(Le * Vy)) + np.power(Le, 4.0) * np.power(Vy, 5.0) * (Oy * Oy) * np.sin(2.0 * Le * Vy) / 2.0) - P * np.sin(Le * Vy) * (by * by) * ((((np.power(Le, 4.0) * np.power(Vy, 5.0) * (Oy * Oy) + np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy) + 2.0 * (Le * Le) * np.power(Vy, 3.0) * Oy) + Le * Le * np.power(Vy, 3.0)) + Vy)) / (Vy * Vy * (cy * cy)) - (((((Emod * Iz * (np.sin(Le * Vy) - np.sin(2.0 * Le * Vy) / 2.0) / Vy - Emod * Iz * Le * (np.cos(Le * Vy) - np.cos(2.0 * Le * Vy))) - Emod * Iz * np.power(Le, 3.0) * (Vy * Vy) * (((5.0 * Oy - 4.0 * Oy * np.cos(Le * Vy)) - Oy * np.cos(2.0 * Le * Vy)) + 1.0) / 2.0) + Emod * Iz * (Le * Le) * Vy * (((4.0 * np.sin(Le * Vy) + np.sin(2.0 * Le * Vy)) + 8.0 * Oy * np.sin(Le * Vy)) - 4.0 * Oy * np.sin(2.0 * Le * Vy)) / 4.0) - Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 4.0) * (Oy * Oy) * ((np.cos(2.0 * Le * Vy) / 2.0 - 3.0 * np.cos(Le * Vy)) + 2.5)) + Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 3.0) * Oy * ((2.0 * np.sin(Le * Vy) + 2.0 * Oy * np.sin(Le * Vy)) - Oy * np.sin(2.0 * Le * Vy)) / 2.0) / (Le * Le * Oy * (cy * cy))) - Emod * Iz * (((((((((Le * (Vy * Vy) / 2.0 - np.cos(2.0 * Le * Vy) * (Oy * np.power(Le, 3.0) * np.power(Vy, 4.0) + Le * (Vy * Vy)) / 2.0) + np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy / 2.0) - Le * (Vy * Vy) * (by * by) / 2.0) + Vy * np.sin(2.0 * Le * Vy) * (((np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy)) + 1.0) / 4.0) + Le * np.sin(Le * Vy) * (Oy * np.power(Le, 3.0) * np.power(Vy, 5.0) + Le * np.power(Vy, 3.0))) - Le * (Vy * Vy) * (((np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) + Le * Le * (Vy * Vy)) + 1.0) / 2.0) - Vy * np.sin(Le * Vy) * (by * by)) + Vy * np.sin(2.0 * Le * Vy) * (by * by) / 4.0) + Le * (Vy * Vy) * np.cos(Le * Vy) * (by * by)) / (cy * cy)) - Iz * P * (((((((((Le * (Vy * Vy) / 2.0 + np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy / 2.0) - Le * (Vy * Vy) * (by * by) / 2.0) + Vy * np.sin(2.0 * Le * Vy) * (((np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy)) + 1.0) / 4.0) + Le * np.sin(Le * Vy) * (Oy * np.power(Le, 3.0) * np.power(Vy, 5.0) + Le * np.power(Vy, 3.0))) - Le * (Vy * Vy) * (((np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) + Le * Le * (Vy * Vy)) + 1.0) / 2.0) - Vy * np.sin(Le * Vy) * (by * by)) + Vy * np.sin(2.0 * Le * Vy) * (by * by) / 4.0) - Le * (Vy * Vy) * np.cos(2.0 * Le * Vy) * (Oy * (Le * Le) * (Vy * Vy) + 1.0) / 2.0) + Le * (Vy * Vy) * np.cos(Le * Vy) * (by * by)) / (A * (cy * cy))) + (Emod * Iz * (by * by) * (((((((((((Vy * np.sin(2.0 * Le * Vy) / 2.0 + np.power(Le, 3.0) * np.power(Vy, 4.0) / 2.0) + 5.0 * np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy / 2.0) - Le * Le * np.power(Vy, 3.0) * np.sin(2.0 * Le * Vy) / 4.0) + np.power(Le, 5.0) * np.power(Vy, 6.0) * (Oy * Oy)) + Le * (Vy * Vy) * np.cos(Le * Vy)) - Le * (Vy * Vy) * np.cos(2.0 * Le * Vy)) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy * np.cos(Le * Vy)) - np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy * np.cos(2.0 * Le * Vy) / 2.0) + Le * Le * np.power(Vy, 3.0) * Oy * np.sin(2.0 * Le * Vy)) - np.power(Le, 5.0) * np.power(Vy, 6.0) * (Oy * Oy) * np.cos(Le * Vy)) + np.power(Le, 4.0) * np.power(Vy, 5.0) * (Oy * Oy) * np.sin(2.0 * Le * Vy) / 2.0) - Emod * Iz * np.sin(Le * Vy) * (by * by) * ((((np.power(Le, 4.0) * np.power(Vy, 5.0) * (Oy * Oy) + np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy) + 2.0 * (Le * Le) * np.power(Vy, 3.0) * Oy) + Le * Le * np.power(Vy, 3.0)) + Vy)) / (Le * Le * (Vy * Vy) * Oy * (cy * cy))) + 2.0 * Emod * Iz * (Oy * (Le * Le) * (Vy * Vy) + 1.0) * ((((((((((((((((2.0 * np.sin(Le * Vy) - np.sin(2.0 * Le * Vy)) - np.power(Le, 3.0) * np.power(Vy, 3.0)) - 5.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy) + 2.0 * (Le * Le) * (Vy * Vy) * np.sin(Le * Vy)) + Le * Le * (Vy * Vy) * np.sin(2.0 * Le * Vy) / 2.0) - 2.0 * Le * Vy * np.cos(Le * Vy)) + 2.0 * Le * Vy * np.cos(2.0 * Le * Vy)) - 2.0 * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(Le * Vy)) + np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(2.0 * Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(Le * Vy)) - 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(2.0 * Le * Vy)) + 2.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * np.sin(Le * Vy)) + 2.0 * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * np.cos(Le * Vy)) + 2.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.sin(Le * Vy)) - np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.sin(2.0 * Le * Vy)) / (Le * Le * Vy * Oy * ((((((((((((((4.0 * np.cos(2.0 * Le * Vy) - 16.0 * np.cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 24.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * np.cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * np.sin(Le * Vy)) + 4.0 * Le * Vy * np.sin(2.0 * Le * Vy)) + 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)) - 32.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(2.0 * Le * Vy)) - 8.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(2.0 * Le * Vy)) - 16.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(2.0 * Le * Vy)) + 12.0))

    et[2 , 3] = et[3 , 2] = Mza / Le
    et[3 , 4] = et[4 , 3] = (Oz * (Le * Le) * (Vz * Vz) + 1.0) * (Mza + Mzb) * ((((((Le * Le * (Vz * Vz) + 8.0 * (sz * sz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * (2.0 * (sz * sz) - 1.0)) - 4.0 * Le * Vz * np.sin(Le * Vz)) - 2.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (2.0 * (sz * sz) - 1.0)) / (Le * Le * (Vz * Vz) * (((4.0 * (sz * sz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Vz * np.sin(Le * Vz)) + 2.0 * (Le * Le) * (Vz * Vz) * Oz * (2.0 * (sz * sz) - 1.0))) - Mza / 2.0
    et[5 , 6] = et[6 , 5] = Mza / Le
    et[8 , 9] = et[9 , 8] = -Mzb / Le
    et[9 , 10] = et[10 , 9] = (Oz * (Le * Le) * (Vz * Vz) + 1.0) * (Mza + Mzb) * ((((((Le * Le * (Vz * Vz) + 8.0 * (sz * sz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * (2.0 * (sz * sz) - 1.0)) - 4.0 * Le * Vz * np.sin(Le * Vz)) - 2.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (2.0 * (sz * sz) - 1.0)) / (Le * Le * (Vz * Vz) * (((4.0 * (sz * sz) + 2.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Vz * np.sin(Le * Vz)) + 2.0 * (Le * Le) * (Vz * Vz) * Oz * (2.0 * (sz * sz) - 1.0))) - Mzb / 2.0

    et[1 , 3] = et[3 , 1] = Mya / Le
    et[2 , 4] = et[4 , 2] = ((P * (bz * bz) * (((4.0 * np.cos(Le * Vz) + Le * Le * (Vz * Vz)) + Le * Vz * np.sin(Le * Vz)) - 4.0) / (2.0 * (((((((((4.0 * np.cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * np.cos(Le * Vz)) + 4.0 * Le * Vz * np.sin(Le * Vz)) - 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) - 4.0)) + Emod * Iy * Le * np.power(Vz, 3.0) * (sz * sz) * (np.sin(Le * Vz) - Le * Vz) / (az * az)) + Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * (np.sin(Le * Vz) + Le * Vz) / (2.0 * (((((((((4.0 * np.cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * np.cos(Le * Vz)) + 4.0 * Le * Vz * np.sin(Le * Vz)) - 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) - 4.0))) + Iy * Le * np.power(Vz, 3.0) * P * (sz * sz) * (np.sin(Le * Vz) - Le * Vz) / (A * (az * az))
    et[3 , 5] = et[5 , 3] = Mya / 2.0 - (Oy * (Le * Le) * (Vy * Vy) + 1.0) * (Mya + Myb) * ((((((Le * Le * (Vy * Vy) + 8.0 * (sy * sy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * (2.0 * (sy * sy) - 1.0)) - 4.0 * Le * Vy * np.sin(Le * Vy)) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (2.0 * (sy * sy) - 1.0)) / (Le * Le * (Vy * Vy) * (((4.0 * (sy * sy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Vy * np.sin(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * (2.0 * (sy * sy) - 1.0)))
    et[4 , 6] = et[6 , 4] = Mya / Le
    et[5 , 7] = et[7 , 5] = ((P * (by * by) * (((4.0 * np.cos(Le * Vy) + Le * Le * (Vy * Vy)) + Le * Vy * np.sin(Le * Vy)) - 4.0) / (2.0 * (((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * np.cos(Le * Vy)) + 4.0 * Le * Vy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) - 4.0)) + Emod * Iz * Le * np.power(Vy, 3.0) * (sy * sy) * (np.sin(Le * Vy) - Le * Vy) / (ay * ay)) + Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * (np.sin(Le * Vy) + Le * Vy) / (2.0 * (((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * np.cos(Le * Vy)) + 4.0 * Le * Vy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) - 4.0))) + Iz * Le * np.power(Vy, 3.0) * P * (sy * sy) * (np.sin(Le * Vy) - Le * Vy) / (A * (ay * ay))
    et[7 , 9] = et[9 , 7] = -Myb / Le
    et[8 , 10] = et[10 , 8] = ((-(P * (bz * bz) * (((4.0 * np.cos(Le * Vz) + Le * Le * (Vz * Vz)) + Le * Vz * np.sin(Le * Vz)) - 4.0)) / (2.0 * (((((((((4.0 * np.cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * np.cos(Le * Vz)) + 4.0 * Le * Vz * np.sin(Le * Vz)) - 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) - 4.0)) - Emod * Iy * Le * np.power(Vz, 3.0) * (sz * sz) * (np.sin(Le * Vz) - Le * Vz) / (az * az)) - Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * (np.sin(Le * Vz) + Le * Vz) / (2.0 * (((((((((4.0 * np.cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * np.cos(Le * Vz)) + 4.0 * Le * Vz * np.sin(Le * Vz)) - 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) - 4.0))) - Iy * Le * np.power(Vz, 3.0) * P * (sz * sz) * (np.sin(Le * Vz) - Le * Vz) / (A * (az * az))
    et[9 , 11] = et[11 , 9] = Myb / 2.0 - (Oy * (Le * Le) * (Vy * Vy) + 1.0) * (Mya + Myb) * ((((((Le * Le * (Vy * Vy) + 8.0 * (sy * sy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * (2.0 * (sy * sy) - 1.0)) - 4.0 * Le * Vy * np.sin(Le * Vy)) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (2.0 * (sy * sy) - 1.0)) / (Le * Le * (Vy * Vy) * (((4.0 * (sy * sy) + 2.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Vy * np.sin(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * (2.0 * (sy * sy) - 1.0)))

    et[1 , 4] = et[4 , 1] = (Vy * Vy * ((4.0 * Vz * Mxb * (np.sin(Le * Vz) - 2.0 * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0)) - 4.0 * Le * (Vz * Vz) * Mxb * (c_z * c_z - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * np.power(Vz, 3.0) * Mxb * (Oz * np.sin(Le * Vz) - 2.0 * Oz * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0))) - np.power(Vy, 3.0) * ((4.0 * Mxb * (np.sin(Le * Vy) - 2.0 * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vz * Vz) * Mxb * (Oz * np.sin(Le * Vy) - 2.0 * Oz * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0))) - Le * Vz * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    et[2 , 5] = et[5 , 2] = -(Vz * Vz * ((4.0 * Vy * Mxb * (np.sin(Le * Vy) - 2.0 * np.cos(Le * Vy / 2.0) * (c_y * c_y) * np.sin(Le * Vy / 2.0)) - 4.0 * Le * (Vy * Vy) * Mxb * (c_y * c_y - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * np.power(Vy, 3.0) * Mxb * (Oy * np.sin(Le * Vy) - 2.0 * Oy * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0))) - np.power(Vz, 3.0) * ((4.0 * Mxb * (np.sin(Le * Vz) - 2.0 * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Mxb * (Oy * np.sin(Le * Vz) - 2.0 * Oy * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0))) - Le * Vy * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))       
    et[4 , 7] = et[7 , 4] = -(Vy * Vy * ((4.0 * Vz * Mxb * (np.sin(Le * Vz) - 2.0 * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0)) - 4.0 * Le * (Vz * Vz) * Mxb * (c_z * c_z - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * np.power(Vz, 3.0) * Mxb * (Oz * np.sin(Le * Vz) - 2.0 * Oz * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0))) - np.power(Vy, 3.0) * ((4.0 * Mxb * (np.sin(Le * Vy) - 2.0 * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vz * Vz) * Mxb * (Oz * np.sin(Le * Vy) - 2.0 * Oz * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0))) - Le * Vz * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    et[5 , 8] = et[8 , 5] = (Vz * Vz * ((4.0 * Vy * Mxb * (np.sin(Le * Vy) - 2.0 * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0)) - 4.0 * Le * (Vy * Vy) * Mxb * (c_y * c_y - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * np.power(Vy, 3.0) * Mxb * (Oy * np.sin(Le * Vy) - 2.0 * Oy * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0))) - np.power(Vz, 3.0) * ((4.0 * Mxb * (np.sin(Le * Vz) - 2.0 * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Mxb * (Oy * np.sin(Le * Vz) - 2.0 * Oy * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0))) - Le * Vy * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    et[7 , 10] = et[10 , 7] = (Vy * Vy * ((4.0 * Vz * Mxb * (np.sin(Le * Vz) - 2.0 * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0)) - 4.0 * Le * (Vz * Vz) * Mxb * (c_z * c_z - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * np.power(Vz, 3.0) * Mxb * (Oz * np.sin(Le * Vz) - 2.0 * Oz * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0))) - np.power(Vy, 3.0) * ((4.0 * Mxb * (np.sin(Le * Vy) - 2.0 * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vz * Vz) * Mxb * (Oz * np.sin(Le * Vy) - 2.0 * Oz * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0))) - Le * Vz * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    et[8 , 11] = et[11 , 8] = -(Vz * Vz * ((4.0 * Vy * Mxb * (np.sin(Le * Vy) - 2.0 * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0)) - 4.0 * Le * (Vy * Vy) * Mxb * (c_y * c_y - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * np.power(Vy, 3.0) * Mxb * (Oy * np.sin (Le * Vy) - 2.0 * Oy * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0))) - np.power(Vz, 3.0) * ((4.0 * Mxb * (np.sin(Le * Vz) - 2.0 * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Mxb * (Oy * np.sin(Le * Vz) - 2.0 * Oy * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0))) - Le * Vy * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))

    et[0 , 4] = et[4 , 0] = -Mya / Le
    et[1 , 5] = et[5 , 1] = ((-(P * (by * by) * (((4.0 * np.cos(Le * Vy) + Le * Le * (Vy * Vy)) + Le * Vy * np.sin(Le * Vy)) - 4.0)) / (2.0 * (((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * np.cos(Le * Vy)) + 4.0 * Le * Vy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) - 4.0)) - Emod * Iz * Le * np.power(Vy, 3.0) * (sy * sy) * (np.sin(Le * Vy) - Le * Vy) / (ay * ay)) - Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * (np.sin(Le * Vy) + Le * Vy) / (2.0 * (((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * np.cos(Le * Vy)) + 4.0 * Le * Vy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) - 4.0))) - Iz * Le * np.power(Vy, 3.0) * P * (sy * sy) * (np.sin(Le * Vy) - Le * Vy) / (A * (ay * ay))
    et[3 , 7] = et[7 , 3] = -Mya / Le
    et[4 , 8] = et[8 , 4] = ((-(P * (bz * bz) * (((4.0 * np.cos(Le * Vz) + Le * Le * (Vz * Vz)) + Le * Vz * np.sin(Le * Vz)) - 4.0)) / (2.0 * (((((((((4.0 * np.cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * np.cos(Le * Vz)) + 4.0 * Le * Vz * np.sin(Le * Vz)) - 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) - 4.0)) - Emod * Iy * Le * np.power(Vz, 3.0) * (sz * sz) * (np.sin(Le * Vz) - Le * Vz) / (az * az)) - Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * (np.sin(Le * Vz) + Le * Vz) / (2.0 * (((((((((4.0 * np.cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * np.cos(Le * Vz)) + 4.0 * Le * Vz * np.sin(Le * Vz)) - 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) - 4.0))) - Iy * Le * np.power(Vz, 3.0) * P * (sz * sz) * (np.sin(Le * Vz) - Le * Vz) / (A * (az * az))
    et[5 , 9] = et[9 , 5] = 2.0 * (Oy * (Le * Le) * (Vy * Vy) + 1.0) * (Mya + Myb) * (((((Le * Le * (Vy * Vy) + 4.0 * (sy * sy)) - 2.0 * Le * Vy * np.sin(Le * Vy)) - Le * Le * (Vy * Vy) * (sy * sy)) - np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) / (Le * Le * (Vy * Vy) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)))
    et[6 , 10] = et[10 , 6] = Myb / Le
    et[7 , 11] = et[11 , 7] = ((P * (by * by) * (((4.0 * np.cos(Le * Vy) + Le * Le * (Vy * Vy)) + Le * Vy * np.sin(Le * Vy)) - 4.0) / (2.0 * (((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * np.cos(Le * Vy)) + 4.0 * Le * Vy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) - 4.0)) + Emod * Iz * Le * np.power(Vy, 3.0) * (sy * sy) * (np.sin(Le * Vy) - Le * Vy) / (ay * ay)) + Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * (np.sin(Le * Vy) + Le * Vy) / (2.0 * (((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * np.cos(Le * Vy)) + 4.0 * Le * Vy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) - 4.0))) + Iz * Le * np.power(Vy, 3.0) * P * (sy * sy) * (np.sin(Le * Vy) - Le * Vy) / (A * (ay * ay))

    et[0 , 5] = et[5 , 0] = -Mza / Le
    et[3 , 8] = et[8 , 3] = -Mza / Le
    et[4 , 9] = et[9 , 4] = -(2.0 * (Oz * (Le * Le) * (Vz * Vz) + 1.0) * (Mza + Mzb) * (((((Le * Le * (Vz * Vz) + 4.0 * (sz * sz)) - 2.0 * Le * Vz * np.sin(Le * Vz)) - Le * Le * (Vz * Vz) * (sz * sz)) - np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz))) / (Le * Le * (Vz * Vz) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    et[5 , 10] = et[10 , 5] = 2.0 * Mxb * (((((((((((((((4.0 * (Vy * Vy) * (sy * sy) * (sz * sz) - 4.0 * (Vz * Vz) * (sy * sy) * (sz * sz)) + Le * Le * (Vy * Vy) * (Vz * Vz) * (sy * sy)) - Le * Le * (Vy * Vy) * (Vz * Vz) * (sz * sz)) + 4.0 * (Le * Le) * np.power(Vy, 4.0) * Oy * (sy * sy) * (sz * sz)) - 4.0 * (Le * Le) * np.power(Vz, 4.0) * Oz * (sy * sy) * (sz * sz)) + 2.0 * Le * Vy * (Vz * Vz) * np.sin(Le * Vy) * (sz * sz)) - 2.0 * Le * (Vy * Vy) * Vz * np.sin(Le * Vz) * (sy * sy)) - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy * (sy * sy) * (sz * sz)) + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz * (sy * sy) * (sz * sz)) - np.power(Le, 3.0) * np.power(Vy, 4.0) * Vz * Oy * np.sin(Le * Vz) * (sy * sy)) + np.power(Le, 3.0) * Vy * np.power(Vz, 4.0) * Oz * np.sin(Le * Vy) * (sz * sz)) + np.power(Le, 3.0) * np.power(Vy, 3.0) * (Vz * Vz) * Oy * np.sin(Le * Vy) * (sz * sz)) - np.power(Le, 3.0) * (Vy * Vy) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz) * (sy * sy)) - 4.0 * np.power(Le, 4.0) * (Vy * Vy) * np.power(Vz, 4.0) * Oy * Oz * (sy * sy) * (sz * sz)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (sy * sy) * (sz * sz)) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    et[6 , 11] = et[11 , 6] = Mzb / Le

    et[0 , 6] = et[6 , 0] = -P / Le - A * Emod / Le  
    et[1 , 7] = et[7 , 1] = (((2.0 * Emod * Iz * np.power(Vy, 3.0) * (sy * sy) * (np.sin(Le * Vy) - Le * Vy) / (ay * ay) - (sy * sy * (2.0 * Le * (Vy * Vy) * P * (((A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) + 2.0 * A * (Le * Le) * (Vy * Vy) * Oy) + Iz * (Vy * Vy)) + A) - Vy * P * np.sin(Le * Vy) * (((-2.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) + 4.0 * A * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Iz * (Vy * Vy)) + 6.0 * A)) + A * Le * (Vy * Vy) * P * (np.sin(Le * Vy) * np.sin(Le * Vy))) / (A * (ay * ay))) - Emod * Iz * Vy * ((((((Le * Vy * (np.sin(Le * Vy) * np.sin(Le * Vy)) - 6.0 * np.sin(Le * Vy) * (sy * sy)) + 2.0 * Le * Vy * (sy * sy)) + 2.0 * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * (sy * sy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * (sy * sy)) - 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(Le * Vy) * (sy * sy)) + 2.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.sin(Le * Vy) * (sy * sy)) / (Le * Le * Oy * (ay * ay))) - 2.0 * Emod * Iz * Vy * ((((2.0 * Le * Vy - 3.0 * np.sin(Le * Vy)) + np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy) + Le * Vy * np.cos(Le * Vy)) - Le * Le * (Vy * Vy) * Oy * np.sin(Le * Vy)) / (Le * Le * Oy * (((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * np.cos(Le * Vy)) + 4.0 * Le * Vy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) - 4.0))) - 2.0 * Emod * Iz * Vy * (sy * sy) * ((2.0 * Le * Vy - 3.0 * np.sin(Le * Vy)) + Le * Vy * np.cos(Le * Vy)) / (Le * Le * Oy * (ay * ay))
    et[2 , 8] = et[8 , 2] = (((2.0 * Emod * Iy * np.power(Vz, 3.0) * (sz * sz) * (np.sin(Le * Vz) - Le * Vz) / (az * az) - (sz * sz * (2.0 * Le * (Vz * Vz) * P * (((A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) + 2.0 * A * (Le * Le) * (Vz * Vz) * Oz) + Iy * (Vz * Vz)) + A) - Vz * P * np.sin(Le * Vz) * (((-2.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) + 4.0 * A * (Le * Le) * (Vz * Vz) * Oz) + 2.0 * Iy * (Vz * Vz)) + 6.0 * A)) + A * Le * (Vz * Vz) * P * (np.sin(Le * Vz) * np.sin(Le * Vz))) / (A * (az * az))) - Emod * Iy * Vz * ((((((Le * Vz * (np.sin(Le * Vz) * np.sin(Le * Vz)) - 6.0 * np.sin(Le * Vz) * (sz * sz)) + 2.0 * Le * Vz * (sz * sz)) + 2.0 * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * (sz * sz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * (sz * sz)) - 4.0 * (Le * Le) * (Vz * Vz) * Oz * np.sin(Le * Vz) * (sz * sz)) + 2.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.sin(Le * Vz) * (sz * sz)) / (Le * Le * Oz * (az * az))) - 2.0 * Emod * Iy * Vz * ((((2.0 * Le * Vz - 3.0 * np.sin(Le * Vz)) + np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz) + Le * Vz * np.cos(Le * Vz)) - Le * Le * (Vz * Vz) * Oz * np.sin(Le * Vz)) / (Le * Le * Oz * (((((((((4.0 * np.cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * np.cos(Le * Vz)) + 4.0 * Le * Vz * np.sin(Le * Vz)) - 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) - 4.0))) - 2.0 * Emod * Iy * Vz * (sz * sz) * ((2.0 * Le * Vz - 3.0 * np.sin(Le * Vz)) + Le * Vz * np.cos(Le * Vz)) / (Le * Le * Oz * (az * az))
    et[3 , 9] = et[9 , 3] = -(K / Le) - (Gmod*Jx / Le)
    et[4 , 10] = et[10 , 4] = (((((((((((((((((((((((((((((((((((((((((((((((((ez - 2.0 * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz) * P) + 2.0 * dz * np.sin(Le * Vz)) - dz * np.sin(2.0 * Le * Vz)) + 8.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * Oz * P * np.sin(Le * Vz)) - 2.0 * fz * P * np.sin(Le * Vz)) + fz * P * np.sin(2.0 * Le * Vz)) + 2.0 * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * Oz * P * np.sin(Le * Vz)) - 7.0 * A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz)) - 2.0 * A * Emod * Iy * np.power(Le, 7.0) * np.power(Vz, 9.0) * np.power(Oz, 3.0)) + 16.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P * np.cos(Le * Vz)) + 2.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P * np.cos(2.0 * Le * Vz)) - A * np.power(Le, 7.0) * np.power(Vz, 7.0) * (Oz * Oz) * P * np.cos(Le * Vz)) + 8.0 * A * np.power(Le, 7.0) * np.power(Vz, 7.0) * np.power(Oz, 3.0) * P * np.cos(Le * Vz)) + A * np.power(Le, 7.0) * np.power(Vz, 7.0) * np.power(Oz, 3.0) * P * np.cos(2.0 * Le * Vz)) + 2.0 * A * np.power(Le, 9.0) * np.power(Vz, 9.0) * np.power(Oz, 4.0) * P * np.cos(Le * Vz)) + 2.0 * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz) * P * np.cos(Le * Vz)) + 12.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * np.sin(Le * Vz)) - 6.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * np.sin(2.0 * Le * Vz)) + 7.0 * A * np.power(Le, 6.0) * np.power(Vz, 6.0) * (Oz * Oz) * P * np.sin(Le * Vz)) + 8.0 * A * np.power(Le, 6.0) * np.power(Vz, 6.0) * np.power(Oz, 3.0) * P * np.sin(Le * Vz)) - 4.0 * A * np.power(Le, 6.0) * np.power(Vz, 6.0) * np.power(Oz, 3.0) * P * np.sin(2.0 * Le * Vz)) + 2.0 * A * np.power(Le, 8.0) * np.power(Vz, 8.0) * np.power(Oz, 3.0) * P * np.sin(Le * Vz)) + 2.0 * A * np.power(Le, 8.0) * np.power(Vz, 8.0) * np.power(Oz, 4.0) * P * np.sin(Le * Vz)) - A * np.power(Le, 8.0) * np.power(Vz, 8.0) * np.power(Oz, 4.0) * P * np.sin(2.0 * Le * Vz)) - 2.0 * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * P * np.sin(Le * Vz)) + Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * P * np.sin(2.0 * Le * Vz)) + 2.0 * A * Emod * Iy * Le * np.power(Vz, 3.0) * np.cos(Le * Vz)) - A * Emod * Iy * Le * np.power(Vz, 3.0) * np.cos(2.0 * Le * Vz)) - 6.0 * A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz) - A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * np.cos(Le * Vz)) + A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * np.sin(Le * Vz)) + 16.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P * np.cos(Le * Vz)) + A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P * np.cos(2.0 * Le * Vz)) - 2.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * Oz * P * np.cos(Le * Vz)) + 4.0 * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * P * np.cos(Le * Vz)) - Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * P * np.cos(2.0 * Le * Vz)) + 6.0 * A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * np.cos(Le * Vz)) - A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * Oz * np.cos(Le * Vz)) - 2.0 * A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * Oz * np.sin(Le * Vz)) + A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * Oz * np.sin(2.0 * Le * Vz)) + 5.0 * A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * Oz * np.sin(Le * Vz)) + 6.0 * A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz) * np.cos(Le * Vz)) + A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz) * np.cos(2.0 * Le * Vz)) + 2.0 * A * Emod * Iy * np.power(Le, 7.0) * np.power(Vz, 9.0) * np.power(Oz, 3.0) * np.cos(Le * Vz)) + 2.0 * A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * np.sin(Le * Vz)) - A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * np.sin(2.0 * Le * Vz)) + 2.0 * A * Emod * Iy * np.power(Le, 6.0) * np.power(Vz, 8.0) * (Oz * Oz) * np.sin(Le * Vz)) + 2.0 * A * Emod * Iy * np.power(Le, 6.0) * np.power(Vz, 8.0) * np.power(Oz, 3.0) * np.sin(Le * Vz)) - A * Emod * Iy * np.power(Le, 6.0) * np.power(Vz, 8.0) * np.power(Oz, 3.0) * np.sin(2.0 * Le * Vz)) / (A * Vz * ((((((((((((((4.0 * np.cos(2.0 * Le * Vz) - 16.0 * np.cos(Le * Vz)) + Le * Le * (Vz * Vz)) + 24.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * np.cos(2.0 * Le * Vz)) - 8.0 * Le * Vz * np.sin(Le * Vz)) + 4.0 * Le * Vz * np.sin(2.0 * Le * Vz)) + 12.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)) - 32.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(2.0 * Le * Vz)) - 8.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(2.0 * Le * Vz)) - 16.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(2.0 * Le * Vz)) + 12.0))    
    et[5 , 11] = et[11 , 5] = (((((((((((((((((((((((((((((((((((((((((((((((((ey - 2.0 * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy) * P) + 2.0 * dy * np.sin(Le * Vy)) - dy * np.sin(2.0 * Le * Vy)) + 8.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * P * np.sin(Le * Vy)) - 2.0 * fy * P * np.sin(Le * Vy)) + fy * P * np.sin(2.0 * Le * Vy)) + 2.0 * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * Oy * P * np.sin(Le * Vy)) - 7.0 * A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy)) - 2.0 * A * Emod * Iz * np.power(Le, 7.0) * np.power(Vy, 9.0) * np.power(Oy, 3.0)) + 16.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P * np.cos(Le * Vy)) + 2.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P * np.cos(2.0 * Le * Vy)) - A * np.power(Le, 7.0) * np.power(Vy, 7.0) * (Oy * Oy) * P * np.cos(Le * Vy)) + 8.0 * A * np.power(Le, 7.0) * np.power(Vy, 7.0) * np.power(Oy, 3.0) * P * np.cos(Le * Vy)) + A * np.power(Le, 7.0) * np.power(Vy, 7.0) * np.power(Oy, 3.0) * P * np.cos(2.0 * Le * Vy)) + 2.0 * A * np.power(Le, 9.0) * np.power(Vy, 9.0) * np.power(Oy, 4.0) * P * np.cos(Le * Vy)) + 2.0 * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy) * P * np.cos(Le * Vy)) + 12.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * np.sin(Le * Vy)) - 6.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * np.sin(2.0 * Le * Vy)) + 7.0 * A * np.power(Le, 6.0) * np.power(Vy, 6.0) * (Oy * Oy) * P * np.sin(Le * Vy)) + 8.0 * A * np.power(Le, 6.0) * np.power(Vy, 6.0) * np.power(Oy, 3.0) * P * np.sin(Le * Vy)) - 4.0 * A * np.power(Le, 6.0) * np.power(Vy, 6.0) * np.power(Oy, 3.0) * P * np.sin(2.0 * Le * Vy)) + 2.0 * A * np.power(Le, 8.0) * np.power(Vy, 8.0) * np.power(Oy, 3.0) * P * np.sin(Le * Vy)) + 2.0 * A * np.power(Le, 8.0) * np.power(Vy, 8.0) * np.power(Oy, 4.0) * P * np.sin(Le * Vy)) - A * np.power(Le, 8.0) * np.power(Vy, 8.0) * np.power(Oy, 4.0) * P * np.sin(2.0 * Le * Vy)) - 2.0 * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * P * np.sin(Le * Vy)) + Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * P * np.sin(2.0 * Le * Vy)) + 2.0 * A * Emod * Iz * Le * np.power(Vy, 3.0) * np.cos(Le * Vy)) - A * Emod * Iz * Le * np.power(Vy, 3.0) * np.cos(2.0 * Le * Vy)) - 6.0 * A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy) - A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * np.cos(Le * Vy)) + A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * np.sin(Le * Vy)) + 16.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P * np.cos(Le * Vy)) + A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P * np.cos(2.0 * Le * Vy)) - 2.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * Oy * P * np.cos(Le * Vy)) + 4.0 * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * P * np.cos(Le * Vy)) - Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * P * np.cos(2.0 * Le * Vy)) + 6.0 * A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * np.cos(Le * Vy)) - A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * Oy * np.cos(Le * Vy)) - 2.0 * A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * Oy * np.sin(Le * Vy)) + A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * Oy * np.sin(2.0 * Le * Vy)) + 5.0 * A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * Oy * np.sin(Le * Vy)) + 6.0 * A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy) * np.cos(Le * Vy)) + A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy) * np.cos(2.0 * Le * Vy)) + 2.0 * A * Emod * Iz * np.power(Le, 7.0) * np.power(Vy, 9.0) * np.power(Oy, 3.0) * np.cos(Le * Vy)) + 2.0 * A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * np.sin(Le * Vy)) - A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * np.sin(2.0 * Le * Vy)) + 2.0 * A * Emod * Iz * np.power(Le, 6.0) * np.power(Vy, 8.0) * (Oy * Oy) * np.sin(Le * Vy)) + 2.0 * A * Emod * Iz * np.power(Le, 6.0) * np.power(Vy, 8.0) * np.power(Oy, 3.0) * np.sin(Le * Vy)) - A * Emod * Iz * np.power(Le, 6.0) * np.power(Vy, 8.0) * np.power(Oy, 3.0) * np.sin(2.0 * Le * Vy)) / (A * Vy * ((((((((((((((4.0 * np.cos(2.0 * Le * Vy) - 16.0 * np.cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 24.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * np.cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * np.sin(Le * Vy)) + 4.0 * Le * Vy * np.sin(2.0 * Le * Vy)) + 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)) - 32.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(2.0 * Le * Vy)) - 8.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(2.0 * Le * Vy)) - 16.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(2.0 * Le * Vy)) + 12.0))
    
    et[2 , 9] = et[9 , 2] = Mzb / Le
    et[3 , 10] = et[10 , 3] = -(2.0 * (Oz * (Le * Le) * (Vz * Vz) + 1.0) * (Mza + Mzb) * (((((Le * Le * (Vz * Vz) + 4.0 * (sz * sz)) - 2.0 * Le * Vz * np.sin(Le * Vz)) - Le * Le * (Vz * Vz) * (sz * sz)) - np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz))) / (Le * Le * (Vz * Vz) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    et[4 , 11] = et[11 , 4] = -(2.0 * Mxb * (((((((((((((((4.0 * (Vy * Vy) * (sy * sy) * (sz * sz) - 4.0 * (Vz * Vz) * (sy * sy) * (sz * sz)) + Le * Le * (Vy * Vy) * (Vz * Vz) * (sy * sy)) - Le * Le * (Vy * Vy) * (Vz * Vz) * (sz * sz)) + 4.0 * (Le * Le) * np.power(Vy, 4.0) * Oy * (sy * sy) * (sz * sz)) - 4.0 * (Le * Le) * np.power(Vz, 4.0) * Oz * (sy * sy) * (sz * sz)) + 2.0 * Le * Vy * (Vz * Vz) * np.sin(Le * Vy) * (sz * sz)) - 2.0 * Le * (Vy * Vy) * Vz * np.sin(Le * Vz) * (sy * sy)) - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy * (sy * sy) * (sz * sz)) + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz * (sy * sy) * (sz * sz)) - np.power(Le, 3.0) * np.power(Vy, 4.0) * Vz * Oy * np.sin(Le * Vz) * (sy * sy)) + np.power(Le, 3.0) * Vy * np.power(Vz, 4.0) * Oz * np.sin(Le * Vy) * (sz * sz)) + np.power(Le, 3.0) * np.power(Vy, 3.0) * (Vz * Vz) * Oy * np.sin(Le * Vy) * (sz * sz)) - np.power(Le, 3.0) * (Vy * Vy) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz) * (sy * sy)) - 4.0 * np.power(Le, 4.0) * (Vy * Vy) * np.power(Vz, 4.0) * Oy * Oz * (sy * sy) * (sz * sz)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (sy * sy) * (sz * sz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))

    et[1 , 9] = et[9 , 1] = Myb / Le
    et[2 , 10] = et[10 , 2] = ((P * (bz * bz) * (((4.0 * np.cos(Le * Vz) + Le * Le * (Vz * Vz)) + Le * Vz * np.sin(Le * Vz)) - 4.0) / (2.0 * (((((((((4.0 * np.cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * np.cos(Le * Vz)) + 4.0 * Le * Vz * np.sin(Le * Vz)) - 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) - 4.0)) + Emod * Iy * Le * np.power(Vz, 3.0) * (sz * sz) * (np.sin(Le * Vz) - Le * Vz) / (az * az)) + Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * (np.sin(Le * Vz) + Le * Vz) / (2.0 * (((((((((4.0 * np.cos(Le * Vz) - Le * Le * (Vz * Vz)) - 8.0 * (Le * Le) * (Vz * Vz) * Oz) - Le * Le * (Vz * Vz) * np.cos(Le * Vz)) + 4.0 * Le * Vz * np.sin(Le * Vz)) - 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz)) + 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cos(Le * Vz)) + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sin(Le * Vz)) + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cos(Le * Vz)) - 4.0))) + Iy * Le * np.power(Vz, 3.0) * P * (sz * sz) * (np.sin(Le * Vz) - Le * Vz) / (A * (az * az))
    et[3 , 11] = et[11 , 3] = 2.0 * (Oy * (Le * Le) * (Vy * Vy) + 1.0) * (Mya + Myb) * (((((Le * Le * (Vy * Vy) + 4.0 * (sy * sy)) - 2.0 * Le * Vy * np.sin(Le * Vy)) - Le * Le * (Vy * Vy) * (sy * sy)) - np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) / (Le * Le * (Vy * Vy) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)))

    et[1 , 10] = et[10 , 1] = -(Vy * Vy * ((4.0 * Vz * Mxb * (np.sin(Le * Vz) - 2.0 * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0)) - 4.0 * Le * (Vz * Vz) * Mxb * (c_z * c_z - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * np.power(Vz, 3.0) * Mxb * (Oz * np.sin(Le * Vz) - 2.0 * Oz * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0))) - np.power(Vy, 3.0) * ((4.0 * Mxb * (np.sin(Le * Vy) - 2.0 * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vz * Vz) * Mxb * (Oz * np.sin(Le * Vy) - 2.0 * Oz * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0))) - Le * Vz * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    et[2 , 11] = et[11 , 2] = (Vz * Vz * ((4.0 * Vy * Mxb * (np.sin(Le * Vy) - 2.0 * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0)) - 4.0 * Le * (Vy * Vy) * Mxb * (c_y * c_y - c_y * c_y * (c_z * c_z))) + 4.0 * (Le * Le) * np.power(Vy, 3.0) * Mxb * (Oy * np.sin(Le * Vy) - 2.0 * Oy * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0))) - np.power(Vz, 3.0) * ((4.0 * Mxb * (np.sin(Le * Vz) - 2.0 * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Mxb * (Oy * np.sin(Le * Vz) - 2.0 * Oy * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0))) - Le * Vy * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz))) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))

    et[0 , 10] = et[10 , 0] = -Myb / Le
    et[1 , 11] = et[11 , 1] = ((-(P * (by * by) * (((4.0 * np.cos(Le * Vy) + Le * Le * (Vy * Vy)) + Le * Vy * np.sin(Le * Vy)) - 4.0)) / (2.0 * (((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * np.cos(Le * Vy)) + 4.0 * Le * Vy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) - 4.0)) - Emod * Iz * Le * np.power(Vy, 3.0) * (sy * sy) * (np.sin(Le * Vy) - Le * Vy) / (ay * ay)) - Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * (np.sin(Le * Vy) + Le * Vy) / (2.0 * (((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy)) - 8.0 * (Le * Le) * (Vy * Vy) * Oy) - Le * Le * (Vy * Vy) * np.cos(Le * Vy)) + 4.0 * Le * Vy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy)) + 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cos(Le * Vy)) - 4.0))) - Iz * Le * np.power(Vy, 3.0) * P * (sy * sy) * (np.sin(Le * Vy) - Le * Vy) / (A * (ay * ay))
    
    et[0 , 11] = et[11 , 0] = -Mzb / Le
    et[4 , 5] = et[5 , 4] = ((np.power(Le, 3.0) * (np.power(Vy, 3.0) * (Vz * Vz) * (Mxb * (4.0 * Oy * np.sin(Le * Vy) - 4.0 * Oz * np.sin(Le * Vy)) / 2.0 - Mxb * np.sin(Le * Vy / 2.0) * (8.0 * Oy * np.cos(Le * Vy / 2.0) * (c_z * c_z) - 8.0 * Oz * np.cos(Le * Vy / 2.0) * (c_z * c_z)) / 2.0) - Vy * Vy * np.power(Vz, 3.0) * (Mxb * (4.0 * Oy * np.sin(Le * Vz) - 4.0 * Oz * np.sin(Le * Vz)) / 2.0 - Mxb * np.sin(Le * Vz / 2.0) * (8.0 * Oy * (c_y * c_y) * np.cos(Le * Vz / 2.0) - 8.0 * Oz * (c_y * c_y) * np.cos(Le * Vz / 2.0)) / 2.0)) - Le * (((np.power(Vy, 3.0) * (2.0 * Mxb * np.sin(Le * Vy) - 4.0 * Mxb * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0)) + np.power(Vz, 3.0) * (2.0 * Mxb * np.sin(Le * Vz) - 4.0 * Mxb * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0))) - Vy * (Vz * Vz) * (2.0 * Mxb * np.sin(Le * Vy) - 4.0 * Mxb * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0))) - Vy * Vy * Vz * (2.0 * Mxb * np.sin(Le * Vz) - 4.0 * Mxb * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0)))) + Le * Le * ((Vy * np.power(Vz, 3.0) * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz) / 2.0 - Vy * Vy * (Vz * Vz) * Mxb * ((4.0 * (c_y * c_y) + 4.0 * (c_z * c_z)) - 8.0 * (c_y * c_y) * (c_z * c_z)) / 2.0) + np.power(Vy, 3.0) * Vz * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz) / 2.0)) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    et[10 , 11] = et[11 , 10] = -((np.power(Le, 3.0) * (np.power(Vy, 3.0) * (Vz * Vz) * (Mxb * (4.0 * Oy * np.sin(Le * Vy) - 4.0 * Oz * np.sin(Le * Vy)) / 2.0 - Mxb * np.sin(Le * Vy / 2.0) * (8.0 * Oy * np.cos(Le * Vy / 2.0) * (c_z * c_z) - 8.0 * Oz * np.cos(Le * Vy / 2.0) * (c_z * c_z)) / 2.0) - Vy * Vy * np.power(Vz, 3.0) * (Mxb * (4.0 * Oy * np.sin(Le * Vz) - 4.0 * Oz * np.sin(Le * Vz)) / 2.0 - Mxb * np.sin(Le * Vz / 2.0) * (8.0 * Oy * (c_y * c_y) * np.cos(Le * Vz / 2.0) - 8.0 * Oz * (c_y * c_y) * np.cos(Le * Vz / 2.0)) / 2.0)) - Le * (((np.power(Vy, 3.0) * (2.0 * Mxb * np.sin(Le * Vy) - 4.0 * Mxb * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0)) + np.power(Vz, 3.0) * (2.0 * Mxb * np.sin(Le * Vz) - 4.0 * Mxb * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0))) - Vy * (Vz * Vz) * (2.0 * Mxb * np.sin(Le * Vy) - 4.0 * Mxb * np.cos(Le * Vy / 2.0) * (c_z * c_z) * np.sin(Le * Vy / 2.0))) - Vy * Vy * Vz * (2.0 * Mxb * np.sin(Le * Vz) - 4.0 * Mxb * (c_y * c_y) * np.cos(Le * Vz / 2.0) * np.sin(Le * Vz / 2.0)))) + Le * Le * ((Vy * np.power(Vz, 3.0) * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz) / 2.0 - Vy * Vy * (Vz * Vz) * Mxb * ((4.0 * (c_y * c_y) + 4.0 * (c_z * c_z)) - 8.0 * (c_y * c_y) * (c_z * c_z)) / 2.0) + np.power(Vy, 3.0) * Vz * Mxb * np.sin(Le * Vy) * np.sin(Le * Vz) / 2.0)) / ((Vy * Vy - Vz * Vz) * ((4.0 * (sy * sy) - Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy)) * ((4.0 * (sz * sz) - Le * Vz * np.sin(Le * Vz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz)))
    
    # Torsion for symmetric sections
    if Iy == Iz:
        et[1 , 4] = et[4 , 1] = Vy * Mxb * (np.sin(Le * Vy) - Le * Vy) / ((((4.0 * np.cos(Le * Vy) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) - 4.0)
        et[2 , 5] = et[5 , 2] = -(Vy * Mxb * (((((((((2.0 * np.sin(Le * Vy) - np.sin(2.0 * Le * Vy)) - 5.0 * Le * Vy / 2.0) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy) + Le * Le * (Vy * Vy) * np.sin(Le * Vy)) + 2.0 * Le * Vy * np.cos(Le * Vy)) + Le * Vy * np.cos(2.0 * Le * Vy) / 2.0) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(Le * Vy)) - Le * Le * (Vy * Vy) * Oy * np.sin(2.0 * Le * Vy))) / (((((((((((((((((((4.0 * np.cos(2.0 * Le * Vy) - 16.0 * np.cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy) + 12.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * np.cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * np.sin(Le * Vy)) + 4.0 * Le * Vy * np.sin(2.0 * Le * Vy)) + 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz) - 16.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) - 16.0 * (Le * Le) * (Vy * Vy) * Oz * np.cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(2.0 * Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.cos(2.0 * Le * Vy)) - 4.0 * np.power(Le,3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sin(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(2.0 * Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sin(2.0 * Le * Vy)) - 16.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cos(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cos(2.0 * Le * Vy)) + 12.0)
        et[4 , 7] = et[7 , 4] = -(Vy * Mxb * (np.sin(Le * Vy) - Le * Vy)) / ((((4.0 * np.cos(Le * Vy) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) - 4.0)
        et[5 , 8] = et[8 , 5] = Vy * Mxb * (((((((((2.0 * np.sin(Le * Vy) - np.sin(2.0 * Le * Vy)) - 5.0 * Le * Vy / 2.0) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy) + Le * Le * (Vy * Vy) * np.sin(Le * Vy)) + 2.0 * Le * Vy * np.cos(Le * Vy)) + Le * Vy * np.cos(2.0 * Le * Vy) / 2.0) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(Le * Vy)) - Le * Le * (Vy * Vy) * Oy * np.sin(2.0 * Le * Vy)) / (((((((((((((((((((4.0 * np.cos(2.0 * Le * Vy) - 16.0 * np.cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy) + 12.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * np.cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * np.sin(Le * Vy)) + 4.0 * Le * Vy * np.sin(2.0 * Le * Vy)) + 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz) - 16.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) - 16.0 * (Le * Le) * (Vy * Vy) * Oz * np.cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(2.0 * Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.cos(2.0 * Le * Vy)) - 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sin(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(2.0 * Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sin(2.0 * Le * Vy)) - 16.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cos(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cos(2.0 * Le * Vy)) + 12.0)
        et[7 , 10] = et[10 , 7] = Vy * Mxb * (np.sin(Le * Vy) - Le * Vy) / ((((4.0 * np.cos(Le * Vy) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) - 4.0)
        et[8 , 11] = et[11 , 8] = -(Vy * Mxb * (((((((((2.0 * np.sin(Le * Vy) - np.sin(2.0 * Le * Vy)) - 5.0 * Le * Vy / 2.0) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy) + Le * Le * (Vy * Vy) * np.sin(Le * Vy)) + 2.0 * Le * Vy * np.cos(Le * Vy)) + Le * Vy * np.cos(2.0 * Le * Vy) / 2.0) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(Le * Vy)) - Le * Le * (Vy * Vy) * Oy * np.sin(2.0 * Le * Vy))) / (((((((((((((((((((4.0 * np.cos(2.0 * Le * Vy) - 16.0 * np.cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy) + 12.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * np.cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * np.sin(Le * Vy)) + 4.0 * Le * Vy * np.sin(2.0 * Le * Vy)) + 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz) - 16.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) - 16.0 * (Le * Le) * (Vy * Vy) * Oz * np.cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(2.0 * Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.cos(2.0 * Le * Vy)) - 4.0 * np.power(Le,3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sin(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(2.0 * Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sin(2.0 * Le * Vy)) - 16.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cos(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cos(2.0 * Le * Vy)) + 12.0)

        et[5 , 10] = et[10 , 5] = -(Mxb * (((((((((((((((((4.0 * np.sin(3.0 * Le * Vy / 2.0) - 12.0 * np.sin(Le * Vy / 2.0)) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * np.cos(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * np.sin(Le * Vy / 2.0)) + 2.0 * Le * Vy * np.cos(Le * Vy / 2.0)) - 2.0 * Le * Vy * np.cos(3.0 * Le * Vy / 2.0)) + np.power(Le,3.0) * np.power(Vy, 3.0) * Oy * np.cos(Le * Vy / 2.0)) + np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.cos(Le * Vy / 2.0)) - np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(3.0 * Le * Vy / 2.0)) - np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.cos(3.0 * Le * Vy / 2.0)) - 12.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(Le * Vy / 2.0)) - 12.0 * (Le * Le) * (Vy * Vy) * Oz * np.sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(3.0 * Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.sin(3.0 * Le * Vy / 2.0)) + 2.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * np.sin(Le * Vy / 2.0)) + 2.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oz * np.sin(Le * Vy / 2.0)) - 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.sin(Le * Vy / 2.0)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.sin(3.0 * Le * Vy / 2.0))) / (2.0 * (((((((((((((((12.0 * np.sin(Le * Vy / 2.0) - 4.0 * np.sin(3.0 * Le * Vy / 2.0)) + Le * Le * (Vy * Vy) * np.sin(Le * Vy / 2.0)) + Le * Le * (Vy * Vy) * np.sin(3.0 * Le * Vy / 2.0)) - 4.0 * Le * Vy * np.cos(Le * Vy / 2.0)) + 4.0 * Le * Vy * np.cos(3.0 * Le * Vy / 2.0)) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(Le * Vy / 2.0)) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.cos(Le * Vy / 2.0)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(3.0 * Le * Vy / 2.0)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.cos(3.0 * Le * Vy / 2.0)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(Le * Vy / 2.0)) + 12.0 * (Le * Le) * (Vy * Vy) * Oz * np.sin(Le * Vy / 2.0)) - 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(3.0 * Le * Vy / 2.0)) - 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.sin(3.0 * Le * Vy / 2.0)) + 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.sin(Le * Vy / 2.0)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.sin(3.0 * Le * Vy / 2.0)))
        et[4 , 11] = et[11 , 4] = Mxb * (((((((((((((((((4.0 * np.sin(3.0 * Le * Vy / 2.0) - 12.0 * np.sin(Le * Vy / 2.0)) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * np.cos(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * np.sin(Le * Vy / 2.0)) + 2.0 * Le * Vy * np.cos(Le * Vy / 2.0)) - 2.0 * Le * Vy * np.cos(3.0 * Le * Vy / 2.0)) + np.power(Le,3.0) * np.power(Vy, 3.0) * Oy * np.cos(Le * Vy / 2.0)) + np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.cos(Le * Vy / 2.0)) - np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(3.0 * Le * Vy / 2.0)) - np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.cos(3.0 * Le * Vy / 2.0)) - 12.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(Le * Vy / 2.0)) - 12.0 * (Le * Le) * (Vy * Vy) * Oz * np.sin(Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(3.0 * Le * Vy / 2.0)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.sin(3.0 * Le * Vy / 2.0)) + 2.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * np.sin(Le * Vy / 2.0)) + 2.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oz * np.sin(Le * Vy / 2.0)) - 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.sin(Le * Vy / 2.0)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.sin(3.0 * Le * Vy / 2.0)) / (2.0 * (((((((((((((((12.0 * np.sin(Le * Vy / 2.0) - 4.0 * np.sin(3.0 * Le * Vy / 2.0)) + Le * Le * (Vy * Vy) * np.sin(Le * Vy / 2.0)) + Le * Le * (Vy * Vy) * np.sin(3.0 * Le * Vy / 2.0)) - 4.0 * Le * Vy * np.cos(Le * Vy / 2.0)) + 4.0 * Le * Vy * np.cos(3.0 * Le * Vy / 2.0)) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(Le * Vy / 2.0)) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.cos(Le * Vy / 2.0)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(3.0 * Le * Vy / 2.0)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.cos(3.0 * Le * Vy / 2.0)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(Le * Vy / 2.0)) + 12.0 * (Le * Le) * (Vy * Vy) * Oz * np.sin(Le * Vy / 2.0)) - 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(3.0 * Le * Vy / 2.0)) - 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.sin(3.0 * Le * Vy / 2.0)) + 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.sin(Le * Vy / 2.0)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.sin(3.0 * Le * Vy / 2.0)))

        et[1 , 10] = et[10 , 1] = -(Vy * Mxb * (np.sin(Le * Vy) - Le * Vy)) / ((((4.0 * np.cos(Le * Vy) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) + 2.0 * Le * Vy * np.sin(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) - 4.0)
        et[2 , 11] = et[11 , 2] = Vy * Mxb * (((((((((2.0 * np.sin(Le * Vy) - np.sin(2.0 * Le * Vy)) - 5.0 * Le * Vy / 2.0) - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy) + Le * Le * (Vy * Vy) * np.sin(Le * Vy)) + 2.0 * Le * Vy * np.cos(Le * Vy)) + Le * Vy * np.cos(2.0 * Le * Vy) / 2.0) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.cos(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.sin(Le * Vy)) - Le * Le * (Vy * Vy) * Oy * np.sin(2.0 * Le * Vy)) / (((((((((((((((((((4.0 * np.cos(2.0 * Le * Vy) - 16.0 * np.cos(Le * Vy)) + Le * Le * (Vy * Vy)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy) + 12.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * np.cos(2.0 * Le * Vy)) - 8.0 * Le * Vy * np.sin(Le * Vy)) + 4.0 * Le * Vy * np.sin(2.0 * Le * Vy)) + 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz) - 16.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) - 16.0 * (Le * Le) * (Vy * Vy) * Oz * np.cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(2.0 * Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.cos(2.0 * Le * Vy)) - 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) - 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sin(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(2.0 * Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sin(2.0 * Le * Vy)) - 16.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cos(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cos(2.0 * Le * Vy)) + 12.0)

        et[4 , 5] = et[5 , 4] = -(np.power(Le, 3.0) * np.power(Vy, 3.0) * Mxb * (Oy - Oz) * (np.sin(Le * Vy) - Le * Vy)) / (2.0 * ((((((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy)) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) - 4.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * np.cos(Le * Vy)) + 4.0 * Le * Vy * np.sin(Le * Vy)) - 4.0 * np.power(Le,4.0) * np.power(Vy, 4.0) * Oy * Oz) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.cos(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sin(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cos(Le * Vy)) - 4.0))
        et[10 , 11] = et[11 , 10] = np.power(Le, 3.0) * np.power(Vy, 3.0) * Mxb * (Oy - Oz) * (np.sin(Le * Vy) - Le * Vy) / (2.0 * ((((((((((((4.0 * np.cos(Le * Vy) - Le * Le * (Vy * Vy)) - 4.0 * (Le * Le) * (Vy * Vy) * Oy) - 4.0 * (Le * Le) * (Vy * Vy) * Oz) - Le * Le * (Vy * Vy) * np.cos(Le * Vy)) + 4.0 * Le * Vy * np.sin(Le * Vy)) - 4.0 * np.power(Le,4.0) * np.power(Vy, 4.0) * Oy * Oz) + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cos(Le * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.cos(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sin(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sin(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cos(Le * Vy)) - 4.0))
    #
    return et
#
#
def Kt_tension(Le:float, Ax:float,
               Jx:float, Iy:float, Iz:float,
               Emod:float, Gmod:float,
               Oy:float, Oz:float,
               Fb:list[float]):
    """Positive axial force"""
    (Fxa, Fya, Fza, Mxa, Mya, Mza,
     Fxb, Fyb, Fzb, Mxb, Myb, Mzb) = Fb
    #
    P = Fxb
    A = Ax
    L2 = Le*Le
    #L3 = L2*Le
    #pi = np.pi
    #
    mu_y = np.sqrt(P / (Emod*Iz))
    mu_z = np.sqrt(P / (Emod*Iy))

    Vy = mu_y / (np.sqrt(1 + Oy*mu_y*mu_y*L2))
    Vz = mu_z / (np.sqrt(1 + Oz*mu_z*mu_z*L2))

    K = P*(Iy + Iz) / Ax
    
    ay = ((Le * Vy * np.cosh(Le * Vy / 2.0) - 2.0 * np.sinh(Le * Vy / 2.0))
          + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.sinh(Le * Vy / 2.0))
    
    az = ((Le * Vz * np.cosh(Le * Vz / 2.0) - 2.0 * np.sinh(Le * Vz / 2.0))
          + 2.0 * (Le * Le) * (Vz * Vz) * Oz * np.sinh(Le * Vz / 2.0))
    
    by = (((Le * Vy - 2.0 * np.sinh(Le * Vy / 2.0))
           + 2.0 * Le * Vy * (np.sinh(Le * Vy / 4.0) * np.sinh(Le * Vy / 4.0)))
          + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.sinh(Le * Vy / 2.0))
    
    bz = (((Le * Vz - 2.0 * np.sinh(Le * Vz / 2.0))
           + 2.0 * Le * Vz * (np.sinh(Le * Vz / 4.0) * np.sinh(Le * Vz / 4.0)))
          + 2.0 * (Le * Le) * (Vz * Vz) * Oz * np.sinh(Le * Vz / 2.0))
    
    cy = 2.0 * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy) * P
    cz = 2.0 * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz) * P
    
    dy = 4.0 * A * (Le * Le) * (Vy * Vy) * Oy * P
    dz = 4.0 * A * (Le * Le) * (Vz * Vz) * Oz * P
    
    ey = A * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * P
    ez = A * np.power(Le, 4.0) * np.power(Vz, 4.0) * Oz * P
    
    fy = 2.0 * Iz * (Le * Le) * np.power(Vy, 4.0) * Oy * P
    fz = 2.0 * Iy * (Le * Le) * np.power(Vz, 4.0) * Oz * P
    
    sy = np.sinh(Le * Vy / 2.0)
    sz = np.sinh(Le * Vz / 2.0)
    
    c_y = np.cosh(Le * Vy / 2.0)
    c_z = np.cosh(Le * Vz / 2.0)
    #
    #
    et = np.zeros((12, 12))
    #
    #
    et[0 , 0] = et[6 , 6] = P / Le + Ax * Emod / Le
    et[1 , 1] = et[7 , 7] = (Vy * ((((((((((((Iz * (Vy * Vy) * P * np.sinh(Le * Vy)
                                              - 3.0 * A * P * np.sinh(Le * Vy)) + 2.0 * A * Le * Vy * P)
                                            - Iz * Le * np.power(Vy, 3.0) * P)
                                           - A * Emod * Iz * Le * np.power(Vy, 3.0))
                                          + A * Emod * Iz * (Vy * Vy) * np.sinh(Le * Vy))
                                         - 2.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P)
                                        + A * Le * Vy * P * np.cosh(Le * Vy))
                                       + A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P)
                                      + 2.0 * A * (Le * Le) * (Vy * Vy) * Oy * P * np.sinh(Le * Vy))
                                     + A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * np.sinh(Le * Vy))
                                    + A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy)
                                   + A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * Oy * np.sinh(Le * Vy))
                             / (2.0 * A * (ay * ay)))
    
    et[2 , 2] = et[8 , 8] = (Vz * ((((((((((((Iy * (Vz * Vz) * P * np.sinh(Le * Vz)
                                              - 3.0 * A * P * np.sinh(Le * Vz)) + 2.0 * A * Le * Vz * P)
                                            - Iy * Le * np.power(Vz, 3.0) * P)
                                           - A * Emod * Iy * Le * np.power(Vz, 3.0))
                                          + A * Emod * Iy * (Vz * Vz) * np.sinh(Le * Vz))
                                         - 2.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P)
                                        + A * Le * Vz * P * np.cosh(Le * Vz))
                                       + A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P)
                                      + 2.0 * A * (Le * Le) * (Vz * Vz) * Oz * P * np.sinh(Le * Vz))
                                     + A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * np.sinh(Le * Vz))
                                    + A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz)
                                   + A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * Oz * np.sinh(Le * Vz))
                             / (2.0 * A * (az * az)))
    
    et[3 , 3] = et[9 , 9] = K / Le + Gmod*Jx / Le
    
    x = (((((((((((((((((((((((((A * P * np.sinh(2.0 * Le * Vz)
                                 - 2.0 * A * P * np.sinh(Le * Vz))
                                - 2.0 * Iy * (Vz * Vz) * P * np.sinh(Le * Vz))
                               + Iy * (Vz * Vz) * P * np.sinh(2.0 * Le * Vz))
                              - A * np.power(Le, 3.0) * np.power(Vz, 3.0) * P)
                             + Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * P)
                            - Iy * Le * np.power(Vz, 3.0) * P)
                           - A * Emod * Iy * Le * np.power(Vz, 3.0))
                          - 2.0 * A * Emod * Iy * (Vz * Vz) * np.sinh(Le * Vz))
                         + A * Emod * Iy * (Vz * Vz) * np.sinh(2.0 * Le * Vz))
                        + 2.0 * Iy * Le * np.power(Vz, 3.0) * P * np.cosh(Le * Vz))
                       - Iy * Le * np.power(Vz, 3.0) * P * np.cosh(2.0 * Le * Vz))
                      + A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0))
                     - 5.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P)
                    + 2.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * Oz * P)
                   + 3.0 * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * P)
                  + 2.0 * A * (Le * Le) * (Vz * Vz) * P * np.sinh(Le * Vz))
                 + A * (Le * Le) * (Vz * Vz) * P * np.sinh(2.0 * Le * Vz) / 2.0)
                - 2.0 * Iy * (Le * Le) * np.power(Vz, 4.0) * P * np.sinh(Le * Vz))
               + Iy * (Le * Le) * np.power(Vz, 4.0) * P * np.sinh(2.0 * Le * Vz) / 2.0)
              + 2.0 * A * Le * Vz * P * np.cosh(Le * Vz))
             - 2.0 * A * Le * Vz * P * np.cosh(2.0 * Le * Vz))
            + 12.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P)
           - A * np.power(Le, 7.0) * np.power(Vz, 7.0) * (Oz * Oz) * P)
          - 9.0 * A * np.power(Le, 7.0) * np.power(Vz, 7.0) * np.power(Oz, 3.0) * P)
         + 2.0 * A * np.power(Le, 9.0) * np.power(Vz, 9.0) * np.power(Oz, 4.0) * P)
    
    et[4 , 4] = et[10 , 10] = ((((((((((((((((((((((((((((((((((((((((((((((((((x - cz)
                                                                               + 2.0 * dz * np.sinh(Le * Vz))
                                                                              - dz * np.sinh(2.0 * Le * Vz))
                                                                             - 6.0 * ez * np.sinh(Le * Vz))
                                                                            - ez * np.sinh(2.0 * Le * Vz))
                                                                           + 2.0 * fz * np.sinh(Le * Vz))
                                                                          - fz * np.sinh(2.0 * Le * Vz))
                                                                         + 2.0 * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * Oz * P * np.sinh(Le * Vz))
                                                                        - 7.0 * A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz))
                                                                       + 2.0 * A * Emod * Iy * np.power(Le, 7.0) * np.power(Vz, 9.0) * np.power(Oz, 3.0))
                                                                      - 8.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P * np.cosh(Le * Vz))
                                                                     - 4.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P * np.cosh(2.0 * Le * Vz))
                                                                    + 8.0 * A * np.power(Le, 7.0) * np.power(Vz, 7.0) * np.power(Oz, 3.0) * P * np.cosh(Le * Vz))
                                                                   + A * np.power(Le, 7.0) * np.power(Vz, 7.0) * np.power(Oz, 3.0) * P * np.cosh(2.0 * Le * Vz))
                                                                  - 2.0 * A * np.power(Le, 9.0) * np.power(Vz, 9.0) * np.power(Oz, 4.0) * P * np.cosh(Le * Vz))
                                                                 + 2.0 * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz) * P * np.cosh(Le * Vz))
                                                                - 12.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * np.sinh(Le * Vz))
                                                               + 6.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * np.sinh(2.0 * Le * Vz))
                                                              + 6.0 * A * np.power(Le, 6.0) * np.power(Vz, 6.0) * (Oz * Oz) * P * np.sinh(Le * Vz))
                                                             + 8.0 * A * np.power(Le, 6.0) * np.power(Vz, 6.0) * np.power(Oz, 3.0) * P * np.sinh(Le * Vz))
                                                            + A * np.power(Le, 6.0) * np.power(Vz, 6.0) * (Oz * Oz) * P * np.sinh(2.0 * Le * Vz) / 2.0)
                                                           - 4.0 * A * np.power(Le, 6.0) * np.power(Vz, 6.0) * np.power(Oz, 3.0) * P * np.sinh(2.0 * Le * Vz))
                                                          - 2.0 * A * np.power(Le, 8.0) * np.power(Vz, 8.0) * np.power(Oz, 3.0) * P * np.sinh(Le * Vz))
                                                         - 2.0 * A * np.power(Le, 8.0) * np.power(Vz, 8.0) * np.power(Oz, 4.0) * P * np.sinh(Le * Vz))
                                                        + A * np.power(Le, 8.0) * np.power(Vz, 8.0) * np.power(Oz, 4.0) * P * np.sinh(2.0 * Le * Vz))
                                                       - 2.0 * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * P * np.sinh(Le * Vz))
                                                      + Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * P * np.sinh(2.0 * Le * Vz)) + 2.0 * A * Emod * Iy * Le * np.power(Vz, 3.0) * np.cosh(Le * Vz))
                                                    - A * Emod * Iy * Le * np.power(Vz, 3.0) * np.cosh(2.0 * Le * Vz))
                                                   + 6.0 * A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz)
                                                  - A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * Oz)
                                                 - 2.0 * A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * np.sinh(Le * Vz))
                                                + A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * np.sinh(2.0 * Le * Vz) / 2.0)
                                               + 5.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P * np.cosh(2.0 * Le * Vz))
                                              - 4.0 * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * P * np.cosh(Le * Vz))
                                             + Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * P * np.cosh(2.0 * Le * Vz))
                                            - 6.0 * A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * np.cosh(Le * Vz))
                                           + 2.0 * A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * Oz * np.sinh(Le * Vz))
                                          - A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * Oz * np.sinh(2.0 * Le * Vz))
                                         + 4.0 * A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * Oz * np.sinh(Le * Vz))
                                        + A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * Oz * np.sinh(2.0 * Le * Vz) / 2.0)
                                       + 6.0 * A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz) * np.cosh(Le * Vz))
                                      + A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz) * np.cosh(2.0 * Le * Vz))
                                     - 2.0 * A * Emod * Iy * np.power(Le, 7.0) * np.power(Vz, 9.0) * np.power(Oz, 3.0) * np.cosh(Le * Vz))
                                    + 2.0 * A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * np.sinh(Le * Vz))
                                   - A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * np.sinh(2.0 * Le * Vz))
                                  - 2.0 * A * Emod * Iy * np.power(Le, 6.0) * np.power(Vz, 8.0) * (Oz * Oz) * np.sinh(Le * Vz))
                                 - 2.0 * A * Emod * Iy * np.power(Le, 6.0) * np.power(Vz, 8.0) * np.power(Oz, 3.0) * np.sinh(Le * Vz))
                                + A * Emod * Iy * np.power(Le, 6.0) * np.power(Vz, 8.0) * np.power(Oz, 3.0) * np.sinh(2.0 * Le * Vz))
                               / (2.0 * A * Vz * (np.cosh(Le * Vz) - 1.0)
                                  * (((((((((4.0 * np.cosh(Le * Vz) + Le * Le * (Vz * Vz))
                                            + 8.0 * (Le * Le) * (Vz * Vz) * Oz) + Le * Le * (Vz * Vz) * np.cosh(Le * Vz))
                                          - 4.0 * Le * Vz * np.sinh(Le * Vz))
                                         - 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz))
                                        - 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cosh(Le * Vz))
                                       + 4.0 * np.power(Le, 3.0)   * np.power(Vz, 3.0) * Oz * np.sinh(Le * Vz))
                                      + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cosh(Le * Vz)) - 4.0)))
    
    x = (((((((((((((((((((((((((A * P * np.sinh(2.0 * Le * Vy)
                                 - 2.0 * A * P * np.sinh(Le * Vy))
                                - 2.0 * Iz * (Vy * Vy) * P * np.sinh(Le * Vy))
                               + Iz * (Vy * Vy) * P * np.sinh(2.0 * Le * Vy))
                              - A * np.power(Le, 3.0) * np.power(Vy, 3.0) * P)
                             + Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * P)
                            - Iz * Le * np.power(Vy, 3.0) * P) - A * Emod * Iz * Le * np.power(Vy, 3.0))
                          - 2.0 * A * Emod * Iz * (Vy * Vy) * np.sinh(Le * Vy))
                         + A * Emod * Iz * (Vy * Vy) * np.sinh(2.0 * Le * Vy))
                        + 2.0 * Iz * Le * np.power(Vy, 3.0) * P * np.cosh(Le * Vy))
                       - Iz * Le * np.power(Vy, 3.0) * P * np.cosh(2.0 * Le * Vy))
                      + A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0))
                     - 5.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P)
                    + 2.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * Oy * P)
                   + 3.0 * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * P)
                  + 2.0 * A * (Le * Le) * (Vy * Vy) * P * np.sinh(Le * Vy))
                 + A * (Le * Le) * (Vy * Vy) * P * np.sinh(2.0 * Le * Vy) / 2.0)
                - 2.0 * Iz * (Le * Le) * np.power(Vy, 4.0) * P * np.sinh(Le * Vy))
               + Iz * (Le * Le) * np.power(Vy, 4.0) * P * np.sinh(2.0 * Le * Vy) / 2.0)
              + 2.0 * A * Le * Vy * P * np.cosh(Le * Vy)) - 2.0 * A * Le * Vy * P * np.cosh(2.0 * Le * Vy))
            + 12.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P)
           - A * np.power(Le, 7.0) * np.power(Vy, 7.0) * (Oy * Oy) * P)
          - 9.0 * A * np.power(Le, 7.0) * np.power(Vy, 7.0) * np.power(Oy, 3.0) * P)
         + 2.0 * A * np.power(Le, 9.0) * np.power(Vy, 9.0) * np.power(Oy, 4.0) * P)
    
    et[5 , 5] = et[11 , 11] = ((((((((((((((((((((((((((((((((((((((((((((((((((x - cy)
                                                                               + 2.0 * dy * np.sinh(Le * Vy))
                                                                              - dy * np.sinh(2.0 * Le * Vy))
                                                                             - 6.0 * ey * np.sinh(Le * Vy))
                                                                            - ey * np.sinh(2.0 * Le * Vy))
                                                                           + 2.0 * fy * np.sinh(Le * Vy))
                                                                          - fy * np.sinh(2.0 * Le * Vy))
                                                                         + 2.0 * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * Oy * P * np.sinh(Le * Vy))
                                                                        - 7.0 * A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy))
                                                                       + 2.0 * A * Emod * Iz * np.power(Le, 7.0) * np.power(Vy, 9.0) * np.power(Oy, 3.0))
                                                                      - 8.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P * np.cosh(Le * Vy))
                                                                     - 4.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P * np.cosh(2.0 * Le * Vy))
                                                                    + 8.0 * A * np.power(Le, 7.0) * np.power(Vy, 7.0) * np.power(Oy, 3.0) * P * np.cosh(Le * Vy))
                                                                   + A * np.power(Le, 7.0) * np.power(Vy, 7.0) * np.power(Oy, 3.0) * P * np.cosh(2.0 * Le * Vy))
                                                                  - 2.0 * A * np.power(Le, 9.0) * np.power(Vy, 9.0) * np.power(Oy, 4.0) * P * np.cosh(Le * Vy))
                                                                 + 2.0 * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy) * P * np.cosh(Le * Vy))
                                                                - 12.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * np.sinh(Le * Vy))
                                                               + 6.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * np.sinh(2.0 * Le * Vy))
                                                              + 6.0 * A * np.power(Le, 6.0) * np.power(Vy, 6.0) * (Oy * Oy) * P * np.sinh(Le * Vy))
                                                             + 8.0 * A * np.power(Le, 6.0) * np.power(Vy, 6.0) * np.power(Oy, 3.0) * P * np.sinh(Le * Vy))
                                                            + A * np.power(Le, 6.0) * np.power(Vy, 6.0) * (Oy * Oy) * P * np.sinh(2.0 * Le * Vy) / 2.0)
                                                           - 4.0 * A * np.power(Le, 6.0) * np.power(Vy, 6.0) * np.power(Oy, 3.0) * P * np.sinh(2.0 * Le * Vy))
                                                          - 2.0 * A * np.power(Le, 8.0) * np.power(Vy, 8.0) * np.power(Oy, 3.0) * P * np.sinh(Le * Vy))
                                                         - 2.0 * A * np.power(Le, 8.0) * np.power(Vy, 8.0) * np.power(Oy, 4.0) * P * np.sinh(Le * Vy))
                                                        + A * np.power(Le, 8.0) * np.power(Vy, 8.0) * np.power(Oy, 4.0) * P * np.sinh(2.0 * Le * Vy))
                                                       - 2.0 * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * P * np.sinh(Le * Vy))
                                                      + Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * P * np.sinh(2.0 * Le * Vy))
                                                     + 2.0 * A * Emod * Iz * Le * np.power(Vy, 3.0) * np.cosh(Le * Vy))
                                                    - A * Emod * Iz * Le * np.power(Vy, 3.0) * np.cosh(2.0 * Le * Vy))
                                                   + 6.0 * A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy)
                                                  - A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * Oy)
                                                 - 2.0 * A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * np.sinh(Le * Vy))
                                                + A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * np.sinh(2.0 * Le * Vy) / 2.0)
                                               + 5.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P * np.cosh(2.0 * Le * Vy))
                                              - 4.0 * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * P * np.cosh(Le * Vy))
                                             + Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * P * np.cosh(2.0 * Le * Vy))
                                            - 6.0 * A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * np.cosh(Le * Vy))
                                           + 2.0 * A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * Oy * np.sinh(Le * Vy))
                                          - A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * Oy * np.sinh(2.0 * Le * Vy))
                                         + 4.0 * A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * Oy * np.sinh(Le * Vy))
                                        + A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * Oy * np.sinh(2.0 * Le * Vy) / 2.0)
                                       + 6.0 * A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy) * np.cosh(Le * Vy))
                                      + A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy) * np.cosh(2.0 * Le * Vy))
                                     - 2.0 * A * Emod * Iz * np.power(Le, 7.0) * np.power(Vy, 9.0) * np.power(Oy, 3.0) * np.cosh(Le * Vy))
                                    + 2.0 * A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * np.sinh(Le * Vy))
                                   - A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * np.sinh(2.0 * Le * Vy))
                                  - 2.0 * A * Emod * Iz * np.power(Le, 6.0) * np.power(Vy, 8.0) * (Oy * Oy) * np.sinh(Le * Vy))
                                 - 2.0 * A * Emod * Iz * np.power(Le, 6.0) * np.power(Vy, 8.0) * np.power(Oy, 3.0) * np.sinh(Le * Vy))
                                + A * Emod * Iz * np.power(Le, 6.0) * np.power(Vy, 8.0) * np.power(Oy, 3.0) * np.sinh(2.0 * Le * Vy))
                               / (2.0 * A * Vy * (np.cosh(Le * Vy) - 1.0) * (((((((((4.0 * np.cosh(Le * Vy) + Le * Le * (Vy * Vy))
                                                                                    + 8.0 * (Le * Le) * (Vy * Vy) * Oy)
                                                                                   + Le * Le * (Vy * Vy) * np.cosh(Le * Vy))
                                                                                  - 4.0 * Le * Vy * np.sinh(Le * Vy))
                                                                                 - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy))
                                                                                - 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy))
                                                                               + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sinh(Le * Vy))
                                                                              + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cosh(Le * Vy)) - 4.0)))

    et[2 , 3] = et[3 , 2] = Mza / Le
    
    et[3 , 4] = et[4 , 3] = (-Mza / 2.0 - (Oz * (Le * Le) * (Vz * Vz) - 1.0) * (Mza + Mzb)
                             * ((((Le * Vz - 2.0 * np.sinh(Le * Vz)) - 2.0 * np.cosh(Le * Vz))
                                 + Le * Vz * (np.cosh(Le * Vz) + np.sinh(Le * Vz))) + 2.0)
                             / (Le * Le * (Vz * Vz) * ((np.cosh(Le * Vz) + np.sinh(Le * Vz)) - 1.0)))
    
    et[5 , 6] = et[6 , 5] = Mza / Le
    et[8 , 9] = et[9 , 8] = -Mzb / Le
    
    et[9 , 10] = et[10 , 9] = (-Mzb / 2.0 - (Oz * (Le * Le) * (Vz * Vz) - 1.0) * (Mza + Mzb)
                               * ((((Le * Vz - 2.0 * np.sinh(Le * Vz)) - 2.0 * np.cosh(Le * Vz))
                                   + Le * Vz * (np.cosh(Le * Vz) + np.sinh(Le * Vz))) + 2.0)
                               / (Le * Le * (Vz * Vz) * ((np.cosh(Le * Vz) + np.sinh(Le * Vz)) - 1.0)))

    et[1 , 3] = et[3 , 1] = Mya / Le
    
    et[2 , 4] = et[4 , 2] = (-((((((((((((((A * (Le * Le) * (Vz * Vz) * P - 8.0 * A * P * (sz * sz))
                                           - Iy * (Le * Le) * np.power(Vz, 4.0) * P)
                                          + Iy * Le * np.power(Vz, 3.0) * P * np.sinh(Le * Vz))
                                         - A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0))
                                        - 2.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * Oz * P)
                                       + A * Le * Vz * P * np.sinh(Le * Vz))
                                      + A * np.power(Le, 6.0) * np.power(Vz, 6.0) * (Oz * Oz) * P)
                                     - 2.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P * np.sinh(Le * Vz))
                                    + 16.0 * A * (Le * Le) * (Vz * Vz) * Oz * P * (sz * sz))
                                   + A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P * np.sinh(Le * Vz))
                                  + A * Emod * Iy * Le * np.power(Vz, 3.0) * np.sinh(Le * Vz))
                                 - 8.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * (sz * sz))
                                + A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * Oz)
                               + A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * np.sinh(Le * Vz))
                             / (4.0 * A * (bz * bz)))
    
    et[3 , 5] = et[5 , 3] =( Mya / 2.0 + (Oy * (Le * Le) * (Vy * Vy) - 1.0) * (Mya + Myb)
                             * ((((Le * Vy - 2.0 * np.sinh(Le * Vy)) - 2.0 * np.cosh(Le * Vy))
                                 + Le * Vy * (np.cosh(Le * Vy) + np.sinh(Le * Vy))) + 2.0)
                             / (Le * Le * (Vy * Vy) * ((np.cosh(Le * Vy) + np.sinh(Le * Vy)) - 1.0)))
    
    et[4 , 6] = et[6 , 4] = Mya / Le
    
    et[5 , 7] = et[7 , 5] = (-((((((((((((((A * (Le * Le) * (Vy * Vy) * P - 8.0 * A * P * (sy * sy))
                                           - Iz * (Le * Le) * np.power(Vy, 4.0) * P)
                                          + Iz * Le * np.power(Vy, 3.0) * P * np.sinh(Le * Vy))
                                         - A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0))
                                        - 2.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * P)
                                       + A * Le * Vy * P * np.sinh(Le * Vy))
                                      + A * np.power(Le, 6.0) * np.power(Vy, 6.0) * (Oy * Oy) * P)
                                     - 2.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P * np.sinh(Le * Vy))
                                    + 16.0 * A * (Le * Le) * (Vy * Vy) * Oy * P * (sy * sy))
                                   + A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P * np.sinh(Le * Vy))
                                  + A * Emod * Iz * Le * np.power(Vy, 3.0) * np.sinh(Le * Vy))
                                 - 8.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * (sy * sy))
                                + A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * Oy)
                               + A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * np.sinh(Le * Vy))
                             / (4.0 * A * (by * by)))
    
    et[7 , 9] = et[9 , 7] = -Myb / Le
    
    et[8 , 10] = et[10 , 8] = (((((((((((((((A * (Le * Le) * (Vz * Vz) * P - 8.0 * A * P * (sz * sz))
                                            - Iy * (Le * Le) * np.power(Vz, 4.0) * P)
                                           + Iy * Le * np.power(Vz, 3.0) * P * np.sinh(Le * Vz))
                                          - A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0))
                                         - 2.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * Oz * P)
                                        + A * Le * Vz * P * np.sinh(Le * Vz))
                                       + A * np.power(Le, 6.0) * np.power(Vz, 6.0) * (Oz * Oz) * P)
                                      - 2.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P * np.sinh(Le * Vz))
                                     + 16.0 * A * (Le * Le) * (Vz * Vz) * Oz * P * (sz * sz))
                                    + A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P * np.sinh(Le * Vz))
                                   + A * Emod * Iy * Le * np.power(Vz, 3.0) * np.sinh(Le * Vz))
                                  - 8.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * (sz * sz))
                                 + A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * Oz)
                                + A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * np.sinh(Le * Vz))
                               / (4.0 * A * (bz * bz)))
    
    et[9 , 11] = et[11 , 9] = Myb / 2.0 + (Oy * (Le * Le) * (Vy * Vy) - 1.0) * (Mya + Myb) * ((((Le * Vy - 2.0 * np.sinh(Le * Vy)) - 2.0 * np.cosh(Le * Vy)) + Le * Vy * (np.cosh(Le * Vy) + np.sinh(Le * Vy))) + 2.0) / (Le * Le * (Vy * Vy) * ((np.cosh(Le * Vy) + np.sinh(Le * Vy)) - 1.0))

    et[1 , 4] = et[4 , 1] = Vy * Vy * Mxb * (Vy * np.cosh(Le * Vy / 2.0) * np.sinh(Le * Vz / 2.0) - Vz * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vy / 2.0)) / (np.sinh(Le * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((Le * Vy * np.cosh(Le * Vy / 2.0) - 2.0 * np.sinh(Le * Vy / 2.0)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.sinh(Le * Vy / 2.0)))
    et[2 , 5] = et[5 , 2] = Vz * Vz * Mxb * (Vy * np.cosh(Le * Vy / 2.0) * np.sinh(Le * Vz / 2.0) - Vz * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vy / 2.0)) / (np.sinh(Le * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((Le * Vz * np.cosh(Le * Vz / 2.0) - 2.0 * np.sinh(Le * Vz / 2.0)) + 2.0 * (Le * Le) * (Vz * Vz) * Oz * np.sinh(Le * Vz / 2.0)))
    et[4 , 7] = et[7 , 4] = -(Vy * Vy * Mxb * (Vy * np.cosh(Le * Vy / 2.0) * np.sinh(Le * Vz / 2.0) - Vz * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vy / 2.0))) / (np.sinh(Le * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((Le * Vy * np.cosh(Le * Vy / 2.0) - 2.0 * np.sinh(Le * Vy / 2.0)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.sinh(Le * Vy / 2.0)))
    et[5 , 8] = et[8 , 5] = -(Vz * Vz * Mxb * (Vy * np.cosh(Le * Vy / 2.0) * np.sinh(Le * Vz / 2.0) - Vz * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vy / 2.0))) / (np.sinh(Le * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((Le * Vz * np.cosh(Le * Vz / 2.0) - 2.0 * np.sinh(Le * Vz / 2.0)) + 2.0 * (Le * Le) * (Vz * Vz) * Oz * np.sinh(Le * Vz / 2.0)))
    et[7 , 10] = et[10 , 7] = Vy * Vy * Mxb * (Vy * np.cosh(Le * Vy / 2.0) * np.sinh(Le * Vz / 2.0) - Vz * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vy / 2.0)) / (np.sinh(Le * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((Le * Vy * np.cosh(Le * Vy / 2.0) - 2.0 * np.sinh(Le * Vy / 2.0)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.sinh(Le * Vy / 2.0)))
    et[8 , 11] = et[11 , 8] = Vz * Vz * Mxb * (Vy * np.cosh(Le * Vy / 2.0) * np.sinh(Le * Vz / 2.0) - Vz * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vy / 2.0)) / (np.sinh(Le * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((Le * Vz * np.cosh(Le * Vz / 2.0) - 2.0 * np.sinh(Le * Vz / 2.0)) + 2.0 * (Le * Le) * (Vz * Vz) * Oz * np.sinh(Le * Vz / 2.0)))

    et[0 , 4] = et[4 , 0] = -Mya / Le
    
    et[1 , 5] = et[5 , 1] = (((((((((((((((A * (Le * Le) * (Vy * Vy) * P - 8.0 * A * P * (sy * sy))
                                          - Iz * (Le * Le) * np.power(Vy, 4.0) * P)
                                         + Iz * Le * np.power(Vy, 3.0) * P * np.sinh(Le * Vy))
                                        - A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0))
                                       - 2.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * P)
                                      + A * Le * Vy * P * np.sinh(Le * Vy))
                                     + A * np.power(Le, 6.0) * np.power(Vy, 6.0) * (Oy * Oy) * P)
                                    - 2.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P * np.sinh(Le * Vy))
                                   + 16.0 * A * (Le * Le) * (Vy * Vy) * Oy * P * (sy * sy))
                                  + A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P * np.sinh(Le * Vy))
                                 + A * Emod * Iz * Le * np.power(Vy, 3.0) * np.sinh(Le * Vy))
                                - 8.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * (sy * sy))
                               + A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * Oy)
                              + A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * np.sinh(Le * Vy))
                             / (4.0 * A * (by * by)))
    
    et[3 , 7] = et[7 , 3] = -Mya / Le
    
    et[4 , 8] = et[8 , 4] = ((((((((((((((A * (Le * Le) * (Vz * Vz) * P - 8.0 * A * P * (sz * sz)) - Iy * (Le * Le) * np.power(Vz, 4.0) * P) + Iy * Le * np.power(Vz, 3.0) * P * np.sinh(Le * Vz)) - A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0)) - 2.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * Oz * P) + A * Le * Vz * P * np.sinh(Le * Vz)) + A * np.power(Le, 6.0) * np.power(Vz, 6.0) * (Oz * Oz) * P) - 2.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P * np.sinh(Le * Vz)) + 16.0 * A * (Le * Le) * (Vz * Vz) * Oz * P * (sz * sz)) + A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P * np.sinh(Le * Vz)) + A * Emod * Iy * Le * np.power(Vz, 3.0) * np.sinh(Le * Vz)) - 8.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * (sz * sz)) + A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * Oz) + A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * np.sinh(Le * Vz)) / (4.0 * A * (bz * bz))
    et[5 , 9] = et[9 , 5] = -((Oy * (Le * Le) * (Vy * Vy) - 1.0) * (Le * Vy * np.cosh(Le * Vy / 2.0) / np.sinh(Le * Vy / 2.0) - 2.0) * (Mya + Myb)) / (Le * Le * (Vy * Vy))
    et[6 , 10] = et[10 , 6] = Myb / Le
    et[7 , 11] = et[11 , 7] = -((((((((((((((A * (Le * Le) * (Vy * Vy) * P - 8.0 * A * P * (sy * sy)) - Iz * (Le * Le) * np.power(Vy, 4.0) * P) + Iz * Le * np.power(Vy, 3.0) * P * np.sinh(Le * Vy)) - A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0)) - 2.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * P) + A * Le * Vy * P * np.sinh(Le * Vy)) + A * np.power(Le, 6.0) * np.power(Vy, 6.0) * (Oy * Oy) * P) - 2.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P * np.sinh(Le * Vy)) + 16.0 * A * (Le * Le) * (Vy * Vy) * Oy * P * (sy * sy)) + A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P * np.sinh(Le * Vy)) + A * Emod * Iz * Le * np.power(Vy, 3.0) * np.sinh(Le * Vy)) - 8.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * (sy * sy)) + A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * Oy) + A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * np.sinh(Le * Vy)) / (4.0 * A * (by * by))

    et[0 , 5] = et[5 , 0] = -Mza / Le
    et[3 , 8] = et[8 , 3] = -Mza / Le
    et[4 , 9] = et[9 , 4] = (Oz * (Le * Le) * (Vz * Vz) - 1.0) * (Le * Vz * np.cosh(Le * Vz / 2.0) / np.sinh(Le * Vz / 2.0) - 2.0) * (Mza + Mzb) / (Le * Le * (Vz * Vz))
    
    et[5 , 10] = et[10 , 5] = (2.0 * Mxb
                               * (((((((((((((((((((((((((((((((((((((((((((((4.0 * Vy * Vy - 4.0 * Vz * Vz)
                                                                             - 4.0 * (Vy * Vy) * c_y * c_y)
                                                                            - 4.0 * Vy * Vy * (c_z * c_z))
                                                                           + 4.0 * Vz * Vz * (c_y * c_y))
                                                                          + 4.0 * (Vz * Vz) * (c_z * c_z))
                                                                         - 4.0 * Le * Le * np.power(Vy, 4.0) * Oy)
                                                                        + 4.0 * (Le * Le) * np.power(Vz, 4.0) * Oz)
                                                                       + 4.0 * (Vy * Vy) * (c_y * c_y) * (c_z * c_z))
                                                                      - 4.0 * (Vz * Vz) * (c_y * c_y) * (c_z * c_z))
                                                                     - 2.0 * Le * Vy * (Vz * Vz) * np.sinh(Le * Vy))
                                                                    + 2.0 * Le * (Vy * Vy) * Vz * np.sinh(Le * Vz))
                                                                   + Le * Le * (Vy * Vy) * (Vz * Vz) * (c_y * c_y))
                                                                  - Le * Le * (Vy * Vy) * (Vz * Vz) * (c_z * c_z))
                                                                 + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy)
                                                                - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz)
                                                               + 4.0 * (Le * Le) * np.power(Vy, 4.0) * Oy * (c_y * c_y))
                                                              + 4.0 * (Le * Le) * np.power(Vy, 4.0) * Oy * (c_z * c_z))
                                                             - 4.0 * (Le * Le) * np.power(Vz, 4.0) * Oz * (c_y * c_y))
                                                            - 4.0 * (Le * Le) * np.power(Vz, 4.0) * Oz * (c_z * c_z))
                                                           - 4.0 * np.power(Le, 4.0) * (Vy * Vy) * np.power(Vz, 4.0) * Oy * Oz)
                                                          + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Vz * Vz) * Oy * Oz)
                                                         + np.power(Le, 3.0) * np.power(Vy, 3.0) * (Vz * Vz) * Oy * np.sinh(Le * Vy))
                                                        - np.power(Le, 3.0) * (Vy * Vy) * np.power(Vz, 3.0) * Oz * np.sinh(Le * Vz))
                                                       - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy * (c_y * c_y))
                                                      - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy * (c_z * c_z))
                                                     + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz * (c_y * c_y))
                                                    + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz * (c_z * c_z))
                                                   - 4.0 * (Le * Le) * np.power(Vy, 4.0) * Oy * (c_y * c_y) * (c_z * c_z))
                                                  + 4.0 * (Le * Le) * np.power(Vz, 4.0) * Oz * (c_y * c_y) * (c_z * c_z))
                                                 + np.power(Le, 3.0) * Vy * np.power(Vz, 4.0) * Oz * np.sinh(Le * Vy))
                                                - np.power(Le, 3.0) * np.power(Vy, 4.0) * Vz * Oy * np.sinh(Le * Vz))
                                               + 4.0 * Le * Vy * (Vz * Vz) * np.cosh(Le * Vy / 2.0) * (c_z * c_z) * np.sinh(Le * Vy / 2.0))
                                              - 4.0 * Le * (Vy * Vy) * Vz * (c_y * c_y) * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vz / 2.0))
                                             + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy * (c_y * c_y) * (c_z * c_z))
                                            - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz * (c_y * c_y) * (c_z * c_z))
                                           + 4.0 * np.power(Le, 4.0) * (Vy * Vy) * np.power(Vz, 4.0) * Oy * Oz * (c_y * c_y))
                                          - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_y * c_y))
                                         + 4.0 * np.power(Le, 4.0) * (Vy * Vy) * np.power(Vz, 4.0) * Oy * Oz * (c_z * c_z))
                                        - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_z * c_z))
                                       - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * (Vz * Vz) * Oy * np.cosh(Le * Vy / 2.0) * (c_z * c_z) * np.sinh(Le * Vy / 2.0))
                                      + 2.0 * np.power(Le, 3.0) * (Vy * Vy) * np.power(Vz, 3.0) * Oz * (c_y * c_y) * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vz / 2.0))
                                     - 4.0 * np.power(Le, 4.0) * (Vy * Vy) * np.power(Vz, 4.0) * Oy * Oz * (c_y * c_y) * (c_z * c_z))
                                    + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_y * c_y) * (c_z * c_z))
                                   - 2.0 * np.power(Le, 3.0) * Vy * np.power(Vz, 4.0) * Oz * np.cosh(Le * Vy / 2.0) * (c_z * c_z) * np.sinh(Le * Vy / 2.0))
                                  + 2.0 * np.power(Le, 3.0) * np.power(Vy, 4.0) * Vz * Oy * (c_y * c_y) * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vz / 2.0))
                               / ((Vy * Vy - Vz * Vz) * ((((Le * Vy * np.sinh(Le * Vy)
                                                            - 2.0 * (Le * Le) * (Vy * Vy) * Oy) - 2.0 * np.cosh(Le * Vy))
                                                          + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy)) + 2.0)
                                  * ((((Le * Vz * np.sinh(Le * Vz) - 2.0 * (Le * Le) * (Vz * Vz) * Oz) - 2.0 * np.cosh(Le * Vz))
                                      + 2.0 * (Le * Le) * (Vz * Vz) * Oz * np.cosh(Le * Vz)) + 2.0)))
    
    et[6 , 11] = et[11 , 6] = (-(2.0 * Mxb
                                 * (((((((((((((((((((((((((((((((((((((((((((((4.0 * (Vy * Vy) - 4.0 * (Vz * Vz))
                                                                               - 4.0 * (Vy * Vy) * c_y * c_y)
                                                                              - 4.0 * (Vy * Vy) * c_z * c_z)
                                                                             + 4.0 * Vz * Vz * (c_y * c_y))
                                                                            + 4.0 * Vz * Vz * (c_z * c_z))
                                                                           - 4.0 * Le * Le * np.power(Vy, 4.0) * Oy)
                                                                          + 4.0 * Le * Le * np.power(Vz, 4.0) * Oz)
                                                                         + 4.0 * (Vy * Vy) * (c_y * c_y) * (c_z * c_z))
                                                                        - 4.0 * (Vz * Vz) * (c_y * c_y) * (c_z * c_z))
                                                                       - 2.0 * Le * Vy * (Vz * Vz) * np.sinh(Le * Vy))
                                                                      + 2.0 * Le * (Vy * Vy) * Vz * np.sinh(Le * Vz))
                                                                     + Le * Le * (Vy * Vy) * (Vz * Vz) * (c_y * c_y))
                                                                    - Le * Le * (Vy * Vy) * (Vz * Vz) * (c_z * c_z))
                                                                   + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy)
                                                                  - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz)
                                                                 + 4.0 * (Le * Le) * np.power(Vy, 4.0) * Oy * (c_y * c_y))
                                                                + 4.0 * (Le * Le) * np.power(Vy, 4.0) * Oy * (c_z * c_z))
                                                               - 4.0 * (Le * Le) * np.power(Vz, 4.0) * Oz * (c_y * c_y))
                                                              - 4.0 * (Le * Le) * np.power(Vz, 4.0) * Oz * (c_z * c_z))
                                                             - 4.0 * np.power(Le, 4.0) * (Vy * Vy) * np.power(Vz, 4.0) * Oy * Oz)
                                                            + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Vz * Vz) * Oy * Oz)
                                                           + np.power(Le, 3.0) * np.power(Vy, 3.0) * (Vz * Vz) * Oy * np.sinh(Le * Vy))
                                                          - np.power(Le, 3.0) * (Vy * Vy) * np.power(Vz, 3.0) * Oz * np.sinh(Le * Vz))
                                                         - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy * (c_y * c_y))
                                                        - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy * (c_z * c_z))
                                                       + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz * (c_y * c_y))
                                                      + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz * (c_z * c_z))
                                                     - 4.0 * (Le * Le) * np.power(Vy, 4.0) * Oy * (c_y * c_y) * (c_z * c_z))
                                                    + 4.0 * (Le * Le) * np.power(Vz, 4.0) * Oz * (c_y * c_y) * (c_z * c_z))
                                                   + np.power(Le, 3.0) * Vy * np.power(Vz, 4.0) * Oz * np.sinh(Le * Vy))
                                                  - np.power(Le, 3.0) * np.power(Vy, 4.0) * Vz * Oy * np.sinh(Le * Vz))
                                                 + 4.0 * Le * Vy * (Vz * Vz) * np.cosh(Le * Vy / 2.0) * (c_z * c_z) * np.sinh(Le * Vy / 2.0))
                                                - 4.0 * Le * (Vy * Vy) * Vz * (c_y * c_y) * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vz / 2.0))
                                               + 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oy * (c_y * c_y) * (c_z * c_z))
                                              - 4.0 * (Le * Le) * (Vy * Vy) * (Vz * Vz) * Oz * (c_y * c_y) * (c_z * c_z))
                                             + 4.0 * np.power(Le, 4.0) * (Vy * Vy) * np.power(Vz, 4.0) * Oy * Oz * (c_y * c_y))
                                            - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_y * c_y))
                                           + 4.0 * np.power(Le, 4.0) * (Vy * Vy) * np.power(Vz, 4.0) * Oy * Oz * (c_z * c_z))
                                          - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_z * c_z))
                                         - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * (Vz * Vz) * Oy * np.cosh(Le * Vy / 2.0) * (c_z * c_z) * np.sinh(Le * Vy / 2.0))
                                        + 2.0 * np.power(Le, 3.0) * (Vy * Vy) * np.power(Vz, 3.0) * Oz * (c_y * c_y) * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vz / 2.0))
                                       - 4.0 * np.power(Le, 4.0) * (Vy * Vy) * np.power(Vz, 4.0) * Oy * Oz * (c_y * c_y) * (c_z * c_z))
                                      + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Vz * Vz) * Oy * Oz * (c_y * c_y) * (c_z * c_z))
                                     - 2.0 * np.power(Le, 3.0) * Vy * np.power(Vz, 4.0) * Oz * np.cosh(Le * Vy / 2.0) * (c_z * c_z) * np.sinh(Le * Vy / 2.0))
                                    + 2.0 * np.power(Le, 3.0) * np.power(Vy, 4.0) * Vz * Oy * (c_y * c_y) * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vz / 2.0)))
                               / ((Vy * Vy - Vz * Vz) * ((((Le * Vy * np.sinh(Le * Vy)
                                                            - 2.0 * (Le * Le) * (Vy * Vy) * Oy) - 2.0 * np.cosh(Le * Vy))
                                                          + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy)) + 2.0)
                                  * ((((Le * Vz * np.sinh(Le * Vz) - 2.0 * (Le * Le) * (Vz * Vz) * Oz)
                                       - 2.0 * np.cosh(Le * Vz)) + 2.0 * (Le * Le) * (Vz * Vz) * Oz * np.cosh(Le * Vz)) + 2.0)))
    
    et[0 , 6] = et[6 , 0] = -P / Le - Ax * Emod / Le
    
    et[1 , 7] = et[7 , 1] = (-(Vy * ((((((((((((Iz * (Vy * Vy) * P * np.sinh(Le * Vy)
                                                - 3.0 * A * P * np.sinh(Le * Vy))
                                               + 2.0 * A * Le * Vy * P) - Iz * Le * np.power(Vy, 3.0) * P)
                                             - A * Emod * Iz * Le * np.power(Vy, 3.0))
                                            + A * Emod * Iz * (Vy * Vy) * np.sinh(Le * Vy))
                                           - 2.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P)
                                          + A * Le * Vy * P * np.cosh(Le * Vy))
                                         + A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P)
                                        + 2.0 * A * (Le * Le) * (Vy * Vy) * Oy * P * np.sinh(Le * Vy))
                                       + A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * np.sinh(Le * Vy))
                                      + A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy)
                                     + A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * Oy * np.sinh(Le * Vy)))
                             / (2.0 * A * (ay * ay)))
    
    et[2 , 8] = et[8 , 2] = (-(Vz * ((((((((((((Iy * (Vz * Vz) * P * np.sinh(Le * Vz)
                                                - 3.0 * A * P * np.sinh(Le * Vz))
                                               + 2.0 * A * Le * Vz * P) - Iy * Le * np.power(Vz, 3.0) * P)
                                             - A * Emod * Iy * Le * np.power(Vz, 3.0))
                                            + A * Emod * Iy * (Vz * Vz) * np.sinh(Le * Vz))
                                           - 2.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P)
                                          + A * Le * Vz * P * np.cosh(Le * Vz))
                                         + A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P)
                                        + 2.0 * A * (Le * Le) * (Vz * Vz) * Oz * P * np.sinh(Le * Vz))
                                       + A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * np.sinh(Le * Vz))
                                      + A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz)
                                     + A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * Oz * np.sinh(Le * Vz)))
                             / (2.0 * A * (az * az)))
    
    et[3 , 9] = et[9 , 3] = -(K / Le) - (Gmod*Jx / Le)
    
    x = ((((((((((((((((((((2.0 * A * P * np.sinh(Le * Vz)
                            - A * P * np.sinh(2.0 * Le * Vz))
                           + 2.0 * Iy * (Vz * Vz) * P * np.sinh(Le * Vz))
                          - Iy * (Vz * Vz) * P * np.sinh(2.0 * Le * Vz))
                         - 6.0 * A * Le * Vz * P) + Iy * Le * np.power(Vz, 3.0) * P)
                       + A * Emod * Iy * Le * np.power(Vz, 3.0))
                      + 2.0 * A * Emod * Iy * (Vz * Vz) * np.sinh(Le * Vz))
                     - A * Emod * Iy * (Vz * Vz) * np.sinh(2.0 * Le * Vz))
                    - 2.0 * Iy * Le * np.power(Vz, 3.0) * P * np.cosh(Le * Vz))
                   + Iy * Le * np.power(Vz, 3.0) * P * np.cosh(2.0 * Le * Vz))
                  + 17.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P)
                 - 3.0 * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * P)
                + A * np.power(Le, 3.0) * np.power(Vz, 3.0) * P * np.cosh(Le * Vz))
               - Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * P * np.cosh(Le * Vz))
              - 3.0 * A * (Le * Le) * (Vz * Vz) * P * np.sinh(Le * Vz))
             + Iy * (Le * Le) * np.power(Vz, 4.0) * P * np.sinh(Le * Vz))
            + 6.0 * A * Le * Vz * P * np.cosh(Le * Vz))
           - 18.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P)
          + 9.0 * A * np.power(Le, 7.0) * np.power(Vz, 7.0) * np.power(Oz, 3.0) * P)
         - 2.0 * A * np.power(Le, 9.0) * np.power(Vz, 9.0) * np.power(Oz, 4.0) * P)
    
    et[4 , 10] = et[10 , 4] = ((((((((((((((((((((((((((((((((((((((((((((((((((x + cz)
                                                                               - 2.0 * dz * np.sinh(Le * Vz))
                                                                              + dz * np.sinh(2.0 * Le * Vz))
                                                                             + 8.0 * ez * np.sinh(Le * Vz))
                                                                            - 2.0 * fz * np.sinh(Le * Vz))
                                                                           + fz * np.sinh(2.0 * Le * Vz))
                                                                          - 2.0 * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * Oz * P * np.sinh(Le * Vz))
                                                                         + 7.0 * A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz))
                                                                        - 2.0 * A * Emod * Iy * np.power(Le, 7.0) * np.power(Vz, 9.0) * np.power(Oz, 3.0))
                                                                       + 16.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P * np.cosh(Le * Vz))
                                                                      + 2.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P * np.cosh(2.0 * Le * Vz))
                                                                     + A * np.power(Le, 7.0) * np.power(Vz, 7.0) * (Oz * Oz) * P * np.cosh(Le * Vz))
                                                                    - 8.0 * A * np.power(Le, 7.0) * np.power(Vz, 7.0) * np.power(Oz, 3.0) * P * np.cosh(Le * Vz))
                                                                   - A * np.power(Le, 7.0) * np.power(Vz, 7.0) * np.power(Oz, 3.0) * P * np.cosh(2.0 * Le * Vz))
                                                                  + 2.0 * A * np.power(Le, 9.0) * np.power(Vz, 9.0) * np.power(Oz, 4.0) * P * np.cosh(Le * Vz))
                                                                 - 2.0 * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz) * P * np.cosh(Le * Vz))
                                                                + 12.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * np.sinh(Le * Vz))
                                                               - 6.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * np.sinh(2.0 * Le * Vz))
                                                              - 7.0 * A * np.power(Le, 6.0) * np.power(Vz, 6.0) * (Oz * Oz) * P * np.sinh(Le * Vz))
                                                             - 8.0 * A * np.power(Le, 6.0) * np.power(Vz, 6.0) * np.power(Oz, 3.0) * P * np.sinh(Le * Vz))
                                                            + 4.0 * A * np.power(Le, 6.0) * np.power(Vz, 6.0) * np.power(Oz, 3.0) * P * np.sinh(2.0 * Le * Vz))
                                                           + 2.0 * A * np.power(Le, 8.0) * np.power(Vz, 8.0) * np.power(Oz, 3.0) * P * np.sinh(Le * Vz))
                                                          + 2.0 * A * np.power(Le, 8.0) * np.power(Vz, 8.0) * np.power(Oz, 4.0) * P * np.sinh(Le * Vz))
                                                         - A * np.power(Le, 8.0) * np.power(Vz, 8.0) * np.power(Oz, 4.0) * P * np.sinh(2.0 * Le * Vz))
                                                        + 2.0 * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * P * np.sinh(Le * Vz))
                                                       - Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * P * np.sinh(2.0 * Le * Vz))
                                                      - 2.0 * A * Emod * Iy * Le * np.power(Vz, 3.0) * np.cosh(Le * Vz))
                                                     + A * Emod * Iy * Le * np.power(Vz, 3.0) * np.cosh(2.0 * Le * Vz))
                                                    - 6.0 * A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz)
                                                   - A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * np.cosh(Le * Vz))
                                                  + A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * np.sinh(Le * Vz))
                                                 - 16.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P * np.cosh(Le * Vz))
                                                - A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P * np.cosh(2.0 * Le * Vz))
                                               - 2.0 * A * np.power(Le, 5.0) * np.power(Vz, 5.0) * Oz * P * np.cosh(Le * Vz))
                                              + 4.0 * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * P * np.cosh(Le * Vz))
                                             - Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * P * np.cosh(2.0 * Le * Vz))
                                            + 6.0 * A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * np.cosh(Le * Vz))
                                           + A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * Oz * np.cosh(Le * Vz))
                                          - 2.0 * A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * Oz * np.sinh(Le * Vz))
                                         + A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0) * Oz * np.sinh(2.0 * Le * Vz))
                                        - 5.0 * A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * Oz * np.sinh(Le * Vz))
                                       - 6.0 * A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz) * np.cosh(Le * Vz))
                                      - A * Emod * Iy * np.power(Le, 5.0) * np.power(Vz, 7.0) * (Oz * Oz) * np.cosh(2.0 * Le * Vz))
                                     + 2.0 * A * Emod * Iy * np.power(Le, 7.0) * np.power(Vz, 9.0) * np.power(Oz, 3.0) * np.cosh(Le * Vz))
                                    - 2.0 * A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * np.sinh(Le * Vz))
                                   + A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * (Oz * Oz) * np.sinh(2.0 * Le * Vz))
                                  + 2.0 * A * Emod * Iy * np.power(Le, 6.0) * np.power(Vz, 8.0) * (Oz * Oz) * np.sinh(Le * Vz))
                                 + 2.0 * A * Emod * Iy * np.power(Le, 6.0) * np.power(Vz, 8.0) * np.power(Oz, 3.0) * np.sinh(Le * Vz))
                                - A * Emod * Iy * np.power(Le, 6.0) * np.power(Vz, 8.0) * np.power(Oz, 3.0) * np.sinh(2.0 * Le * Vz))
                               / (A * Vz * ((((((((((((((4.0 * np.cosh(2.0 * Le * Vz)
                                                         - 16.0 * np.cosh(Le * Vz)) - Le * Le * (Vz * Vz))
                                                       - 24.0 * (Le * Le) * (Vz * Vz) * Oz)
                                                      + Le * Le * (Vz * Vz) * np.cosh(2.0 * Le * Vz))
                                                     + 8.0 * Le * Vz * np.sinh(Le * Vz))
                                                    - 4.0 * Le * Vz * np.sinh(2.0 * Le * Vz))
                                                   + 12.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz))
                                                  + 32.0 * (Le * Le) * (Vz * Vz) * Oz * np.cosh(Le * Vz))
                                                 - 8.0 * (Le * Le) * (Vz * Vz) * Oz * np.cosh(2.0 * Le * Vz))
                                                - 8.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sinh(Le * Vz))
                                               + 4.0 * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * np.sinh(2.0 * Le * Vz))
                                              - 16.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cosh(Le * Vz))
                                             + 4.0 * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * np.cosh(2.0 * Le * Vz)) + 12.0)))
    
    x = ((((((((((((((((((((2.0 * A * P * np.sinh(Le * Vy)
                            - A * P * np.sinh(2.0 * Le * Vy))
                           + 2.0 * Iz * (Vy * Vy) * P * np.sinh(Le * Vy))
                          - Iz * (Vy * Vy) * P * np.sinh(2.0 * Le * Vy))
                         - 6.0 * A * Le * Vy * P) + Iz * Le * np.power(Vy, 3.0) * P)
                       + A * Emod * Iz * Le * np.power(Vy, 3.0))
                      + 2.0 * A * Emod * Iz * (Vy * Vy) * np.sinh(Le * Vy))
                     - A * Emod * Iz * (Vy * Vy) * np.sinh(2.0 * Le * Vy))
                    - 2.0 * Iz * Le * np.power(Vy, 3.0) * P * np.cosh(Le * Vy))
                   + Iz * Le * np.power(Vy, 3.0) * P * np.cosh(2.0 * Le * Vy))
                  + 17.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P)
                 - 3.0 * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * P)
                + A * np.power(Le, 3.0) * np.power(Vy, 3.0) * P * np.cosh(Le * Vy))
               - Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * P * np.cosh(Le * Vy))
              - 3.0 * A * (Le * Le) * (Vy * Vy) * P * np.sinh(Le * Vy))
             + Iz * (Le * Le) * np.power(Vy, 4.0) * P * np.sinh(Le * Vy))
            + 6.0 * A * Le * Vy * P * np.cosh(Le * Vy))
           - 18.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P)
          + 9.0 * A * np.power(Le, 7.0) * np.power(Vy, 7.0) * np.power(Oy, 3.0) * P)
         - 2.0 * A * np.power(Le, 9.0) * np.power(Vy, 9.0) * np.power(Oy, 4.0) * P)
    
    et[5 , 11] = et[11 , 5] = ((((((((((((((((((((((((((((((((((((((((((((((((((x + cy) - 2.0 * dy * np.sinh(Le * Vy))
                                                                              + dy * np.sinh(2.0 * Le * Vy))
                                                                             + 8.0 * ey * np.sinh(Le * Vy))
                                                                            - 2.0 * fy * np.sinh(Le * Vy))
                                                                           + fy * np.sinh(2.0 * Le * Vy))
                                                                          - 2.0 * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * Oy * P * np.sinh(Le * Vy))
                                                                         + 7.0 * A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy))
                                                                        - 2.0 * A * Emod * Iz * np.power(Le, 7.0) * np.power(Vy, 9.0) * np.power(Oy, 3.0))
                                                                       + 16.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P * np.cosh(Le * Vy))
                                                                      + 2.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P * np.cosh(2.0 * Le * Vy))
                                                                     + A * np.power(Le, 7.0) * np.power(Vy, 7.0) * (Oy * Oy) * P * np.cosh(Le * Vy))
                                                                    - 8.0 * A * np.power(Le, 7.0) * np.power(Vy, 7.0) * np.power(Oy, 3.0) * P * np.cosh(Le * Vy))
                                                                   - A * np.power(Le, 7.0) * np.power(Vy, 7.0) * np.power(Oy, 3.0) * P * np.cosh(2.0 * Le * Vy))
                                                                  + 2.0 * A * np.power(Le, 9.0) * np.power(Vy, 9.0) * np.power(Oy, 4.0) * P * np.cosh(Le * Vy))
                                                                 - 2.0 * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy) * P * np.cosh(Le * Vy))
                                                                + 12.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * np.sinh(Le * Vy))
                                                               - 6.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * np.sinh(2.0 * Le * Vy))
                                                              - 7.0 * A * np.power(Le, 6.0) * np.power(Vy, 6.0) * (Oy * Oy) * P * np.sinh(Le * Vy))
                                                             - 8.0 * A * np.power(Le, 6.0) * np.power(Vy, 6.0) * np.power(Oy, 3.0) * P * np.sinh(Le * Vy))
                                                            + 4.0 * A * np.power(Le, 6.0) * np.power(Vy, 6.0) * np.power(Oy, 3.0) * P * np.sinh(2.0 * Le * Vy))
                                                           + 2.0 * A * np.power(Le, 8.0) * np.power(Vy, 8.0) * np.power(Oy, 3.0) * P * np.sinh(Le * Vy))
                                                          + 2.0 * A * np.power(Le, 8.0) * np.power(Vy, 8.0) * np.power(Oy, 4.0) * P * np.sinh(Le * Vy))
                                                         - A * np.power(Le, 8.0) * np.power(Vy, 8.0) * np.power(Oy, 4.0) * P * np.sinh(2.0 * Le * Vy))
                                                        + 2.0 * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * P * np.sinh(Le * Vy))
                                                       - Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * P * np.sinh(2.0 * Le * Vy))
                                                      - 2.0 * A * Emod * Iz * Le * np.power(Vy, 3.0) * np.cosh(Le * Vy))
                                                     + A * Emod * Iz * Le * np.power(Vy, 3.0) * np.cosh(2.0 * Le * Vy))
                                                    - 6.0 * A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy)
                                                   - A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * np.cosh(Le * Vy))
                                                  + A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * np.sinh(Le * Vy))
                                                 - 16.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P * np.cosh(Le * Vy))
                                                - A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P * np.cosh(2.0 * Le * Vy))
                                               - 2.0 * A * np.power(Le, 5.0) * np.power(Vy, 5.0) * Oy * P * np.cosh(Le * Vy))
                                              + 4.0 * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * P * np.cosh(Le * Vy))
                                             - Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * P * np.cosh(2.0 * Le * Vy))
                                            + 6.0 * A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * np.cosh(Le * Vy))
                                           + A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * Oy * np.cosh(Le * Vy))
                                          - 2.0 * A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * Oy * np.sinh(Le * Vy))
                                         + A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0) * Oy * np.sinh(2.0 * Le * Vy))
                                        - 5.0 * A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * Oy * np.sinh(Le * Vy))
                                       - 6.0 * A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy) * np.cosh(Le * Vy))
                                      - A * Emod * Iz * np.power(Le, 5.0) * np.power(Vy, 7.0) * (Oy * Oy) * np.cosh(2.0 * Le * Vy))
                                     + 2.0 * A * Emod * Iz * np.power(Le, 7.0) * np.power(Vy, 9.0) * np.power(Oy, 3.0) * np.cosh(Le * Vy))
                                    - 2.0 * A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * np.sinh(Le * Vy))
                                   + A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * (Oy * Oy) * np.sinh(2.0 * Le * Vy))
                                  + 2.0 * A * Emod * Iz * np.power(Le, 6.0) * np.power(Vy, 8.0) * (Oy * Oy) * np.sinh(Le * Vy))
                                 + 2.0 * A * Emod * Iz * np.power(Le, 6.0) * np.power(Vy, 8.0) * np.power(Oy, 3.0) * np.sinh(Le * Vy))
                                - A * Emod * Iz * np.power(Le, 6.0) * np.power(Vy, 8.0) * np.power(Oy, 3.0) * np.sinh(2.0 * Le * Vy))
                               / (A * Vy * ((((((((((((((4.0 * np.cosh(2.0 * Le * Vy) - 16.0 * np.cosh(Le * Vy))
                                                        - Le * Le * (Vy * Vy)) - 24.0 * (Le * Le) * (Vy * Vy) * Oy)
                                                      + Le * Le * (Vy * Vy) * np.cosh(2.0 * Le * Vy)) + 8.0 * Le * Vy * np.sinh(Le * Vy))
                                                    - 4.0 * Le * Vy * np.sinh(2.0 * Le * Vy))
                                                   + 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy))
                                                  + 32.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy))
                                                 - 8.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(2.0 * Le * Vy))
                                                - 8.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sinh(Le * Vy))
                                               + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sinh(2.0 * Le * Vy))
                                              - 16.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cosh(Le * Vy))
                                             + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * np.cosh(2.0 * Le * Vy)) + 12.0)))

    et[2 , 9] = et[9 , 2] = Mzb / Le
    
    et[3 , 10] = et[10 , 3] = (Oz * (Le * Le) * (Vz * Vz) - 1.0) * (Le * Vz * np.cosh(Le * Vz / 2.0) / np.sinh(Le * Vz / 2.0) - 2.0) * (Mza + Mzb) / (Le * Le * (Vz * Vz))
    
    et[4 , 11] = et[11 , 4] = 0

    et[1 , 9] = et[9 , 1] = Myb / Le
    
    et[2 , 10] = et[10 , 2] = (-((((((((((((((A * (Le * Le) * (Vz * Vz) * P - 8.0 * A * P * (sz * sz))
                                             - Iy * (Le * Le) * np.power(Vz, 4.0) * P)
                                            + Iy * Le * np.power(Vz, 3.0) * P * np.sinh(Le * Vz))
                                           - A * Emod * Iy * (Le * Le) * np.power(Vz, 4.0))
                                          - 2.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * Oz * P)
                                         + A * Le * Vz * P * np.sinh(Le * Vz))
                                        + A * np.power(Le, 6.0) * np.power(Vz, 6.0) * (Oz * Oz) * P)
                                       - 2.0 * A * np.power(Le, 3.0) * np.power(Vz, 3.0) * Oz * P * np.sinh(Le * Vz))
                                      + 16.0 * A * (Le * Le) * (Vz * Vz) * Oz * P * (sz * sz))
                                     + A * np.power(Le, 5.0) * np.power(Vz, 5.0) * (Oz * Oz) * P * np.sinh(Le * Vz))
                                    + A * Emod * Iy * Le * np.power(Vz, 3.0) * np.sinh(Le * Vz))
                                   - 8.0 * A * np.power(Le, 4.0) * np.power(Vz, 4.0) * (Oz * Oz) * P * (sz * sz))
                                  + A * Emod * Iy * np.power(Le, 4.0) * np.power(Vz, 6.0) * Oz)
                                 + A * Emod * Iy * np.power(Le, 3.0) * np.power(Vz, 5.0) * Oz * np.sinh(Le * Vz))
                               / (4.0 * A * (bz * bz)))
    
    et[3 , 11] = et[11 , 3] = -((Oy * (Le * Le) * (Vy * Vy) - 1.0) * (Le * Vy * np.cosh(Le * Vy / 2.0) / np.sinh(Le * Vy / 2.0) - 2.0) * (Mya + Myb)) / (Le * Le * (Vy * Vy))

    et[1 , 10] = et[10 , 1] = -(Vy * Vy * Mxb * (Vy * np.cosh(Le * Vy / 2.0) * np.sinh(Le * Vz / 2.0) - Vz * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vy / 2.0))) / (np.sinh(Le * Vz / 2.0) * (Vy * Vy - Vz * Vz) * ((Le * Vy * np.cosh(Le * Vy / 2.0) - 2.0 * np.sinh(Le * Vy / 2.0)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.sinh(Le * Vy / 2.0)))
    
    et[2 , 11] = et[11 , 2] = -(Vz * Vz * Mxb * (Vy * np.cosh(Le * Vy / 2.0) * np.sinh(Le * Vz / 2.0) - Vz * np.cosh(Le * Vz / 2.0) * np.sinh(Le * Vy / 2.0))) / (np.sinh(Le * Vy / 2.0) * (Vy * Vy - Vz * Vz) * ((Le * Vz * np.cosh(Le * Vz / 2.0) - 2.0 * np.sinh(Le * Vz / 2.0)) + 2.0 * (Le * Le) * (Vz * Vz) * Oz * np.sinh(Le * Vz / 2.0)))

    et[0 , 10] = et[10 , 0] = -Myb / Le
    
    et[1 , 11] = et[11 , 1] = (((((((((((((((A * (Le * Le) * (Vy * Vy) * P
                                             - 8.0 * A * P * (sy * sy))
                                            - Iz * (Le * Le) * np.power(Vy, 4.0) * P)
                                           + Iz * Le * np.power(Vy, 3.0) * P * np.sinh(Le * Vy))
                                          - A * Emod * Iz * (Le * Le) * np.power(Vy, 4.0))
                                         - 2.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * P)
                                        + A * Le * Vy * P * np.sinh(Le * Vy))
                                       + A * np.power(Le, 6.0) * np.power(Vy, 6.0) * (Oy * Oy) * P)
                                      - 2.0 * A * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * P * np.sinh(Le * Vy))
                                     + 16.0 * A * (Le * Le) * (Vy * Vy) * Oy * P * (sy * sy))
                                    + A * np.power(Le, 5.0) * np.power(Vy, 5.0) * (Oy * Oy) * P * np.sinh(Le * Vy))
                                   + A * Emod * Iz * Le * np.power(Vy, 3.0) * np.sinh(Le * Vy))
                                  - 8.0 * A * np.power(Le, 4.0) * np.power(Vy, 4.0) * (Oy * Oy) * P * (sy * sy))
                                 + A * Emod * Iz * np.power(Le, 4.0) * np.power(Vy, 6.0) * Oy)
                                + A * Emod * Iz * np.power(Le, 3.0) * np.power(Vy, 5.0) * Oy * np.sinh(Le * Vy))
                               / (4.0 * A * (by * by)))

    et[0 , 11] = et[11 , 0] = -Mzb / Le

    s_y = np.sinh(Le * Vy / 4.0)
    s_z = np.sinh(Le * Vz / 4.0)
    
    et[4 , 5] = et[5 , 4] = (((Le * Le * ((2.0 * (Vy * Vy) * (Vz * Vz) * Mxb *
                                           (((sy * sy - 2.0 * (sy * sy + 1.0) * (sz * sz + 1.0)) + sz * sz) + 2.0)
                                           + Vy * np.power(Vz, 3.0) * Mxb * np.sinh(Le * Vy) * np.sinh(Le * Vz) / 2.0)
                                          + np.power(Vy, 3.0) * Vz * Mxb * np.sinh(Le * Vy) * np.sinh(Le * Vz) / 2.0)
                               + Le * (((2.0 * np.power(Vy, 3.0) * Mxb * (np.sinh(Le * Vy) - 2.0 * np.sinh(Le * Vy / 2.0)
                                                                          * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0))
                                         + 2.0 * np.power(Vz, 3.0) * Mxb * (np.sinh(Le * Vz) - 2.0 * np.sinh(Le * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)))
                                        - 2.0 * Vy * (Vz * Vz) * Mxb * (np.sinh(Le * Vy) - 2.0 * np.sinh(Le * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)))
                                       - 2.0 * (Vy * Vy) * Vz * Mxb * (np.sinh(Le * Vz) - 2.0 * np.sinh(Le * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0))))
                              + np.power(Le, 3.0) * (2.0 * np.power(Vy, 3.0) * (Vz * Vz) * Mxb * (((Oy * np.sinh(Le * Vy) - Oz * np.sinh(Le * Vy)) - 2.0 * Oy * np.sinh(Le * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0))
                                                                                                  + 2.0 * Oz * np.sinh(Le * Vy / 2.0) * (sy * sy + 1.0) * (2.0 * (s_y * s_y) + 1.0))
                                                     - 2.0 * (Vy * Vy) * np.power(Vz, 3.0) * Mxb * (((Oy * np.sinh(Le * Vz) - Oz * np.sinh(Le * Vz))
                                                                                                     - 2.0 * Oy * np.sinh(Le * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0))
                                                                                                    + 2.0 * Oz * np.sinh(Le * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0))))
                             / ((Vy * Vy - Vz * Vz) * ((Le * Vy * np.sinh(Le * Vy) - 4.0 * (sy * sy))
                                                       + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy))
                                * ((Le * Vz * np.sinh(Le * Vz) - 4.0 * (sz * sz)) + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz))))
    
    et[10 , 11] = et[11 , 10] = (-((Le * Le * ((2.0 * (Vy * Vy) * (Vz * Vz) * Mxb
                                                * (((sy * sy - 2.0 * (sy * sy + 1.0) * (sz * sz + 1.0)) + sz * sz) + 2.0)
                                                + Vy * np.power(Vz, 3.0) * Mxb * np.sinh(Le * Vy) * np.sinh(Le * Vz) / 2.0)
                                               + np.power(Vy, 3.0) * Vz * Mxb * np.sinh(Le * Vy) * np.sinh(Le * Vz) / 2.0)
                                    + Le * (((2.0 * np.power(Vy, 3.0) * Mxb * (np.sinh(Le * Vy) - 2.0 * np.sinh(Le * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0))
                                              + 2.0 * np.power(Vz, 3.0) * Mxb * (np.sinh(Le * Vz) - 2.0 * np.sinh(Le * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0)))
                                             - 2.0 * Vy * (Vz * Vz) * Mxb * (np.sinh(Le * Vy) - 2.0 * np.sinh(Le * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0)))
                                            - 2.0 * (Vy * Vy) * Vz * Mxb * (np.sinh(Le * Vz) - 2.0 * np.sinh(Le * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0))))
                                   + np.power(Le, 3.0) * (2.0 * np.power(Vy, 3.0) * (Vz * Vz) * Mxb * (((Oy * np.sinh(Le * Vy) - Oz * np.sinh(Le * Vy))
                                                                                                        - 2.0 * Oy * np.sinh(Le * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0))
                                                                                                       + 2.0 * Oz * np.sinh(Le * Vy / 2.0) * (sz * sz + 1.0) * (2.0 * (s_y * s_y) + 1.0))
                                                          - 2.0 * (Vy * Vy) * np.power(Vz, 3.0) * Mxb * (((Oy * np.sinh(Le * Vz) - Oz * np.sinh(Le * Vz))
                                                                                                          - 2.0 * Oy * np.sinh(Le * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0))
                                                                                                         + 2.0 * Oz * np.sinh(Le * Vz / 2.0) * (sy * sy + 1.0) * (2.0 * (s_z * s_z) + 1.0))))
                                 / ((Vy * Vy - Vz * Vz) * ((Le * Vy * np.sinh(Le * Vy)
                                                            - 4.0 * (sy * sy))
                                                           + 4.0 * (Le * Le) * (Vy * Vy) * Oy * (sy * sy))
                                    * ((Le * Vz * np.sinh(Le * Vz) - 4.0 * (sz * sz))
                                       + 4.0 * (Le * Le) * (Vz * Vz) * Oz * (sz * sz))))
    #
    # Torsion for symmetric sections
    if Iy == Iz:
        et[1 , 4] = et[4 , 1] = Vy * Mxb * (np.sinh(Le * Vy) - Le * Vy) / (2.0 * ((((Le * Vy * np.sinh(Le * Vy) - 2.0 * (Le * Le) * (Vy * Vy) * Oy) - 2.0 * np.cosh(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy)) + 2.0))
        et[2 , 5] = et[5 , 2] = Vy * Mxb * (np.sinh(Le * Vy) - Le * Vy) / (2.0 * ((((Le * Vy * np.sinh(Le * Vy) - 2.0 * (Le * Le) * (Vy * Vy) * Oz) - 2.0 * np.cosh(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oz * np.cosh(Le * Vy)) + 2.0))
        et[4 , 7] = et[7 , 4] = -(Vy * Mxb * (np.sinh(Le * Vy) - Le * Vy)) / (2.0 * ((((Le * Vy * np.sinh(Le * Vy) - 2.0 * (Le * Le) * (Vy * Vy) * Oy) - 2.0 * np.cosh(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy)) + 2.0))
        et[5 , 8] = et[8 , 5] = -(Vy * Mxb * (np.sinh(Le * Vy) - Le * Vy)) / (2.0 * ((((Le * Vy * np.sinh(Le * Vy) - 2.0 * (Le * Le) * (Vy * Vy) * Oz) - 2.0 * np.cosh(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oz * np.cosh(Le * Vy)) + 2.0))
        et[7 , 10] = et[10 , 7] = Vy * Mxb * (np.sinh(Le * Vy) - Le * Vy) / (2.0 * ((((Le * Vy * np.sinh(Le * Vy) - 2.0 * (Le * Le) * (Vy * Vy) * Oy) - 2.0 * np.cosh(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy)) + 2.0))
        et[8 , 11] = et[11 , 8] = Vy * Mxb * (np.sinh(Le * Vy) - Le * Vy) / (2.0 * ((((Le * Vy * np.sinh(Le * Vy) - 2.0 * (Le * Le) * (Vy * Vy) * Oz) - 2.0 * np.cosh(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oz * np.cosh(Le * Vy)) + 2.0))

        et[5 , 10] = et[10 , 5] = (-(Mxb * ((((((((((((((((((((((((6.0 * Vy - 8.0 * Vy * np.cosh(Le * Vy))
                                                                  + 2.0 * Vy * np.cosh(2.0 * Le * Vy))
                                                                 + 2.0 * (Le * Le) * np.power(Vy, 3.0))
                                                                - 6.0 * (Le * Le) * np.power(Vy, 3.0) * Oy)
                                                               - 6.0 * (Le * Le) * np.power(Vy, 3.0) * Oz)
                                                              - np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy)
                                                             - np.power(Le,4.0) * np.power(Vy, 5.0) * Oz)
                                                            - 2.0 * (Le * Le) * np.power(Vy, 3.0) * np.cosh(Le * Vy))
                                                           + np.power(Le, 3.0) * np.power(Vy, 4.0) * np.sinh(Le * Vy))
                                                          + 2.0 * Le * (Vy * Vy) * np.sinh(Le * Vy))
                                                         - Le * (Vy * Vy) * np.sinh(2.0 * Le * Vy))
                                                        + 6.0 * np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy * Oz)
                                                       + 8.0 * (Le * Le) * np.power(Vy, 3.0) * Oy * np.cosh(Le * Vy))
                                                      + 8.0 * (Le * Le) * np.power(Vy,3.0) * Oz * np.cosh(Le * Vy))
                                                     - 2.0 * (Le * Le) * np.power(Vy, 3.0) * Oy * np.cosh(2.0 * Le * Vy))
                                                    - 2.0 * (Le * Le) * np.power(Vy, 3.0) * Oz * np.cosh(2.0 * Le * Vy))
                                                   + np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy * np.cosh(Le * Vy))
                                                  + np.power(Le, 4.0) * np.power(Vy, 5.0) * Oz * np.cosh(Le * Vy))
                                                 - np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy * np.sinh(Le * Vy))
                                                - np.power(Le, 3.0) * np.power(Vy, 4.0) * Oz * np.sinh(Le * Vy))
                                               + np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy * np.sinh(2.0 * Le * Vy) / 2.0)
                                              + np.power(Le, 3.0) * np.power(Vy, 4.0) * Oz * np.sinh(2.0 * Le * Vy) / 2.0)
                                             - 8.0 * np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy * Oz * np.cosh(Le * Vy))
                                            + 2.0 * np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy * Oz * np.cosh(2.0 * Le * Vy)))
                                   / (Vy * (((((((((((((((((((16.0 * np.cosh(Le * Vy) - 4.0 * np.cosh(2.0 * Le * Vy))
                                                             + Le * Le * (Vy * Vy)) + 12.0 * (Le * Le) * (Vy * Vy) * Oy)
                                                           + 12.0 * (Le * Le) * (Vy * Vy) * Oz)
                                                          - Le * Le * (Vy * Vy) * np.cosh(2.0 * Le * Vy))
                                                         - 8.0 * Le * Vy * np.sinh(Le * Vy))
                                                        + 4.0 * Le * Vy * np.sinh(2.0 * Le * Vy))
                                                       - 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz)
                                                      - 16.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy))
                                                     - 16.0 * (Le * Le) * (Vy * Vy) * Oz * np.cosh(Le * Vy))
                                                    + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(2.0 * Le * Vy))
                                                   + 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.cosh(2.0 * Le * Vy))
                                                  + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sinh(Le * Vy))
                                                 + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sinh(Le * Vy))
                                                - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sinh(2.0 * Le * Vy))
                                               - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sinh(2.0 * Le * Vy))
                                              + 16.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cosh(Le * Vy))
                                             - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cosh(2.0 * Le * Vy)) - 12.0)))
        
        et[4 , 11] = et[11 , 4] = (Mxb * ((((((((((((((((((((((((6.0 * Vy - 8.0 * Vy * np.cosh(Le * Vy))
                                                                + 2.0 * Vy * np.cosh(2.0 * Le * Vy))
                                                               + 2.0 * (Le * Le) * np.power(Vy, 3.0))
                                                              - 6.0 * (Le * Le) * np.power(Vy, 3.0) * Oy)
                                                             - 6.0 * (Le * Le) * np.power(Vy, 3.0) * Oz)
                                                            - np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy)
                                                           - np.power(Le,4.0) * np.power(Vy, 5.0) * Oz)
                                                          - 2.0 * (Le * Le) * np.power(Vy, 3.0) * np.cosh(Le * Vy))
                                                         + np.power(Le, 3.0) * np.power(Vy, 4.0) * np.sinh(Le * Vy))
                                                        + 2.0 * Le * (Vy * Vy) * np.sinh(Le * Vy))
                                                       - Le * (Vy * Vy) * np.sinh(2.0 * Le * Vy))
                                                      + 6.0 * np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy * Oz)
                                                     + 8.0 * (Le * Le) * np.power(Vy, 3.0) * Oy * np.cosh(Le * Vy))
                                                    + 8.0 * (Le * Le) * np.power(Vy,3.0) * Oz * np.cosh(Le * Vy))
                                                   - 2.0 * (Le * Le) * np.power(Vy, 3.0) * Oy * np.cosh(2.0 * Le * Vy))
                                                  - 2.0 * (Le * Le) * np.power(Vy, 3.0) * Oz * np.cosh(2.0 * Le * Vy))
                                                 + np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy * np.cosh(Le * Vy))
                                                + np.power(Le, 4.0) * np.power(Vy, 5.0) * Oz * np.cosh(Le * Vy))
                                               - np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy * np.sinh(Le * Vy))
                                              - np.power(Le, 3.0) * np.power(Vy, 4.0) * Oz * np.sinh(Le * Vy))
                                             + np.power(Le, 3.0) * np.power(Vy, 4.0) * Oy * np.sinh(2.0 * Le * Vy) / 2.0)
                                            + np.power(Le, 3.0) * np.power(Vy, 4.0) * Oz * np.sinh(2.0 * Le * Vy) / 2.0)
                                           - 8.0 * np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy * Oz * np.cosh(Le * Vy))
                                          + 2.0 * np.power(Le, 4.0) * np.power(Vy, 5.0) * Oy * Oz * np.cosh(2.0 * Le * Vy))
                                   / (Vy * (((((((((((((((((((16.0 * np.cosh(Le * Vy)
                                                              - 4.0 * np.cosh(2.0 * Le * Vy)) + Le * Le * (Vy * Vy))
                                                            + 12.0 * (Le * Le) * (Vy * Vy) * Oy)
                                                           + 12.0 * (Le * Le) * (Vy * Vy) * Oz)
                                                          - Le * Le * (Vy * Vy) * np.cosh(2.0 * Le * Vy))
                                                         - 8.0 * Le * Vy * np.sinh(Le * Vy))
                                                        + 4.0 * Le * Vy * np.sinh(2.0 * Le * Vy))
                                                       - 12.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz)
                                                      - 16.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy))
                                                     - 16.0 * (Le * Le) * (Vy * Vy) * Oz * np.cosh(Le * Vy))
                                                    + 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(2.0 * Le * Vy))
                                                   + 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.cosh(2.0 * Le * Vy))
                                                  + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sinh(Le * Vy))
                                                 + 4.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sinh(Le * Vy))
                                                - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sinh(2.0 * Le * Vy))
                                               - 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sinh(2.0 * Le * Vy))
                                              + 16.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cosh(Le * Vy))
                                             - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cosh(2.0 * Le * Vy)) - 12.0)))

        et[1 , 10] = et[10 , 1] = -(Vy * Mxb * (np.sinh(Le * Vy) - Le * Vy)) / (2.0 * ((((Le * Vy * np.sinh(Le * Vy) - 2.0 * (Le * Le) * (Vy * Vy) * Oy) - 2.0 * np.cosh(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy)) + 2.0))
        et[2 , 11] = et[11 , 2] = -(Vy * Mxb * (np.sinh(Le * Vy) - Le * Vy)) / (2.0 * ((((Le * Vy * np.sinh(Le * Vy) - 2.0 * (Le * Le) * (Vy * Vy) * Oz) - 2.0 * np.cosh(Le * Vy)) + 2.0 * (Le * Le) * (Vy * Vy) * Oz * np.cosh(Le * Vy)) + 2.0))

        et[4 , 5] = et[5 , 4] = -(np.power(Le, 3.0) * np.power(Vy, 3.0) * Mxb * (Oy - Oz) * (np.sinh(Le * Vy) - Le * Vy)) / (2.0 * ((((((((((((4.0 * np.cosh(Le * Vy) + Le * Le * (Vy * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy) + 4.0 * (Le * Le) * (Vy * Vy) * Oz) + Le * Le * (Vy * Vy) * np.cosh(Le * Vy)) - 4.0 * Le * Vy * np.sinh(Le * Vy)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz) - 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy)) - 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.cosh(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sinh(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sinh(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cosh(Le * Vy)) - 4.0))
        et[10 , 11] = et[11 , 10] = np.power(Le, 3.0) * np.power(Vy, 3.0) * Mxb * (Oy - Oz) * (np.sinh(Le * Vy) - Le * Vy) / (2.0 * ((((((((((((4.0 * np.cosh(Le * Vy) + Le * Le * (Vy * Vy)) + 4.0 * (Le * Le) * (Vy * Vy) * Oy) + 4.0 * (Le * Le) * (Vy * Vy) * Oz) + Le * Le * (Vy * Vy) * np.cosh(Le * Vy)) - 4.0 * Le * Vy * np.sinh(Le * Vy)) - 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz) - 4.0 * (Le * Le) * (Vy * Vy) * Oy * np.cosh(Le * Vy)) - 4.0 * (Le * Le) * (Vy * Vy) * Oz * np.cosh(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oy * np.sinh(Le * Vy)) + 2.0 * np.power(Le, 3.0) * np.power(Vy, 3.0) * Oz * np.sinh(Le * Vy)) + 4.0 * np.power(Le, 4.0) * np.power(Vy, 4.0) * Oy * Oz * np.cosh(Le * Vy)) - 4.0))
    #
    return et
 #
 #
 #
    