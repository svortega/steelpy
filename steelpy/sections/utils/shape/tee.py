# 
# Copyright (c) 2009 steelpy
#
#
# Python stdlib imports
from __future__ import annotations
#from array import array
from collections import namedtuple
from dataclasses import dataclass
import math
#import re
#
#
# package imports
from steelpy.sections.utils.shape.utils import ShapeProperty
from steelpy.sections.utils.shape.stress import ShapeStressBasic
#
#
#
# ----------------------------------------
#
points = namedtuple('Points', ['y', 'z'])
axis = namedtuple('Axis', ['y', 'z'])
#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------
#
#
@dataclass
class TeeBasic(ShapeStressBasic):
    """
    Calculate the section properties of a T section

        *  b  *
    +   +-----+
           |         
    d      |       Z
           |       ^
    +      +       + > Y

    Parameters
    ----------
    d  : Section Heigh
    tw : Web thickness
    b  : Base
    tf : Flange thickness

    Returns
    ----------
    area: Section area
    Zc  : Elastic neutral centre
    Yc  : Elastic neutral centre
    Iy  : Second moment of area about mayor axis
    Zey : Elastic modulus about mayor axis
    Zpy : Plastic modulus about mayor axis
    SFy : Shape factor mayor axis
    ry  : Radius of gyration about mayor Axis
    Iz  : Second moment of area about minor axis
    Zez : Elastic modulus about minor axis
    Zpz : Plastic modulus about minor axis
    SFz : Shape factor minor axis
    rz  : Radius of gyration about minor Axis
    SC  : Shear centre
    Cw  : Warping constant

    Notes
    ----------
    Uses formulas from:
    1.- Formulas for stress, strain and strucutral matrices [W.D. Pilkey]
    2.- Roark's formulas for stress and strain [7th Edition]
    3.- Wikipedia

    Examples
    ----------

    """
    d: float
    tw: float
    b: float
    tb: float
    shape:str = 'T Section'
    #
    #
    def _stressX(self, actions, stress=None, stress_type: str='average'):
        """
        """
        # get section's coordinates
        coord =  self.section_coordinates()
        prop = self.properties()
        #
        # ----------------------------------------------
        # FIXME: torsion
        tau_x = [actions.Mx * 0
                 for item in coord.y]
        #
        # In Plane
        tau_y = [actions.Fy / prop.area
                 for item in coord.y]
        # Out Plane
        tau_z = [actions.Fz / prop.area
                 for item in coord.z]
        #
        # get bending stress
        sigma_x = [actions.Fx / prop.area for item in coord.y]
        sigma_y = [actions.My * item / prop.Iy for item in coord.z]
        sigma_z = [actions.Mz * item / prop.Iz for item in coord.y]
        #
        stress_out = BeamStress(sigma_x, sigma_y, sigma_z, 
                                tau_x, tau_y, tau_z, coord)
        #
        if stress:
            stress_out = self.add_stress(stress=stress, other=stress_out)
        #
        return stress_out         
    #  
    #
    def _properties(self, poisson):
        """
        """
        #
        _C = self.b - self.tw
        _h = self.d - self.tb / 2.0
        _D2 = self.d - self.tb
        
        #-------------------------------------------------
        #   Cross-Sectional Area
        area = self.b * self.tb + self.tw * _D2
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        #Zc = (((self.d**2 * self.tw) + (_C * self.tb**2)) /
        #           (2 * (self.b*self.tb + _D2*self.tw)))
        #Yc = 0
        Yc, Zc = self.centroid
        #-------------------------------------------------
        #   Shear Centre 
        SCz = self.tb / 2.0
        SCy = 0
        #-------------------------------------------------
        #   Warping Constant Cw
        Cw = ((self.tb**3 * self.b**3 / 144.0)
                   + (self.tw**3 * _h**3 / 36.0))
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        Iy = ((self.tw * (self.d - Zc)**3 + self.b * Zc**3 -
                    _C * (Zc - self.tb)**3) / 3.0)
        #   Elastic Modulus about Mayor Axis
        Zey = (Iy / (self.d - Zc))
        #   Plastic Modulus about Mayor Axis
        if self.tw * _D2 > self.b * self.tb :
            Zpy = ((_D2**2 * self.tw / 4.0) -
                        (self.b**2 * self.tb**2 / (4.0 * self.tw)) + 
                        (self.b * self.tb * self.d / 2.0))
        else:
            Zpy = ((self.tb**2 * self.b / 4.0) +
                        0.50 * self.tw 
                        * _D2*(self.d - self.tw*_D2 / (2 * self.b)))
        #   Shape Factor
        SFy = (Zpy * (self.d - Zc) / Iy)
        #   Radius of gyration about Mayor Axis
        ry = (Iy / area)**0.50
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        Iz = (self.b**3 * self.tb + _D2 * self.tw**3) / 12.0
        #   Elastic Modulus about Minor Axis
        Zez = 2 * Iz / self.b
        #   Plastic Modulus about Minor Axis
        Zpz = (self.b**2 * self.tb + self.tw**2 * _D2) / 4.0
        #   Shape Factor
        SFz = Zpz * self.b / (2 * Iz)
        #   Radius of gyration about Minor Axis 
        rz = ( Iz / area)**0.50
        #-------------------------------------------------
        #   Torsional Constant
        J = (self.b * self.tb**3 + _h * self.tw**3) / 3.0
        #   Product of inertia
        _Iyz = 0
        Jx = Iy + Iz
        rp = (Jx / area)**0.50
        #
        #
        alpha_sy, alpha_sz = self.alpha_s(poisson=poisson)
        #
        return ShapeProperty(area=area, Zc=Zc, Yc=Yc,
                             Iy=Iy, Sy=Zey, Zy=Zpy, ry=ry,
                             Iz=Iz, Sz=Zez, Zz=Zpz, rz=rz,
                             J=J, Cw=Cw,
                             alpha_sy=alpha_sy,
                             alpha_sz=alpha_sz)
    #
    #
    def curved(self, R):
        """
        ---------
        R = Radio
        """
        b = self.b
        b1 = self.tw
        t = self.tb
        d = self.d
    
        # shear area
        warea = self.area
    
        # extreme fibre distances c
        c = (d * (((b1 / b) + (1.0 - (b1 / b))*(t / d)**2) /
                    (2.0*((b1 / b) + (1.0 - (b1 / b))*(t / d)))))
    
        c1 = c * ((d / c) - 1.0)
    
        # centroidal radius
        #_R = R
        #_R = orad - c
    
        # Shift of neutral axis from neutral axis
        e = (c * ((R/c) - (((d/c)*(b1/b + (1.0 - b1/b)*(t/d))) / 
                               (((b1/b) * math.log((d/c + R/c -1.0) / 
                                                     ((d/c)*(t/d) + R/c - 1.0))) +
                                math.log(((d/c)*(t/d) + R/c - 1.0) /
                                         (R/c - 1.0))))))
    
        # where
        Ic = ((warea * c**2) * (((((d/c)**2 *((b1/b + (1.0 - b1/b)*(t/d)**3)
                                                  / (b1/b + (1.0 - b1/b)*(t/d))))) 
                                    / 3.0) - 1.0))
    
        # stress factors Ki
        self.ki = ((Ic / (warea * c**2 * (R/c - 1.0))) 
                   * ((1.0 - e / c) / (e / c)))
    
        # stress factors Ko
        self.ko = ((Ic / (warea * c**2 * (e/c))) 
                   * ((d/c + e/c -1.0) / (R/c + d/c - 1.0)) 
                   * (1.0 / (d/c - 1.0)))
        
        # Modulus of rigidity factor (section 8.10)
        F = 1.0
        
        # Shear factor (section 8.1 equ 8.1-13)
        #    
    #
    #
    def taux_max(self, Mt, alpha: float = 1.12):
        """
        Mt : Torsional moment
        
        Roark Torsion chapter
        """
        ti = self.tb + self.tw
        h = self.d - 0.50 * self.tb
        J = alpha / 3.0 * (self.b * self.tb**3
                           + h * self.tw ** 3)
        return Mt * ti / J
    #    
    #
    # ----------------------------------------
    #
    @property
    def centroid(self):
        """ """
        C = self.b - self.tw
        D2 = self.d - self.tb        
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        Zc = ((self.d**2 * self.tw + C * self.tb**2)
              / (2 * (self.b*self.tb + D2*self.tw)))
        Yc = 0
        #
        return axis(Yc, Zc)
    #
    @property
    def I(self):
        """Seconf modemnt of inertia"""
        Yc, Zc = self.centroid
        C = self.b - self.tw
        D2 = self.d - self.tb
        #
        #   Second Moment of Area about Mayor Axis
        Iy = ((self.tw * (self.d - Zc)**3 + self.b * Zc**3
               - C * (Zc - self.tb)**3) / 3.0)
        #   Second Moment of Area about Minor Axis
        Iz = (self.b**3 * self.tb + D2 * self.tw**3) / 12.0
        #
        return axis(Iy, Iz)
    #
    def Qb(self):
        """
        Returns
        ----------
        Qb: Q/bI
        Where: 
        Q : Ax - Shear's first moment of area
        b : cross section width
        """
        #
        coord = self.section_coordinates()
        Yc, Zc = self.centroid
        #
        # -----------------------------------------------------
        #
        bmin = self.b * 0.50
        #
        Qy = []
        for item in coord.y:
            Hi = abs(item)
            di = bmin - Hi
            area = di * self.tb
            Zi = 0.5 * di + Hi
            Qy.append(area * Zi)
        #
        # -----------------------------------------------------
        #
        
        bf = self.b
        tf = self.tb
        Qz = []
        for item in coord.z:
            Hi = abs(item)
            if item < 0: # Bottom section
                d = self.d - Zc
            else: # Top flange
                d = Zc
            #
            di = d - Hi
            if Hi < (d - tf): # T section
                area = bf * tf + self.tw * di
                C = bf - self.tw
                D = di - tf
                Zi = (d - ((di**2 * self.tw + C * tf**2)
                           / (2 * (bf * tf + D * self.tw))))
                #print('T', item, area, Zi)
            else: # Flange
                area = bf * di
                Zi = 0.5 * di + Hi
                #print('flange', item, area, Zi)
            #
            Qz.append(area * Zi)
        #
        #1 / 0
        return axis(Qy, Qz)
    #
    #
    def alpha_s(self, poisson: float):
        """Shear correction factor"""
        h = self.d - self.tb * 0.50
        j = (self.b * self.tb) / (h * self.tw)
        k = self.b / h
        alpha_sz = (((12 + 96 * j + 276 * j**2 + 192 * j**3)
                    + poisson * (11 + 88 * j + 248 * j**2 + 216 * j**3)
                    + 30 * k**2 * (j + j**2)
                    + 10 * poisson * k**2 * (4 * j + 5 * j**2 + j**3))
                    / (10 * (1 + poisson) * (1 + 4 * j)**2))
        alpha_sy = self.b * self.tb
        return alpha_sy, alpha_sz
    #
    # ----------------------------------------
    #
    def section_coordinates(self):
        """
        1    2     3
        +----------+
        |____+_____|      ^ z
             +      4     |
             |            |
             + 5          +--> y
             |
             +  6
        """
        CoG = self.centroid
        # horizontal
        h1 = self.b * 0.50
        coord_y = [-1 * h1, 0, h1, # 1-3
                   0, 0, 0]        # 4-6
        #
        # vertical
        Zc = CoG.z
        Zcb = Zc - self.d
        coord_z = [Zc, Zc - self.tb * 0.50, Zc, # 1-3
                   Zc - self.tb, # 4
                   0 , Zcb]      # 5, 6
        #
        return points(coord_y, coord_z)
    #
    def _dimension(self) -> str:
        """ Print section dimensions"""
        out = "{:<32s}{:1.4e} {:1.4e}\n"\
               .format(self.shape, self.d, self.b)
        out += "{:<48s}{:1.4e} {:1.4e}\n"\
               .format("", self.tw, self.tb)         
        return out
    #
    #
    @property
    def Dh(self):
        """Hydrodynamic diametre"""
        return self.d   
#
