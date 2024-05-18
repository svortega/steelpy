# 
# Copyright (c) 2009 steelpy
#
#
# Python stdlib imports
from __future__ import annotations
#from array import array
from dataclasses import dataclass
from collections import namedtuple
import math
#

# package imports
from steelpy.sections.utils.shape.utils import ShapeProperty
from steelpy.sections.utils.shape.stress import ShapeStressBasic

#
#-------------------------------------------------
#
points = namedtuple('Points', ['y', 'z'])
axis = namedtuple('Axis', ['y', 'z'])
#
#
#-------------------------------------------------
#      Standard Section Profiles
#-------------------------------------------------
#
@dataclass
class ChannelBasic(ShapeStressBasic):
    """
    Calculate the section properties of a channel section

    +   +-----+
        |         
    D   |         Z
        |         ^
    +   +-----+   + > Y
        *  B  *

    Parameters
    ----------
    D  : Section Heigh
    Tw : Web thickness
    B  : Base
    Tf : Flange thickness

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
    #
    d: float
    tw: float
    b: float
    tb: float
    shape:str = 'Channel'
    #
    # --------------------------------------------
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
    def _properties(self):
        """ """
        # self.a *= factors[0]
        # self.ta *= factors[0]
        _b = self.b - 0.50 * self.tw
        _C = self.b - self.tw
        _h = self.d - self.tb
        _D2 = self.d - 2 * self.tb
        # -------------------------------------------------
        #   Cross-Sectional Area
        area = (_h * self.tw) + (2 * _b * self.tb)
        # -------------------------------------------------
        #   Elastic Neutral Centre
        Yc, Zc = self.centroid
        #Zc = 0.50 * self.d
        #Yc = ((2 * self.b ** 2 * self.tb + _D2 * self.tw ** 2) /
        #      (2 * self.b * self.d - 2 * _D2 * _C))
        # -------------------------------------------------
        #   Shear Centre 
        SCz = 0.50 * self.d
        SCy = ((3 * self.tb * self.b ** 2)
               / (6 * _b * self.tb + _h * self.tw))
        # -------------------------------------------------
        #   Warping Constant Cw
        Cw = ((_b ** 3 * _h ** 2 * self.tb / 12.0) *
              ((2 * _h * self.tw + 3 * _b * self.tb) /
               (_h * self.tw + 6 * _b * self.tb)))
        # -------------------------------------------------
        #               Section Properties
        # -------------------------------------------------
        Iy, Iz = self.I
        #   Second Moment of Area about Mayor Axis
        #Iy = (self.b * self.d ** 3 - _C * _D2 ** 3) / 12.0
        #   Elastic Modulus about Mayor Axis
        Zey = 2 * Iy / self.d
        #   Plastic Modulus about Mayor Axis
        Zpy = (((self.d ** 2 * self.tw) / 4.0) +
               (self.tb * _C * (self.d - self.tb)))
        #   Shape Factor
        SFy = Zpy * self.d / (2 * Iy)
        #   Radius of gyration about Mayor Axis
        ry = math.sqrt(Iy / area)
        # -------------------------------------------------
        #   Second Moment of Area about Minor Axis
        #Iz = ((self.d * self.b ** 3 / 3.0) - (_C ** 3 * _D2 / 3.0)
        #      - (area * (self.b - Yc) ** 2))
        #   Elastic Modulus about Minor Axis
        Zez = Iz / (self.b - Yc)
        #   Plastic Modulus about Minor Axis
        if (2 * self.tb * _C) > (self.d * self.tw):
            Zpz = ((_C ** 2 * self.tb / 2.0)
                   - (self.d ** 2 * self.tw ** 2 / (8 * self.tb)) +
                   (self.d * self.tw * self.b / 2.0))
        else:
            Zpz = (((self.tw ** 2 * self.d) / 4.0) +
                   (self.tb * _C * (self.b - (self.tb * _C / self.d))))
        #   Shape Factor
        SFz = (Zpz * (self.b - Yc) / Iz)
        #   Radius of gyration about Minor Axis 
        rz = math.sqrt(Iz / area)
        # -------------------------------------------------
        #   Torsional Constant
        J = (2 * _b * self.tb ** 3 + _h * self.tw ** 3) / 3.0
        #   Product of inertia
        _Iyz = 0.0
        Jx = Iy + Iz
        rp = math.sqrt(Jx / area)
        #
        #
        return ShapeProperty(area=area, Zc=Zc, Yc=Yc,
                             Iy=Iy, Sy=Zey, Zy=Zpy, ry=ry,
                             Iz=Iz, Sz=Zez, Zz=Zpz, rz=rz,
                             J=J, Cw=Cw)

    #
    def curved(self, R):
        """
        ---------
        R = Radio
        """
        _b = self.b
        _b1 = 2 * self.tw
        _t = self.tb
        _d = self.d

        # shear area
        _warea = self.area

        # extreme fibre distances c
        _c = (_d * (((_b1 / _b) + (1.0 - (_b1 / _b)) * (_t / _d) ** 2) /
                    (2.0 * ((_b1 / _b) + (1.0 - (_b1 / _b)) * (_t / _d)))))

        _c1 = _c * ((_d / _c) - 1.0)

        # centroidal radius
        _R = R
        # _R = orad - _c

        # Shift of neutral axis from neutral axis
        _e = (_c * ((_R / _c) - (((_d / _c) * (_b1 / _b + (1.0 - _b1 / _b) * (_t / _d))) /
                                 (((_b1 / _b) * math.log((_d / _c + _R / _c - 1.0) /
                                                         ((_d / _c) * (_t / _d) + _R / _c - 1.0))) +
                                  math.log(((_d / _c) * (_t / _d) + _R / _c - 1.0) /
                                           (_R / _c - 1.0))))))

        # where
        _Ic = ((_warea * _c ** 2) * (((((_d / _c) ** 2 * ((_b1 / _b + (1.0 - _b1 / _b) * (_t / _d) ** 3)
                                                          / (_b1 / _b + (1.0 - _b1 / _b) * (_t / _d)))))
                                      / 3.0) - 1.0))

        # stress factors Ki
        self.ki = ((_Ic / (_warea * _c ** 2 * (_R / _c - 1.0)))
                   * ((1.0 - _e / _c) / (_e / _c)))

        # stress factors Ko
        self.ko = ((_Ic / (_warea * _c ** 2 * (_e / _c)))
                   * ((_d / _c + _e / _c - 1.0) / (_R / _c + _d / _c - 1.0))
                   * (1.0 / (_d / _c - 1.0)))

        # Modulus of rigidity factor (section 8.10)
        _F = 1.0

        # Shear factor (section 8.1 equ 8.1-13)
        #

    #    
    # --------------------------------------------
    #
    def tau_t(self, phi, G: float, alpha: float = 1.31):
        """
        Torsional stress due to pure torsion
        
        psi : theta' (Rate of angle of twist theta with respect to lengh)
        G   : Shear modulus of elasticity 
        
        sigma_t : 
        """
        #ti = self.tw + self.tb
        h = self.d - 0.50 * self.tb
        J = alpha / 3.0 * (self.b * self.tb**3
                           + h * self.tw**3)
        return phi * G * J
    #  
    def warping_stress(self, psi, Tw, B, G: float):
        """
        psi : theta' (Rate of angle of twist theta with respect to lengh)
        Tw  : Twisting moment
        B   : Bimoment warping moment
        G   : Shear modulus of elasticity 

        sigma_t : 
        sigma_w : Normal stress due to warping
        tau_w   : Shear stress due to warping 
        """
        1 / 0
    #
    def warping_properties(self):
        """
        Returns:
        Cw, J, K, w, Qw
        """
        1 / 0
    #
    # --------------------------------------------
    #
    @property
    def I(self):
        """Second moment of inertia"""
        b = self.b - 0.50 * self.tw
        h = self.d - self.tb
        area = h * self.tw + 2 * b * self.tb
        #C = self.b - self.tw
        #D = h - self.tb
        #
        Yc, Zc = self.centroid
        #
        #   Second Moment of Area about Mayor Axis
        #Iy = (self.b * self.d**3 - C * D**3) / 12.0
        #
        Iy = (self.b  * self.d**3
              - (self.b - self.tw)
              * (self.d - 2 * self.tb)**3) / 12.0
        #
        #   Second Moment of Area about Minor Axis
        #Iz = ((2 * self.tb * self.b**3 + D * self.tw**3) / 3.0
        #      - area * (self.b - Yc)**2)
        #
        Iz = (self.d / 3.0 * self.b**3
              - (self.b - self.tw)**3 / 3.0 * (self.d - 2 * self.tb)
              - area * (self.b - Yc)**2)
        #
        return axis(Iy, Iz)
    #
    def Qb(self):
        """
        Returns
        ----------
        Qb: Q/b
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
        #h = self.d - self.tb
        #D = self.d - 2 * self.tb
        #
        Qy = []
        for item in coord.y:
            Hi = abs(item)
            if item < 0: # left
                b = self.b - Yc
                #bmin = b
                #di = b - Hi
            else: # right 
                b = Yc
                #bmin = b - self.tw
                #di = Hi
            #
            #if Hi > bmin: # C area
            #    area = di * self.tw + 2 * bmin * self.tb
            #    Zi = (b - ((2 * Hi**2 * self.tb + D * self.tw**2)
            #               / (2 * Hi * self.d - 2 * D * bmin)))
            #    print('C', item, area, Zi)
            #else: # two flanges
            di = b - Hi
            area = 2 * (di * self.tb)
            Zi = 0.5 * di
            #print('flange', item, area, Zi)
            #
            Qy.append(area * Zi)
        #
        # -----------------------------------------------------
        #
        Qz = []
        for item in coord.z:
            Hi = abs(item)
            #if item < 0: # bottom
            #    d = self.d - Zc
            #    #hmin = d
            #    #di = d - Hi
            #else: # top 
            #    d = Zc
            #    #hmin = d - self.tb
            #    #di = Hi            
            #
            d = Zc
            #hmin = self.tb
            di = d - Hi
            if Hi < (d - self.tb): # L area
                Zi = (d - ((self.b**2 + di * self.tw)
                           / (2 * (self.b + di))))
                #
                b1 = self.b - self.tw
                area = di * self.tw + b1 * self.tb
                #print('L', item, di, area, Zi)
            else: # one area
                #di = hmin - Hi
                Zi = 0.5 * di + Hi
                area =  di * self.tb
                #print('web', item, di, area, Zi)            
            #
            Qz.append(area * Zi)
        #
        #1 / 0
        return axis(Qy, Qz)
    #
    # --------------------------------------------
    #
    @property
    def centroid(self):
        """ """
        #h = self.d - self.tb
        C = self.b - self.tw
        D = self.d - 2 * self.tb        
        # -------------------------------------------------
        #   Elastic Neutral Centre 
        Zc = 0.50 * self.d
        Yc = ((2 * self.b ** 2 * self.tb + D * self.tw ** 2)
              / (2 * self.b * self.d - 2 * D * C))
        #
        return axis(Yc, Zc)
    #    
    @property
    def shear_center(self):
        """ """
        #   Shear Centre
        b = self.b - 0.50 * self.tw
        h = self.d - self.tb        
        Zc = 0 # 0.50 * self.d
        Yc = ((3 * self.tb * self.b ** 2)
               / (6 * b * self.tb + h * self.tw))
        return axis(Yc, Zc)
    #
    def section_coordinates(self):
        """
        1   2    3 4
        +---+----+-+
        |           
        +5  ^ z
        |   |      
        +6  +-->    
        |       y
        + 7
        | 
        +---+----+-+
        8   9   10 11
        """
        #
        CoG = self.centroid
        #CoG = self.shear_center        
        # horizontal
        h1 = CoG.y - self.b + self.tw * 0.50
        h2 = CoG.y - self.tw
        h3 = CoG.y # - self.tw * 0.50
        #
        coord_y = [h1, 0.0, h2, h3, # 1,2,3,4
                   h1, h1, h1,     # 5,6,7
                   h1, 0.0, h2, h3] # 8,9,10,11
        #
        # vertical
        #d = self.d * 0.50
        top1 = CoG.z - self.tb * 0.50
        top2 = CoG.z - self.tb
        top3 = -1 * top1
        coord_z = [top1, top1, top1, top1, # 1-4
                   top2, 0.0,  -1 * top2,    # 5,6,7
                   top3, top3, top3, top3] # 8-11
        
        return points(coord_y, coord_z)
    #
    #
    def _dimension(self) -> str:
        """ Print section dimensions"""
        out = "{:<32s}{:1.4e} {:1.4e} {:1.4e}\n" \
            .format(self.shape, self.d, self.b, self.b)
        out += "{:<48s}{:1.4e} {:1.4e} {:1.4e}\n" \
            .format("", self.tw, self.tb, self.tb)
        return out
    #
    #
    @property
    def Dh(self):
        """Hydrodynamic diametre"""
        return math.hypot(self.d, self.b)    
#    
#
