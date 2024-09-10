# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from array import array
from collections import namedtuple
from dataclasses import dataclass
import math
#

# package imports
from steelpy.sections.utils.shape.utils import ShapeProperty
from steelpy.sections.utils.shape.stress import BeamStress, ShapeStressBasic
#from steelpy.sections.process.shape import SectionBasic, ShapeBasic

#
#-------------------------------------------------
#
points = namedtuple('Points', ['y', 'z'])
axis = namedtuple('Axis', ['y', 'z'])
#
#-------------------------------------------------
#      Standard Section Profiles
#-------------------------------------------------
#
@dataclass(slots=True)
class AngleBasic(ShapeStressBasic):
    """
    Calculate the section properties of an angle section

         +   +
             |         
        d    |         Z
             |         ^
         +   +-----+   + > Y
             *  b  *

    Parameters
    ----------
    d  : Section Heigh
    tw : Web/flange thickness
    b  : Base
    r  : Root Radious

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
    1.- Full plastic capacity of equal angle sections under biaxial
        bending and normal force [A.E. Charalampakis]
    2.- Formulas for stress, strain and strucutral matrices [W.D. Pilkey]
    3.- Roark's formulas for stress and strain [7th Edition]
    4.- Wikipedia

    Examples
    ----------

    """
    #name:str|int
    d:float
    tw:float
    b:float
    r:float
    #shape:str = 'Angle'
    #
    # --------------------------------------------
    @property
    def shape(self):
        return 'Angle'
    #    
    def _properties(self, poisson: float):
        # -------------------------------------------------
        #
        Yc, Zc = self.centroid
        Iy, Iz = self.I
        #
        # Simmetric angle
        if math.isclose(self.d, self.b):
            #
            # -------------------------------------------------
            #   Cross-Sectional Area (A.1)
            area = ((2 * self.d * self.tw)
                            - self.tw ** 2 + (2 - math.pi / 2.0) * self.r ** 2)
        
            # Distance v1 of elastic centroid from heel (A.2)
            _v1 = (((6 * math.sqrt(2)) * (((50 / 3.0 - 5 * math.pi) * self.r ** 3)
                                          - ((2 - math.pi / 2.0) * (self.d - 3 * self.tw) * self.r ** 2)
                                          + ((self.d ** 2 + self.d * self.tw - self.tw * 2) * self.tw)))
                   / (((24 - 6 * math.pi) * self.r ** 2) + (24 * self.d * self.tw) - (12 * self.tw ** 2)))
        
            #  Distance v2 (A.3)
            _v2 = ((self.d / math.sqrt(2)) + (1 - math.sqrt(2)) +
                           (self.tw / math.sqrt(2)) - _v1)
        
            # Distance v3 (A.4)
            _v3 = self.d / math.sqrt(2)
        
            # Moment of inertia around u-u (A.5)
            _Iu = (((70 - 21 * math.pi) * self.r ** 4 / 24.0) +
                   ((self.d - self.tw) ** 2 * (math.pi - 4) * self.r ** 2 / 4.0)
                   - (self.tw ** 4 / 12.0) + (self.d ** 3 * self.tw / 3.0) + (self.d * self.tw ** 3 / 3.0)
                   - (self.d ** 2 * self.tw ** 2 / 2.0))
        
            # Elastic section modulus around u-u (A.6)
            _Wu = _Iu / _v3
        
            # Moment of inertia around v-v (A.7)
            _Iv = ((1.0 / ((288 - 72 * math.pi) * self.r ** 2 + 288 * (self.d - (self.tw / 2.0)) * self.tw))
                   * ((((7926 * math.pi) - (1233 * math.pi ** 2) - 12776) * self.r ** 6)
                      + (432 * (math.pi - (10.0 / 3.0)) * (self.d - self.tw) * (math.pi - 4) * self.r ** 5)
                      + (((1422 * math.pi) - (36 * math.pi ** 2) - 4188) * self.tw ** 2)
                      + (((72 * ((-79 * math.pi / 2.0) + math.pi ** 2 + (349 / 3.0)) * self.d * self.tw)
                          - (36 * self.d ** 2 * (math.pi - 4) ** 2) * self.r ** 4))
                      + (432 * (math.pi - (10.0 / 3.0)) * (4 * self.tw ** 2 / 3.0 + self.d ** 2 -
                                                           4 * self.d * self.tw) * self.tw * self.r ** 3)
                      - (24 * (self.d ** 3 - (9 * self.tw * self.d ** 2) + (13 * self.tw ** 2 * self.d)
                               - (13 * self.tw ** 3 / 4.0))
                         * ((math.pi - 4) * self.tw * self.r ** 2))
                      + (-(72 * self.d * self.tw ** 5) + (12 * self.tw ** 6) + (24 * self.tw ** 2 * self.d ** 4)
                         - (48 * self.tw ** 3 * self.d ** 3) + (96 * self.tw ** 4 * self.d ** 2))))
        
            # Elastic section modulus around v-v (A.8)
            _Wv = _Iv / _v1
        
            # Distance e (A.9)
            _e = _v1 / math.sqrt(2)
            #Zc = _e
            #Yc = _e
        
            # Moment of inertia around y-y or z-z (A.10)
            #Iy = ((1.0 / (((288 - 72 * math.pi) * self.r ** 2)
            #              + (288 * (self.d - self.tw / 2.0) * self.tw)))
            #      * (((3732 * math.pi - 5968 - 585 * math.pi ** 2) * self.r ** 6)
            #         + ((216 * (math.pi - 4) * (math.pi - 10.0 / 3.0) * (self.d - self.tw) * self.r ** 5)
            #            + (60 * self.d ** 4 * self.tw ** 2))
            #         + (-(120 * self.d ** 3 * self.tw ** 3) + (132 * self.d ** 2 * self.tw ** 4)
            #            - (72 * self.d * self.tw ** 5) +
            #            (12 * self.tw ** 6) + (((216 * (math.pi - 10.0 / 3.0))) *
            #                                   ((self.d ** 2 - 4 * self.d * self.tw + 4 * self.tw ** 2 / 3.0)
            #                                    * self.r ** 3 * self.tw)))
            #         + (((-27 * self.d ** 2 * (math.pi - 4) ** 2)
            #             + (54 * (272 / 3.0 - 94 * math.pi / 3.0 + math.pi ** 2) * self.d * self.tw)
            #             + ((846 * math.pi - 2448 - 27 * math.pi ** 2) * self.tw ** 2)) * self.r ** 4)
            #         + (12 * (math.pi - 4) * (self.d ** 3 + 3 * self.d ** 2 * self.tw
            #                                  - 8 * self.d * self.tw ** 2 + 2 * self.tw ** 3) * self.r ** 2 * self.tw)))
            #
            #Iz = Iy
        
            # Elastic section modulus araund y-y or z-z (A.11)
            Zey = Iy / (self.d - _e)
            Zez = Zey
        
            # Radii of gyration around any axis (A.12)
            ry = (Iy / area) ** 0.50
            rz = ry
        
            # Plastic Properties
            # Distance v1p of the plastic centroid from heel (A.13)
            _v1p = ((3 * (math.pi - 4) * self.r ** 2
                             + 4 * self.d * self.tw + 6 * self.tw ** 2)
                            / (8 * math.sqrt(2) * self.tw))
        
            # Major plastic section modulus around u'-u' (A.14)
            Zpy = ((((48 * self.r ** 3
                      + 3 * (math.pi - 4) * (self.d - self.tw) * self.r ** 2
                      - 6 * self.d * self.tw ** 2
                      + 6 * self.d ** 2 * self.tw
                      + 2 * self.tw ** 3) * math.sqrt(2))
                    / 12.0) - (16 * self.r ** 3 / 3.0))
        
            # Minor Plastic section modulus around v'-v' (A.15)
            Zpz = ((math.sqrt(2) / (192 * self.tw))
                   * ((-27 * (math.pi - 4) ** 2 * self.r ** 4)
                      + (96 * (3 * math.pi - 10) * self.r ** 3 * self.tw)
                      - ((12 * (math.pi - 4) * self.r ** 2) * ((2 * self.d - 11 * self.tw) * self.tw))
                      + (4 * self.tw ** 2 * (12 * self.d ** 2
                                             - 12 * self.d * self.tw + 13 * self.tw ** 2))))
        
            _phi = 1.0
        
            # -------------------------------------------------
            #   Plastic Modulus about Minor Axis
            #   error, needs fix
        
            # _Zpy = 0
            # _Zpz = 0
            #            
            ####
            ta = self.tw / self.d
            if float(ta) > 0.40:
                yp = 0.3536 * (self.d + 1.5 * self.tw)
                Zpy = area * (Yc - 0.6667 * yp)
            else:
                yp = self.d * math.sqrt(ta - ta**2 / 2.0)
                Zpy = (area * Yc - 2.8284 * yp**2 * self.tw
                       + 1.8856 * self.tw**3)
            Zpz = Zpy
        
        else:
            #
            _b1 = self.b - self.tw
            _h1 = self.d - self.tw
        
            # -------------------------------------------------
            #   Cross-Sectional Area
            area = (self.d + _b1) * self.tw
        
            # -------------------------------------------------
            #   Elastic Neutral Centre
            #Zc = (self.d ** 2 + _b1 * self.tw) / (2 * (self.d + _b1))
            #Yc = (self.b ** 2 + _h1 * self.tw) / (2 * (self.b + _h1))
        
            # -------------------------------------------------
            #   Shear Centre 
            _SCz = self.tw / 2.0
            _SCy = self.tw / 2.0
        
            # -------------------------------------------------
            #               Section Properties
            # -------------------------------------------------
            #   Second Moment of Area about Mayor Axis
            #Iy = ((self.tw * (self.d - Zc) ** 3 + self.b * Zc ** 3 -
            #               _b1 * (Zc - self.tw) ** 3) / 3.0)
        
            Zey = Iy / (self.d - Zc)
            #
            ta = self.tw / self.d
            if float(ta) > 0.40:
                yp = 0.3536 * (self.d + 1.5 * self.tw)
                Zpy = area * (Yc - 0.6667 * yp)
            else:
                yp = self.d * math.sqrt(ta - ta**2 / 2.0)
                Zpy = (area * Yc - 2.8284 * yp**2 * self.tw
                       + 1.8856 * self.tw**3)
        
            #   Radius of gyration about Mayor Axis
            ry = (Iy / area) ** 0.50
            _SFy = 0
            # -------------------------------------------------
            #   Second Moment of Area about Minor Axis
            #Iz = ((self.tw * (self.b - Yc) ** 3 + self.d * Yc ** 3 -
            #               _h1 * (Yc - self.tw) ** 3) / 3.0)
        
            Zez = Iz / (self.b - Yc)
            Zpz = Zpy
        
            #   Radius of gyration about Minor Axis 
            rz = (Iz / area) ** 0.50
            _SFz = 0
        
            # -------------------------------------------------
            #   Product of inertia
            _Izy = (self.b * _b1 * self.d * _h1 * self.tw) / (4 * (self.b + _h1))
        
            _Iu = (0.50 * (Iz + Iy) +
                           0.50 * ((Iz - Iy) ** 2 + 4 * _Izy ** 2) ** 0.50)
        
            _Iv = (0.50 * (Iz + Iy) -
                           0.50 * ((Iz - Iy) ** 2 + 4 * _Izy ** 2) ** 0.50)
        
            #_phi = (math.atan(2 * _Izy / (Iy - Iz))) / 2.0
            #_Wu = 0
            #_Wv = 0
            #_Wup = 0
            #_Wvp = 0
            #_v1p = 0            
            # Minor Plastic section modulus around v'-v' (A.15)
            Zpy = ((((48 * self.r ** 3
                      + 3 * (math.pi - 4) * (self.d - self.tw) * self.r ** 2
                      - 6 * self.d * self.tw ** 2
                      + 6 * self.d ** 2 * self.tw
                      + 2 * self.tw ** 3) * math.sqrt(2))
                    / 12.0) - (16 * self.r ** 3 / 3.0))
            #
            # Major plastic section modulus around u'-u' (A.14)
            Zpzz = ((math.sqrt(2) / (192 * self.tw))
                   * ((-27 * (math.pi - 4) ** 2 * self.r ** 4)
                      + (96 * (3 * math.pi - 10) * self.r ** 3 * self.tw)
                      - ((12 * (math.pi - 4) * self.r ** 2) * ((2 * self.d - 11 * self.tw) * self.tw))
                      + (4 * self.tw ** 2 * (12 * self.d ** 2
                                             - 12 * self.d * self.tw + 13 * self.tw ** 2))))           
        #
        b1 = self.d - 0.50 * self.tw
        b2 = self.b - 0.50 * self.tw
        #_c1 = _b1 - 0.50 * self.tw
        #_c2 = _b2 - 0.50 * self.tw
        #
        # -------------------------------------------------
        #   Warping Constant Cw (Bleich 1952, Picard and Beaulie 1991)
        Cw = (((b1 ** 3 + b2 ** 3) * (self.tw ** 3)) / 36.0)

        # -------------------------------------------------
        #   Torsional Constant (fillets neglected)
        J = (((b1 + b2) * (self.tw ** 3)) / 3.0)

        # -------------------------------------------------
        #   Product of inertia
        Iyz = (self.tw * (b1 * Yc * (b1 - 2 * Zc)) +
               self.tw * (b2 * Zc * (b2 - 2 * Yc)))

        Jx = Iy + Iz
        rp = (Jx / area) ** 0.50
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
    def alpha_s(self, poisson: float):
        """Shear correction factor"""
        alpha_s = (4 + 3 * poisson) / (2 * (1 + poisson))
        print('fix angle shear correction factor')
        return alpha_s, alpha_s    
    #
    def _stressX(self, actions, stress=None):
        """
        """
        #
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
    # -----------------------------------------------------
    #
    #
    def taux_max(self, Mt, alpha: float = 1.0):
        """
        Mt : Torsional moment
        
        Roark Torsion chapter
        """
        ti = 2 * self.tw 
        h = self.d - 0.50 * self.tw
        J = alpha / 3.0 * (self.b * self.tw**3
                           + h * self.tw**3)
        return Mt * ti / J    
    #
    # -----------------------------------------------------
    #
    @property
    def centroid(self):
        """ """
        # Simmetric angle
        #
        Yc = ((self.b**2 + self.d * self.tw - self.tw**2)
               / (2 * (self.b + self.d - self.tw)))
        
        Zc = ((self.d**2 + self.b * self.tw - self.tw**2)
               / (2 * (self.b + self.d - self.tw)))
        #
        return axis(Yc, Zc)
    #
    @property
    def shear_center(self):
        """ """
        #   Shear Centre
        1 / 0
    #
    @property
    def I(self):
        """Second moment of inertia"""
        # Simmetric angle
        if math.isclose(self.d, self.b):
            # Moment of inertia around y-y or z-z (A.10)
            Iy = ((1.0 / (((288 - 72 * math.pi) * self.r**2)
                          + (288 * (self.d - self.tw / 2.0) * self.tw)))
                  * (((3732 * math.pi - 5968 - 585 * math.pi**2) * self.r**6)
                     + ((216 * (math.pi - 4) * (math.pi - 10.0 / 3.0) * (self.d - self.tw) * self.r**5)
                        + (60 * self.d**4 * self.tw**2))
                     + (-(120 * self.d**3 * self.tw**3) + (132 * self.d**2 * self.tw**4)
                        - (72 * self.d * self.tw**5) +
                        (12 * self.tw**6) + (((216 * (math.pi - 10.0 / 3.0))) *
                                             ((self.d**2 - 4 * self.d * self.tw + 4 * self.tw**2 / 3.0)
                                              * self.r**3 * self.tw)))
                     + (((-27 * self.d**2 * (math.pi - 4)**2)
                         + (54 * (272 / 3.0 - 94 * math.pi / 3.0 + math.pi**2) * self.d * self.tw)
                         + ((846 * math.pi - 2448 - 27 * math.pi**2) * self.tw**2)) * self.r**4)
                     + (12 * (math.pi - 4) * (self.d**3 + 3 * self.d**2 * self.tw
                                              - 8 * self.d * self.tw**2 + 2 * self.tw**3) * self.r**2 * self.tw)))
            Iz = Iy
        else:
            Yc, Zc = self.centroid
            b1 = self.b - self.tw
            h1 = self.d - self.tw
            # -------------------------------------------------
            #   Second Moment of Area about Mayor Axis
            Iy = ((self.tw * (self.d - Zc)**3 + self.b * Zc**3
                   - b1 * (Zc - self.tw)**3) / 3.0)
            #
            #   Second Moment of Area about Minor Axis
            Iz = ((self.tw * (self.b - Yc)**3 + self.d * Yc**3
                   - h1 * (Yc - self.tw)**3) / 3.0)            
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
        Qy = []
        for item in coord.y:
            Hi = abs(item)
            if item < 0: # left
                b = Yc
            else: # right 
                b = self.b - Yc
            #
            di = b - Hi
            area = di * self.tw
            Zi = 0.5 * di + Hi
            #print('flange', item, area, Zi)
            #
            Qy.append(area * Zi)
        #
        # -----------------------------------------------------
        #
        Qz = []
        for item in coord.z:
            Hi = abs(item)
            if item < 0: # bottom
                d = Zc
            else: # top
                d = self.d - Zc
            #
            di = d - Hi
            area = di * self.tw
            Zi = 0.5 * di + Hi
            #print('web', item, area, Zi)
            #
            Qz.append(area * Zi)            
        #
        #1 / 0
        return axis(Qy, Qz)        
    #
    # -----------------------------------------------------
    #
    def section_coordinates(self):
        """
        1 
        +
        |    ^ z    
        + 2  |      
        |    |
        + 3  +--> y  
        |       
        + 4
        | 
        +-+--+---+-+
        5 6  7   8 9
        """
        # horizontal
        
        CoG = self.centroid
        #
        h1 = -1 * CoG.y + self.tw * 0.50
        h2 = -1 * CoG.y - self.tw
        h3 = CoG.y - self.tw
        coord_y = [h1, h1, h1, h1, # 1,2,3,4
                   h1, h2 , 0.0,   # 5,6,7
                   h3, -1 * CoG.y]    # 8,9
        #
        # vertical
        D = self.d - CoG.z
        top1 = D
        top2 = D - self.tw
        top3 = -1 * ( CoG.z - self.tw)
        top4 = -1 * (CoG.z - self.tw * 0.50)
        coord_z = [top1, top2, # 1,2
                   0,          # 3
                   top3,       # 4
                   top4, top4, top4, top4, top4]  # 5-9
        #
        return points(coord_y, coord_z)      
    #
    # -----------------------------------------------------
    #
    def print_file(self, file_name):

        check_out = print_header()

        check_out.append("{:23s} {:>19} {:1.4E} {:1.4E} {:1.4E} {:1.4E}\n"
                         .format(self.shape, "", self.d, self.tw, self.b, self.tw))

        check_out.extend(print_properties(self))

        # file_checkout = split_file_name(file_name)
        # file_checkout = str(file_checkout[0]) +'_check_me.txt'
        file_checkout = str(file_name) + '.txt'
        add_out = open(file_checkout, 'w')
        add_out.write("".join(check_out))
        add_out.close()
        print('ok')
    #
    def _dimension(self) -> str:
        """ Print section dimensions"""
        out = "{:<32s}{:1.4e} {:1.4e}\n"\
               .format(self.shape, self.d, self.b)
        out += "{:<48s}{:1.4e} {:1.4e}\n"\
               .format("", self.tw, self.tw)         
        #return ("{:23s} {:>19} {:1.4E} {:1.4E} {:1.4E} {:1.4E}\n"
        #        .format(self.shape, "", self.d, self.tw, self.b, self.tw))
        return out

    #
    #
    @property
    def Dh(self):
        """Hydrodynamic diametre"""
        return math.hypot(self.d, self.b)