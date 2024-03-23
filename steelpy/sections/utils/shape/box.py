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
#
@dataclass
class BoxBasic(ShapeStressBasic):
    """
    Calculate the section properties of a box section

    +   +------+
        | +--+ |
        | |  | |    
    d   | |  | |  
        | |  | |   Z
        | +--+ |   ^
    +   +------+   + > Y
        *  b   *

    Parameters
    ----------
    d  : Section Heigh
    tw : Web thickness
    b  : Base
    tb : Flange thickness

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
    #name:str | int
    d: float
    tw: float
    b: float
    tb: float
    type:str = 'box'
    #
    # --------------------------------------------
    #      
    #
    #
    def _properties(self, poisson: float):
        """ """
        # self.units_in = _units_output

        # self.d *= factors[0]
        # self.tw *= factors[0]

        # self.a *= factors[0]
        # self.ta *= factors[0]

        # self.b *= factors[0]
        # self.tb *= factors[0]
        # -------------------------------------------------
        #
        _hi = self.d - 2 * self.tb
        _bi = self.b - 2 * self.tw
        # -------------------------------------------------
        #   Cross-Sectional Area
        area = self.b * self.d - _bi * _hi
        # -------------------------------------------------
        #   Elastic Neutral Centre 
        Zc = 0
        Yc = 0
        # -------------------------------------------------
        #   Shear Centre 
        SCz = 0
        SCy = 0
        # -------------------------------------------------
        #   Warping Constant Cw
        Cw = 0
        # -------------------------------------------------
        #               Section Properties
        # -------------------------------------------------
        Iy, Iz = self.I
        #   Second Moment of Area about Mayor Axis
        #Iy = ((self.b * self.d ** 3 - _bi * _hi ** 3) / 12.0)
        #   Elastic Modulus about Mayor Axis
        Zey = ((self.b * self.d ** 3 - _bi * _hi ** 3) / (6.0 * self.d))
        #   Plastic Modulus about Mayor Axis
        Zpy = ((self.b * self.d ** 2 - _bi * _hi ** 2) / 4.0)
        #   Shape Factor
        SFy = (1.50 * (self.d * (self.b * self.d ** 2 - _bi * _hi ** 2)) /
               (self.b * self.d ** 3 - _bi * _hi ** 3))
        #   Radius of gyration about Mayor Axis
        ry = math.sqrt(Iy / area)
        # -------------------------------------------------
        #   Second Moment of Area about Minor Axis
        #Iz = ((self.d * self.b ** 3 - _hi * _bi ** 3) / 12.0)
        #   Elastic Modulus about Minor Axis
        Zez = ((self.d * self.b ** 3 - _hi * _bi ** 3) / (6.0 * self.b))
        #   Plastic Modulus about Minor Axis
        Zpz = ((self.d * self.b ** 2 - _hi * _bi ** 2) / 4.0)
        #   Shape Factor
        SFz = (1.50 * (self.b * (self.d * self.b ** 2 - _hi * _bi ** 2)) /
               (self.d * self.b ** 3 - _hi * _bi ** 3))
        #   Radius of gyration about Minor Axis 
        rz = math.sqrt(Iz / area)
        # -------------------------------------------------
        #   Torsional Constant
        # Mean corner radious
        _Rc = 1.50 * (self.tw + self.tb) / 2.0
        # Mid countour length
        _p = (2 * ((self.d - self.tb) + (self.b - self.tw))
              - 2 * _Rc ** 2 * (4 - math.pi))
        # Enclosed Area
        _Ap = ((self.d - self.tb) * (self.b - self.tw)
               - _Rc ** 2 * (4 - math.pi))
        # for thin walled sections b/t >= 10
        J = ((4 * _Ap ** 2 * ((self.tw + self.tb) / 2.0)) / _p)
        # -------------------------------------------------
        #   Product of inertia
        _Iyz = 0.0
        Jx = Iy + Iz
        rp = math.sqrt(Jx / area)
        #
        #
        alpha_sy = self.alpha_s(poisson=poisson)
        #
        return ShapeProperty(area=area, Zc=Zc, Yc=Yc,
                             Iy=Iy, Sy=Zey, Zy=Zpy, ry=ry,
                             Iz=Iz, Sz=Zez, Zz=Zpz, rz=rz,
                             J=J, Cw=Cw,
                             alpha_sy=alpha_sy,
                             alpha_sz=alpha_sy)
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
        _c = _D / 2.0

        _c1 = _d - _c

        # centroidal radius
        _R = R
        # _R = orad - _c1

        # Shift of neutral axis from neutral axis
        _e = (_c * ((_R / _c) - ((2.0 * (_t / _c + (1 - _t / _c) * (_b1 / _b)))
                                 / ((math.log(((_R / _c) ** 2 + (_R / _c + 1) * (_t / _c) - 1.0)
                                              / ((_R / _c) ** 2 - (_R / _c - 1.0) * (_t / _c) - 1.0)))
                                    + ((_b1 / _b) * math.log((_R / _c - _t / _c + 1.0)
                                                             / (_R / _c + _t / _c - 1.0)))))))

        # where
        _Ic = self.Iy

        # stress factors Ki
        self.ki = ((_Ic / (_warea * _c ** 2 * (_R / _c - 1.0)))
                   * ((1.0 - _e / _c) / (_e / _c)))

        # stress factors Ko
        self.ko = ((_Ic / (_warea * _c ** 2 * (_R / _c + 1.0)))
                   * ((1.0 + _e / _c) / (_e / _c)))

        # Modulus of rigidity factor (section 8.10)
        _nai = _c - _e  # neautral axis inner fiber
        _nao = _c1 + _e  # neautral axis outer fiber

        _D1 = _nai - _Tfb
        _D2 = _nai
        _t1 = 2 * _Tw
        _t2 = _Bfb
        _r = _rip

        self.F = ((1 + (((3 * (_D2 ** 2 - _D1 ** 2) * _D1) / (2.0 * _D2 ** 3)) * (_t2 / _t1 - 1.0)))
                  * (4 * _D2 ** 2 / (10 * _r ** 2)))
        #
        # Shear factor (section 8.1 equ 8.1-13)
    #
    def taux_max(self, Mt):
        """
        """
        #K = ((2 * self.tw * self.tb
        #      * (self.d * self.tw)**2 * (self.b * self.tb)**2)
        #     / (self.d * self.tw + self.b * self.tb - self.tw**2 - self.tb**2))
        #
        #tau_short = Mt / (2 * self.tb * (self.d - self.tb) * (self.b - self.tw))
        tau_long = Mt / (2 * self.tw * (self.d - self.tb) * (self.b - self.tw))
        #
        #tau_max = max(tau_short, tau_long)
        return tau_long    
    #
    #
    def alpha_s(self, poisson: float):
        """Shear correction factor"""
        j = self.b * self.tb / (self.d * self.tw)
        k = self.b / self.d
        alpha_sy = (((12 + 72 * j + 150 * j**2 + 90 * j**3)
                    + poisson * (11 + 66 * j + 135 * j**2 + 90 * j**3)
                    + 10 * k**2 * ((3 + poisson) * j + 3 * j**2))
                    / (10 * (1 + poisson) * (1 + 3 * j)**2))
        return alpha_sy    
    #
    # --------------------------------------------
    #
    @property
    def I(self):
        """Second moment of inertia"""
        hi = self.d - 2 * self.tb
        bi = self.b - 2 * self.tw        
        #   Second Moment of Area about Mayor Axis
        Iy = ((self.b * self.d ** 3 - bi * hi ** 3) / 12.0)
        #   Second Moment of Area about Minor Axis
        Iz = ((self.d * self.b ** 3 - hi * bi ** 3) / 12.0)
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
        #
        # -----------------------------------------------------
        #
        def _qi(coord: list, d: float, td: float,
                b: float, tb: float):
            """qi = Ax/Ib"""
            H = d - tb
            D = b - 2 * td
            Qout = []
            for item in coord:
                Hi = abs(item)
                di = d - Hi
                if Hi < H: # C section
                    area = tb * D + 2 * di * td
                    B = di
                    C = di - tb * 0.50
                    Zi = (d - ((2 * B**2 * td + D * tb**2)
                               / (2 * B * b - 2 * D * C)))
                    #print('web', Hi, area, Zi)
                else: # Flange
                    area = di * b
                    Zi = 0.5 * di + Hi
                    #print('flange', Hi, area, Zi)
                #
                Qout.append(area * Zi)
            return Qout
        #
        # -----------------------------------------------------
        #
        d = self.b * 0.50
        Qy = _qi(coord=coord.y, d=d, td=self.tb,
                 b=self.d, tb=self.tw)
        #
        #
        # -----------------------------------------------------
        #
        d = self.d * 0.50
        Qz = _qi(coord=coord.z, d=d, td=self.tw,
                 b=self.b, tb=self.tb)
        #
        return axis(Qy, Qz)
    #
    # --------------------------------------------
    #
    def section_coordinates(self):
        """
        1 2  3   4 5
        +-+--+---+-+
        |    :     |       ^ z
      6 +    :     + 7     |
        |    :     |       |
      8 +    :     + 9     +--> y
        |    :     |
     10 +    :     + 11
        |    :     | 
        +-+--+---+-+      
       12 13 14 15 16
        """
        # horizontal
        Yc = self.b * 0.50
        h1 = Yc - self.tw * 0.50
        h2 = Yc - self.tw
        coord_y = [-1 * h1, -1 * h2, 0, h2,  h1,   # 1-5
                   -1 * h1, h1, # 6-7
                   -1 * h1, h1, # 8-9
                   -1 * h1, h1, # 10-11
                   -1 * h1, -1 * h2, 0, h2,  h1,] # 12-16
        # vertical
        Zc = self.d * 0.50
        top1 = (Zc - self.tb * 0.50)
        top2 = (Zc - self.tb)
        top3 = -1 * top1
        coord_z = [top1, top1, top1, top1, top1, # 1-5
                   top2, top2,           # 6,7
                   0, 0,                 # 8,9
                   -1 * top2, -1 * top2, # 10,11
                   top3, top3, top3, top3, top3] # 12-16
        #
        return points(coord_y, coord_z)
    #
    def _dimension(self) -> str:
        """ Print section dimensions"""
        out = "{:<32s}{:1.4e} {:1.4e} {:1.4e}\n" \
            .format(self.type, self.d, self.b, self.b)
        out += "{:<48s}{:1.4e} {:1.4e} {:1.4e}\n" \
            .format("", self.tw, self.tb, self.tb)
        return out
    #    