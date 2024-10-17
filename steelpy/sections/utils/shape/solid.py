# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from array import array
from collections import namedtuple
from dataclasses import dataclass
import math
import re
from typing import NamedTuple

# package imports
from steelpy.sections.utils.shape.utils import (ShapeProperty, get_sect_list,
                                                get_prop_dict)
from steelpy.sections.utils.shape.stress import ShapeStressBasic


# ----------------------------------------
#      Basic Solid Shapes
# ----------------------------------------
#
points = namedtuple('Points', ['y', 'z'])
axis = namedtuple('Axis', ['y', 'z'])
#
# ----------------------------------------
#
#
@dataclass
class TrapezoidBasic(ShapeStressBasic):
    """ """
    
    def Qb(self):
        """
        Returns
        ----------
        Qb: Q/b
        Where: 
        Q : Shear's first moment of area
        b : cross section width
        """
        #
        coord = self.section_coordinates()
        Yc, Zc = self.centroid
        c = self.c
        #I = self.I
        #
        # -----------------------------------------------------
        #
        def _area(d: float, a: float, b: float):
            """ section area
            d: depth
            a: top width
            b: bottom width
            """
            return 0.50 * d * (a + b)        
        #
        def _base(dx: float, d: float,
                  b: float, c: float):
            """ """
            theta = math.atan(c/d)
            c1 = dx * math.tan(theta)
            b1 = b - 2 * c1
            return b1        
        #
        def _Cz(d: float, a: float, b: float,
                c: float):
            """geometric centre"""
            return (d / 3.0 * ((2 * b + a) / (a + b)))
        #
        # -----------------------------------------------------
        #
        #
        diz = [(Yc + item) if item < 0
               else ((self.width - Yc) - item)
               for item in coord.y]
        #
        minbase = min(self.a, self.width)
        #
        #
        Yy = [item / 3.0 if abs(coord.y[x]) > minbase
              else item / 2.0
              for x, item in enumerate(diz)]
        #
        Qay = [_area(d=self.depth, a=self.a, b=item)
               if abs(coord.y[x]) > minbase
               else _area(d=self.depth, a=item, b=item)
               for x, item in enumerate(diz)]
        #
        Qy = [item * Yy[x] / self.depth
              for x, item in enumerate(Qay)]
        #
        # -----------------------------------------------------
        #
        diy = [(Zc + item) if item < 0
               else ((self.depth - Zc) - item)
               for item in coord.z]
        #
        #
        b1 = _base(dx=Zc, d=self.depth,
                   b=self.width, c=c)
        theta =  math.atan(c/self.depth)
        c1 = Zc * math.tan(theta)
        #
        Zbase = [_base(dx=item, d=Zc, b=b1, c=c1)
                 for item in diy]
        #
        Yz = [_Cz(d=item, a=self.width,
                  b=Zbase[x], c=c) if coord.z[x] < 0
              else _Cz(d=item, a=self.a,
                       b=Zbase[x], c=c)
              for x, item in enumerate(diy)]
        #
        Qaz = [_area(d=item, a=self.width, b=Zbase[x])
               if coord.z[x] < 0
               else _area(d=item, a=self.a, b=Zbase[x]) 
               for x, item in enumerate(diy)]
        #
        Qz = [item * Yz[x] / Zbase[x] 
              for x, item in enumerate(Qaz)]
        #
        # -----------------------------------------------------
        #
        return Qy, Qz
    #
    @property
    def centroid(self):
        """geometric centre"""
        #
        Zc = (self.depth / 3.0 * ((2 * self.width + self.a)
                                  / (self.a + self.width)))

        Yc = ((2 * self.width**2
               + 2 * self.a * self.width - self.c * self.width
               - 2 * self.a * self.c - self.a**2)
              / (3 * (self.a + self.width)))
        
        return axis(Yc, Zc)
    #
    def alpha_s(self, poisson: float):
        """Shear correction factor"""
        alpha_s = (12 + 11 * poisson) / (10 * (1 + poisson))
        return alpha_s, alpha_s
    #
    @property
    def d(self):
        return self.depth
    #
    @property
    def Dh(self):
        """Hydrodynamic diametre"""
        return math.hypot(self.d, self.width)
#
#
@dataclass(slots=True)
class RectangleSolid(TrapezoidBasic):
    """
    Calculate the section properties of a rectangular solid section\n

+   +-----+
    |     |
h   |     |   Z
    |     |   ^
+   +-----+   + > Y
    *  b  *

    Parameters
    ----------
    h : Height
    b : width|base

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
    Cw  : Warping constant = None

    Notes
    ----------
    Uses formulas from:
    1.- Formulas for stress, strain and strucutral matrices [W.D. Pilkey]
    2.- Roark's formulas for stress and strain [7th Edition]
    3.- Wikipedia

    Examples
    ----------

    """
    #name: str | int
    depth:float
    width:float
    #
    # --------------------------------------------
    #
    @property
    def shape(self):
        return 'rectangle'    
    #
    @property
    def a(self):
        """ """
        return self.width
    #
    @property
    def c(self):
        """ """
        return abs(self.a - self.width) / 2.0 
    #
    #
    def _properties(self, poisson: float):
        """
        """
        #-------------------------------------------------
        #   Cross-Sectional Area
        area = self.width * self.depth
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        Zc = 0 # self.depth / 2.0
        Yc = 0 # self.width / 2.0
        #-------------------------------------------------
        #   Plastic Neutral Centre 
        Zp = 0
        Yp = 0
        #-------------------------------------------------
        #   Shear Centre 
        SCz = 0 * self.depth
        SCy = 0 * self.width
        #-------------------------------------------------
        #   Warping Constant Cw
        Cw = 0
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        Iy = self.width*self.depth**3 / 12.0
        #   Elastic Modulus about Mayor Axis
        Zey = self.width*self.depth**2 / 6.0
        #   Plastic Modulus about Mayor Axis
        Zpy = self.width*self.depth**2 / 4.0
        #   Shape Factor
        SFy = 1.50
        #   Radius of gyration about Mayor Axis
        ry = self.depth / 12**0.50
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        Iz = self.width**3 *self.depth / 12.0
        #   Elastic Modulus about Minor Axis
        Zez = self.width**2 *self.depth / 6.0
        #   Plastic Modulus about Minor Axis
        Zpz = self.width**2 *self.depth / 4.0
        #   Shape Factor
        SFz = 1.50
        #   Radius of gyration about Minor Axis 
        rz = self.width / 12**0.50
        #-------------------------------------------------
        #   Torsional Constant
        if self.depth == self.width:
            J = 0.1406 * self.depth**4
            # Polar area section module
            Zej = 0.208 * self.depth**3
        else:
            J = ((self.depth**3 * self.width / 3.0) 
                 * (1 - (0.630 *  self.depth / self.width) 
                    + (0.052 * self.depth**5 / self.width**5)))
            # Polar area section module
            Zej = (self.depth**2 * self.width
                   / (3 + 1.8 * self.depth / self.width))
            if self.depth > self.width:
                J = ((self.depth * self.width**3 / 3.0) 
                     * (1 - (0.630 * self.width / self.depth) 
                        + (0.052 * self.width**5 / self.depth**5)))
                # Polar area section module
                Zej = (self.depth * self.width**2
                       / (3 + 1.8 * self.width / self.depth))
        #
        #   Product of inertia
        Iyz = 0.0
        Jx = self.width*self.depth*(self.width**2 + self.depth**2) / 12.0
        rp = ((self.width**2 + self.depth**2) / 12.0)**0.50
        #
        #-------------------------------------------------
        #self._get_section_coordinates()        
        #
        #-------------------------------------------------
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
    def taux_max(self, Mt):
        """
        Mt : Torsional moment
        
        Roark Torsion chapter
        tau_max = Mt*r/J
        """
        d =  self.depth
        w =  self.width
        
        if d == w:
            tau_max = 4.81 * Mt / d**3
            #k = 2.25 * (d / 2.0) ** 4
            tau_max2 = 0.601 ** Mt / (d / 2.0)**3
        else:
            if d > w:
                tau_max = (3 * Mt / (d * w**2 * (1 - 0.630 * w / d + 0.250 * w**2 / d**2)))
                a = d / 2
                b = w / 2
            else:
                tau_max = (3 * Mt / (w * d**2 * (1 - 0.630 * d / w + 0.250 * d**2 / w**2)))
                a = w / 2
                b = d / 2
            #
            #k = (a * b**3 * (16/3.0 - 3.36 * b / a 
            #                 * (1.0 - b**4 / (12.0 * a**4))))
            #
            tau_max2 = ((3 * Mt * (1.0 + 0.6095 * b / a 
                                  + 0.8865 * (b / a)**3
                                  + 0.91 * (b / a)**4))
                        / (8 * a * b**2))            
        #
        #
        #tau_x = [tau_max for item in coord.y]
        return tau_max
    #
    #
    @property
    def I(self):
        """Moments of inertia"""
        #   Second Moment of Area about Mayor Axis
        Iy = self.width*self.depth**3 / 12.0
        #   Second Moment of Area about Minor Axis
        Iz = self.width**3 *self.depth / 12.0
        #
        return axis(Iy, Iz)
    #
    def curved(self, R):
        """
        ---------
        R = Radio
        """
        # shear area
        #warea = self.area
        #
        d = self.depth
        b = self.width
        # extreme fibre distances c
        c = d/2.0
        # 
        c1 = d - c
        #self.c1 = c1
        #
        # centroidal radius
        #_R = R
        #_R = R - c1
        #self.R = R
        # Shift of neutral axis from neutral axis
        e = c * (R/c - 2.0 / math.log((R/c + 1.0) / (R/c - 1.0)))
        #self.e = e
        # where
        #
        # Correction factors
        # stress factors Ki
        ki = (1.0 / (3.0*e / c) 
              * ((1.0 - (e / c)) / ((R / c) - 1.0)))
        # stress factors Ko
        ko = (1.0 / (3.0*e / c) 
              * ((1.0 + (e / c)) / ((R / c) + 1.0)))
    
        # Modulus of rigidity factor (section 8.10)
        F = 3/2
        
        #self.tau_y, self.tau_z = self.shear_stress(stress_type='true')
        return c, c1, e, ki, ko, F
    #
    def print_file(self, file_name):
        
        check_out = print_header()
        
        check_out.append("{:23s} {:>19} {:1.4E} {:>9} {:1.4E}\n"
                         .format(self.type, "", self.depth, "", self.width))
        
        check_out.extend(print_properties(self))
        
        #file_checkout = split_file_name(file_name)
        #file_checkout = str(file_checkout[0]) +'_check_me.txt'
        file_checkout = str(file_name) + '.txt'
        add_out = open(file_checkout,'w')
        add_out.write("".join(check_out))
        add_out.close()        
        print('ok')
    #
    #
    def _shape(self):
        """
        """
        _section = []
        _section.append("+   +-----+{:35s}{:1.3E} {:1.3E}\n"
                        .format("", self.d.convert('millimetre').value, 
                                self.w.convert('millimetre').value))
        _section.append("    |     |\n")
        _section.append("d   |     |   Z\n")
        _section.append("    |     |   ^\n")
        _section.append("+   +-----+   + > Y\n")
        _section.append("    +  w  +\n")
        return _section
    #
    #
    def section_coordinates(self):
        """
        1    2     3
        +----+-----+
        |    :     |       ^ z
        |    :     |       |
      4 +    + 5   + 6     +--> y
        |    :     |
        |    :     | 
        +----+-----+      
        7    8     9
        """
        # horizontal
        width = self.width * 0.50
        coord_y = [-1 * width, 0 * width, width,
                   -1 * width, 0 * width, width,
                   -1 * width, 0 * width, width]
        # vertical
        h = self.depth * 0.50
        coord_z = [h , h , h,
                   #0.5 * h, 0.5 * h, 0.5 * h,
                   0 * h, 0 * h, 0 * h,
                   #-0.5 * h, -0.5 * h, -0.5 * h,
                   -1 * h, -1 * h, -1 * h]
        #
        #return RectanglePoint(coord_y, coord_z,
        #                      depth=self.depth,
        #                      width=self.width,
        #                      a=self.width)
        #print('ok')
        return points(coord_y, coord_z)
    #
    def _dimension(self) -> str:
        """ """
        return  ("{:32s}{:1.4E} {:1.4E} {:1.4E}\n"
                 .format(self.shape, self.depth, self.width, self.width))
    #
    #   
    #
#
#
@dataclass(slots=True)
class TrapezoidSolid(TrapezoidBasic):
    """
    Calculate the section properties of a trapezoidal solid section\n  
    
        | c |  a  |
    +   +   +------+
           *        *
    h     *          *     Z
         *            *    ^
    +   +--------------+   + > Y
        |      b      |
    
    Parameters
    ----------
    h  : Section height
    b  : Width bottom
    a  : Width top (default Wb)
    c  : (default [wb-wt]/2)

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
    """
    #name: str | int
    depth:float
    width:float
    a:float
    c:float
    #
    # --------------------------------------------
    #
    @property
    def shape(self):
        return 'trapezoid'   
    #
    def _properties(self, poisson: float):
        """
        :return:
        area : Geometric Area (m^2)
        Zc : Z Distance to Centroid (m)
        Yc : Y Distance to Centroid (m)
        Iy : Second moment of area (m^4)
        Zey : Elastic Section Modulus (m^3)
        Zpy : Plastic Section Modulus (m^3)
        ry : Radius of gyration (m)
        Iz : Second moment of area (m^4)
        Zez : Elastic Section Modulus (m^3)
        Zpz : Plastic Section Modulus (m^3)
        rz : Radius of gyration (m)
        J : Torsional Constant (m^4)
        Cw : Warping Constant ()
        """
        # get geometry
        depth = self.depth
        width = self.width
        a = self.a
        c = self.c
        #-------------------------------------------------
        #   Cross-Sectional Area
        area = 0.50 * depth * (a + width)
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        Zc = (depth / 3.0 * ((2 * a + width) / (a + width)))
        #
        #
        #
        #Yc =  ((2*width**2 + 2*a*width - c*width
        #        - 2*c*a - a**2)/(3*(a + width)))
        #
        Yc = max(width, a) * 0.50
        #
        #-------------------------------------------------
        #   Plastic Neutral Centre 
        #_Zp = 'N/A'
        #_Yp = 'N/A'
        
        #-------------------------------------------------
        #   Shear Centre 
        #_SCz = 'N/A'
        #_SCy = 'N/A'
        
        #-------------------------------------------------
        #   Warping Constant Cw
        Cw = 0 * width
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        #
        Iy = (depth**3 / (36.0 * (a + width))
              * (a**2 + 4 * a * width + width**2))
        #   Elastic Modulus about Mayor Axis
        Zey = Iy / Zc
        #   Plastic Modulus about Mayor Axis
        #Zpy = (depth**2/(12*(a + width))
        #       * (a**2 + 4* a*width + width**2))
        factor = 2
        if c != 0:
            factor = 1
        Zpy = (depth**2/(48*factor * (a + width))
               * (11* a**2 + 26*a * width + 11*width**2))
        #   Shape Factor
        #_SFy = 'N/A'
        #   Radius of gyration about Mayor Axis
        ry = (Iy / area)**0.50
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        Iz = (depth/(36*(a + width))
              * (4*a*width*c**2 + 3*a**2 *width*c
                 - 3*a*width**2 *c + a**4 + width**4
                 + 2*a**3 *width + a**2 * c**2 + a**3 * c
                 + 2*a*width**3 - c*width**3 + width**2 * c**2))
       #
       #self.Iz = (depth / (36 * (a + width))
       #            * (width**4 + a**4 + 2*width*a * (width**2 + a**2))
       #            - c * (width**3 + 3*width**2 * a
       #                        - 3*width*a**2 - a**3)
       #            + c**2 * (width**2 + 4*width*a + a**2))
       #
        #   Elastic Modulus about Minor Axis
        Zez = Iz / Yc
        #   Plastic Modulus about Minor Axis
        Zpz = (depth/12 * (2*a**2 - a*width + 2*width**2))
        #   Shape Factor
        #_SFz = 'N/A'
        #   Radius of gyration about Minor Axis 
        rz = (Iz / area)**0.50
        #
        #-------------------------------------------------
        #   Torsional Constant
        if c == 0:
            if depth == width:
                J = 0.1406 * depth ** 4
                # Polar area section module
                Zej = 0.208 * depth ** 3
            else:
                J = ((depth ** 3 * width / 3.0) *
                     (1 - (0.630 * depth / width) +
                      (0.052 * depth ** 5 / width ** 5)))
                # Polar area section module
                Zej = (depth ** 2 * width
                       / (3 + 1.8 * depth / width))
                if depth > width:
                    J = ((depth * width ** 3 / 3.0) *
                         (1 - (0.630 * width / depth) +
                          (0.052 * width ** 5 / depth ** 5)))
                    # Polar area section module
                    Zej = (depth * width ** 2
                           / (3 + 1.8 * width / depth))
        else:
            s = abs(width-a) / depth
            Vl = 0.10504 - 0.10 * s + 0.0848 * s**2 - 0.06746 * s**3 + 0.0515 * s**4
            Vs = 0.10504 + 0.10 * s + 0.0848 * s**2 + 0.06746 * s**3 + 0.0515 * s**4
            J = (width/12.0 * (a + width) * (a**2 + width**2)
                 - Vl * width**4 - Vs * a**4)
        #
        #-------------------------------------------------
        #   Product of inertia
        Jx = Iy + Iz
        rp = (Jx / area)**0.50
        #
        # reset Yc to the centre of the section
        #self.Yc =  Yc # - max(a, width) * 0.50
        #self.Zc =  Zc
        #
        Yc = 0.50 * max(a, width) - Yc
        Zc = 0.50 * self.depth - Zc
        #
        #
        alpha_sy, alpha_sz = self.alpha_s(poisson=poisson)
        #
        #
        return ShapeProperty(area=area, Zc=Zc, Yc=Yc,
                             Iy=Iy, Sy=Zey, Zy=Zpy, ry=ry,
                             Iz=Iz, Sz=Zez, Zz=Zpz, rz=rz,
                             J=J, Cw=Cw,
                             alpha_sy=alpha_sy,
                             alpha_sz=alpha_sz)
    #
    #
    def taux_max(self, Mt):
        """
        Mt : Torsional moment
        
        Roark Torsion chapter
        """
        #coord =  self.section_coordinates()
        b =  self.depth
        m = self.width 
        n = self.a
        if self.width < self.a:
            m = self.a
            n = self.width
        #
        s = (m - n) / b
        A = (m + n) * b / 2.0
        r = (m - n) / 2.0
        D = n + r
        #
        Vl = (0.10504 - 0.10 * s + 0.0848 * s ** 2
              - 0.06746 * s ** 3 + 0.0515 * s ** 4)
        
        Vs = (0.10504 + 0.10 * s + 0.0848 * s ** 2
              + 0.06746 * s ** 3 + 0.0515 * s ** 4)
        
        K = (b / 12 * (m + n) * (m ** 2 + n ** 2)
             - Vl * m ** 4 - Vs * n ** 4)
        
        try:
            rstep = D / (2 * r)
        except ZeroDivisionError:
            rstep = 0
        
        C = (D / (1 + (math.pi ** 2 * D ** 4) / (16 * A ** 2))
             * (1 + 0.15 * ((math.pi ** 2 * D ** 4) / (16 * A ** 2)
                            -  rstep)))
        
        tau_max = Mt / K * C
        return tau_max    
    #
    #
    def curved(self, R: float):
        """
        ---------
        R: radius of curvature measured to centroid of section
        c : distance from centroidal axis to extreme fiber on concave side of beam

        :return:
        e: distance from centroidal axis to neutral axis measured toward center of curvature
        ki: sigma_i/sigma
        k0: sigma_0/sigma
        sigma_i = actual stress in extreme fiber on concave side
        sigma_0 = actual stress in extreme fiber on convex side
        sigma = fictitious unit stress in corresponding fiber as computed by ordinary flexure formula for a straight beam
        """
        b1 = self.a
        b = self.width
        d = self.depth
        # Distance from centroidal axis to extre fiber on concave side of beam
        c = d / (3*(1+b1/b)/(1+2*b1/b))
        # distance from centroidal axis to extreme fiber on convex side of beam
        c1 = c * (d/c - 1.0)
        # Conditions
        # if R/d >= 8 then the beam should be considered thin
        # if R/d < 8 then the beam should be considered thick
        thin = 1
        thick = 0
        if R / d < 8:
            thin = 0
            thick = 1.0
        #
        area = d/2 * (b + b1)
        Ic = d**3/36 * (b**2 + 4*b*b1 + b1**2)/ (b+b1)
        #
        # distance from centroidal axis to neutral axis measured towards centre of curvature
        h = ((R - c * (0.50*(1+b1/b)*(d/c)**2)
                 /((R/c + c1/c - b1/b*(R/c - 1.0)))*math.log((R/c + c1/c)/(R/c - 1.0))
                 - (1.0 - b1/b)*d/c) * thick + Ic/(R*area) * thin)
        # stress factors Ki
        ki = (1.0/(2*h/c) * (1.0 - h/c)/(R/c - 1.0)
              * (1.0 + 4*b1/b + (b1/b)**2)/(1.0+b1/b)**2)
        # stress factors Ko
        k0 = ((c1/c)/(2*h/c) * (c1/c + h/c)/(R/c + c1/c)
              * (1 + 4*b1/b + (b1/b)**2)/(2 + b1/b)**2)
        #
        # Modulus of rigidity factor (section 8.1)
        F = 3/2
        return ki, k0

    #
    #
    @property
    def I(self):
        """Moments of inertia"""
        depth = self.depth
        width = self.width
        a = self.a
        c = self.c
        #   Second Moment of Area about Mayor Axis
        Iy = (depth**3 / (36.0 * (a + width))
              * (a**2 + 4 * a * width + width**2))
        #   Second Moment of Area about Minor Axis
        Iz = (depth/(36*(a + width))
              * (4*a*width*c**2 + 3*a**2 *width*c
                 - 3*a*width**2 *c + a**4 + width**4
                 + 2*a**3 *width + a**2 * c**2 + a**3 * c
                 + 2*a*width**3 - c*width**3 + width**2 * c**2))
        #
        return axis(Iy, Iz)    
    #
    #
    def _get_section_table(self) -> tuple:
        """
        """
        project = (self.name, None, self.type,
                   None, None,
                   self.height, None,
                   self.width, None,
                   self.a, None,)
        return project
    #
    #def __str__(self) -> str:
    #    """ """
    #    output = "{:<12s} {:>19} {:1.4E}\n".format("Bar", "", self.depth)
    #    output += "{:12s}{:<12s}{:40s} {:1.4E}\n".format("", self.type.capitalize(), "", self.a)
    #    output += "{:>64} {:1.4E}\n".format("", self.width)       
    #    #return  ("{:23s} {:1.4E} {:1.4E}"
    #    #         .format(self.type, self.d, self.t))
    #    return output
    #
    def _dimension(self) -> str:
        """ Print section dimensions"""
        output = "{:<8s}{:>24}{:1.4e} {:1.4e} {:1.4e}\n"\
            .format("BarTrapz", " ", self.depth, self.a, self.width)
        #output += "{:15s}{:<8s}\n".format("", "Trapzd")
        #output += "{:67} {:1.4e}\n".format("", self.width)
        return output
    #
    #def _properties(self) -> str:
    #    """ Print section properties"""
    #    output = "{:<12s} {:>19} {:1.4E}\n".format("Bar", "", self.depth)
    #    output += "{:12s}{:<12s}{:40s} {:1.4E}\n".format("", self.type.capitalize(), "", self.a)
    #    output += "{:>64} {:1.4E}\n".format("", self.width)
    #    return output
    #
    #
    def section_coordinates(self):
        """
        1    2     3
        +----+-----+
        |    :     |       ^ z
        |    :     |       |
      4 +    + 5   + 6     +--> y
        |    :     |
        |    :     | 
        +----+-----+      
        7    8     9
        """
        # horizontal
        width = self.width * 0.50
        coord_y = [-1 * width, 0 * width, width, 
                   -1 * width, 0 * width, width, 
                   -1 * width, 0 * width, width]
        # vertical
        h = self.depth * 0.50
        coord_z = [h , h , h, 
                   0 * h, 0 * h, 0 * h, 
                   -1 * h, -1 * h, -1 * h]
        
        #return RectanglePoint(coord_y, coord_z,
        #                      depth=self.depth,
        #                      width=self.width,
        #                      a=self.a, c=self.c)
        return points(coord_y, coord_z)
#
# ------------------------------
#
@dataclass(slots=True)
class CircleSolid(ShapeStressBasic):
    """
    Calculate the section properties of a circular solid section\n
    
    Parameters
    ----------
    d : Diameter

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
    """
    #name: str | int
    d:float
    #shape: str = 'Circular'
    #
    # --------------------------------------------
    #
    @property
    def shape(self):
        return 'Circular'
    #
    # --------------------------------------------
    #
    def _properties(self, poisson: float):
        """
        """
        #
        #self.depth *= factors[0]
        R = 0.50 * self.d
        #-------------------------------------------------
        #   Cross-Sectional Area
        area = math.pi * R**2
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        Zc = 0 # self.d / 2.0
        Yc = 0
        #-------------------------------------------------
        #   Shear Centre 
        _SCz = 0
        _SCy = 0
        #-------------------------------------------------
        #   Warping Constant Cw
        Cw = 0
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        #Iy = math.pi * R**4 / 4.0
        Iy, Iz = self.I
        #   Elastic Modulus about Mayor Axis
        Zey = math.pi * R**3 / 4.0
        #   Plastic Modulus about Mayor Axis
        Zpy = 4 * R**3 / 3.0
        #   Shape Factor
        SFy = 1.698
        #   Radius of gyration about Mayor Axis
        ry = R / 2.0
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        #Iz = math.pi * R**4 / 4.0
        #   Elastic Modulus about Minor Axis
        Zez = math.pi * R**3 / 4.0
        #   Plastic Modulus about Minor Axis
        Zpz = 4 * R**3 / 3.0
        #   Shape Factor
        SFz = 1.698
        #   Radius of gyration about Minor Axis 
        rz = R / 2.0
        #-------------------------------------------------
        #   Torsional Constant
        J = math.pi * self.d**4 / 32.0
        #   Product of inertia
        _Iyz = 0.0
        Jx = Iy + Iz
        rp = self.d / math.sqrt(8.0)
        #
        #
        alpha_sy,alpha_sz = self.alpha_s(poisson=poisson)
        #
        return ShapeProperty(area=area, Zc=Zc, Yc=Yc,
                             Iy=Iy, Sy=Zey, Zy=Zpy, ry=ry,
                             Iz=Iz, Sz=Zez, Zz=Zpz, rz=rz,
                             J=J, Cw=Cw,
                             alpha_sy=alpha_sy,
                             alpha_sz=alpha_sz)
    #
    def taux_max(self, Mt):
        """
        Max shear stress due to pure torsion
        tau_max = Mt*ro/J
        """
        r = self.d * 0.50
        #K = math.pi * r**4
        #J = math.pi * r ** 4 / 2.0
        tau_max = 2 * Mt / (math.pi * r**3)
        return tau_max
    #
    def curved(self, R):
        """
        ---------
        R = Radio
        """
        # shear area
        warea = self.area
        D = self.depth
    
        # extreme fibre distances c
        c = D/2.0
    
        c1 = D - c
    
        # centroidal radius
        #_R = R
        #_R = R - c1
    
        # Shift of neutral axis from neutral axis
        e = c*(((R/c) - math.sqrt((R/c)**2 - 1.0))/2.0)
    
        # where
        _Ic = self.Iy
    
        # stress factors Ki
        self.ki = ((1.0 / (4.0*e / c)) * 
                   ((1.0 - (e / c)) / ((R / c) - 1.0)))
    
        # stress factors Ko
        self.ko = ((1.0 / (4.0*e / c)) * 
                   ((1.0 + (e / c)) / ((R / c) + 1.0)))
    
        # Modulus of rigidity factor (section 8.10)
        self.F = 4/3
    #
    # --------------------------------------------
    #
    @property
    def r(self):
        """radious"""
        return self.d * 0.50
    #
    @property
    def I(self):
        """Moments of inertia"""
        R = self.d * 0.50
        Iy = math.pi * R**4 / 4.0
        return axis(Iy, Iy)        
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
        R = self.d * 0.50
        coord = self.section_coordinates()
        #
        # -----------------------------------------------------
        #
        def _base(R: float, d: float):
            """area"""
            return 2 * math.sqrt(R**2 - d**2)        
        #
        def _qi(items: list, R: float):
            """qi = Ax/Ib"""
            Qout = []
            for item in items:
                try:
                    Qout.append((2 / 3.0 * (R ** 2 - item ** 2) ** (3 / 2))
                                / _base(R=R, d=item))
                except ZeroDivisionError:
                    Qout.append(0)
            return Qout
        #
        # -----------------------------------------------------
        #
        Qy = _qi(items=coord.y, R=R)    
        #
        # -----------------------------------------------------
        #
        Qz = _qi(items=coord.z, R=R)
        #
        # -----------------------------------------------------
        #
        return Qy, Qz
    #
    #
    #
    def alpha_s(self, poisson: float):
        """Shear correction factor"""
        alpha_s = (7 + 6 * poisson) / (6 * (1 + poisson))
        return alpha_s, alpha_s
    #
    # --------------------------------------------
    #
    def section_coordinates(self, theta: float = 90, steps: int = 2):
        """
                 1   
              .  +   . 
          2 +    :      + 3    ^ z
           .     :       .     |
        4 +      + 5      + 6  +--> y
           .     :       .
          7 +    :      + 8
              .  +   .
                 9
        """
        # horizontal
        radius = self.d * 0.50
        arc_length = radius * math.tau * theta / 360
        sinc = arc_length / steps
        r_theta = 360 * sinc / (radius * math.tau)
        #
        coord_1 = []
        for i in range(steps + 1):
            rad = math.radians(i * r_theta)
            _x, _z = self._circunference_line(x=rad, r=radius)
            coord_1.append(_x)
            #coord_2.append(_z)
        #
        #coordy = list(reversed(coord_1))
        #coordy.extend([-item for item in coord_1[1:]])
        #coordy.extend(list(reversed(coordy[1:-1])))
        coordy = [coord_1[0],
                  -1 * coord_1[1], coord_1[1],
                  -1 * coord_1[2], coord_1[0], coord_1[2],
                  -1 * coord_1[1], coord_1[1],
                  coord_1[0]]
        
        #
        #coordz = coord_1.copy()
        #coordz.extend(list(reversed(coord_1[:steps])))
        #coordz.extend([-item for item in coordz[1:-1]])
        #
        coordz = [coord_1[2],
                  coord_1[1], coord_1[1],
                  coord_1[0], coord_1[0], coord_1[0],
                  -1 * coord_1[1], -1 * coord_1[1],
                  -1 * coord_1[2]]
        #
        return points(coordy, coordz)
        #return CircularlePoint(coordy, coordz,
        #                      d=self.d)        
    #
    def _circunference_line(self, x: float, r: float, xp1: float = 0, yp1: float = 0):
        """
        Calculating the coordinates of a point on a circles
        circumference from the radius, an origin and the
        arc between the points
        """
        xp2 = xp1 + r * math.sin(x)
        yp2 = yp1 - r * (1 - math.cos(x))
        return xp2, yp2    
    #
    def _dimension(self) -> str:
        """ """
        return  ("{:9s} {:1.4E}\n"
                 .format(self.shape, self.d))
    #
    #
    @property
    def Dh(self):
        """Hydrodynamic diametre"""
        return self.d     
    #
    #
#
#
#
class SolidDim(NamedTuple):
    """ """
    shape:str
    h:float
    b:float|None
    a:float|None
    #
    FAvy:float
    FAvz:float
    shear_stress:str
    build:str
    compactness:str|None
    title:str|None
    #
    @property
    def d(self):
        """"""
        return self.h
#
#
def get_solid_section(shape:str, parameters: list|tuple|dict)->list:
    """Return : [diameter/height, base, a,
                 FAvy, FAvz, shear_stress,
                 build, compactness, title]"""
    if isinstance(parameters,(list,tuple)):
        prop = get_sect_list(parameters, number= 9, step=3)
    elif isinstance(parameters, dict):
        prop, sect = get_SolidSect_dict(shape, parameters)
        prop = [*prop[:3], *prop[4:]]
    else:
        raise IOError('Section data not valid')
    #
    match shape:
        case 'Rectangle bar':
            # check if symmetry
            if not prop[1]:
                prop[1] = prop[0]
            section = [shape,       # shape type
                        prop[0],    # h - height
                        prop[1],    # b - base bottom
                        prop[2],    # None
                        *prop[3:]]  # FAvy, FAvz, shear_stress, build, compactness, title
        case 'Trapezoid bar':
            #if not prop[4]:
            #    prop[4] = abs(prop[3] - prop[2]) / 2.0
            section = [shape,       # shape type
                        prop[0],    # h - height
                        prop[1],    # b - base bottom
                        prop[2],    # a - base top
                        *prop[3:]]  # FAvy, FAvz, shear_stress, build, compactness, title
        case 'Circular bar':
            section = [shape,       # shape type
                       prop[0],     # diameter
                       prop[1],     # None
                       prop[2],     # None                     
                       *prop[3:]]   # FAvy, FAvz, shear_stress, build, compactness, title
    #
    return SolidDim(*section)
#
#
def get_SolidSect_dict(shape:str, parameters: dict,
                  number:int = 9, step:int=4)->list:
    """Return : [diameter/height, base, a, c,
                 FAvy, FAvz, shear_stress,
                 build, compactness, title]"""
    # basic information
    section = [None] * step
    name = 'Solid'
    #
    for key, item in parameters.items():
        if re.match(r"\b(d(iamet(re|er)|epth)?|h(eight)?)\b", key, re.IGNORECASE):
            try:
                section[0] = item.value
            except AttributeError:
                raise IOError("units required")

        elif re.match(r"\b(b(ase)?)\b", key, re.IGNORECASE):
            try:
                section[1] = item.value
            except AttributeError:
                raise IOError("units required")
        
        elif re.match(r"\b(a)\b", key, re.IGNORECASE):
            try:
                section[2] = item.value
            except AttributeError:
                raise IOError("units required")
        
        elif re.match(r"\b(c)\b", key, re.IGNORECASE):
            try:
                section[3] = item.value
            except AttributeError:
                raise IOError("units required")        

        elif re.match(r"\b((section|shape)?(_|-|\s*)?(name|id))\b", key, re.IGNORECASE):
            name = item
    #
    #
    match shape:
        case 'Rectangle bar':
            if not section[1]:
                section[1] = section[0]
            properties = RectangleSolid(name=name,
                                        depth=section[0],
                                        width=section[1])
        case 'Trapezoid bar':
            if not section[1]:
                section[1] = section[0]
                #section[2] = section[0]
                properties = RectangleSolid(name=name,
                                            depth=section[0],
                                            width=section[1])
            else:
                if not section[2]:
                    section[2] = section[1]
                if not section[3]:
                    section[3] = 0
                properties = TrapezoidSolid(name=name,
                                            depth=section[0],
                                            width=section[1],
                                            a=section[2],
                                            c=section[3])
        case 'Circular bar':
            properties = CircleSolid(name=name,
                                     d=section[0])
    #
    section.extend(get_prop_dict(parameters))
    return section, properties._properties(poisson=0.30)