# 
# Copyright (c) 2009 steelpy
#
#
# Python stdlib imports
from __future__ import annotations
#from array import array
from dataclasses import dataclass
from collections import namedtuple
from typing import NamedTuple
import math
#import re
#
#
# package imports
#
from steelpy.sections.utils.shape.utils import ShapeProperty
from steelpy.sections.utils.shape.stress import BeamStress, ShapeStressBasic
#
#from steelpy.utils.dataframe.main import DBframework
#
#
#def find_tubular_dimensions(line_in: str) -> str:
#    """
#    """
#    _key = {"diameter": r"\b(d(iamet(ro|er|re))?(y)?)\b",
#            "thickness": r"\b((w(all)?(\_)?)?t(hickness|hk)?)\b"}
#
#    keyWord, line_out, _match = search_line(line_in, _key)
#    return keyWord
#
#
#
#
#def get_compactness(diameter: float, thickness: float) -> str:
#    """
#    """
#    dt = diameter / thickness
#    if dt > 80.0:
#        compactness = 'slender'
#        
#    elif dt > 60.0:
#        compactness = 'noncompact'
#        
#    else:
#        compactness = 'compact'
#        
#    return compactness
#
#
#-------------------------------------------------
#
points = namedtuple('Points', ['y', 'z', 'alpha_y', 'alpha_z'])
axis = namedtuple('Axis', ['y', 'z'])
#
#-------------------------------------------------
#
class TubularPoint(NamedTuple):
    """
                 1   
              .  +   . 
          2 +    :      + 3   ^ z
           .     :       .    |
        4 +      :       + 5  +--> y
           .     :       .
          6 +    :      + 7
              .  +   .
                 8
    """
    y: list
    z: list
    #
    d: float
    tw: float
    #
    alpha_y: list[float]
    alpha_z: list[float]
    #
    def Qarea(self):
        """
        returns:  
        Qy,z : Shear first moment of area (m^3)
        """
        #
        R = self.d * 0.50
        #alpha = self.alpha       
        #
        # -----------------------------------------------------
        #        
        Qay = [self._area(R=R, tw=self.tw, alpha=item)
               for item in self.alpha_y]
        #
        Yy = [self._Cz(R=R, tw=self.tw, alpha=item)
              for item in self.alpha_y]
        #
        #Qy = [Qay[x] * Yy[x] / self._base(R=R, alpha=item)
        #      for x, item in enumerate(self.alpha)]
        #
        Qy = []
        for x, item in enumerate(self.alpha_y):
            try:
                Qy.append(Qay[x] * Yy[x] / (2*self.tw))
            except ZeroDivisionError:
                Qy.append(0)
        #
        # -----------------------------------------------------
        #
        Qaz = [self._area(R=R, tw=self.tw, alpha=item)
               for item in self.alpha_z]        
        #
        Yz = [self._Cz(R=R, tw=self.tw, alpha=item)
              for item in self.alpha_z]
        #
        Qz = []
        for x, item in enumerate(self.alpha_z):
            try:
                Qz.append(Qaz[x] * Yz[x] / (2*self.tw))
            except ZeroDivisionError:
                Qz.append(0)        
        #
        # -----------------------------------------------------
        #
        return Qy, Qz
    #
    def _area(self, R: float, tw: float,
              alpha: float):
        """Sector of hollow circle"""
        return alpha * tw * (2 * R - tw)
    #
    def _base(self, R: float, alpha: float):
        """base """
        Xc = R * math.sin(alpha)
        return 2 * Xc
    #
    def _Cz(self, R: float, tw: float, alpha: float):
        """ """
        try:
            yc = (R * (1 - (2 * math.sin(alpha) / (3 * alpha))
                       * (1 - tw / R + 1 / (2 - tw / R))))
        except ZeroDivisionError:
            yc = tw * 0.50
        return R - yc

#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------
#
#
#
@dataclass
class TubularBasic(ShapeStressBasic):
    #name:str | int
    diameter:float
    thickness:float
    shape:str = 'tubular'
    #
    # --------------------------------------------
    #
    def _properties(self, poisson: float):
        """
        """
        # get geometry
        diameter, thickness = self.get_geometry()
        #
        # -------------------------------------------------
        #   Cross-Sectional Area
        area = (math.pi / 4.
                * (diameter ** 2 - (diameter - 2 * thickness) ** 2))
        # Centroid
        #Zc = diameter / 2.0
        #Yc = diameter / 2.0
        Yc, Zc = self.centroid
        # Shear centre
        SCz = diameter * 0
        SCy = diameter * 0
        # -------------------------------------------------
        #               Section Properties
        # -------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        #   --------------------------------------
        #Iy = (math.pi / 64.0
        #      * (diameter ** 4 - (diameter - 2 * thickness) ** 4))
        #Iz = Iy
        Iy, Iz = self.I
        #
        #   Elastic Modulus about Mayor Axis
        #   --------------------------------------
        Zey = (2 * Iy / diameter)
        Zez = Zey
        # -------------------------------------------------
        #   Plastic Modulus about Mayor Axis
        Zpy = ((diameter ** 3
                - (diameter - 2 * thickness) ** 3) / 6.0)
        Zpz = Zpy
        # -------------------------------------------------
        #   Radius of gyration about Mayor Axis
        ry = (diameter ** 2 + (diameter - 2 * thickness) ** 2) ** 0.50 / 4.0
        rz = ry
        # -------------------------------------------------
        # Shear Factor
        SFy = Zpy / Zey
        SFz = SFy
        # -------------------------------------------------
        #   Warping Constant Cw
        Cw = 0 * Iy
        # -------------------------------------------------
        #   Torsional Constant
        J = 2 * Iy
        # -------------------------------------------------
        #   Polar Moment of Inertia
        Ip = (math.pi / 32.0
              * (diameter ** 4 - (diameter - 2 * thickness) ** 4))
        #   Product of inertia
        _Iyz = 0
        Jx = 2 * Iy
        rp = (Jx / area) ** 0.50
        #
        # -------------------------------------------------
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
        """
        ro = self.d * 0.50
        #ri = ro - self.tw
        #di = self.d - 2 * self.tw
        #alpha = di / self.d 
        #
        thin = Mt / (2 * math.pi * ro**2 * self.tw)
        #
        #hollow = 2 * Mt / (math.pi * ro ** 3 * (1 - alpha ** 4))
        #roar = 2 * Mt * ro/ (math.pi * (ro ** 4 - ri ** 4) )
        #
        return thin
    #
    def curved(self, R: float):
        """
        ---------
        R = Radio
        """
        # shear area
        warea = self.area
        D = self.diameter
        Tw = self.thickness

        # extreme fibre distances c
        c = D / 2.0

        c1 = c - Tw

        # centroidal radius
        # _R = R
        # _R = orad - c1

        # Shift of neutral axis from neutral axis
        e = (c * ((2.0 * R / c) - math.sqrt((R / c) ** 2 - 1.0) -
                  math.sqrt((R / c) ** 2 - (c1 / c) ** 2)) / 2.0)

        # where
        _Ic = self.Iy

        # stress factors Ki
        ki = ((1.0 / (4.0 * e / c))
              * ((1.0 - (e / c)) / ((R / c) - 1.0))
              * (1.0 + (c1 / c) ** 2))

        # stress factors Ko
        k0 = ((1.0 / (4.0 * e / c))
              * ((1.0 + (e / c)) / ((R / c) + 1.0))
              * (1.0 + (c1 / c) ** 2))

        # Modulus of rigidity factor (section 8.10)
        F = 2.0
        # Shear factor (section 8.1 equ 8.1-13)
        return ki, k0

    #
    def __str__(self, units: str = "si") -> str:
        """ """
        unit_sec = " m"
        unit_mas = "kg/m"
        space = " "
        output = "\n"
        # output += "\n"
        output += "{:}\n".format(80 * "_")
        output += "\n"
        output += f"{30 * space}SECTION PROPERTIES [{unit_sec}]\n"
        output += "\n"
        output += "Member ID      Type      Diametre   Thickness\n"
        output += "{:}\n".format("." * 80)
        output += "{:<14s} ".format(str(self.name))
        output += self._dimension()
        output += "\n"
        output += "{:}\n".format(80 * "_")
        output += "\n"
        output += (f"{15 * space}Area[{unit_sec}^2] Ixx [{unit_sec}^4] Iyy [{unit_sec}^4]"
                   f" Yp    [{unit_sec}] rx    [{unit_sec}] J   [{unit_sec}^4]\n")
        output += (f"{26 * space}Sxx [{unit_sec}^3] Syy [{unit_sec}^3] SCeny [{unit_sec}]"
                   f" ry    [{unit_sec}] Cw  [{unit_sec}^6]\n")
        output += f"{26 * space}Zxx [{unit_sec}^3] Zyy [{unit_sec}^3] SCenx [{unit_sec}] Mass[{unit_mas}]\n"
        output += "{:}\n".format(80 * ".")
        # output += "\n"
        output += "{:<14s} ".format("")
        output += self.properties.__str__()
        return output

    #
    #
    # --------------------------------------------
    #
    @property
    def I(self):
        """Moments of inertia"""
        #   Second Moment of Area about Mayor Axis
        #   --------------------------------------
        Iy = (math.pi / 64.0
              * (self.diameter ** 4 - (self.diameter - 2 * self.thickness)**4))
        #
        return axis(Iy, Iy)        
        
    #
    @property
    def centroid(self):
        """ Elastic Neutral Centre """
        Zc = self.diameter / 2.0
        return axis(Zc, Zc)
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
        R = self.diameter * 0.50 
        coord = self.section_coordinates()
        #I = self.I
        #
        # -----------------------------------------------------
        #
        def _area(R: float, tw: float,
                  alpha: float):
            """Sector of hollow circle"""
            return alpha * tw * (2 * R - tw)
        #
        def _Cz(R: float, tw: float, alpha: float):
            """ """
            try:
                yc = (R * (1 - (2 * math.sin(alpha) / (3 * alpha))
                           * (1 - tw / R + 1 / (2 - tw / R))))
            except ZeroDivisionError:
                yc = tw * 0.50
            return R - yc
        #
        def _qi(alpha: list, Qa: list,
                Yi: list, tw: float):
            """
            qi = Ax/Ib
            """
            Qout = []
            for x, item in enumerate(alpha):
                try:
                    Qout.append(Qa[x] * Yi[x] / (2*tw))
                except ZeroDivisionError:
                    Qout.append(0)
            return Qout
        #
        # -----------------------------------------------------
        #        
        Qay = [_area(R=R, tw=self.tw, alpha=item)
               for item in coord.alpha_y]
        #
        Yy = [_Cz(R=R, tw=self.tw, alpha=item)
              for item in coord.alpha_y]
        #
        Qy = _qi(alpha=coord.alpha_y, Qa=Qay,
                 Yi=Yy, tw=self.tw)
        #
        # -----------------------------------------------------
        #
        Qaz = [_area(R=R, tw=self.tw, alpha=item)
               for item in coord.alpha_z]        
        #
        Yz = [_Cz(R=R, tw=self.tw, alpha=item)
              for item in coord.alpha_z]
        #
        Qz = _qi(alpha=coord.alpha_z, Qa=Qaz,
                 Yi=Yz, tw=self.tw)
        #
        return axis(Qy, Qz)
    #
    #
    def alpha_s(self, poisson: float):
        """Shear correction factor"""
        alpha_s = (4 + 3 * poisson) / (2 * (1 + poisson))
        return alpha_s, alpha_s
    #
    # --------------------------------------------
    #
    def section_coordinates(self, theta: float = 90, steps: int = 2):
        """
        theta : Arch internal angle
        steps : arch division
        :return:
        arch coordinates: list[y, z]
        

                 1   
              .  +  . 
          2 +    :     + 3     ^ z
           .     :       .     |
        4 +      :        + 5  +--> y
           .     :       .
          6 +    :     + 7
              .  +  .
                 8
        """
        diameter, thickness = self.get_geometry()
        radius = diameter * 0.50
        arc_length = radius * math.tau * theta / 360
        sinc = arc_length / steps
        r_theta = 360 * sinc / (radius * math.tau)
        coord_1 = []
        #coord_2 = []
        alpha = []
        for i in range(steps + 1):
            rad = math.radians(i * r_theta)
            _x, _z = self._circunference_line(x=rad, r=radius)
            coord_1.append(_x)
            #coord_2.append(_z)
            alpha.append(rad)
        #
        alpha_z = [alpha[0],
                   alpha[1], alpha[1],
                   alpha[2], alpha[2],
                   alpha[1], alpha[1],
                   alpha[0]]
        #
        alpha_y = [alpha[2],
                   alpha[1], alpha[1],
                   alpha[0], alpha[0],
                   alpha[1], alpha[1],
                   alpha[2]]
        #
        #coordx = list(reversed(coord_1))
        #coordx.extend([-item for item in coord_1[1:]])
        #coordx.extend(list(reversed(coordx[1:-1])))
        coordy = [coord_1[0],
                  -1 * coord_1[1], coord_1[1],
                  -1 * coord_1[2], coord_1[2],
                  -1 * coord_1[1], coord_1[1],
                  coord_1[0]]
        #
        #coordz = coord_1.copy()
        #coordz.extend(list(reversed(coord_1[:steps])))
        #coordz.extend([-item for item in coordz[1:-1]])
        #coordz.extend([-item for item in coord_1[1:]])
        #
        coordz = [coord_1[2],
                  coord_1[1], coord_1[1],
                  coord_1[0], coord_1[0],
                  -1 * coord_1[1], -1 * coord_1[1],
                  -1 * coord_1[2]]
        #
        #coord_1.reverse()
        #coord = coord_1 + [-item for item in coord_1[1:]]
        # return [coord, coord]
        #self.section_coordinates = points(coordx, coordz)
        return points(coordy, coordz, alpha_y, alpha_z)
        #
        #
        #return TubularPoint(y=coordy, z=coordz,
        #                    d=self.diameter,
        #                    tw=self.thickness,
        #                    alpha_y=alpha_y, 
        #                    alpha_z=alpha_z)

    #
    def _circunference_line(self, x: float, r: float, xp1: float = 0, yp1: float = 0):
        """
        Calculating the coordinates of a point on a circles
        circumference from the radius, an origin and the
        arc between the points
        """
        xp2 = xp1 + r * math.sin(x)
        yp2 = yp1 - r * (1 - math.cos(x))
        # try:
        #    xp2 = xp1 + r * math.cos(math.tau/x)
        #    yp2 = yp1 + r * math.sin(math.tau/x)
        # except ZeroDivisionError:
        #    xp2 = x
        #    yp2 = yp1
        return xp2, yp2

    #   
    #
    def _dimension(self) -> str:
        """ Print section dimensions"""
        return "{:<9s} {:1.4e} {:1.4e}\n" \
            .format(self.shape, self.d, self.t)

    #
    #
    def get_geometry(self):
        """ """
        return self.diameter, self.thickness
    #
    # --------------------------------------------
    #
    #def __getattr__(self, attr):
    #    """
    #    Getter for myattr
    #    :param attr:
    #    :return:
    #    """
    #    # if attr in self.__slots__:
    #    #    return self[attr]
    #    if re.search(r"\bd\b", attr, re.IGNORECASE):
    #        return self.diameter
    #    elif re.search(r"\bt(w)?\b", attr, re.IGNORECASE):
    #        return self.thickness
    #    else:
    #        try:
    #            return self.__dict__[attr]
    #        except KeyError:
    #            raise AttributeError(f"Variable {attr} not found")
    #
    ##
    #def __setattr__(self, attr, value):
    #    """
    #    Setter for myattr
    #    :param attr:
    #    :return:
    #    """
    #
    #    if re.search(r"\bd\b", attr, re.IGNORECASE):
    #        value = get_sect_properties([value])
    #        self.diameter = value[0]
    #    elif re.search(r"\bt(w)?\b", attr, re.IGNORECASE):
    #        value = get_sect_properties([value])
    #        self.thickness = value[0]
    #    else:
    #        try:
    #            # super().__setattr__(attr, value)
    #            self.__dict__[attr] = value
    #        except KeyError:
    #            raise AttributeError(f"Variable {attr} not found")
    #
    # --------------------------------------------
    #
    @property
    def d(self):
        """
        d : Diametre
        """
        return self.diameter
    @d.setter
    def d(self, value):
        """
        d : Diametre
        """
        self.diameter = value
    #
    @property
    def t(self):
        """
        t : wall thickness
        """
        return self.thickness
    @t.setter
    def t(self, value):
        """
        t : wall thickness
        """
        self.thickness = value
    #
    @property
    def tw(self):
        """
        tw : wall thickness
        """
        return self.thickness
    @t.setter
    def tw(self, value):
        """
        tw : wall thickness
        """
        self.thickness = value    
    #     
    #
    @property
    def Dh(self):
        """Hydrodynamic diametre"""
        return self.diameter
    #
    #
    def _data_df(self):
        """ """
        return {'type': self.shape,
                'diameter': self.d,
                'wall_thickness': self.t}
#
