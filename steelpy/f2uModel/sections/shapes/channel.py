# 
# Copyright (c) 2019-2021 steelpy
#
#
# Python stdlib imports
import math
from typing import Union
import sys
#

# package imports
from steelpy.process.units.main import Units
from steelpy.f2uModel.sections.shapes.sqlite.main import SectionSQLite
from steelpy.f2uModel.sections.process.io_sections import SectionProperty, PropertyOut, get_sect_properties

#
#
#
#
# ----------------------------------------
#      Standard Sections Profiles
# ----------------------------------------
#
class ChannelBasic:
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
    def __init__(self):
        """ """
        # Build [WELDED / ROLLED]
        self.build = 'welded'
        self.type = 'Channel'
        self._units = Units ()
    #
    def _get_properties(self):
        """ """
        #self.a *= factors[0]
        #self.ta *= factors[0]
        _b = self.b - 0.50 * self.tw
        _C = self.b - self.tw
        _h = self.d - self.tb
        _D2 = self.d - 2 * self.tb
        #-------------------------------------------------
        #   Cross-Sectional Area
        area = (_h * self.tw) + (2 * _b * self.tb)
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        Zc = 0.50 * self.d
        Yc = ((2 * self.b**2 * self.tb + _D2 * self.tw**2) /
                   (2 * self.b * self.d - 2 * _D2 * _C))
        #-------------------------------------------------
        #   Shear Centre 
        SCz = 0.50 * self.d
        SCy = ((3 * self.tb * self.b**2)
                    / (6 * _b * self.tb + _h * self.tw))
        #-------------------------------------------------
        #   Warping Constant Cw
        Cw = ((_b**3 * _h**2 * self.tb / 12.0) *
                   ((2 * _h * self.tw + 3 * _b * self.tb) /
                    (_h * self.tw + 6 * _b * self.tb)))
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        Iy = (self.b * self.d**3 - _C * _D2**3) / 12.0
        #   Elastic Modulus about Mayor Axis
        Zey = 2 * Iy / self.d
        #   Plastic Modulus about Mayor Axis
        Zpy = (((self.d**2 * self.tw) / 4.0) +
                    (self.tb * _C * (self.d - self.tb)))
        #   Shape Factor
        SFy = Zpy * self.d / (2 * Iy)
        #   Radius of gyration about Mayor Axis
        ry = math.sqrt(Iy / area)
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        Iz = ((self.d * self.b**3 / 3.0) - (_C**3 * _D2 / 3.0)
                   - (area * (self.b - Yc)**2))
        #   Elastic Modulus about Minor Axis
        Zez = Iz / (self.b - Yc)
        #   Plastic Modulus about Minor Axis
        if (2*self.tb*_C) > (self.d*self.tw):
            Zpz = ((_C**2 *self.tb / 2.0)
                        - (self.d**2 * self.tw**2 / (8 * self.tb)) +
                        (self.d * self.tw * self.b / 2.0))
        else:
            Zpz = (((self.tw**2 *self.d) / 4.0) +
                        (self.tb * _C * (self.b - (self.tb * _C / self.d))))
        #   Shape Factor
        SFz = (Zpz * (self.b - Yc) / Iz)
        #   Radius of gyration about Minor Axis 
        rz = math.sqrt(Iz / area)
        #-------------------------------------------------
        #   Torsional Constant
        J = (2 * _b * self.tb**3 + _h * self.tw**3) / 3.0
        #   Product of inertia
        _Iyz = 0.0
        Jx = Iy + Iz
        rp = math.sqrt(Jx / area)
        #
        #
        return PropertyOut( area=area, Zc=Zc, Yc=Yc,
                             Iy=Iy, Zey=Zey, Zpy=Zpy, ry=ry,
                             Iz=Iz, Zez=Zez, Zpz=Zpz, rz=rz,
                             J=J, Cw=Cw)
    #
    @property
    def properties(self):
        """
        """
        if not self._properties:
            self._properties = self._get_properties ()
        return self._properties

    @properties.setter
    def properties(self, values):
        """
        --------------------------
        General Beam Element Data
        --------------------------

        Parameters
        ----------
        area: Section area
        Zc  : Elastic neutral centre
        Yc  : Elastic neutral centre

        Iy  : Second moment of area about mayor axis
        Zy : Elastic modulus about mayor axis
        Sy : Plastic modulus about mayor axis
        Avy : Shear area mayor axis
        ry  : Radius of gyration about mayor Axis

        Iz  : Second moment of area about minor axis
        Zz : Elastic modulus about minor axis
        Sz : Plastic modulus about minor axis
        Avz : Shear area minor axis
        rz  : Radius of gyration about minor Axis

        SCz  : Shear centre about z axis
        SCy  : Shear centre about y axis

        Cw  : Warping constant
        """
        self._properties = SectionProperty( *values )
    #
    def rolled(self):
        """
        """
        self.build='rolled'
    #
    def curved(self, R):
        """
        ---------
        R = Radio
        """
        _b = self.b
        _b1 = 2*self.tw
        _t = self.tb
        _d = self.d
    
        # shear area
        _warea = self.area
    
        # extreme fibre distances c
        _c = (_d * (((_b1 / _b) + (1.0 - (_b1 / _b))*(_t / _d)**2) /
                    (2.0*((_b1 / _b) + (1.0 - (_b1 / _b))*(_t / _d)))))
    
        _c1 = _c * ((_d / _c) - 1.0)
    
        # centroidal radius
        _R = R
        #_R = orad - _c
    
        # Shift of neutral axis from neutral axis
        _e = (_c * ((_R/_c) - (((_d/_c)*(_b1/_b + (1.0 - _b1/_b)*(_t/_d))) / 
                               (((_b1/_b) * math.log((_d/_c + _R/_c -1.0) / 
                                                     ((_d/_c)*(_t/_d) + _R/_c - 1.0))) +
                                math.log(((_d/_c)*(_t/_d) + _R/_c - 1.0) /
                                         (_R/_c - 1.0))))))
    
        # where
        _Ic = ((_warea * _c**2) * (((((_d/_c)**2 *((_b1/_b + (1.0 - _b1/_b)*(_t/_d)**3)
                                                   / (_b1/_b + (1.0 - _b1/_b)*(_t/_d))))) 
                                    / 3.0) - 1.0))
    
        # stress factors Ki
        self.ki = ((_Ic / (_warea * _c**2 * (_R/_c - 1.0))) 
                   * ((1.0 - _e / _c) / (_e / _c)))
    
        # stress factors Ko
        self.ko = ((_Ic / (_warea * _c**2 * (_e/_c))) 
                   * ((_d/_c + _e/_c -1.0) / (_R/_c + _d/_c - 1.0)) 
                   * (1.0 / (_d/_c - 1.0)))
    
        # Modulus of rigidity factor (section 8.10)
        _F = 1.0
    
        # Shear factor (section 8.1 equ 8.1-13)
        #
    #
    #
    def _dimension(self) -> str:
        """ Print section dimensions"""
        out = "{:<32s}{:1.4e} {:1.4e} {:1.4e}\n"\
               .format(self.type, self.d, self.b, self.b)
        out += "{:<48s}{:1.4e} {:1.4e} {:1.4e}\n"\
               .format("", self.tw, self.tb, self.tb)
        return out
#    
#
class ChannelSQLite(ChannelBasic, SectionSQLite):
    __slots__ = ['name', 'number', 'db_file', '_properties']

    def __init__(self, name: Union[ str, int ],
                 d: Union[ float, None ], tw: Union[ float, None ],
                 b: Union[ float, None ], tb: Union[ float, None ],
                 db_file: str,
                 build: str = 'welded',
                 shear_stress: str = 'maximum',
                 FAvy: float = 1.0, FAvz: float = 1.0):
        """ """
        ChannelBasic.__init__(self)
        self.name = name
        self._properties = None
        self.db_file = db_file
        compactness = None
        section = (self.name,
                   None,       # title
                   "Channel",  # shape type
                   None, None, # diameter, wall_thickess
                   d, tw,      # height, web_thickness
                   b, tb,      # top_flange_width, top_flange_thickness
                   b, tb,      # bottom_flange_width, bottom_flange_thickness
                   FAvy, FAvz,
                   shear_stress, build,
                   compactness,)
        # push data to sqlite table
        SectionSQLite.__init__(self, db_file=self.db_file,
                               section=section)
    #
    #
    @property
    def d(self):
        """ """
        return self.get_item ( item="height" )

    @d.setter
    def d(self, diameter: Union[ Units, float ]):
        """
        """
        diameter = get_sect_properties ( [ diameter ] )
        self.update_item ( item='height', value=diameter[ 0 ] )
        self.push_property ()

    #
    @property
    def tw(self):
        """
        """
        return self.get_item ( item="web_thickness" )

    @tw.setter
    def tw(self, thickness: Union[ Units, float ]):
        """
        """
        thickness = get_sect_properties ( [ thickness ] )
        self.update_item ( item='web_thickness', value=thickness[ 0 ] )
        self.push_property ()
    #
    #
    #
    @property
    def b(self):
        """
        D: diameter
        """
        return self.get_item ( item="top_flange_width" )

    @b.setter
    def b(self, diameter: Union[ Units, float ]):
        """
        """
        diameter = get_sect_properties ( [ diameter ] )
        self.update_item ( item='top_flange_width', value=diameter[ 0 ] )
        self.push_property ()

    #
    @property
    def tb(self):
        """
        """
        return self.get_item ( item="top_flange_thickness" )

    @tb.setter
    def tb(self, thickness: Union[ Units, float ]):
        """
        """
        thickness = get_sect_properties ( [ thickness ] )
        self.update_item ( item='top_flange_thickness', value=thickness[ 0 ] )
        self.push_property ()
    #
#
#
class ChannelInMemory(ChannelBasic):
    __slots__ = ['name', 'build', 'shear_stress', 'compactness'
                 '_properties', 'FAvy', 'FAvz',
                 '_d', '_tw', '_b', '_tb']

    def __init__(self, name: Union[ str, int ],
                 d: Union[ float, None ], tw: Union[ float, None ],
                 b: Union[ float, None ], tb: Union[ float, None ],
                 build: str = 'welded',
                 shear_stress: str = 'maximum',
                 FAvy: float = 1.0, FAvz: float = 1.0):
        """ """
        ChannelBasic.__init__ ( self )
        self.name = name
        self.build = build
        self.shear_stress = shear_stress
        # Shear factor
        self.FAvy = FAvy
        self.FAvz = FAvz
        #
        self.compactness: Union[ str, None ] = None
        self._properties = None
        self._d = d
        self._tw = tw
        self._b = b
        self._tb = tb
    #
    #
    @property
    def d(self):
        return self._d

    @d.setter
    def d(self, value:Union[Units,float]):
        """ """
        value = get_sect_properties([value])
        self._d = value[0]
    #
    @property
    def tw(self):
        return self._tw

    @tw.setter
    def tw(self, value:Union[Units,float]):
        """ """
        value = get_sect_properties([value])
        self._tw = value[0]
    #
    #
    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, value:Union[Units,float]):
        """ """
        value = get_sect_properties([value])
        self._b = value[0]
    #
    @property
    def tb(self):
        return self._tb

    @tb.setter
    def tb(self, value:Union[Units,float]):
        """ """
        value = get_sect_properties([value])
        self._tb = value[0]
    #
    #
    def push_property(self):
        """ """
        self.properties