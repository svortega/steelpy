# 
# Copyright (c) 2019-2022 steelpy
#
#
# Python stdlib imports
import math
from typing import Union
#

# package imports
from steelpy.process.units.main import Units
from steelpy.sections.shapes.sqlite.main import SectionSQLite
from steelpy.sections.process.io_sections import SectionProperty, PropertyOut, get_sect_properties

#
#
#
#
# ----------------------------------------
#      Standard Sections Profiles
# ----------------------------------------
#
class BoxBasic:
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
    #
    def __init__(self):
        """ """
        self.type = 'Box'
        self._units = Units()
    
    #def geometry(self, **kwargs):
    #    #
    #    self.type = 'Symmetrical Box Section'
    #
    #    for key, value in kwargs.items():
    #        _dim = find_section_dimensions(key)
    #
    #        get_dimension(self, _dim, value)
    #
    #    # check data
    #    try:
    #        self.b
    #
    #    except AttributeError:
    #        try:
    #            self.b = self.a
    #
    #        except AttributeError:
    #            self.b = self.d
    #
    #    try:
    #        self.tb
    #
    #    except AttributeError:
    #        try:
    #            self.tb = self.ta
    #
    #        except AttributeError:
    #            self.tb = self.tw
    #
    #    #
    #
    #
    def _get_properties(self):
        """ """
        #self.units_in = _units_output
        
        #self.d *= factors[0]
        #self.tw *= factors[0]
        
        #self.a *= factors[0]
        #self.ta *= factors[0]
        
        #self.b *= factors[0]
        #self.tb *= factors[0]
        #-------------------------------------------------
        #
        _hi = self.d - 2 * self.tb
        _bi = self.b - 2 * self.tw
        #-------------------------------------------------
        #   Cross-Sectional Area
        area = self.b * self.d - _bi * _hi
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        Zc = self.b / 2.0
        Yc = 0
        #-------------------------------------------------
        #   Shear Centre 
        SCz = 0
        SCy = 0
        #-------------------------------------------------
        #   Warping Constant Cw
        Cw = 0
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        Iy = ((self.b * self.d**3 - _bi * _hi**3) / 12.0)
        #   Elastic Modulus about Mayor Axis
        Zey = ((self.b * self.d**3 - _bi * _hi**3) / (6.0 * self.d))
        #   Plastic Modulus about Mayor Axis
        Zpy = ((self.b * self.d**2 - _bi * _hi**2) / 4.0)
        #   Shape Factor
        SFy = (1.50 * (self.d * (self.b*self.d**2 - _bi * _hi**2)) /
                    (self.b * self.d**3 - _bi * _hi**3))
        #   Radius of gyration about Mayor Axis
        ry = math.sqrt(Iy / area)
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        Iz = ((self.d *self.b**3 - _hi * _bi**3) / 12.0)
        #   Elastic Modulus about Minor Axis
        Zez = ((self.d * self.b**3 - _hi * _bi**3) / (6.0 * self.b))
        #   Plastic Modulus about Minor Axis
        Zpz = ((self.d *self.b**2 - _hi * _bi**2) / 4.0)
        #   Shape Factor
        SFz = (1.50 * (self.b * (self.d  * self.b**2 - _hi * _bi**2)) /
                    ( self.d * self.b**3 -  _hi * _bi**3))
        #   Radius of gyration about Minor Axis 
        rz = math.sqrt(Iz / area)
        #-------------------------------------------------
        #   Torsional Constant
        # Mean corner radious
        _Rc = 1.50 * (self.tw + self.tb) / 2.0
        # Mid countour length
        _p = (2*((self.d - self.tb) + (self.b - self.tw)) 
              - 2 * _Rc**2 * (4 - math.pi))
        # Enclosed Area
        _Ap = ((self.d - self.tb) * (self.b - self.tw) 
               - _Rc**2 * (4 - math.pi))
        # for thin walled sections b/t >= 10
        J = ((4*_Ap**2 * ((self.tw + self.tb) / 2.0)) / _p)
        #-------------------------------------------------
        #   Product of inertia
        _Iyz = 0.0
        Jx = Iy + Iz
        rp = math.sqrt(Jx / area)
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
        _b1 = 2* self.tw
        _t = self.tb
        _d = self.d
    
        # shear area
        _warea = self.area
    
        # extreme fibre distances c
        _c = _D/2.0
    
        _c1 = _d - _c
    
        # centroidal radius
        _R = R
        # _R = orad - _c1
    
        # Shift of neutral axis from neutral axis
        _e = (_c *((_R/_c)- ((2.0*(_t/_c + (1 - _t/_c)*(_b1/_b))) 
                             / ((math.log(((_R/_c)**2 + (_R/_c + 1)*(_t/_c) - 1.0) 
                                          / ((_R/_c)**2 - (_R/_c - 1.0)*(_t/_c) - 1.0))) 
                                + ((_b1/_b)*math.log((_R/_c - _t/_c + 1.0) 
                                                     /(_R/_c + _t/_c - 1.0)))))))
    
        # where
        _Ic = self.Iy
    
        # stress factors Ki
        self.ki = ((_Ic / (_warea * _c**2 * (_R/_c - 1.0))) 
                   * ((1.0 - _e / _c) / (_e / _c)))
    
        # stress factors Ko
        self.ko = ((_Ic / (_warea * _c**2 * (_R/_c + 1.0))) 
                   * ((1.0 + _e / _c) / (_e / _c)))
    
        # Modulus of rigidity factor (section 8.10)
        _nai = _c - _e    # neautral axis inner fiber
        _nao = _c1 + _e   # neautral axis outer fiber
    
        _D1 = _nai - _Tfb
        _D2 = _nai 
        _t1 = 2*_Tw
        _t2 = _Bfb
        _r = _rip
    
        self.F = ((1 + (((3*(_D2**2 - _D1**2)*_D1)/(2.0*_D2**3)) * (_t2/_t1 - 1.0)))
                  * (4*_D2**2 / (10*_r**2)))
        #
        # Shear factor (section 8.1 equ 8.1-13)
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
class BoxSQLite(BoxBasic, SectionSQLite):
    __slots__ = [ '_properties',
                  'name', 'number', 'db_file' ]

    def __init__(self, name: Union[ str, int ],
                 d: Union[ float, None ], tw: Union[ float, None ],
                 b: Union[ float, None ], tb: Union[ float, None ],
                 db_file: str,
                 build: str = 'welded',
                 shear_stress: str = 'maximum',
                 FAvy: float = 1.0, FAvz: float = 1.0):
        """ """
        BoxBasic.__init__(self)
        self.name = name
        self._properties = None
        self.db_file = db_file
        compactness = None
        section = (self.name,
                   None,       # title
                   "Box",  # shape type
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
    @property
    def d(self):
        """
        D: diameter
        """
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
#
class BoxInMemory(BoxBasic):
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
        BoxBasic.__init__ ( self )
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
    def push_property(self):
        """ """
        self.properties