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
class TeeBasic:
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
    #
    def __init__(self):
        """
        """
        self._units = Units()
        self.type = 'T Section'
        # Build [WELDED / ROLLED]
        self.build = 'welded'
    #  
    #
    def _get_properties(self):
        """
        """
        #
        #self.d *= factors[0]
        #self.tw *= factors[0]
        #self.a *= factors[0]
        #self.ta *= factors[0]
        #self.b *= factors[0]
        #self.tb *= factors[0]
        #
        _C = self.b - self.tw
        _h = self.d - self.tb / 2.0
        _D2 = self.d - self.tb
        
        #-------------------------------------------------
        #   Cross-Sectional Area
        area = self.b * self.tb + self.tw * _D2
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        Zc = (((self.d**2 * self.tw) + (_C * self.tb**2)) /
                   (2 * (self.b*self.tb + _D2*self.tw)))
        Yc = 0
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
        return PropertyOut(area=area, Zc=Zc, Yc=Yc,
                             Iy=Iy, Zey=Zey, Zpy=Zpy, ry=ry,
                             Iz=Iz, Zez=Zez, Zpz=Zpz, rz=rz,
                             J=J, Cw=Cw)
    #
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
    def _dimension(self) -> str:
        """ Print section dimensions"""
        out = "{:<32s}{:1.4e} {:1.4e}\n"\
               .format(self.type, self.d, self.b)
        out += "{:<48s}{:1.4e} {:1.4e}\n"\
               .format("", self.tw, self.tb)         
        return out
    #
    #def set_default(self):
    #    """ """
    #    self.cls._default = self.name    
    #
#
#
class TeeSQLite(TeeBasic, SectionSQLite):
    __slots__ = ['name', 'number', 'db_file', '_properties']

    def __init__(self, name: Union[ str, int ],
                 d: Union[ float, None ], tw: Union[ float, None ],
                 b: Union[ float, None ], tb: Union[ float, None ],
                 db_file: str,
                 build: str = 'welded',
                 shear_stress: str = 'maximum',
                 FAvy: float = 1.0, FAvz: float = 1.0):
        """ """
        TeeBasic.__init__(self)
        self.name = name
        self._properties = None
        self.db_file = db_file
        compactness = None
        section = (self.name,
                   None,       # title
                   "Tee",      # shape type
                   None, None, # diameter, wall_thickess
                   d, tw,      # height, web_thickness
                   b, tb,      # top_flange_width, top_flange_thickness
                   None, None,      # bottom_flange_width, bottom_flange_thickness
                   FAvy, FAvz,
                   shear_stress, build,
                   compactness,)
        # push data to sqlite table
        SectionSQLite.__init__(self, db_file=self.db_file,
                               section=section)
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
class TeeInMemory(TeeBasic):
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
        TeeBasic.__init__( self )
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
#