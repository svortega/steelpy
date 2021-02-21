# 
# Copyright (c) 2019-2021 steelpy
#

# Python stdlib imports
#from dataclasses import dataclass
import math
from typing import Union
#

# package imports
from steelpy.process.units.main import Units
from steelpy.process.io_module.text import search_line
from steelpy.f2uModel.sections.process.io_sections import SectionProperty, PropertyOut, get_sect_properties
import steelpy.f2uModel.sections.process.print_report as print_report
from steelpy.f2uModel.sections.shapes.sqlite.main import SectionSQLite
#
#
#
def find_tubular_dimensions(line_in: str) -> str:
    """
    """
    _key = {"diameter":r"\b(d(iamet(ro|er|re))?(y)?)\b",
            "thickness":r"\b((w(all)?(\_)?)?t(hickness|hk)?)\b"}

    keyWord, line_out, _match = search_line(line_in, _key)
    return keyWord
#
def get_dimension(self, dim:str, value:float):
    """
    """
    # Section Definition
    if dim == 'diameter':
        self.diameter = value
    elif dim == 'thickness':
        self.thickness = value
    else:
        raise IOError('   *** error Tubular geometry {:} not found'.format(dim))
#
def get_compactness(self) -> str:
    """
    """
    _dt = self.diameter/self.thickness
    if _dt > 80.0:
        _compactness = 'slender'
    elif _dt > 60.0:
        _compactness = 'noncompact'
    else:
        _compactness = 'compact'
    return _compactness 
#
#
# ----------------------------------------
#      Standard Sections Profiles
# ----------------------------------------
#
#
class TubularBasic:
    #__slots__ =  ['type', 'diameter', 'thickness',
    #              'compactness']    
    #
    #def __init__(self):
    #    """  
    #    """
    #    TubularSection.__init__(self)
    #    self.type = 'Tubular'
    #
    #
    def _get_properties(self):
        """
        """
        # get geometry
        diameter, thickness = self.get_geometry()
        #
        #-------------------------------------------------
        #   Cross-Sectional Area
        area = (math.pi / 4.
                * (diameter**2 -(diameter - 2*thickness)**2))
        # Centroid
        Zc = diameter / 2.0
        Yc = diameter / 2.0
        # Shear centre
        SCz = diameter * 0
        SCy = diameter * 0
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        #   --------------------------------------
        Iy = (math.pi / 64.0
              * (diameter**4 - (diameter - 2 * thickness)**4))
        Iz = Iy
        #   Elastic Modulus about Mayor Axis
        #   --------------------------------------
        Zey = (2 * Iy / diameter)
        Zez = Zey
        #-------------------------------------------------
        #   Plastic Modulus about Mayor Axis
        Zpy = ((diameter**3
                - (diameter - 2*thickness)**3) / 6.0)
        Zpz = Zpy
        #-------------------------------------------------
        #   Radius of gyration about Mayor Axis
        ry = (diameter**2 + (diameter - 2*thickness)**2)**0.50 / 4.0
        rz = ry
        #-------------------------------------------------
        # Shear Factor
        SFy = Zpy / Zey
        SFz = SFy
        #-------------------------------------------------
        #   Warping Constant Cw
        Cw = 0 * Iy
        #-------------------------------------------------
        #   Torsional Constant
        J = 2 * Iy
        #-------------------------------------------------
        #   Polar Moment of Inertia
        Ip = (math.pi/32.0
              * (diameter**4 - (diameter - 2*thickness)**4))
        #   Product of inertia
        _Iyz = 0
        Jx = 2 * Iy
        rp = (Jx / area)**0.50
        #
        #return _Area, _Zc, _Yc, _Iy, _Zey, _Zpy, _ry, _Iz, _Zez, _Zpz, _rz
        return PropertyOut(area=area, Zc=Zc, Yc=Yc,
                           Iy=Iy, Zey=Zey, Zpy=Zpy, ry=ry,
                           Iz=Iz, Zez=Zez, Zpz=Zpz, rz=rz,
                           J=J, Cw=Cw)
    #
    def shear_stress(self, stress='average', Vz=1.0, Vy=1.0, alpha = 2.0):
        """
        alpha: Shape factor (section 8.10 roakrs 7ed) 
        """
        #-------------------------------------------------        
        #            Shear Stress Calculation
        #
        # Area of Web
        # The overall depth times the web thickness
        self.Aw = self.area       
        #
        # Area of Flange
        self.Af = self.area
        #
        self.tau_z = Vz / self.Aw
        self.tau_y = Vy / self.Af
        #
        if stress != 'average':
            # Shape factor (section 8.10 roakrs 7ed)        
            self.tau_z = self.tau_z * alpha
            self.tau_y = self.tau_y * alpha
    #
    def curved(self, R:Union[Units, float]):
        """
        ---------
        R = Radio
        """
        # shear area
        warea = self.area
        D = self.diameter
        Tw = self.thickness
    
        # extreme fibre distances c
        c = D/2.0
    
        c1 = c - Tw
    
        # centroidal radius
        #_R = R
        #_R = orad - c1
    
        # Shift of neutral axis from neutral axis
        e = (c * ((2.0*R/c)- math.sqrt((R/c)**2 - 1.0) -
                    math.sqrt((R/c)**2 - (c1/c)**2))/2.0)
    
        # where
        _Ic = self.Iy
    
        # stress factors Ki
        ki = ((1.0 / (4.0*e / c))
              * ((1.0 - (e / c)) / ((R / c) - 1.0))
              * (1.0 + (c1/c)**2))
    
        # stress factors Ko
        k0 = ((1.0 / (4.0*e / c))
              * ((1.0 + (e / c)) / ((R / c) + 1.0))
              * (1.0 + (c1/c)**2))
    
        # Modulus of rigidity factor (section 8.10)
        F = 2.0
        # Shear factor (section 8.1 equ 8.1-13)
        return ki, k0
    #
    def print_properties(self, file_name:Union[str,bool]=False):
        """
        """
        self.properties
        check_out = print_report.print_header()       
        check_out.append("{:23s} {:1.4E} {:1.4E}"
                         .format(self.type, 
                                 self.diameter.convert("millimetre").value, 
                                 self.thickness.convert("millimetre").value))
        check_out.extend(print_report.print_properties(self))
        if file_name:
            #file_checkout = split_file_name(file_name)
            #file_checkout = str(file_checkout[0]) +'_check_me.txt'
            file_checkout = str(file_name) + '.txt'
            add_out = open(file_checkout,'w')
            add_out.write("".join(check_out))
            add_out.close()
        else:
            #for line in check_out:
            #    print(line.rstrip())
            return check_out
        #print('ok')
    #
    def __str__(self) -> str:
        return  ("{:23s} {:1.4E} {:1.4E}"
                 .format("Tubular", 
                         self.d.convert("millimetre").value, 
                         self.t.convert("millimetre").value))
    #
    @property
    def properties(self):
        """
        """
        if not self._properties:
            self._properties = self._get_properties()
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
        self._properties = SectionProperty(*values)
    #     
    #
    def set_default(self):
        """ """
        self.cls._default = self.name    
    #    
#
#
#@dataclass
class TubularInmemory(TubularBasic):
    """ """
    __slots__ = ['diameter', 'thickness', 'build', 
                 'shear_stress', 'compactness', '_properties',
                 'FAvy', 'FAvz', 'name', 'number']
    
    def __init__(self, name:Union[str, int],
                 d:Union[float,None], t:Union[float,None], 
                 build:str = 'welded', 
                 shear_stress:str = 'average',
                 FAvy:float = 1.0, FAvz:float = 1.0):
        """
        Shear Stress: MAXIMUM / AVERAGE
        """
        self.name = name
        self.build = build
        self.shear_stress = shear_stress
        # Shear factor
        self.FAvy = FAvy
        self.FAvz = FAvz
        #
        self.compactness:Union[str,None] = None
        self._properties = None
        #
        if d:
            self.diameter = d
        if t:
            self.thickness = t            
    #
    def geometry(self, d: Union[Units, float], 
                 t: Union[Units, float]):
        """
        Parameters  
        ----------
        diameter : Diameter    
        thickness : Wall Thickness
        """
        # Geometry
        self.diameter = d
        self.thickness = t
        #for key, value in kwargs.items():
        #    _dimensions = find_tubular_dimensions(key)
        #    get_dimension(self, _dimensions, value)
        # Geometry
        self.compactness = get_compactness(self)
    #
    def get_geometry(self):
        """ """
        try:
            diameter = self.diameter.value
            thickness = self.thickness.value
        except AttributeError:
            diameter = self.diameter
            thickness = self.thickness
        return diameter, thickness
    #
    @property
    def D(self):
        """
        D: diametre
        """
        return self.diameter
    
    @D.setter
    def D(self, value):
        """
        """
        self.diameter = value
    #
    #
    @property
    def d(self):
        """
        D: diametre
        """
        return self.diameter

    @d.setter
    def d(self, diameter: Union[Units, float]):
        """
        """
        self.diameter = diameter
    #
    @property
    def t(self):
        """
        d: diametre
        """
        return self.thickness
    
    @t.setter
    def t(self, thickness: Union[Units, float]):
        """
        """
        self.thickness = thickness  
    #
    def rolled(self):
        """
        """
        self.build='rolled'
    #
#
#
class TubularSQLite(TubularBasic, SectionSQLite):
    """ """
    __slots__ = ['_properties', # 'diameter', 'thickness', 
                 'name', 'number', 'db_file']
    
    def __init__(self, name:Union[str, int],
                 d:Union[float,None], t:Union[float,None], 
                 db_file:str,
                 build:str = 'welded', 
                 shear_stress:str = 'maximum',
                 FAvy:float = 1.0, FAvz:float = 1.0):
        """
        Shear Stress: MAXIMUM / AVERAGE
        """
        self.name = name
        #self.compactness:Union[str,None] = None
        self._properties = None
        #
        self.db_file = db_file
        compactness = None
        section = (self.name, 
                   None,      # title
                   "Tubular", # type
                   d, t,
                   None, None,
                   None, None,
                   None, None,
                   FAvy, FAvz,
                   shear_stress, build,
                   compactness,)        
        #
        SectionSQLite.__init__(self, db_file=self.db_file,
                               section=section)
    #
    #
    def geometry(self, d: Union[Units, float], 
                 t: Union[Units, float]):
        """
        Parameters  
        ----------
        diameter : Diameter    
        thickness : Wall Thickness
        """
        # Geometry
        self.diameter = d
        self.thickness = t
        #for key, value in kwargs.items():
        #    _dimensions = find_tubular_dimensions(key)
        #    get_dimension(self, _dimensions, value)
        # Geometry
        self.compactness = get_compactness(self)
    #
    def get_geometry(self):
        """ """
        return self.d, self.t
    #
    @property
    def d(self):
        """
        D: diameter
        """
        return self.get_item(item="diameter")

    @d.setter
    def d(self, diameter: Union[Units, float]):
        """
        """
        diameter = get_sect_properties([diameter])
        self.update_item(item='diameter', value=diameter[0])
        self.push_property()
    #
    @property
    def t(self):
        """
        d: diametre
        """
        return self.get_item( item="wall_thickess" )
    
    @t.setter
    def t(self, thickness: Union[Units, float]):
        """
        """
        thickness = get_sect_properties([ thickness])
        self.update_item(item='wall_thickess', value=thickness[0])
        self.push_property()
    #
    #
    def _get_section_table(self) -> tuple:
        """
        """
        project = (self.name, None, "Tubular",
                   self.diameter, self.thickness,
                   None, None,
                   None, None,
                   None, None,)
        return project
#
#
#
#@dataclass
class HollowSemicircle:
    """
    Calculate the section properties of a Hollow Semicircle

    Parameters
    ----------
    diameter  : Diameter
    thickness : Thickness wall

    Returns
    ----------
    area: Section area
    Zc  : Elastic neutral centre
    Yc  : Elastic neutral centre
    Iy  : Second moment of area about mayor axis
    Wey : Elastic modulus about mayor axis
    Wpy : Plastic modulus about mayor axis
    SFy : Shape factor mayor axis
    ry  : Radius of gyration about mayor Axis
    Iz  : Second moment of area about minor axis
    Wez : Elastic modulus about minor axis
    Wpz : Plastic modulus about minor axis
    SFz : Shape factor minor axis
    rz  : Radius of gyration about minor Axis
    SC  : Shear centre
    Cw  : Warping constant

    Notes
    ----------
    Uses formulas from:
    1.- Roark's formulas for stress and strain [7th Edition]
    2.- Wikipedia

    Examples
    ----------

    """
    #
    def __init__(self, cls):
        """ """
        #TubularSection.__init__(self, cls) 
        self.type = 'Symmetrical Hollow Semicircle Section'
    #
    def _get_property(self):
        """
        """
        #
        try:
            diameter = self.diameter.value
            thickness = self.thickness.value
        except AttributeError:
            diameter = self.diameter
            thickness = self.thickness
        #
        R = diameter / 2.0
        _Ri = R  - thickness
        _b = (R + _Ri) / 2.0
        
        #-------------------------------------------------
        #   Cross-Sectional Area
        area = (math.pi * (R**2 - _Ri**2)) / 2.0
        
        # Centroid
        Zc = (2 * _b / math.pi) * (1 + (thickness / _b)**2 / 12.0)
        
        Zc1 = R - Zc
        
        Yc = 0
        _Yc2 = R
        
        # Shear Centre
        _SCz = 2 * R
        _SCy = 0
        
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        
        #   Second Moment of Area about Mayor Axis
        #   --------------------------------------
        Iy = (((math.pi / 8.0) * (R**4 - _Ri**4)) 
              - ((8.0 / (9 * math.pi)) 
                 * ((R**3 - _Ri**3)**2) / (R**2 - _Ri**2)))
        
        Iz = ((math.pi / 8.0) * (R**4 - _Ri**4))
        
        #   Elastic Modulus about Mayor Axis
        #   --------------------------------------
        Zey = min(Iy / Zc, Iy / Zc)
        
        Zez = Iz / _Yc2
        
        #-------------------------------------------------
        #   Plastic Modulus about Mayor Axis
        _C = thickness / float(R)
        
        # Plastic neutral axis
        Zp = (R * (0.7071 - 0.2716 * _C 
                   - 0.4299 * _C**2 + 0.3983 * _C**3))
        
        Yp = 0
        
        # Plastic section moduli mayor axis
        Zpy = (R**2 * thickness * (0.8284 - 0.9140 * _C +
                                   0.7245 * _C**2 - 0.2850 * _C**3))
        
        # Plastic section moduli minor axis
        Zpz = 0.6667 * (R**3 - _Ri**3)
        
        #-------------------------------------------------
        #   Radius of gyration
        ry = math.sqrt(Iy / area)
        rz = math.sqrt(Iz / area)
        
        #-------------------------------------------------
        #   Torsional Constant
        J = ((thickness * R**5 / 3.0) 
             * (2 * math.pi**3 
                - (12 * math.pi**2 / math.pi)))
        #
        #
        return PropertyOut(area=area, Zc=Zc, Yc=Yc,
                           Iy=Iy, Zey=Zey, Zpy=Zpy, ry=ry,
                           Iz=Iz, Zez=Zez, Zpz=Zpz, rz=rz,
                           J=J, Cw=Cw)
    #
    #def print_file(self, file_name):
    #
    #    check_out = print_header()
    #
    #    check_out.append("{:23s} {:1.4E} {:1.4E}"
    #                     .format(self.type, self.diameter, self.thickness))
    #
    #    check_out.extend(print_properties(self))
    #
    #    #file_checkout = split_file_name(file_name)
    #    #file_checkout = str(file_checkout[0]) +'_check_me.txt'
    #    file_checkout = str(file_name) + '.txt'
    #    add_out = open(file_checkout,'w')
    #    add_out.write("".join(check_out))
    #    add_out.close()
    #    print('ok')
#
#
