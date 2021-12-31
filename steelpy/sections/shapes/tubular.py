# 
# Copyright (c) 2019-2022 steelpy
#

# Python stdlib imports
#from dataclasses import dataclass
from collections import namedtuple
import math
from typing import Union
#

# package imports
from steelpy.process.units.main import Units
from steelpy.process.io_module.text import search_line
from steelpy.sections.process.io_sections import SectionProperty, PropertyOut, get_sect_properties
#import steelpy.sections.process.print_report as print_report
from steelpy.sections.shapes.sqlite.main import SectionSQLite
from steelpy.sections.process.stress import BeamStress
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
points = namedtuple('Points', ['y', 'z'])
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
    def __init__(self):
        """
        """
        self._units = Units()
        self.type = 'Tubular'
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
        # -------------------------------------------------
        self._get_section_coordinates()
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
    #def print_properties(self, file_name:Union[str,bool]=False):
    #    """
    #    """
    #    self.properties
    #    check_out = print_report.print_header()       
    #    check_out.append("{:23s} {:1.4E} {:1.4E}"
    #                     .format(self.type, 
    #                             self.diameter.convert("millimetre").value, 
    #                             self.thickness.convert("millimetre").value))
    #    check_out.extend(print_report.print_properties(self))
    #    if file_name:
    #        #file_checkout = split_file_name(file_name)
    #        #file_checkout = str(file_checkout[0]) +'_check_me.txt'
    #        file_checkout = str(file_name) + '.txt'
    #        add_out = open(file_checkout,'w')
    #        add_out.write("".join(check_out))
    #        add_out.close()
    #    else:
    #        #for line in check_out:
    #        #    print(line.rstrip())
    #        return check_out
    #    #print('ok')
    #
    def __str__(self, units:str="si") -> str:
        """ """
        unit_sec = " m"
        unit_mas = "kg/m"
        space = " "
        output = "\n"
        #output += "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += f"{30*space}SECTION PROPERTIES [{unit_sec}]\n"
        output += "\n"
        output += "Member ID      Type      Diametre   Thickness\n"
        output += "{:}\n".format("."*80)
        output += "{:<14s} ".format(str(self.name))
        output += self._dimension()        
        output += "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += (f"{15*space}Area[{unit_sec}^2] Ixx [{unit_sec}^4] Iyy [{unit_sec}^4]"
                   f" Yp    [{unit_sec}] rx    [{unit_sec}] J   [{unit_sec}^4]\n")
        output += (f"{26*space}Sxx [{unit_sec}^3] Syy [{unit_sec}^3] SCeny [{unit_sec}]"
                   f" ry    [{unit_sec}] Cw  [{unit_sec}^6]\n")
        output += f"{26*space}Zxx [{unit_sec}^3] Zyy [{unit_sec}^3] SCenx [{unit_sec}] Mass[{unit_mas}]\n"
        output += "{:}\n".format(80*".")
        #output += "\n"
        output += "{:<14s} ".format("")
        output += self.properties.__str__()        
        return output
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
    def stress(self, actions, stress=None):
        """
        """
        #
        coord_y = self.section_coordinates.y # lateral
        coord_z = self.section_coordinates.z # vertical
        #
        D, t = self.get_geometry()
        prop = self.properties
        #
        #sigma_x = actions.Fx / prop.area
        # In Plane
        #tau_y = 2*actions.Fy / prop.area
        tau_y = [2 * actions.Fy / prop.area for _ in coord_y]
        # Out Plane
        #tau_z = 2*actions.Fy / prop.area
        tau_z = [2 * actions.Fz / prop.area for _ in coord_z]
        #
        # FIXME: what is Ip?
        #tau_x = 0 # actions.Mx * D / (2 * prop.Ip)
        tau_x = [0 for _ in coord_y]
        #sigma_z = actions.My/ prop.Zey
        #sigma_y = actions.Mz / prop.Zez
        #
        # get bending stress
        sigma_x = [(actions.Fx / prop.area) for _ in coord_y]
        sigma_y = [(actions.My * _coord / prop.Zey) for _coord in coord_z]
        sigma_z = [(actions.Mz * _coord / prop.Zez) for _coord in coord_y]
        #
        return BeamStress(sigma_x, sigma_y, sigma_z, 
                          tau_x, tau_y, tau_z, actions.x)
    #
    def _get_section_coordinates(self, theta:float=90, steps:int=6):
        """
        theta : Arch internal angle
        steps : arch division
        :return:
        arch coordinates: list[y, z]
        """
        diameter, thickness = self.get_geometry()
        radius = diameter * 0.50
        arc_length = radius * math.tau * theta / 360
        sinc = arc_length / steps
        r_theta = 360 * sinc/(radius * math.tau)
        coord_1 = [ ]
        coord_2 = [ ]
        for i in range(steps+1):
            rad = math.radians(i * r_theta)
            _x, _z = self._circunference_line(x=rad, r=radius)
            coord_1.append(_x)
            coord_2.append(_x)
        #
        coord_1.reverse()
        coord = coord_1 + [-item for item in coord_2[1:]]
        #return [coord, coord]
        self.section_coordinates = points (coord, coord)
    #
    def _circunference_line(self, x:float, r:float, xp1:float=0, yp1:float=0):
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
    @property
    def D(self):
        """
        D: diametre
        """
        return self.d

    @D.setter
    def D(self, value:Union[Units, float]):
        """
        """
        self.d = value
    #
    #
    @property
    def tw(self):
        """
        d: diametre
        """
        return self.t

    @tw.setter
    def tw(self, value: Union[ Units, float ]):
        """
        """
        self.t = value
    #
    #def set_default(self):
    #    """ """
    #    self._default = self.name    
    #
    def _dimension(self) -> str:
        """ Print section dimensions"""
        return "{:<9s} {:1.4e} {:1.4e}\n"\
               .format(self.type, self.d, self.t)
#
#
#@dataclass
class TubularInMemory(TubularBasic):
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
        TubularBasic.__init__(self)
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
    #def geometry(self, d: Union[Units, float],
    #             t: Union[Units, float]):
    #    """
    #    Parameters
    #    ----------
    #    diameter : Diameter
    #    thickness : Wall Thickness
    #    """
    #    # Geometry
    #    self.diameter = d
    #    self.thickness = t
    #    #for key, value in kwargs.items():
    #    #    _dimensions = find_tubular_dimensions(key)
    #    #    get_dimension(self, _dimensions, value)
    #    # Geometry
    #    self.compactness = get_compactness(self)
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
    #
    #
    @property
    def d(self):
        """
        D: diametre
        """
        return self.diameter

    @d.setter
    def d(self, value: Union[Units, float]):
        """
        """
        value = get_sect_properties([value])
        self.diameter = value[0]
    #
    @property
    def t(self):
        """
        d: diametre
        """
        return self.thickness
    
    @t.setter
    def t(self, value: Union[Units, float]):
        """
        """
        value = get_sect_properties([ value ])
        self.thickness = value[0]
    #
    def rolled(self):
        """
        """
        self.build='rolled'
    #
    def push_property(self):
        """ """
        self.properties
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
        Parameters
        ----------
        d : diametre
        t : wall Thickness
        Shear Stress: MAXIMUM / AVERAGE
        """
        TubularBasic.__init__(self)
        self.name = name
        self._properties = None
        self.db_file = db_file
        compactness = None
        section = (self.name, 
                   None,       # title
                   "Tubular",  # shape type
                   d, t,       # diameter, wall_thickess
                   None, None, # height, web_thickness
                   None, None, # top_flange_width, top_flange_thickness
                   None, None, # bottom_flange_width, bottom_flange_thickness
                   FAvy, FAvz,
                   shear_stress, build,
                   compactness,)        
        # push data to sqlite table
        SectionSQLite.__init__(self, db_file=self.db_file,
                               section=section)
    #
    #
    #def geometry(self, d: Union[Units, float],
    #             t: Union[Units, float]):
    #    """
    #    Parameters
    #    ----------
    #    diameter : Diameter
    #    thickness : Wall Thickness
    #    """
    #    # Geometry
    #    self.diameter = d
    #    self.thickness = t
    #    #for key, value in kwargs.items():
    #    #    _dimensions = find_tubular_dimensions(key)
    #    #    get_dimension(self, _dimensions, value)
    #    # Geometry
    #    self.compactness = get_compactness(self)
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
    #def _get_section_table(self) -> tuple:
    #    """
    #    """
    #    project = (self.name, None, "Tubular",
    #               self.diameter, self.thickness,
    #               None, None,
    #               None, None,
    #               None, None,)
    #    return project
#
#
