# 
# Copyright (c) 2019-2021 steelpy
#

# Python stdlib imports
import math
from collections import namedtuple
from dataclasses import dataclass
import re
from typing import NamedTuple, List, Union

# package imports
from steelpy.process.units.main import Units
from steelpy.sections.process.io_sections import SectionProperty, PropertyOut, get_sect_properties
from steelpy.sections.shapes.sqlite.main import SectionSQLite

# ----------------------------------------
#      Basic Solid Shapes
# ----------------------------------------
#
points = namedtuple('Points', ['y', 'z'])
#
#
@dataclass
class SolidSection:
    """ """
    __slots__ = ['shear_stress', 'compactness', 'build',
                 'FAvy', 'FAvz', 'name', 'number', 'cls',
                 '_properties', 'units', 'depth', 'width']
    
    def __init__(self, cls):
        """ """
        self.units = Units()
        self.cls = cls
        self.build:str = 'welded'
        self.shear_stress:str = 'average'
        self.compactness:Union[str,None] = None
        # Shear factor
        self.FAvy:float = 1.0
        self.FAvz:float = 1.0
        #
        #self.name = None
    #
    #
    #def set_default(self):
    #    """ """
    #    self.cls._default = self.name
    #
    @property
    def properties(self):
        """
        """
        try:
            self._properties
        except AttributeError:
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
    @property
    def height(self):
        """
        d : height of rectangular bar
        """
        return self.depth    
    #
    @property
    def d(self):
        """
        d : Section height of rectangular bar
        """
        return self.depth
    @d.setter
    def d(self, value):
        """
        d : Section height of rectangular bar
        """
        self.depth = value
    #
    @property
    def w(self):
        """
        w : width of rectangular bar
        """
        return self.width
    @w.setter
    def w(self, value):
        """
        w : width of rectangular bar
        """
        self.width = value
    #
    @property
    def wb(self):
        """
        wb : width bottom of rectangular bar
        """
        return self.width
    
    @wb.setter
    def wb(self, value):
        """
        wb : width bottom of rectangular bar
        """
        self.width = value
    #
    @property
    def wt(self):
        """
        wt : width top of rectangular bar
        """
        return self.a
    
    @wt.setter
    def wt(self, value):
        """
        wt : width top of rectangular bar
        """
        self.a = value
    #
#
#
class RectangleBasic:
    """
    Calculate the section properties of a rectangular solid section\n

+   +-----+
    |     |
d   |     |   Z
    |     |   ^
+   +-----+   + > Y
    *  w  *

    Parameters
    ----------
    d : Height
    w : Width

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
    #
    def __init__(self):
        """
        """
        self.type = 'SquareBar'
    #    SolidSection.__init__(self, cls)
    #    #self.name = 'rectangular bar'
    #    self.type = 'rectangular'
    #
    #
    @property
    def depth(self):
        return self.d

    @ depth.setter
    def depth(self, value:Union[Units,float]):
        """ """
        self.d = value
    #
    #
    @property
    def width(self):
        return self.w

    @width.setter
    def width(self, value:Union[Units,float]):
        """ """
        self.w = value
    #
    #def geometry(self, d:Union[float,Units], w:Union[float,Units]):
    #    """
    #    """
    #    self.depth = float(d) * self.units.m
    #    self.width = float(w) * self.units.m
    #    #self.type = 'rectangular bar'
    #
    #
    def _get_properties(self):
        """
        """
        #-------------------------------------------------
        #   Cross-Sectional Area
        area = self.width * self.depth
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        Zc = self.depth / 2.0
        Yc = 0 * self.depth
        #-------------------------------------------------
        #   Plastic Neutral Centre 
        _Zp = 0
        _Yp = 0
        #-------------------------------------------------
        #   Shear Centre 
        SCz = 0 * self.depth
        SCy = 0 * self.width
        #-------------------------------------------------
        #   Warping Constant Cw
        Cw = 0 * self.width
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
        _Iyz = 0.0
        Jx = self.width*self.depth*(self.width**2 + self.depth**2) / 12.0
        rp = ((self.width**2 + self.depth**2) / 12.0)**0.50
        #
        #-------------------------------------------------
        self._get_section_coordinates()        
        #
        #-------------------------------------------------
        #
        return PropertyOut(area=area, Zc=Zc, Yc=Yc,
                           Iy=Iy, Zey=Zey, Zpy=Zpy, ry=ry,
                           Iz=Iz, Zez=Zez, Zpz=Zpz, rz=rz,
                           J=J, Cw=Cw)
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
    def shear_stress(self, Vz=1.0, Vy=1.0, stress_type='average'):
        """
        """
        #-------------------------------------------------        
        #            Shear Stress Calculation
        #
        # get section's coordinates
        coord_y = self.section_coordinates.y # lateral
        coord_z = self.section_coordinates.z # vertical        
        # Area of Web
        # The overall depth times the web thickness
        #self.Aw = self.area       
        #
        # Area of Flange
        #self.Af = self.area
        #
        tau_z = Vz / self.area
        tau_y = Vy / self.area
        #
        #
        if stress_type != 'average':
            # Shape factor (section 8.10 roakrs 7ed)
            #_alpha = 3.0 / 2.0  
            #tau_z *= _alpha
            #tau_y *= _alpha
            #
            qz = [0.50 * (self.depth**2 / 4 - _z**2) * self.width 
                  for _z in coord_z]
            qy = [0.50 * (self.width**2 / 4 - _y**2) * self.depth 
                  for _y in coord_y]
            #
            tau_y = [_qy * Vy / (self.Iz * self.depth) for _qy in qy]
            tau_z = [_qz * Vz / (self.Iy * self.width) for _qz in qz]
            #
        else:
            tau_y = [tau_y for _ in coord_y]
            tau_z = [tau_z for _ in coord_z]

        return tau_y, tau_z
    #
    def torsional_stress(self, T=1.0):
        """
        Roark Torsion chapter
        """
        d, w = self.get_geometry()
        a = w / 2.0
        b = d / 2.0
        if a == b:
            K = 2.25 * a**4
            tau_max = 0.601 * T / a**3
        elif a > b :
            K = (a * b**3 * (16/3.0 - 3.36 * b / a 
                               * (1.0 - b**4 / (12.0 * a**4))))
            
            tau_max = ((3 * T * (1.0 + 0.6095 * b / a 
                                 + 0.8865 * (b / a)**3
                                 + 0.91 * (b / a)**4))
                       / (8 * a * b**2))
        else:
            raise ValueError(' section not applicable')
        return tau_max
    #
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
    def _print_section_properties(self):
        """
        """
        file = shape_io.print_header()
        file.extend(self._shape())
        file.extend(shape_io.print_properties(self))
        return file
    # 
    def stress(self, actions, stress=None):
        """
        """
        # get section's coordinates
        coord_y = self.section_coordinates.y # lateral
        coord_z = self.section_coordinates.z # vertical
        #
        tau_y, tau_z = self.shear_stress(actions.Fz, actions.Fy, 
                                           stress_type=self.shear_stress_type)
        #
        sigma_x = [actions.Fx / self.area for _ in coord_y]
        sigma_y = [actions.My * _coord / self.Iy for _coord in coord_z]
        sigma_z = [actions.Mz * _coord / self.Iz for _coord in coord_y]
        tau_x = [tau_y[x] * 0 for x in range(len(tau_y))]

        
        if stress:
            if isinstance(stress.tau_x, list):
                stress.tau_x = self._combine_stress(tau_x, stress.tau_x)
                stress.tau_y = self._combine_stress(tau_y, stress.tau_y)
                stress.tau_z = self._combine_stress(tau_z, stress.tau_z)
                #
                stress.sigma_x = self._combine_stress(sigma_x, stress.sigma_x)
                stress.sigma_y = self._combine_stress(sigma_y, stress.sigma_y)
                stress.sigma_z = self._combine_stress(sigma_z, stress.sigma_z)
            else:
                # Assuming global stress
                stress.tau_x = self._add_global_stress(tau_x, stress.tau_x)
                stress.tau_y = self._add_global_stress(tau_y, stress.tau_y)
                stress.tau_z = self._add_global_stress(tau_z, stress.tau_z)
                #
                stress.sigma_x = self._add_global_stress(sigma_x, stress.sigma_x)
                stress.sigma_y = self._add_global_stress(sigma_y, stress.sigma_y)
                stress.sigma_z = self._add_global_stress(sigma_z, stress.sigma_z)
        else:
            stress = PlateStress(sigma_x, sigma_y, sigma_z, 
                                 tau_x, tau_y, tau_z)
        #
        return stress
    #
    def _add_global_stress(self, stress_local, stress_global):
        """
        """  
        _new_stress = [ _item + math.copysign(1, _item.value) * abs(stress_global)  
                        if _item.value != 0  else stress_global 
                        for _item in stress_local] #aldh6850
        
        return _new_stress
    #
    def _combine_stress(self, stress_1, stress_2):
        """
        """
        _steps = len(stress_2)
        _new_stress = [stress_1[x] + stress_2[x] for x in range(_steps)]
        return _new_stress    
    #
    def _get_section_coordinates(self):
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
        _width = self.width * 0.50
        coord_y = [-1 * _width, 0 * _width, _width, 
                   -1 * _width, 0 * _width, _width, 
                   -1 * _width, 0 * _width, _width]
        # vertical
        _h = self.depth * 0.50
        coord_z = [_h , _h , _h, 
                   0 * _h, 0 * _h, 0 * _h, 
                   -1 * _h, -1 * _h, -1 * _h]
        
        self.section_coordinates = points(coord_y, coord_z)
        #print('ok')
    #
    def _dimension(self) -> str:
        """ """
        return  ("{:32s}{:1.4E} {:1.4E} {:1.4E}\n"
                 .format(self.type, self.d, self.w, self.w))
#
#
class RectangleSQLite(RectangleBasic, SectionSQLite):
    __slots__ = ['_properties', 'name', 'number', 'db_file']

    def __init__(self, name:Union[str, int],
                 d: Union[float, Units], w: Union[float, Units],
                 db_file:str,
                 build:str = 'welded',
                 shear_stress:str = 'maximum',
                 FAvy:float = 1.0, FAvz:float = 1.0):
        """
        Parameters
        ----------
        d : Height
        w : Width
        """
        RectangleBasic.__init__(self)
        self.name = name
        self._properties = None
        self.db_file = db_file
        compactness = None
        section = (self.name,
                   None,  # title
                   "Rectangle",   # shape type
                   None, None,    # diameter, wall_thickess
                   d, None,       # height, web_thickness
                   w, None,       # top_flange_width, top_flange_thickness
                   w, None,       # bottom_flange_width, bottom_flange_thickness
                   FAvy, FAvz,
                   shear_stress, build,
                   compactness,)
        # push data to sqlite table
        SectionSQLite.__init__(self, db_file=self.db_file, section=section)
    #
    #
    @property
    def d(self):
        return self.get_item(item="height")

    @d.setter
    def d(self, value:Union[Units,float]):
        """ """
        value = get_sect_properties([value])
        self.update_item(item='height', value=value[0])
        self.push_property()
    #
    #
    @property
    def w(self):
        return self.get_item(item="top_flange_width")

    @w.setter
    def w(self, value:Union[Units,float]):
        """ """
        value = get_sect_properties([value])
        self.update_item(item='top_flange_width', value=value[0])
        self.push_property()
    #
    #def __setattr__(self, shape_type:str, value:Union[Units,float]):
    #    """ """
    #    value = get_sect_properties([value])
    #    if re.match (r"\b(d(epth)?|h(eight)?)\b", shape_type, re.IGNORECASE):
    #        self.update_item(item='height', value=value[0])
    #        self.push_property()
    #    elif re.match (r"\b(w(idth)?)\b", shape_type, re.IGNORECASE):
    #        self.update_item(item='top_flange_width', value=value[0])
    #        self.push_property()
    #    else:
    #        # in python3+ you can omit the arguments to super:
    #        super().__setattr__(shape_type, value[0])
    #
    #def __getattr__(self, shape_type:str):
    #    """ """
    #    if re.match (r"\b(d(epth)?|h(eight)?)\b", shape_type, re.IGNORECASE):
    #        return self.get_item(item="height")
    #    elif re.match (r"\b(w(idth)?)\b", shape_type, re.IGNORECASE):
    #        return self.get_item(item="top_flange_width")
    #    else:
    #        raise AttributeError(shape_type)
#
class RectangleInMemory(RectangleBasic):

    def __init__(self, name:Union[str, int],
                 d:Union[float,None], w:Union[float,None],
                 build:str = 'welded',
                 shear_stress:str = 'average',
                 FAvy:float = 1.0, FAvz:float = 1.0):
        """
        Parameters
        ----------
        name : section name/id
        d    : Height
        w    : Width
        Shear Stress: MAXIMUM / AVERAGE
        """
        RectangleBasic.__init__(self)
        self.name = name
        self.build = build
        self.shear_stress = shear_stress
        # Shear factor
        self.FAvy = FAvy
        self.FAvz = FAvz
        #
        self.compactness:Union[str, None] = None
        self._properties = None
        self.type = 'rectangular'
        #
        if d:
            self.d = d
        if w:
            self.w = w
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
    def w(self):
        return self._w

    @w.setter
    def w(self, value:Union[Units,float]):
        """ """
        value = get_sect_properties([value])
        self._w = value[0]
    #
    #def __setattr__(self, shape_type:str, value:Union[Units,float]):
    #    """ """
    #    value = get_sect_properties([value])
    #    if re.match (r"\b(d(epth)?|h(eight)?)\b", shape_type, re.IGNORECASE):
    #        self.depth = value[0]
    #    elif re.match (r"\b(w(idth)?)\b", shape_type, re.IGNORECASE):
    #        self.width = value[0]
    #    else:
    #        # in python3+ you can omit the arguments to super:
    #        super().__setattr__(shape_type, value[0])
    ##
    #def __getattr__(self, shape_type:str):
    #    """ """
    #    if re.match (r"\b(d(epth)?|h(eight)?)\b", shape_type, re.IGNORECASE):
    #        return self.depth
    #    elif re.match (r"\b(w(idth)?)\b", shape_type, re.IGNORECASE):
    #        return self.width
    #    else:
    #        raise AttributeError(shape_type)
    #
    def push_property(self):
        """ """
        self.properties
#
#
class CircleBasic:
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
    #
    def __init__(self):
        """
        """
        self.type = 'RoundBar'
    #
    def _get_properties(self):
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
        Zc = self.d / 2.0
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
        Iy = math.pi * R**4 / 4.0
        #   Elastic Modulus about Mayor Axis
        Zey = math.pi * R**3 / 4.0
        #   Plastic Modulus about Mayor Axis
        Zpy = 4 * math.pi * R**3 / 3.0
        #   Shape Factor
        SFy = 1.698
        #   Radius of gyration about Mayor Axis
        ry = R / 2.0
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        Iz = math.pi * R**4 / 4.0
        #   Elastic Modulus about Minor Axis
        Zez = math.pi * R**3 / 4.0
        #   Plastic Modulus about Minor Axis
        Zpz = 4 * math.pi * R**3 / 3.0
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
        return PropertyOut(area=area, Zc=Zc, Yc=Yc,
                           Iy=Iy, Zey=Zey, Zpy=Zpy, ry=ry,
                           Iz=Iz, Zez=Zez, Zpz=Zpz, rz=rz,
                           J=J, Cw=Cw)
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
    def shear_stress(self, stress='average', Vz=1.0, Vy=1.0):
        """
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
            _alpha = 4.0 / 3.0            
            self.tau_z = self.tau_z * _alpha
            self.tau_y = self.tau_y * _alpha
        
    #
    def torsional_stress(self, T):
        """
        """
        _r = self.depth / 2.0
        K = math.pi * _r**4
        tau_max = 2 * T / (math.pi * _r**3)
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
    def _dimension(self) -> str:
        """ """
        return  ("{:9s} {:1.4E}\n"
                 .format(self.type, self.d))
#
#
class CircleSQLite(CircleBasic, SectionSQLite):
    __slots__ = ['_properties', 'name', 'number', 'db_file']
    
    def __init__(self, name:Union[str, int],
                 d:Union[float,None], 
                 db_file:str,
                 build:str = 'welded', 
                 shear_stress:str = 'maximum',
                 FAvy:float = 1.0, FAvz:float = 1.0):
        """
        Parameters
        ----------
        d : diametre
        Shear Stress: MAXIMUM / AVERAGE
        """
        CircleBasic.__init__(self)
        self.name = name
        self._properties = None
        self.db_file = db_file
        compactness = None
        section = (self.name, 
                   None,            # title
                   "Circular Bar",  # shape type
                   d, None,         # diameter, wall_thickess
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
#
#
class CircleInMemory(CircleBasic):
    __slots__ = ['diameter', 'build', 
                 'shear_stress', 'compactness', '_properties',
                 'FAvy', 'FAvz', 'name', 'number']
    
    def __init__(self, name:Union[str, int],
                 d:Union[float,None], 
                 build:str = 'welded', 
                 shear_stress:str = 'average',
                 FAvy:float = 1.0, FAvz:float = 1.0):
        """
        Shear Stress: MAXIMUM / AVERAGE
        """
        CircleBasic.__init__(self)
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
        if d: self.diameter = d          
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
        value = get_sect_properties ( [ value ] )
        self.diameter = value[0]
    #
    #
    def push_property(self):
        """ """
        self.properties
#
#
class Trapeziod(SolidSection):
    """
    Calculate the section properties of a trapezoidal solid section\n  
    
        | c |  wt  |
    +   +   +------+
           *        *
    d     *          *     Z
         *            *    ^
    +   +--------------+   + > Y
        |      wb      |
    
    Parameters
    ----------
    d  : Section height
    wb : Width bottom
    wt : Width top (default Wb)
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
    
    def __init__(self, cls):
        """ """
        SolidSection.__init__(self, cls)
        self.type = 'trapeziodal'

    #
    def geometry(self, d:Union[float,Units], wb:Union[float,Units],
                 wt:Union[float,Units]=0, c:Union[float,Units]=0) -> None:
        """
        d  : Section height
        wb : Width bottom
        wt : Width top (default Wb)
        c  : (default [wb-wt]/2)
        :return:
        """
        self.depth = float(d)  #* self.units.m
        self.width = float(wb) #* self.units.m
        #
        self.a = self.width
        if wt != 0:
            self.a = float(wt) #* self.units.m
        #else:
        #    self.a = self.width
        #    #self.width += 0.000001 * self.units.m
        #
        self.c = abs(self.a - self.width) / 2.0
        if c != 0:
            self.c = float(c) #* self.units.m
        #else:
        #    self.c = abs(self.a - self.width) / 2.0
        #    #self.c = max(abs(self.a - self.width) / 2.0, 0.000001*self.units.m)
        #
        #print('--')
    #
    def _get_properties(self):
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
        try:
            depth = self.depth.value
            width = self.width.value
            a = self.a.value
            c = self.c.value
            try:
                a = self.a.value
            except AttributeError:
                a = width
            try:
                c = self.c.value
            except AttributeError:
                c = abs(a - width) / 2.0            
        except AttributeError:
            depth = self.depth
            width = self.width
            try:
                a = self.a
            except AttributeError:
                a = width
            try:
                c = self.c
            except AttributeError:
                c = abs(a - width) / 2.0
        #-------------------------------------------------
        #   Cross-Sectional Area
        area = depth * (a + width) / 2.0
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        Zc = (depth / 3.0 * ((2 * a + width) / (a + width)))

        Yc = ((2*width**2 + 2*a*width - c*width
               - 2*c*a - a**2)/(3*(a + width)))
        #
        #self.Yc = ((2*width**2 + 2*a*width - c*width
        #            - 2*c*a - a**2)/(3*(a + width)))
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
        return shape_io.PropertyOut(area=area, Zc=Zc, Yc=Yc,
                                    Iy=Iy, Zey=Zey, Zpy=Zpy, ry=ry,
                                    Iz=Iz, Zez=Zez, Zpz=Zpz, rz=rz,
                                    J=J, Cw=Cw)
    #
    def shear_stress(self, Vz=1.0, Vy=1.0, stress_type='average'):
        """
        """
        #-------------------------------------------------
        #            Shear Stress Calculation
        #
        # get section's coordinates
        coord_y = self.section_coordinates.y # lateral
        coord_z = self.section_coordinates.z # vertical
        # Area of Web
        # The overall depth times the web thickness
        #self.Aw = self.area
        #
        # Area of Flange
        #self.Af = self.area
        #
        tau_z = Vz / self.area
        tau_y = Vy / self.area
        #
        #
        if stress_type != 'average':
            # Shape factor (section 8.10 roakrs 7ed)
            #_alpha = 3.0 / 2.0
            #tau_z *= _alpha
            #tau_y *= _alpha
            #
            qz = [0.50 * (self.depth**2 / 4 - _z**2) * self.width
                  for _z in coord_z]
            qy = [0.50 * (self.width**2 / 4 - _y**2) * self.depth
                  for _y in coord_y]
            #
            tau_y = [_qy * Vy / (self.Iz * self.depth) for _qy in qy]
            tau_z = [_qz * Vz / (self.Iy * self.width) for _qz in qz]
            #
        else:
            tau_y = [tau_y for _ in coord_y]
            tau_z = [tau_z for _ in coord_z]

        return tau_y, tau_z
    #
    def curved(self, R:Union[Units, float]):
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
    def print_file(self, file_name):
        """ """
        check_out = print_header ()
        check_out.append ( "{:23s} {:>19} {:1.4E}\n"
                           .format ( self.type, "", self.depth ) )
        check_out.append ( "{:>64} {:1.4E}\n"
                           .format ( "", self.a ) )
        check_out.append ( "{:>64} {:1.4E}\n"
                           .format ( "", self.width ) )
        check_out.extend ( print_properties ( self ) )
        # file_checkout = split_file_name(file_name)
        # file_checkout = str(file_checkout[0]) +'_check_me.txt'
        file_checkout = str ( file_name ) + '.txt'
        add_out = open ( file_checkout, 'w' )
        add_out.write ( "".join ( check_out ) )
        add_out.close ()
        print ( 'ok' )
        #
    #
    #
    def _get_section_table(self) -> tuple:
        """
        """
        project = (self.name, None, self.type,
                   None, None,
                   self.height.value, None,
                   self.width.value, None,
                   self.a.value, None,)
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
#
#