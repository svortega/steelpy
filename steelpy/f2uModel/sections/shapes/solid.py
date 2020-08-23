# 
# Copyright (c) 2019 iLift
#

# Python stdlib imports
import math
from collections import namedtuple
from dataclasses import dataclass
from typing import NamedTuple, List, Union

# package imports
from steelpy.process.units.main import Units
from steelpy.f2uModel.material.main import Materials
import steelpy.f2uModel.sections.process.io_sections as shape_io
#from steelpy.sections.process.stress import BeamStress
#from iLift.load.process.actions import Actions

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
                 '_properties']
    
    def __init__(self, cls):
        """ """
        self.cls = cls
        self.build:str = 'welded'
        self.shear_stress:str = 'average'
        self.compactness:Union[str,None] = None
        # Shear factor
        self.FAvy:float = 1.0
        self.FAvz:float = 1.0
        #
        self._properties = None
    #
    #
    def set_default(self):
        """ """
        self.cls._default = self.name    
#
#
class Rectangle(SolidSection):
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
    d : Depth
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
    def __init__(self, cls):
        """
        """
        SolidSection.__init__(self, cls)
        #self.name = 'rectangular bar'
        self.type = 'rectangular bar'
    #
    #@property
    #def material(self):
    #    """
    #    """
    #    return self._material
    #
    #@material.setter
    #def material(self, value):
    #    """
    #    """
    #    self._material = value
    #
    #@property
    #def section_mass(self):
    #    """
    #    section mass in g/m
    #    """
    #    return self.area * self._material.rho   
    #
    #
    def geometry(self, d, w):
        """
        """
        self.depth = float(d)
        self.width = float(w)
        #self.type = 'rectangular bar'
       
    #
    def _get_properties(self):
        """
        """
        #-------------------------------------------------
        #   Cross-Sectional Area
        self.area = self.width * self.depth
        
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        self.Zc = self.depth / 2.0
        
        self.Yc = 0 * self.depth
        
        #-------------------------------------------------
        #   Plastic Neutral Centre 
        _Zp = 0
        
        _Yp = 0
        
        #-------------------------------------------------
        #   Shear Centre 
        self.SCz = 0 * self.depth
        
        self.SCy = 0 * self.width
        
        #-------------------------------------------------
        #   Warping Constant Cw
        self.Cw = 0 * self.width
        
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        
        #   Second Moment of Area about Mayor Axis
        self.Iy = self.width*self.depth**3 / 12.0
        
        #   Elastic Modulus about Mayor Axis
        self.Zey = self.width*self.depth**2 / 6.0
        
        #   Plastic Modulus about Mayor Axis
        self.Zpy = self.width*self.depth**2 / 4.0
        
        #   Shape Factor
        self.SFy = 1.50
        
        #   Radius of gyration about Mayor Axis
        self.ry = self.depth / 12**0.50
        
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        self.Iz = self.width**3 *self.depth / 12.0
        
        #   Elastic Modulus about Minor Axis
        self.Zez = self.width**2 *self.depth / 6.0
        
        #   Plastic Modulus about Minor Axis
        self.Zpz = self.width**2 *self.depth / 4.0
        
        #   Shape Factor
        self.SFz = 1.50
    
        #   Radius of gyration about Minor Axis 
        self.rz = self.width / 12**0.50
        
        #-------------------------------------------------
        #   Torsional Constant
        if self.depth == self.width:
            self.J = 0.1406 * self.depth**4
            # Polar area section module
            self.Zej = 0.208 * self.depth**3
        else:
            self.J = ((self.depth**3 * self.width / 3.0) * 
                      (1 - (0.630 *  self.depth.value / self.width.value) +
                       (0.052 * self.depth.value**5 / self.width.value**5)))
            # Polar area section module
            self.Zej = (self.depth**2 * self.width 
                        / (3 + 1.8 * self.depth.value / self.width.value))
            if self.depth > self.width:
                self.J = ((self.depth * self.width**3 / 3.0) * 
                          (1 - (0.630 * self.width.value / self.depth.value) +
                           (0.052 * self.width.value**5 / self.depth.value**5)))
                # Polar area section module
                self.Zej = (self.depth * self.width**2 
                            / (3 + 1.8 * self.width.value / self.depth.value))
        #
        #   Product of inertia
        _Iyz = 0.0
        
        self.Jx = self.width*self.depth*(self.width**2 + self.depth**2) / 12.0
        
        self.rp = ((self.width**2 + self.depth**2) / 12.0)**0.50
        #
        #-------------------------------------------------
        self._get_section_coordinates()        
        #
        #-------------------------------------------------
        #file = shape_io.print_header()
        #file.extend(self._shape())
        #file.extend(shape_io.print_properties(self))
        #for row in file:
        #    print(row.rstrip())
        #
        #return _Area, _Zc, _Yc, _Iy, _Zey, _Zpy, _ry, _Iz, _Zez, _Zpz, _rz
        return shape_io.PropertyOut(area=self.area.value, Zc=self.Zc.value, Yc=self.Yc.value,
                                   Iy=self.Iy.value, Zey=self.Zey.value, Zpy=self.Zpy.value, ry=self.ry.value,
                                   Iz=self.Iz.value, Zez=self.Zez.value, Zpz=self.Zpz.value, rz=self.rz.value,
                                   J=self.J.value, Cw=self.Cw.value)
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
        self._properties = shape_io.SectionProperty(*values)
    #    
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
            _tau_y = [_qy * Vy / (self.Iz * self.depth) for _qy in qy]
            _tau_z = [_qz * Vz / (self.Iy * self.width) for _qz in qz]
            #
        else:
            _tau_y = [tau_y for _ in coord_y]
            _tau_z = [tau_z for _ in coord_z]

        return _tau_y, _tau_z
    #
    def torsional_stress(self, T=1.0):
        """
        Roark Torsion chapter
        """
        _a = self.w / 2.0
        _b = self.d / 2.0
        if _a == _b:
            K = 2.25 * _a**4
            tau_max = 0.601 * T / _a**3
        elif _a > _b :
            K = (_a * _b**3 * (16/3.0 - 3.36 * _b / _a 
                               * (1.0 - _b**4 / (12.0 * _a**4))))
            
            tau_max = ((3 * T * (1.0 + 0.6095 * _b / _a 
                                 + 0.8865 * (_b / _a)**3
                                 + 0.91 * (_b / _a)**4))
                       / (8 * _a * _b**2))
        else:
            raise ValueError(' section not applicable')
        #
        return tau_max
    #
    #
    def curved(self, R):
        """
        ---------
        R = Radio
        """
        # shear area
        _warea = self.area
        _D = self.depth
    
        # extreme fibre distances c
        _c = _D/2.0
        self.c = _c
        
        _c1 = _D - _c
        self.c1 = _c1
    
        # centroidal radius
        _R = R
        #_R = R - _c1
        self.R = _R
    
        # Shift of neutral axis from neutral axis
        _e = _c*(_R.value/_c.value 
                 - (2.0 / (math.log((_R.value/_c.value + 1.0) 
                                    / (_R.value/_c.value - 1.0)))))
        self.e = _e
        # where
        _Ic = self.Iy
    
        # stress factors Ki
        self.ki = ((1.0 / (3.0*_e / _c)) * 
                   ((1.0 - (_e / _c)) / ((_R / _c) - 1.0)))
        
    
        # stress factors Ko
        self.ko = ((1.0 / (3.0*_e / _c)) * 
                   ((1.0 + (_e / _c)) / ((_R / _c) + 1.0)))
    
        # Modulus of rigidity factor (section 8.10)
        self.F = 6.0/5.0
        
        self.tau_y, self.tau_z = self.shear_stress(stress_type='true')
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
    @property
    def height(self):
        """
        t : height of rectangular bar
        """
        return self.depth    
    #
    @property
    def d(self):
        """
        t : height of rectangular bar
        """
        return self.depth
    @d.setter
    def d(self, value):
        """
        t : width of rectangular bar
        """
        self.depth = value
    #
    @property
    def t(self):
        """
        t : width of rectangular bar
        """
        return self.width
    @t.setter
    def t(self, value):
        """
        t : width of rectangular bar
        """
        self.width = value
    #
    @property
    def w(self):
        """
        t : width of rectangular bar
        """
        return self.width
    @w.setter
    def w(self, value):
        """
        t : width of rectangular bar
        """
        self.width = value
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
        _tau_y, _tau_z = self.shear_stress(actions.Fz, actions.Fy, 
                                           stress_type=self.shear_stress_type)
        #
        _sigma_x = [actions.Fx / self.area for _ in coord_y]
        _sigma_y = [actions.My * _coord / self.Iy for _coord in coord_z]
        _sigma_z = [actions.Mz * _coord / self.Iz for _coord in coord_y]
        _tau_x = [_tau_y[x] * 0 for x in range(len(_tau_y))]

        
        if stress:
            if isinstance(stress.tau_x, list):
                stress.tau_x = self._combine_stress(_tau_x, stress.tau_x)
                stress.tau_y = self._combine_stress(_tau_y, stress.tau_y)
                stress.tau_z = self._combine_stress(_tau_z, stress.tau_z)
                #
                stress.sigma_x = self._combine_stress(_sigma_x, stress.sigma_x)
                stress.sigma_y = self._combine_stress(_sigma_y, stress.sigma_y)
                stress.sigma_z = self._combine_stress(_sigma_z, stress.sigma_z)
            else:
                # Assuming global stress
                stress.tau_x = self._add_global_stress(_tau_x, stress.tau_x)
                stress.tau_y = self._add_global_stress(_tau_y, stress.tau_y)
                stress.tau_z = self._add_global_stress(_tau_z, stress.tau_z)
                #
                stress.sigma_x = self._add_global_stress(_sigma_x, stress.sigma_x)
                stress.sigma_y = self._add_global_stress(_sigma_y, stress.sigma_y)
                stress.sigma_z = self._add_global_stress(_sigma_z, stress.sigma_z)
        else:
            stress = PlateStress(_sigma_x, _sigma_y, _sigma_z, 
                                 _tau_x, _tau_y, _tau_z)
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
    #
 
#
class Circle(SolidSection):
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
    def __init__(self, cls):
        """
        """
        SolidSection.__init__(self, cls)
        self.type = 'circular bar'
    #
    def geometry(self, D):
        """
        """
        self.depth = float(D)
        #self.type = 'circular bar'
    #
    def get_property(self):
        
        #
        if self.units_in[0]:
            _units_input = self.units_in
            
        else:
            print('  **  error input units not provided')
            print('      process terminated')
            sys.exit()            
        
        
        # units          
        try:
            _units_output = self.units_out
        
        except AttributeError:
            _units_output = self.units_in
            self.units_out = self.units_in
        
        factors = units.get_length_mass(_units_input, 
                                        _units_output)
        
        self.units_in = _units_output
        
        self.depth *= factors[0]
        
        _R = 0.50 * self.depth
        
        #-------------------------------------------------
        #   Cross-Sectional Area
        self.area = math.pi * _R**2
        
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        self.Zc = self.depth / 2.0
        
        self.Yc = 0
        
        #-------------------------------------------------
        #   Shear Centre 
        _SCz = 0
        
        _SCy = 0
        
        #-------------------------------------------------
        #   Warping Constant Cw
        _Cw = 0
        
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        
        #   Second Moment of Area about Mayor Axis
        self.Iy = math.pi * _R**4 / 4.0
        
        #   Elastic Modulus about Mayor Axis
        self.Zey = math.pi * _R**3 / 4.0
        
        #   Plastic Modulus about Mayor Axis
        self.Zpy = 4 * math.pi * _R**3 / 3.0
        
        #   Shape Factor
        self.SFy = 1.698
        
        #   Radius of gyration about Mayor Axis
        self.ry = _R / 2.0
        
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        
        self.Iz = math.pi * _R**4 / 4.0
        
        #   Elastic Modulus about Minor Axis
        self.Zez = math.pi * _R**3 / 4.0
        
        #   Plastic Modulus about Minor Axis
        self.Zpz = 4 * math.pi * _R**3 / 3.0
        
        #   Shape Factor
        self.SFz = 1.698
    
        #   Radius of gyration about Minor Axis 
        self.rz = _R / 2.0
        
        #-------------------------------------------------
        #   Torsional Constant
        self.J = math.pi * self.depth**4 / 32.0
        
        #   Product of inertia
        _Iyz = 0.0
        
        self.Jx = self.Iy + self.Iz
        
        self.rp = self.depth / math.sqrt(8.0)
        #
        #
        #return _Area, _Zc, _Yc, _Iy, _Zey, _Zpy, _ry, _Iz, _Zez, _Zpz, _rz
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
        _warea = self.area
        _D = self.depth
    
        # extreme fibre distances c
        _c = _D/2.0
    
        _c1 = _D - _c
    
        # centroidal radius
        _R = R
        #_R = R - _c1
    
        # Shift of neutral axis from neutral axis
        _e = _c*(((_R/_c) - math.sqrt((_R/_c)**2 - 1.0))/2.0)
    
        # where
        _Ic = self.Iy
    
        # stress factors Ki
        self.ki = ((1.0 / (4.0*_e / _c)) * 
                   ((1.0 - (_e / _c)) / ((_R / _c) - 1.0)))
    
        # stress factors Ko
        self.ko = ((1.0 / (4.0*_e / _c)) * 
                   ((1.0 + (_e / _c)) / ((_R / _c) + 1.0)))
    
        # Modulus of rigidity factor (section 8.10)
        self.F = 10.0/9.0
    #    
    def print_file(self, file_name):

        check_out = print_header()

        check_out.append("{:23s} {:1.4E}\n"
                         .format(self.type, self.depth))

        check_out.extend(print_properties(self))

        #file_checkout = split_file_name(file_name)
        #file_checkout = str(file_checkout[0]) +'_check_me.txt'
        file_checkout = str(file_name) + '.txt'
        add_out = open(file_checkout,'w')
        add_out.write("".join(check_out))
        add_out.close()        
        print('ok')
#
class SemiCircle(SolidSection):
    """
    Calculate the section properties of a semi-circular solid section\n
    
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
    def __init__(self, cls):
        """ """
        SolidSection.__init__(self, cls)
        self.type = 'semicircular bar'
    #
    def geometry(self, D):
        """ """
        self.depth = float(D)
    #
    def get_property(self):
        
        #
        if self.units_in[0]:
            _units_input = self.units_in
            
        else:
            print('  **  error input units not provided')
            print('      process terminated')
            sys.exit()            
        
        
        # units          
        try:
            _units_output = self.units_out
        
        except AttributeError:
            _units_output = self.units_in
            self.units_out = self.units_in
        
        factors = units.get_length_mass(_units_input, 
                                        _units_output)
        
        self.units_in = _units_output
        
        self.depth *= factors[0]
        
        _R = 0.50 * self.depth
        
        #-------------------------------------------------
        #   Cross-Sectional Area
        self.area = math.pi * _R**2 / 2.0
        
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        self.Zc = 4 * _R / (3 * math.pi)
        self.Yc = 0
        
        #-------------------------------------------------
        #   Plastic Neutral Centre 
        self.Zp = 0.4040 * _R
        self.Yp = 0
        
        #-------------------------------------------------
        #   Shear Centre 
        _SCz = 0 # (8.0 / (15.0 * math.pi)) * _R * (3 + 4 * v) / (1 + v)
        _SCy = 0
        
        #-------------------------------------------------
        #   Warping Constant Cw
        _Cw = 'N/A'
        
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        
        #-------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        self.Iy = math.pi / 8.0 * _R**4
        #   Elastic Modulus about Mayor Axis
        self.Zey = math.pi / 8.0 * _R**3
        #   Plastic Modulus about Mayor Axis
        self.Zpy = 2 / 3.0 * _R**3
        #   Shape Factor
        self.SFy = 1.698
        #   Radius of gyration about Mayor Axis
        self.ry = _R / 2.0
        
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        self.Iz = 0.1098 * _R**4 
        #   Elastic Modulus about Minor Axis
        self.Zez = 0.1908 * _R**3 
        #   Plastic Modulus about Minor Axis
        self.Zpz = 0.354 * _R**3 
        #   Shape Factor
        self.SFz = 1.856 
        #   Radius of gyration about Minor Axis 
        self.rz = 0.2643 * _R 
        
        #-------------------------------------------------
        #   Torsional Constant
        _J = 'N/A'
        
        #-------------------------------------------------
        #   Product of inertia
        _Iyz = 0.0
        self.Jx = self.Iy + self.Iz
        self.rp = math.sqrt(self.Jx / self.area)
        #
        
        # return _Area, _Zc, _Yc, _Iy, _Zey, _Zpy, _ry, _Iz, _Zez, _Zpz, _rz
    #
    def print_file(self, file_name):
    
        check_out = print_header()
    
        check_out.append("{:23s} {:1.4E}\n"
                         .format(self.type, self.depth))
    
        check_out.extend(print_properties(self))
    
        #file_checkout = split_file_name(file_name)
        #file_checkout = str(file_checkout[0]) +'_check_me.txt'
        file_checkout = str(file_name) + '.txt'
        add_out = open(file_checkout,'w')
        add_out.write("".join(check_out))
        add_out.close()
        print('ok')
#
class Trapeziod:
    """
    Calculate the section properties of a trapezoidal solid section\n
    
    Parameters
    ----------
    a : top base
    b : bottom base
    d : Section height
    c :

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
    
    def __init__(self, csl):
        """ """
        SolidSection.__init__(self, cls)
        self.type = 'trapeziodal bar'
    
    def units_input(self, **kwargs):
        """
        Input:
        ======
        length : [mandatory]  
        force  :   
        temperature : 
        gravity     : [default : 9.81ms^2]

        ------
        units [length, mass, time, temperature, force, pressure/stress]
        """

        for key, value in kwargs.items():
            _unit = units.find_unit_case(key)
            self.units_in = units.units_module(_unit, value, 
                                               self.units_in)
        
        if self.units_in[0]:
            pass
        
        else:
            print('error length unit must be provided')
            print('      program aborted')
            sys.exit()
    #
    def geometry(self, a, b, d, c=0):
        
        self.a = float(a)
        self.depth = float(d)
        self.width = float(b)
        
        if c != 0:
            self.c = float(c)
        
        else:
            self.c = (self.a - self.width) / 2.0
        
        self.type = 'trapezoidal bar'
    #
    def get_property(self):
        
        #
        if self.units_in[0]:
            _units_input = self.units_in
            
        else:
            print('  **  error input units not provided')
            print('      process terminated')
            sys.exit()            
        
        
        # units          
        try:
            _units_output = self.units_out
        
        except AttributeError:
            _units_output = self.units_in
            self.units_out = self.units_in
        
        factors = units.get_length_mass(_units_input, 
                                        _units_output)
        
        self.units_in = _units_output
        
        self.a *= factors[0]
        self.depth *= factors[0]
        self.width *= factors[0]
        self.c *= factors[0]
        
        #-------------------------------------------------
        #   Cross-Sectional Area
        self.area = self.depth * (self.a + self.width) / 2.0
        
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        self.Zc = (self.depth/2.0 - 
                   (self.depth / 3.0 * ((2 * self.a + self.width) / (self.a + self.width))))
        
        self.Yc = (self.a/2.0 - 
                   (2 * self.a**2 + 2 * self.a * self.width - self.a * self.c 
                    - 2 * self.width * self.c - self.width**2) / (3 * (self.a + self.width)))
        #-------------------------------------------------
        #   Plastic Neutral Centre 
        _Zp = 'N/A'
        _Yp = 'N/A'
        
        #-------------------------------------------------
        #   Shear Centre 
        _SCz = 'N/A'
        _SCy = 'N/A'
        
        #-------------------------------------------------
        #   Warping Constant Cw
        _Cw = 'N/A'
        
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        self.Iy = (self.depth / (36 * (self.a + self.width)) 
                   * (self.a**4 + self.width**4 + 2 * self.a * self.width * (self.a**2 + self.width**2)
                      - self.c * (self.a**3 + 3 * self.a**2 * self.width 
                                  - 3 * self.a * self.width**2 - self.width**3)
                      + self.c**2 * (self.a**2 + 4 * self.a * self.width + self.width**2)))
        #   Elastic Modulus about Mayor Axis
        self.Zey = (self.Iy / ((2 * self.a**2 + 2 * self.a * self.width - self.a * self.c 
                                - 2 * self.width * self.c - self.width**2) / (3 * (self.a + self.width))))
        #   Plastic Modulus about Mayor Axis
        _Zpy = 'N/A'
        #   Shape Factor
        _SFy = 'N/A'
        #   Radius of gyration about Mayor Axis
        self.ry = math.sqrt((self.a**2 + self.width**2) / 24.0)
        
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        self.Iz = ((self.depth**3 / 36.0) 
                   * ((self.a**2 + 4 * self.a * self.width + self.width**2) / (self.a + self.width)))
        #   Elastic Modulus about Minor Axis
        self.Zez = self.Iz / (self.depth / 3.0 * ((2*self.a + self.width) / (self.a + self.width))) 
        #   Plastic Modulus about Minor Axis
        _Zpz = 'N/A'
        #   Shape Factor
        _SFz = 'N/A'
        #   Radius of gyration about Minor Axis 
        self.rz = (self.depth * (math.sqrt(self.a**2 + 4 * self.a * self.width + self.width**2) 
                             / (math.sqrt(18) * (self.a + self.width))))
        
        #-------------------------------------------------
        #   Torsional Constant
        _J = 'N/A'
        
        #-------------------------------------------------
        #   Product of inertia
        self.Iyz = (self.depth**2 / (72 * (self.a + self.width)) 
                    * (self.width * (3 * self.a**2 - 3 * self.a * self.width - self.width**2)
                       + self.a**3 - self.c * (2 * self.a**2 
                                               + 8 * self.a * self.width + 2 * self.width**2)))
        self.Jx = self.Iy + self.Iz
        self.rp = math.sqrt(self.Jx / self.area)
        
        #
        #
        #return _Area, _Zc, _Yc, _Iy, _Zey, _Zpy, _ry, _Iz, _Zez, _Zpz, _rz 
    #
    def print_file(self, file_name):

        check_out = print_header()       

        check_out.append("{:23s} {:>19} {:1.4E}\n"
                         .format(self.type, "", self.depth))
        
        check_out.append("{:>64} {:1.4E}\n"
                         .format("",  self.a))
        
        check_out.append("{:>64} {:1.4E}\n"
                         .format("", self.width))        

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