# 
# Copyright (c) 2019-2020 steelpy
#

# Python stdlib imports
import math
from typing import Union
#

# package imports
from steelpy.process.units.main import Units
from steelpy.process.io_module.text import search_line
from steelpy.f2uModel.sections.process.io_sections import SectionProperty
import steelpy.f2uModel.sections.process.print_report as print_report
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
class Tubular:
    """
    =====================================================   
    Calculate the section properties of a Tubular section   
    =====================================================   

    Parameters  
    ----------
    diameter  : Diameter   
    thickness : Thickness wall
    ----------
    FAvy = Shear factor y (1.0 default)
    FAvz = Shear factor z (1.0 default)
    
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
        self.name = "tubular"
        self._units = Units()
        self.build = 'welded'
        # Shear Stress [MAXIMUM / AVERAGE]
        self.shear_stress = 'average'
        self.compactness = None
        self.type = 'Tubular Section'
        #       
        # Shear factor
        self.FAvy = 1.0
        self.FAvz = 1.0
    #
    @property
    def units(self):
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
        return self._units
    #
    #
    def geometry(self, **kwargs):
        """
        Parameters  
        ----------
        diameter : Diameter    
        thickness : Wall Thickness
        """
        # Geometry
        for key, value in kwargs.items():
            _dimensions = find_tubular_dimensions(key)
            get_dimension(self, _dimensions, value)

        # Geometry
        self.compactness = get_compactness(self)
    #
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
    @property
    def t(self):
        """
        d: diametre
        """
        return self.thickness
    @t.setter
    def t(self, value):
        """
        """
        self.thickness = value    
    #
    #
    @property
    def properties(self):
        """
        """
        #
        #-------------------------------------------------
        #   Cross-Sectional Area
        self.area = (math.pi / 4. *
                    (self.diameter**2 -(self.diameter - 2*self.thickness)**2))
        # Centroid
        self.Zc = self.diameter / 2.0
        self.Yc = self.diameter * 0
        # Shear centre
        self.SCz = self.diameter * 0
        self.SCy = self.diameter * 0
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        #   --------------------------------------
        self.Iy = (math.pi / 64.0
                   * (self.diameter**4 - (self.diameter - 2 * self.thickness)**4))
        self.Iz = self.Iy
        #   Elastic Modulus about Mayor Axis
        #   --------------------------------------
        self.Zey = (2 * self.Iy / self.diameter)
        self.Zez = self.Zey
        #-------------------------------------------------
        #   Plastic Modulus about Mayor Axis
        self.Zpy = ((self.diameter**3 -
                     (self.diameter - 2*self.thickness)**3) / 6.0)
        self.Zpz = self.Zpy
        #-------------------------------------------------
        #   Radius of gyration about Mayor Axis
        self.ry = (self.diameter**2 + (self.diameter - 2*self.thickness)**2)**0.50 / 4.0
        self.rz = self.ry
        #-------------------------------------------------
        # Shear Factor
        self.SFy = self.Zpy.value / self.Zey.value
        self.SFz = self.SFy
        #-------------------------------------------------
        #   Warping Constant Cw
        _Cw = 0
        #-------------------------------------------------
        #   Torsional Constant
        self.J = (2 * self.Iy)
        #-------------------------------------------------
        #   Polar Moment of Inertia
        self.Ip = ((math.pi/32.0)*
                   (self.diameter**4 - 
                    (self.diameter - 2*self.thickness)**4))
        #   Product of inertia
        _Iyz = 0
        self.Jx = 2 * self.Iy 
        self.rp = (self.Jx / self.area)**0.50
        #
        #return _Area, _Zc, _Yc, _Iy, _Zey, _Zpy, _ry, _Iz, _Zez, _Zpz, _rz
    
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
            _alpha = 2.0
            self.tau_z = self.tau_z * _alpha
            self.tau_y = self.tau_y * _alpha
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
        # shear area
        _warea = self.area
        _D = self.diameter
        _Tw = self.thickness
    
        # extreme fibre distances c
        _c = _D/2.0
    
        _c1 = _c - _Tw
    
        # centroidal radius
        _R = R
        #_R = orad - _c1
    
        # Shift of neutral axis from neutral axis
        _e = (_c * ((2.0*_R/_c)- math.sqrt((_R/_c)**2 - 1.0) -
                    math.sqrt((_R/_c)**2 - (_c1/_c)**2))/2.0)
    
        # where
        _Ic = self.Iy
    
        # stress factors Ki
        self.ki = ((1.0 / (4.0*_e / _c)) * 
                   ((1.0 - (_e / _c)) / ((_R / _c) - 1.0)) *
                   (1.0 + (_c1/_c)**2))
    
        # stress factors Ko
        self.ko = ((1.0 / (4.0*_e / _c)) * 
                   ((1.0 + (_e / _c)) / ((_R / _c) + 1.0)) *
                   (1.0 + (_c1/_c)**2))
    
        # Modulus of rigidity factor (section 8.10)
        self.F = 2.0
    
        # Shear factor (section 8.1 equ 8.1-13)
        
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
    def __init__(self):
        
        # Build [WELDED / ROLLED]
        self.build = 'welded'
        # Shear Stress [MAXIMUM / AVERAGE]
        self.shear_stress = 'average'
        self.compactness = None        
        self.units_in = ["", "", "second", "", "", ""] 
        self.type = 'Symmetrical Hollow Semicircle Section'
    
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
    def geometry(self, **kwargs):
        """
        Parameters  
        ----------
        diameter : Diameter    
        thickness : Wall Thickness
        """
        # Geometry
        for key, value in kwargs.items():
            _dimensions = find_tubular_dimensions(key)
            get_dimension(self, _dimensions, value)

        # Geometry
        self.compactness = get_compactness(self)
    #
    #
    def get_property(self):
        """
        """
        #self.diameter *= factors[0]
        #self.thickness *= factors[0]
        #
        _R = self.diameter / 2.0
        _Ri = _R  - self.thickness
        _b = (_R + _Ri) / 2.0
           
        #-------------------------------------------------
        #   Cross-Sectional Area
        self.area = (math.pi * (_R**2 - _Ri**2)) / 2.0
        
        # Centroid
        self.Zc = (2 * _b / math.pi) * (1 + (self.thickness / _b)**2 / 12.0)
        
        self.Zc1 = _R - self.Zc
        
        self.Yc = 0
        _Yc2 = _R
        
        # Shear Centre
        _SCz = 2 * _R
        _SCy = 0
        
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        
        #   Second Moment of Area about Mayor Axis
        #   --------------------------------------
        self.Iy = (((math.pi / 8.0) * (_R**4 - _Ri**4)) 
                   - ((8.0 / (9 * math.pi)) 
                      * ((_R**3 - _Ri**3)**2) / (_R**2 - _Ri**2)))
        
        self.Iz = ((math.pi / 8.0) * (_R**4 - _Ri**4))
        
        #   Elastic Modulus about Mayor Axis
        #   --------------------------------------
        self.Zey = min(self.Iy / self.Zc, self.Iy / self.Zc)
        
        self.Zez = self.Iz / _Yc2
        
        #-------------------------------------------------
        #   Plastic Modulus about Mayor Axis
        _C = self.thickness / float(_R)
        
        # Plastic neutral axis
        self.Zp = (_R * (0.7071 - 0.2716 * _C 
                         - 0.4299 * _C**2 + 0.3983 * _C**3))
        
        self.Yp = 0
        
        # Plastic section moduli mayor axis
        self.Zpy = (_R**2 * self.thickness * (0.8284 - 0.9140 * _C +
                                       0.7245 * _C**2 - 0.2850 * _C**3))
        
        # Plastic section moduli minor axis
        self.Zpz = 0.6667 * (_R**3 - _Ri**3)
        
        #-------------------------------------------------
        #   Radius of gyration
        self.ry = math.sqrt(self.Iy / self.area)
        self.rz = math.sqrt(self.Iz / self.area)
        
        #-------------------------------------------------
        #   Torsional Constant
        self.J = ((self.thickness * _R**5 / 3.0) * (2 * math.pi**3 
                                             - (12 * math.pi**2 / math.pi)))
        #
        #
        #return self.area, _Zc, _Yc, self.Iy, _Zey, _Zpy, _ry, _Iz, _Zez, _Zpz, _rz, _Zp
    #
    def rolled(self):
        """
        """
        self.build='rolled'
    #
    def print_file(self, file_name):

        check_out = print_header()       

        check_out.append("{:23s} {:1.4E} {:1.4E}"
                         .format(self.type, self.diameter, self.thickness))

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
#
#if __name__ == "__main__":
    #
    #tubo = Tubular()
    #tubo.geometry(d=1, tw=0.025)
    #print('hello')