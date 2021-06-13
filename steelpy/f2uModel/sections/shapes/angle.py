# 
# Copyright (c) 2015 steelpy


# Python stdlib imports
import math
import sys
#

# package imports
import steelpy.units.control as units
from steelpy.sectionproperty.shapes.iomodule import (find_section_dimensions,
                                                     get_dimension)

#
#
#
#
# ----------------------------------------
#      Standard Sections Profiles
# ----------------------------------------
#
class Angle:
    """
    Calculate the section properties of an angle section

         +   +
             |         
        d    |         Z
             |         ^
         +   +-----+   + > Y
             *  b  *

    Parameters
    ----------
    d  : Section Heigh
    tw : Web/flange thickness
    b  : Base
    r  : Root Radious

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
    1.- Full plastic capacity of equal angle sections under biaxial
        bending and normal force [A.E. Charalampakis]
    2.- Formulas for stress, strain and strucutral matrices [W.D. Pilkey]
    3.- Roark's formulas for stress and strain [7th Edition]
    4.- Wikipedia

    Examples
    ----------

    """
    #
    def __init__(self):
        # Build [WELDED / ROLLED]
        self.build = 'welded'
        # Shear Stress [MAXIMUM / AVERAGE]
        self.shear_stress = 'average'
        self.compactness = 'N/A'        
        self.units_in = ["", "", "second", "", "", ""]
    #
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
        #
        #  Beam Section Definition
        for key, value in kwargs.items():
            _dim = find_section_dimensions(key)
            
            get_dimension(self, _dim, value)
        
        # check data
        try:
            self.b
        
        except AttributeError:
            try:
                self.b = self.a
            
            except AttributeError:
                self.b = self.d
        
        try:
            self.tb
        
        except AttributeError:
            try:
                self.tb = self.ta
            
            except AttributeError:
                self.tb = self.tw
        
        #
        try:
            self.r = 0.50 * float(self.r)
        
        except AttributeError:
            self.r = 0.50 * (self.tw / math.sqrt(2.0))
        #
        # Simmetric angle
        self.type = 'Asymmetrical Angle Section'
        
        if self.d == self.b :
            self.type = 'Symmetrical Angle Section'

        #
        #
    #
    def units_output(self, **kwargs):
        """
        Input:\n
        length : [mandatory]\n
        force  : [mandatory]\n
        temperature : \n
        gravity     : [default : 9.81ms^2]\n

        ------
        units [length, mass, time, temperature, force, pressure/stress]/n
        """
        _units_in = ["", "", "second", "", "", ""]
        for key, value in kwargs.items():
            _unit = units.find_unit_case(key)
            self.units_out = units.units_module(_unit, value, 
                                                _units_in)
        #
        #    
    #    
    def get_property(self):
        #-------------------------------------------------
        #
        # Simmetric angle
        if self.type == 'Symmetrical Angle Section':
            
            #-------------------------------------------------
            #   Cross-Sectional Area (A.1)
            self.area = ((2 * self.d * self.tw) 
                         - self.tw**2 + (2 - math.pi/2.0) * self.r**2)
            
            # Distance v1 of elastic centroid from heel (A.2)
            _v1 = (((6 * math.sqrt(2))*(((50/3.0 - 5*math.pi)*self.r**3)
                                        - ((2 - math.pi / 2.0)*(self.d - 3*self.tw)*self.r**2)
                                        + ((self.d**2 + self.d*self.tw - self.tw*2)*self.tw))) /
                   (((24 - 6*math.pi)*self.r**2) + (24*self.d*self.tw) - (12*self.tw**2)))
            
            #  Distance v2 (A.3)
            _v2 = ((self.d/math.sqrt(2)) + (1 - math.sqrt(2)) + 
                   (self.tw / math.sqrt(2)) - _v1)
            
            # Distance v3 (A.4)
            _v3 = self.d/math.sqrt(2)
            
            # Moment of inertia around u-u (A.5)
            _Iu = (((70 - 21*math.pi)*self.r**4 / 24.0) + 
                   ((self.d - self.tw)**2 * (math.pi - 4) * self.r**2 / 4.0) -
                   (self.tw**4 / 12.0) + (self.d**3 * self.tw / 3.0) + (self.d * self.tw**3 / 3.0) -
                   (self.d**2 * self.tw**2 / 2.0))
            
            # Elastic section modulus around u-u (A.6)
            _Wu = _Iu / _v3
            
            # Moment of inertia around v-v (A.7)
            _Iv = ((1.0 / ((288 - 72 * math.pi)*self.r**2 + 288*(self.d - (self.tw / 2.0))*self.tw)) *
                   ((((7926 * math.pi) - (1233 * math.pi**2) - 12776)*self.r**6)
                    + (432*( math.pi - (10.0 / 3.0))*(self.d - self.tw)*( math.pi - 4)*self.r**5)
                    + (((1422 * math.pi) - (36 * math.pi**2) - 4188)*self.tw**2)
                    + (((72*((-79 * math.pi/2.0) +  math.pi**2 + (349/3.0)) * self.d*self.tw) 
                       - (36 * self.d**2 * ( math.pi - 4)**2)*self.r**4))
                    + (432*( math.pi - (10.0/3.0))*(4*self.tw**2 / 3.0 +  self.d**2 -
                                                    4 * self.d*self.tw)*self.tw*self.r**3)
                    - (24*(self.d**3 - (9*self.tw * self.d**2) + (13*self.tw**2 *self.d) - (13*self.tw**3 / 4.0)) *
                       ((math.pi - 4)*self.tw*self.r**2))
                    + (-(72 * self.d*self.tw**5) + (12*self.tw**6) + (24*self.tw**2 *self.d**4) - 
                       (48*self.tw**3 *self.d**3) + (96*self.tw**4 *self.d**2))))
            
            # Elastic section modulus around v-v (A.8)
            _Wv = _Iv / _v1
            
            # Distance e (A.9)
            _e = _v1 / math.sqrt(2)
            self.Zc = _e
            self.Yc = _e
            
            # Moment of inertia around y-y or z-z (A.10)
            self.Iy = ((1.0 / (((288 - 72 * math.pi) * self.r**2) 
                               + (288 * (self.d - self.tw / 2.0) * self.tw)))
                       * (((3732 * math.pi - 5968 - 585 * math.pi**2) * self.r**6)
                          + ((216*(math.pi - 4) * (math.pi - 10.0/3.0) * (self.d - self.tw) * self.r**5)
                             + (60 * self.d**4 * self.tw**2))
                          + (-(120 * self.d**3 * self.tw**3) + (132 * self.d**2 * self.tw**4) 
                             - (72 * self.d * self.tw**5) +
                             (12 * self.tw**6) + (((216 * (math.pi - 10.0/3.0))) * 
                                            ((self.d**2 - 4 * self.d * self.tw + 4 * self.tw**2 / 3.0)
                                             * self.r**3 * self.tw))) 
                          + (((-27*self.d**2 *(math.pi - 4)**2)
                              + (54*(272/3.0 - 94 * math.pi/3.0 + math.pi**2) * self.d * self.tw)
                              + ((846 * math.pi - 2448 - 27 * math.pi**2) * self.tw**2)) * self.r**4)
                          + (12*(math.pi - 4)*(self.d**3 + 3*self.d**2 *self.tw 
                                               - 8 * self.d * self.tw**2 + 2 * self.tw**3) * self.r**2 * self.tw)))
            
            self.Iz = self.Iy
            
            # Elastic section modulus araund y-y or z-z (A.11)
            self.Zey = self.Iy / (self.d - _e)
            self.Zez = self.Zey
            
            # Radii of gyration around any axis (A.12)
            self.ry = math.sqrt(self.Iy / self.area)
            self.rz = self.ry
            
            # Plastic Properties
            # Distance v1p of the plastic centroid from heel (A.13)
            _v1p = ((3 * (math.pi - 4) * self.r**2 
                     + 4 * self.d * self.tw + 6 * self.tw**2) 
                    / (8 * math.sqrt(2) * self.tw))
            
            # Major plastic section modulus around u'-u' (A.14)
            self.Zpy = ((((48 * self.r**3 
                           + 3 * (math.pi - 4) * (self.d - self.tw) * self.r**2 
                           - 6 * self.d * self.tw**2 
                           + 6 * self.d**2 * self.tw 
                           + 2 * self.tw**3) * math.sqrt(2))
                         / 12.0) - (16 * self.r**3 / 3.0))
            
            # Minor Plastic section modulus around v'-v' (A.15)
            self.Zpz = ((math.sqrt(2) / (192 * self.tw)) 
                        * ((-27 * (math.pi - 4)**2 * self.r**4)
                           + (96 * (3 * math.pi - 10) * self.r**3 * self.tw) 
                           - ((12 * (math.pi - 4) * self.r**2) * ((2*self.d - 11*self.tw)*self.tw))
                           + (4 * self.tw**2 * (12 * self.d**2 
                                                - 12 * self.d * self.tw + 13 * self.tw**2))))

            _phi = 1.0
            
            #-------------------------------------------------
            #   Plastic Modulus about Minor Axis
            #   error, needs fix
            
            #_Zpy = 0
            #_Zpz = 0        
            #
        #
        else:
            
            _b1 = self.b - self.tw
            _h1 = self.d - self.tw
            
            #-------------------------------------------------
            #   Cross-Sectional Area
            self.area = (self.d + _b1)  * self.tw
            
            #-------------------------------------------------
            #   Elastic Neutral Centre
            self.Zc = (self.d**2 + _b1 * self.tw) / (2 * (self.d + _b1))
            self.Yc = (self.b**2 + _h1 * self.tw) / (2 * (self.b + _h1))
            
            #-------------------------------------------------
            #   Shear Centre 
            _SCz = self.tw / 2.0
            _SCy = self.tw / 2.0
            
            #-------------------------------------------------
            #               Section Properties
            #-------------------------------------------------
            
            #   Second Moment of Area about Mayor Axis
            
            self.Iy = ((self.tw*(self.d - self.Zc)**3 + self.b*self.Zc**3 -
                    _b1*(self.Zc - self.tw)**3) / 3.0)
            
            self.Zey = self.Iy / (self.d - self.Zc)
            
            _Zpy = 0
            
            #   Radius of gyration about Mayor Axis
            self.ry = math.sqrt(self.Iy / self.area)
            
            _SFy = 0
            
            #-------------------------------------------------
            #   Second Moment of Area about Minor Axis
            
            self.Iz = ((self.tw*(self.b - self.Yc)**3 + self.d*self.Yc**3 -
                    _h1*(self.Yc - self.tw)**3) / 3.0)
            
            self.Zez = self.Iz / (self.b - self.Yc)
            _Zpz = 0
            
            #   Radius of gyration about Minor Axis 
            self.rz = math.sqrt(self.Iz / self.area)
            
            _SFz = 0
            
            #-------------------------------------------------
            #   Product of inertia
            _Izy = (self.b*_b1*self.d*_h1*self.tw)/(4*(self.b + _h1))
            
            _Iu = (0.50 * (self.Iz + self.Iy) +
                   0.50 * math.sqrt((self.Iz - self.Iy)**2 + 4 * _Izy**2))
            
            _Iv = (0.50 * (self.Iz + self.Iy) -
                   0.50 * math.sqrt((self.Iz - self.Iy)**2 + 4 * _Izy**2))
            
            _phi = (math.atan(2 * _Izy / (self.Iy - self.Iz)))/2.0
            
            _Wu = 0
            _Wv = 0
            _Wup = 0
            _Wvp = 0
            _v1p = 0
        #
        _b1 = self.d - 0.50 * self.tw
        _b2 = self.b - 0.50 * self.tw
        _c1 = _b1 - 0.50 * self.tw
        _c2 = _b2 - 0.50 * self.tw
        #
        #-------------------------------------------------
        #   Warping Constant Cw (Bleich 1952, Picard and Beaulie 1991)
        self.Cw = (((_b1**3 + _b2**3) * (self.tw **3)) / 36.0)
        
        #-------------------------------------------------
        #   Torsional Constant (fillets neglected)
        self.J = (((_b1 + _b2) * (self.tw**3)) / 3.0)
        
        #-------------------------------------------------
        #   Product of inertia
        self.Iyz = (self.tw * (_b1 * self.Yc * (_b1 - 2 * self.Zc)) +
                    self.tw * (_b2 * self.Zc * (_b2 - 2 * self.Yc)))
        
        self.Jx = self.Iy + self.Iz
        
        self.rp = math.sqrt(self.Jx / self.area)
        #
        #return _Area, _Zc, _Yc, _Iy, _Zey, _Zpy, _ry, _Iz, _Zez, _Zpz, _rz
    #
    #    
    def print_file(self, file_name):

        check_out = print_header()       

        check_out.append("{:23s} {:>19} {:1.4E} {:1.4E} {:1.4E} {:1.4E}\n"
                         .format(self.type, "", self.d, self.tw, self.b, self.tw))     

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