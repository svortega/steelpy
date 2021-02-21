# 
# Copyright (c) 2015 steelpy


# Python stdlib imports
import math
import sys
#

# package imports
from steelpy.process.units.main import Units
from steelpy.f2uModel.material.main import Materials

#
#
#
#
# ----------------------------------------
#      Standard Sections Profiles
# ----------------------------------------
#    
class Tee:
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
    def __init__(self, cls):
        """
        """
        self.cls = cls
        self._units = Units()
        #_material = Materials()
        #_material[1] = 'plastic'
        #self._material = _material[1]
        #
        self.type = 'Symmetrical Tee Section'
        # Build [WELDED / ROLLED]
        self.build = 'welded'
        # Shear Stress [MAXIMUM / AVERAGE]
        self.shear_stress = 'average'
        self.compactness = 'N/A'        
        #self.units_in = ["", "", "second", "", "", ""]
        #    
    #
    #
    @property
    def material(self):
        """
        """
        return self._material
    
    @material.setter
    def material(self, value):
        """
        """
        self._material = value
    #
    def geometry(self, **kwargs):
        """
        """
        for key, value in kwargs.items():
            _dim = shape_io.find_section_dimensions(key)
            shape_io.get_dimension(self, _dim, value)
        
        # check data
        #try:
        #    self.b
        #except AttributeError:
        #    try:
        #        self.b = self.a    
        #    except AttributeError:
        #        self.b = self.d
        
        #try:
        #    self.tb
        #except AttributeError:
        #    try:
        #        self.tb = self.ta
        #    except AttributeError:
        #        self.tb = self.tw
        #
        #
        #
        
    #
    #
    @property
    def H(self):
        """
        d : Section Heigh
        """
        return self.d
    @H.setter
    def H(self, value):
        """
        """
        self.d = value
    #
    @property
    def t_web(self):
        """
        t_web: Web thickness
        """
        return self.tw
    @t_web.setter
    def t_web(self, value):
        """
        """
        self.tw = value
    #
    #
    @property
    def base(self):
        """
        b : Base
        """
        return self.b
    @base.setter
    def base(self, value):
        """
        """
        self.b = value
    #
    @property
    def t_flange(self):
        """
        t_flange: Flange thickness
        """
        return self.tb
    @t_flange.setter
    def t_flange(self, value):
        """
        """
        self.tb = value       
    #
    @property
    def properties(self):
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
        self.area = self.b * self.tb + self.tw * _D2
        
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        self.Zc = (((self.d**2 * self.tw) + (_C * self.tb**2)) /
                   (2 * (self.b*self.tb + _D2*self.tw)))
        
        self.Yc = 0
        
        #-------------------------------------------------
        #   Shear Centre 
        self.SCz = self.tb / 2.0
        
        self.SCy = 0
        
        #-------------------------------------------------
        #   Warping Constant Cw
        self.Cw = ((self.tb**3 * self.b**3 / 144.0) 
                   + (self.tw**3 * _h**3 / 36.0))
        
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        
        #   Second Moment of Area about Mayor Axis
        self.Iy = ((self.tw * (self.d - self.Zc)**3 + self.b * self.Zc**3 -
                    _C * (self.Zc - self.tb)**3) / 3.0)
        
        #   Elastic Modulus about Mayor Axis
        self.Zey = (self.Iy / (self.d - self.Zc))
        
        #   Plastic Modulus about Mayor Axis
        if self.tw * _D2 > self.b * self.tb :
            self.Zpy = ((_D2**2 * self.tw / 4.0) - 
                        (self.b**2 * self.tb**2 / (4.0 * self.tw)) + 
                        (self.b * self.tb * self.d / 2.0))
        
        else:
            self.Zpy = ((self.tb**2 * self.b / 4.0) +
                        0.50 * self.tw 
                        * _D2*(self.d - self.tw*_D2 / (2 * self.b)))
        
        #   Shape Factor
        self.SFy = (self.Zpy * (self.d - self.Zc) / self.Iy)
        
        #   Radius of gyration about Mayor Axis
        self.ry = (self.Iy / self.area)**0.50
        
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        self.Iz = (self.b**3 * self.tb + _D2 * self.tw**3) / 12.0
        
        #   Elastic Modulus about Minor Axis
        self.Zez = 2 * self.Iz / self.b
        
        #   Plastic Modulus about Minor Axis
        self.Zpz = (self.b**2 * self.tb + self.tw**2 * _D2) / 4.0
        
        #   Shape Factor
        self.SFz = self.Zpz * self.b / (2 * self.Iz)
        
        #   Radius of gyration about Minor Axis 
        self.rz = ( self.Iz / self.area)**0.50
        
        #-------------------------------------------------
        #   Torsional Constant
        self.J = (self.b * self.tb**3 + _h * self.tw**3) / 3.0
        
        #   Product of inertia
        _Iyz = 0
        
        self.Jx = self.Iy + self.Iz
        
        self.rp = (self.Jx / self.area)**0.50
        #
        
        #return _Area, _Zc, _Yc, _Iy, _Zey, _Zpy, _ry, _Iz, _Zez, _Zpz, _rz
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
    def print_file(self, file_name):

        check_out = print_header()       

        check_out.append("{:23s} {:>19} {:1.4E} {:1.4E} {:1.4E} {:1.4E}\n"
                         .format(self.type, "", self.d, self.tw, self.b, self.tb))
      

        check_out.extend(print_properties(self))

        #file_checkout = split_file_name(file_name)
        #file_checkout = str(file_checkout[0]) +'_check_me.txt'
        file_checkout = str(file_name) + '.txt'
        add_out = open(file_checkout,'w')
        add_out.write("".join(check_out))
        add_out.close()
        print('ok')
    #
    def set_default(self):
        """ """
        self.cls._default = self.name    
    #
#
#