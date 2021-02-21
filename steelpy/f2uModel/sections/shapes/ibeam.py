# 
# Copyright (c) 2020 steelpy
#

# Python stdlib imports
#from dataclasses import dataclass
import math
#import sys
from collections import namedtuple
#from typing import NamedTuple, Tuple, List # Iterator, Dict, Iterable
#

# package imports
from steelpy.process.units.main import Units
from steelpy.f2uModel.material.main import Materials
#import iLift.beam.section.process.io_module as shape_io
from steelpy.f2uModel.sections.process.stress import BeamStress
#from steelpy.process.load.actions import Actions

#
#
points = namedtuple('Points', ['y', 'z'])
#
def radial_shear_factor(D: float, Tw: float, Tft: float, 
                        Tfb: float, c: float, c1: float, 
                        R: float, e: float) -> float:
    '''
    Radial/horizontal web shear stress factor calculation
    Ref Roark 7 ed chapter 9 : Shear stress due to the radial 
                               shear force V
    
    tr : thickness of the section normal to the plane of curvature
         at the radial position r
    
    '''  
    # total area of the web 
    Aw = (D - Tft - Tfb) * Tw  
    #
    # Average Shear
    tau_average = 1.0/Aw
    print(' ')
    print('shear average : {: 1.4E}'.format(tau_average))
    #
    # --------------------------------------------
    # Shear stress due to the radial shear force V
    rn = R -e             # radial position to neautral axis
    b = R - c - Tft      # radial position to inner extreme fiber
    tr = Tw   # width of cross section where stresses are calculated
    # fraction of the web dep
    dwo = (c1 - Tft - e) / 3.0 
    dwi = (-c + Tfb + e) / 3.0
    CoorZ = [(c1 - Tft), 2 * dwo + e, dwo + e, e, 0,
              -e, dwi -e, 2*dwi -e, (-c + Tfb)]
    #
    tau_radial = []
    for i in range(len(CoorZ)):
        r = R + CoorZ[i]      # radial position to point of interest
        _cr =  (r - b) / 2.0  # distance from centroide to extreme fibre
        _r1 = _cr + b           # point of cross section from radial centre 
        _Ar = (tr * math.log(r /b) )* r  # area of part portion below r
        #_Ar = (tr * math.log((_r1 / _cr + 1) / (_r1 / _cr - 1))) * r
        _Qr = _Ar * _r1
        # roark Ed 7 equ (9.1-4)
        tau_radial.append( (rn / (tr*Aw*e*r**2)) * (R*_Ar - _Qr))
    #
    #print(' ')
    #for i in range(len(tau_radial)):
    #    print('tau : {: 1.4E} {: 1.4E}'
    #          .format(tau_radial[i], CoorZ[i]))
    #print(' ')
    print('shear radial  : {: 1.4E} '.format(max(tau_radial)))
    
    tau_y = max(tau_average, max(tau_radial))
    #print(' ')
    print('-----------------------------')
    print('Max Shear (tau) : {: 1.4E}'. format(tau_y))
    return tau_y
    #
    #-------------------------------------------------
#
#   
class Ibeam:
    """
    ============================================  
    Calculate the section properties of a I beam   
    ============================================   
         
         *  a  *
    +    +-----+  
            |  
    d       |         Z  
            |         ^  
    +  +---------+    + > Y  
       *    b    *  

    Parameters
    ----------
    d  : Section Height   
    tw : Web thickness   
    a  : Top compression flange base   
    ta : Top flange thickness   
    b  : Bottom tension flange base   
    tb : Bottom flange thickness   
    r  : root radious
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
    Sxc : Elastic section modulus referred to compression flange
    Sxt : Elastic section modulus referred to tension flange
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
        _material = Materials()
        _material[1] = 'plastic'
        self._material = _material[1]
        # Build [WELDED / ROLLED]
        self.build = 'welded'
        # Shear Stress [MAXIMUM / AVERAGE]
        self.shear_stress_type = 'average'
        self.compactness = None
        #
        self._properties = None
        # Shear factor
        self.FAvy = 1.0
        self.FAvz = 1.0
        #self.r = 0
        #
        self.root_radius = 0 * self._units.m
        self.type = 'I section'
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
    
    @property
    def section_mass(self):
        """
        section mass in g/m
        """
        return self.area * self._material.rho    
    #
    #
    #
    def geometry(self, **kwargs):
        """
        """
        for key, value in kwargs.items():
            _dim = shape_io.find_section_dimensions(key)
            shape_io.get_dimension(self, _dim, value)
        
        self.type = 'Symmetrical I section'
        #
        try:
            if self.top_flange_width != self.bottom_flange_width:
                self.type = 'Asymmetrical I section'
        except AttributeError:
            try:
                self.bottom_flange_width = self.top_flange_width
            except AttributeError:
                self.top_flange_width = self.bottom_flange_width
        #
        try:
            if self.top_flange_thickness != self.bottom_flange_thickness:
                self.type = 'Asymmetrical I section'
        except AttributeError:
            try:
                self.bottom_flange_thickness = self.top_flange_thickness
            except AttributeError:
                self.top_flange_thickness = self.bottom_flange_thickness
        #
        self._get_properties()
        self._get_section_coordinates()
    #   
    #
    def shear_stress(self, Vy, Vz,
                     stress_type:str ='average'):
        """
        Vy : horizontal force
        Vz : vertical force
        stress_type: average/true
        -------------------------
        tau_y : stress horizontal
        tau_z : stress vertical
        """
        #
        #-------------------------------------------------        
        #            Shear Stress Calculation
        # vertical section coordinates
        coord_z = self.section_coordinates.z
        #        
        if 'average' in stress_type.lower():
            # Area of Web
            # The overall depth times the web thickness
            self.Aw = self.height * self.web_thickness
            # Area of Flange
            self.Af = (self.bottom_flange_width * self.bottom_flange_thickness 
                       + self.top_flange_width * self.top_flange_thickness)
            #
            tau_z = [Vz / self.Aw for _ in coord_z] # vertical
            _index = [0, 2, 6, 8]
            for x in _index:
                tau_z[x] *= 0
            #
            tau_y = [Vy / self.Af for _ in coord_z] # horizontal
            _index = [4] # [0, 2, 4, 6, 8]
            for x in _index:
                tau_y[x] *= 0
        else:
            # True Shear Stress
            _Zcb = self.d - self.Zc
            # q5a = (((self.Zc - self.tft) / 2.0) * (self.tw * (self.Zc - self.tft/ 2.0))) #aldh6850
            # q5b = (((_Zcb - self.tfb) / 2.0) * (self.tw * (_Zcb - self.tfb/ 2.0))) #aldh6850
            q5a = (((self.Zc - self.tft) / 2.0) * (self.tw * (self.Zc - self.tft))) #aldh6850 - centroid of top half of web from Neutral-Axis times area of top half of web
            q5b = (((_Zcb - self.tfb) / 2.0) * (self.tw * (_Zcb - self.tfb))) #aldh6850 - centroid of bottom half of web from Neutral-Axis times area of bottom half of web
            #
            # q = [0 * self.Zc,  (self.Zc - self.tft/2.0) * self.tft * self.bft / 2.0, #aldh6850
            #      0 * self.Zc,  (self.Zc - self.tft/2.0) * self.tft * self.bft, #aldh6850
            #      None, (_Zcb - self.tfb/2.0) * self.tfb * self.bfb, 0 * self.bfb, #aldh6850
            #      (_Zcb - self.tfb/2.0) * self.tfb * self.bfb / 2.0, 0 * self.bfb] #aldh6850
            #
            # aldh6850 - section points 2 and 8 are on the extreme edges of the flange
            # where there is no shear stress
            q = [0 * (self.Zc)**3,  0 * (self.Zc)**3,
                 0 * (self.Zc)**3,  (self.Zc - self.tft/2.0) * self.tft * self.bft,
                 None, (_Zcb - self.tfb/2.0) * self.tfb * self.bfb, 0 * (_Zcb)**3,
                 0 * (_Zcb)**3, 0 * (_Zcb)**3]
            #aldh6850 ---------------------------------------------------------------------
            #
            q[4] = max(q5a + q[3], q5b + q[5])
            q = [_q * Vz / self.Iy for _q in q]
            # vertical
            # tau_z = [Vz / (2*self.Iy) * (self.d**2 / 4.0 - coord_z[0]**2), #aldh6850
            #          q[1]/self.bft,  #aldh6850
            #          Vz / (2*self.Iy) * (self.d**2 / 4.0 - coord_z[2]**2), #aldh6850
            #          q[3]/self.tw, q[4]/self.tw,  #aldh6850
            #          q[5]/self.tw, #aldh6850
            #          Vz / (2*self.Iy) * (self.d**2 / 4.0 - coord_z[6]**2), #aldh6850
            #          q[7]/self.bfb, #aldh6850
            #          Vz / (2*self.Iy) * (self.d**2 /4.0 - coord_z[8]**2)] #aldh6850
            #aldh6850 ---------------------------------------------------------------------
            tau_z = [q[0]/self.bft,
                     q[1]/self.bft, 
                     q[2]/self.bft,
                     q[3]/self.tw,
                     q[4]/self.tw, 
                     q[5]/self.tw,
                     q[6]/self.bfb,
                     q[7]/self.bfb,
                     q[8]/self.bfb]
            #aldh6850 ---------------------------------------------------------------------
            # get load proportion
            Vtop = Vy * (_Zcb / self.d)
            Vbot = Vy * (self.Zc / self.d)
            # get area flange
            bft = (2 * self.bft * self.tft)
            bfb = (2 * self.bfb * self.tfb)
            # horizontal
            tau_y = [0 * Vtop / bft,
                     3.0 * Vtop / bft, 
                     0 * Vtop / bft, 
                     3.0 * Vtop / bft, 
                     0 * Vtop / bfb, 
                     3.0 * Vbot / bfb,
                     0 * Vtop / bfb, 
                     3.0 * Vbot / bfb, 
                     0 * Vtop / bfb]
        # 
        return tau_y, tau_z
    #
    def torsional_stress(self, theta):
        """
        Roark's Torsion chapter
        """
        #if not G:
        G = self._material.E / (2 * (1.0 + self._material.poisson))
        #
        #if not theta:
        #    theta = self._get_rotation(To, E, G, l)
        #
        #
        if 'symmetrical' in self.type.lower():
            tau_1 = self.tf * G * theta[0]
            tau_2 = - self.ho * self.bf**2 * E * theta[2] / 16.0
            # bending
            sigma_y = self.ho * self.bf * E * theta[1] / 4.0
        else:
            # shear
            tau_1 = max(self.tw, self.tft, self.tfb) * G * theta[0]
            #
            if self.tfb * self.bfb > self.tft * self.bft:
                tau_2 = ((self.ho * self.tfb * self.bfb**3 * self.bft**2)
                         / (8 * self.tft * self.bft**3 + self.tfb * self.bfb**3))
            else:
                tau_2 = ((self.ho * self.tft * self.bft**3 * self.bfb**2)
                         / (8 * self.tft * self.bft**3 + self.tfb * self.bfb**3))                
            #
            tau_2 *= E * theta[2]
            #
            # bending 
            #
            e = (self.tft * self.bft**3 * self.ho
                 / (self.tft * self.bft**3 + self.tfb * self.bfb**3))
            #
            if self.tfb * self.bfb**2 > self.tft * self.bft**2:
                sigma_y = ((self.ho * self.bft / 2.0) 
                           * self.tfb * self.bfb**3 
                           / (self.tft * self.bft**3 
                              + self.tfb * self.bfb**3))
            else:
                sigma_y = ((self.ho * self.bfb / 2.0) 
                           * self.tft * self.bft**3 
                           / (self.tft * self.bft**3 
                              + self.tfb * self.bfb**3))
            #
            sigma_y *= E * theta[1]
        #
        # bending stress
        self.sigma_y = [_sigma + math.copysign(sigma_y, _sigma)
                        if _sigma != 0 else 0 for _sigma in self.sigma_y]
        #
        # shear stress
        self.tau_z = [_tau + math.copysign(tau_1, _tau) - math.copysign(tau_2, _tau)
                      if _tau != 0 else 0 for _tau in self.tau_z ]
        #
        #print('ok')
    #
    #
    #
    def rolled(self):
        """
        """
        self.build='rolled'
    #
    def curved(self, R:float):
        """
        ---------
        R = Radio
        """
        if 'symmetrical' in self.type.lower():
            b = self.bottom_flange_width
            b1 = self.web_thickness
            t = self.bottom_flange_thickness
            d = self.height
            ry = self.ry
        
            # shear area
            warea = self.area
        
            # extreme fibre distances c
            c = d/2.0
            self.c = c
        
            c1 = d - c
            self.c1 = c1
        
            # centroidal radius
            #R = R
            # R = orad - c1
            self.R = R
        
            # Shift of neutral axis from neutral axis
            e = (c *((R/c)- ((2.0*(t/c + (1 - t/c)*(b1/b))) 
                                 / ((math.log(((R/c)**2 + (R/c + 1)*(t/c) - 1.0) 
                                              / ((R/c)**2 - (R/c - 1.0)*(t/c) - 1.0))) 
                                    + ((b1/b)*math.log((R/c - t/c + 1.0) 
                                                         /(R/c + t/c - 1.0)))))))
            self.e = e
            # where
            Ic = self.Iy
        
            # stress factors Ki
            self.ki = ((Ic / (warea * c**2 * (R/c - 1.0))) 
                       * ((1.0 - e / c) / (e / c)))
        
            # stress factors Ko
            self.ko = ((Ic / (warea * c**2 * (R/c + 1.0))) 
                       * ((1.0 + e / c) / (e / c)))
        
            # Modulus of rigidity factor (section 8.10)
            nai = c - e    # neautral axis inner fiber
            nao = c1 + e   # neautral axis outer fiber
        
            D1 = nai - t
            D2 = nai 
            t1 = b1
            t2 = b
            r = ry
        
            self.F = ((1 + (((3*(D2**2 - D1**2)*D1)/(2.0*D2**3)) * (t2/t1 - 1.0)))
                      * (4*D2**2 / (10*r**2)))
            #
            # Shear factor (section 8.1 equ 8.1-13)
            self.tau_y = radial_shear_factor(d, b1, t, t, c, c1, R, e)
        
        else:
            b = self.bottom_flange_width
            t = self.bottom_flange_thickness
            b1 = self.top_flange_width
            t1 = self.top_flange_thickness
            d = self.height
            b2 = self.web_thickness
        
            # shear area
            warea = d * b2
        
            # extreme inner fibre distance c
            c = (d * (((b1/b - b2/b)*(2.0 - t1/d)*(t1/d) 
                         + (1.0 - b2/b)*(t/d)**2 + (b2/b))
                        / (2*self.area / (b*d))))
            self.c = c
            # extreme outer fibre distance c
            c1 = d - c
            self.c1 = c1
            # centroidal radius
            #_R = R
            #R = R - c1
            self.R = R
        
            # Shift of neutral axis from neutral axis
            e = (c * ((R/c)-(((self.area/(b*d))*(d/c)) 
                                 / (math.log((R/c + t/c - 1)/(R/c - 1)) 
                                    + ((b2/b)*math.log((R/c + c1/c - t1/c )
                                                         / (R/c + t/c - 1)))
                                    + ((b1/b)*math.log((R/c + c1/c)
                                                         / (R/c + c1/c - t1/c)))))))
            self.e = e
            # where
            Ic = self.Iy
        
            # stress factors Ki
            self.ki = ((Ic / (self.area * c**2 * (R/c - 1.0))) 
                       * ((1.0 - e / c) / (e / c)))
        
            # stress factors Ko
            self.ko = ((Ic / (self.area * c**2 * (e/c ))) 
                       * ((d/c + e/c - 1.0) / (R/c  + d/c - 1.0))
                       * (1.0  / (d / c - 1.0)))
        
            # Modulus of rigidity factor (section 8.10)
            nai = c - e    # neautral axis inner fiber
            nao = c1 + e   # neautral axis outer fiber
        
            if nai <= nao:
                D1 = nai - t1
                D2 = nai 
                t1 = b2
                t2 = b1
                r = self.ry
                print ('inner fiber ', nai)
            else:
                D1 = nao - t
                D2 = nao 
                t1 = b2
                t2 = b
                r = self.ry
                print ('outer fiber ', nao)
            
            
            self.F = ((1 + (((3*(D2**2 - D1**2)*D1)/(2.0*D2**3)) * (t2/t1 - 1.0)))
                      * (4*D2**2 / (10*r**2)))
            
            #
            # Shear factor (section 8.1 equ 8.1-13)
            #_shearFactor = shear_factor(c, c1, R, e)
            #
            #
            #R, _F, e, c, c1, _ki, _ko, _shearFactor
            self.tau_y = radial_shear_factor(d, b1, t, t1, c, c1, R, e)
        #
        #
    #
    def print_file(self, file_name):
        """
        """
        check_out = print_header()       

        check_out.append("{:23s} {:>19} {:1.4E} {:1.4E} {:1.4E} {:1.4E}\n"
                         .format(self.type, "", self.height, self.web_thickness, self.top_flange_width, self.top_flange_thickness))
        
        check_out.append("{:>65} {:1.4E} {:1.4E}\n"
                         .format("", self.bottom_flange_width, self.bottom_flange_thickness))        

        check_out.extend(print_properties(self))

        #file_checkout = split_file_name(file_name)
        #file_checkout = str(file_checkout[0]) +'_check_me.txt'
        file_checkout = str(file_name) + '.txt'
        add_out = open(file_checkout,'w')
        add_out.write("".join(check_out))
        add_out.close()
        print('ok')
    #
    def get_compactness(self, material: str):
        """
        """
        _class_B41a, _class_B41b = open_section_compactness(self, material)
    
        return _class_B41a, _class_B41b
    #
    #@property
    def _get_properties(self):
        """
        """
        #
        self.type = 'Symmetrical I section'
        #
        try:
            if self.top_flange_width != self.bottom_flange_width:
                self.type = 'Asymmetrical I section'
        except AttributeError:
            try:
                self.bottom_flange_width = self.top_flange_width
            except AttributeError:
                self.top_flange_width = self.bottom_flange_width
        #
        try:
            if self.top_flange_thickness != self.bottom_flange_thickness:
                self.type = 'Asymmetrical I section'
        except AttributeError:
            try:
                self.bottom_flange_thickness = self.top_flange_thickness
            except AttributeError:
                self.top_flange_thickness = self.bottom_flange_thickness        
        #
        #
        #-------------------------------------------------   
        #
        _hw = (self.height - self.top_flange_thickness - self.bottom_flange_thickness) # - 2 * self.root_radius
        _ho = (self.height - 0.5 * self.top_flange_thickness - 0.5 * self.bottom_flange_thickness)
        self.ho = _ho
        self.hw = _hw
        #-------------------------------------------------
        #   Cross-Sectional Area
        self.area = (self.top_flange_width*self.top_flange_thickness 
                     + self.bottom_flange_width*self.bottom_flange_thickness 
                     + _hw*self.web_thickness)
        
        #-------------------------------------------------
        #   Elastic Neutral Centre 
        self.Zc = ((self.top_flange_width * self.top_flange_thickness**2 / 2.0 
                    + self.bottom_flange_width * self.bottom_flange_thickness 
                    * (_hw + self.top_flange_thickness + self.bottom_flange_thickness / 2.0) 
                    + _hw * self.web_thickness * (_hw / 2.0 + self.top_flange_thickness)) 
                   / (self.top_flange_width * self.top_flange_thickness 
                    + _hw*self.web_thickness 
                    + self.bottom_flange_width * self.bottom_flange_thickness))
        
        self.Yc = 0 * self.Zc
        
        #   Plastic Neutral Centre    # @hami2230 - added
        
        if (self.bottom_flange_width * self.bottom_flange_thickness >
                (self.top_flange_width * self.top_flange_thickness
                + self.hw * self.web_thickness)):
            self.Zp = (self.height 
                       - (0.5 * self.area / self.bottom_flange_width)
                       - self.top_flange_thickness)
        elif (self.top_flange_width * self.top_flange_thickness > 
                 (self.bottom_flange_width * self.bottom_flange_thickness 
                 + self.hw * self.web_thickness)):
            self.Zp = (self.height - 
                      (self.bottom_flange_thickness + self.hw 
                      + ((0.5 * self.area - self.bottom_flange_width 
                      * self.bottom_flange_thickness - self.hw * self.web_thickness)
                      / self.top_flange_width)) - self.top_flange_thickness) 
        else:
            self.Zp = (self.height - 
                      ((0.5 * self.area - self.bottom_flange_width 
                      * self.bottom_flange_thickness)
                      / self.web_thickness + self.bottom_flange_thickness) 
                      - self.top_flange_thickness)                           
       
        # Warping Constant Cw
        # Picard and Beaulieu 1991
        d = self.height - (self.top_flange_thickness + self.bottom_flange_thickness) / 2.0
        
        _alpha = (1.0 / (1 + (self.top_flange_width / self.bottom_flange_width)**3 
                         * (self.top_flange_thickness / self.bottom_flange_thickness)))
        
        self.Cw = (d**2 * self.top_flange_width**3 * self.top_flange_thickness * _alpha) / 12.0
        
        #-------------------------------------------------
        #   Torsional Constant
        if self.top_flange_width == self.bottom_flange_width and self.top_flange_thickness == self.bottom_flange_thickness :
            self.J = ((2 * self.top_flange_width * self.top_flange_thickness**3 / 3.0) 
                      + (self.height * self.web_thickness**3 / 3.0))
            #
            self.K = ((2 * self.tf**3 * self.bf +
                       self.tw**3 * self.ho) / 3.0)
        else:
            self.J = ((self.top_flange_width * self.top_flange_thickness**3 
                       + self.bottom_flange_width * self.bottom_flange_thickness**3 
                       + d * self.web_thickness**3) / 3.0)
            #
            self.K = (self.tft**3 * self.bft + self.tfb**3 * self.bfb 
                      + self.tw**3 * self.ho)            
        #
        #-------------------------------------------------
        #   Shear Centre
        self.SCz = (((self.height - self.top_flange_thickness / 2.0 - self.bottom_flange_thickness / 2.0) 
                     * self.top_flange_thickness * self.top_flange_width**3) 
                    / (self.top_flange_thickness * self.top_flange_width**3 
                       + self.bottom_flange_thickness * self.bottom_flange_width**3))
        
        self.SCy = 0 * self.SCz
        
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        
        #   Second Moment of Area about Mayor Axis
        self.Iy = (self.top_flange_width * self.top_flange_thickness**3 / 12.0 
                   + (self.top_flange_width * self.top_flange_thickness 
                      * (self.Zc - self.top_flange_thickness / 2.0)**2 )
                   + self.bottom_flange_width * self.bottom_flange_thickness**3 / 12.0 
                   + (self.bottom_flange_width * self.bottom_flange_thickness 
                      * (_hw + self.bottom_flange_thickness / 2.0 + self.top_flange_thickness - self.Zc)**2) 
                   + self.web_thickness * _hw**3 / 12.0 
                   + self.web_thickness * _hw * (_hw / 2.0 + self.top_flange_thickness - self.Zc)**2)
                      
        #   Second Moment of Area of Compression Flange about Mayor Axis
        
        self._Iy_ft = (self.top_flange_width * self.top_flange_thickness**3 / 12.0 + # @hami2230 - added
                  (self.top_flange_width * self.top_flange_thickness
                   * (self.Zc - self.top_flange_thickness / 2.0)**2 ))
                  
        #   Elastic Modulus about Mayor Axis
        if self.Zc >= (self.height - self.Zc):
            self.Zey = self.Iy / self.Zc
        else:
            self.Zey = self.Iy / (self.height - self.Zc)
        
        #   Plastic Modulus about Mayor Axis
        self.Zpy = ((self.web_thickness * _hw**2 / 4.0) 
                    + (self.top_flange_width * self.top_flange_thickness 
                       * (self.Zc - self.top_flange_thickness / 2.0)) 
                    + (self.bottom_flange_width * self.bottom_flange_thickness 
                       * (self.height - self.Zc - self.bottom_flange_thickness / 2.0)))
        
        #   Radius of gyration about Mayor Axis
        self.ry = (self.Iy / self.area)**0.5
        
        self.SFy = self.Zpy / self.Zey
        
        # Elastic modulus of compression flange about major axis               # @hami2230 - Sxc and Sxt added
        
        self.Sxc = (self.top_flange_width * self.top_flange_thickness**3 / 12.0 
                   + (self.top_flange_width * self.top_flange_thickness 
                      * (self.Zc - self.top_flange_thickness / 2.0)**2 )) / self.Zc

        # Elastic modulus of tension flange about major axis                    

        self.Sxt = (self.bottom_flange_width * self.bottom_flange_thickness**3 / 12.0 
                   + (self.bottom_flange_width * self.bottom_flange_thickness 
                      * ((self.height-self.Zc) - self.bottom_flange_thickness / 2.0)**2 )) / (self.height-self.Zc)
        #
        #-------------------------------------------------
        #   Second Moment of Area about Minor Axis
        self.Iz = (self.top_flange_width**3 * self.top_flange_thickness / 12.0 +
                   self.bottom_flange_width**3 * self.bottom_flange_thickness / 12.0 +
                   self.web_thickness**3 * _hw / 12.0)
        
        #   Elastic Modulus about Minor Axis
        if self.top_flange_width >= self.bottom_flange_width:
            self.Zez = 2 * self.Iz / self.top_flange_width
        else:
            self.Zez = 2 * self.Iz / self.bottom_flange_width
        
        #   Plastic Modulus about Minor Axis  
        self.Zpz = ((self.top_flange_thickness * self.top_flange_width**2 
                     + self.bottom_flange_thickness * self.bottom_flange_width**2 
                     + _hw * self.web_thickness**2) / 4.0)
        
        #   Radius of gyration about Minor Axis  
        self.rz = (self.Iz / self.area)**0.5
        
        self.SFz = self.Zpz / self.Zez
        
        #-------------------------------------------------
        #   Product of inertia
        _Iyz = 0
        self.Jx = self.Iy + self.Iz
        self.rp = (self.Jx / self.area)**0.50
        
        # warping statical moment at point s on cross-section
        #self.Sws = self.ho * self.bf**2 * self.tw / 16.0
        self.Sws = self.ho * self.bf**2 * self.tf / 16.0 #aldh6850
        # normalized warping function at point s on cross-section
        self.Wns = self.ho * self.bf / 4.0
        #
        #-------------------------------------------------
        self._get_section_coordinates()
        #
        #file = shape_io.print_header()
        #file.extend(self._shape())
        #file.extend(shape_io.print_properties(self))
        #for row in file:
        #    print(row.rstrip())
        #print('ok')
        #file_control = 'testXXX.txt'
        #_control = open(file_control, 'w+')
        #_control.write("".join(file))
        #_control.close()
        #return _Area, _Zc, _Yc, _Iy, _Zey, _Zpy, ry, _Iz, _Zez, _Zpz, _rz
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
    # Double symetrical section
    #
    @property
    def d(self):
        """
        d : section height
        """
        return self.height
    @d.setter
    def d(self, value):
        """
        d : section height
        """
        self.height = value
    #
    @property
    def tw(self):
        """
        tw : section web_thickness 
        """
        return self.web_thickness
    @tw.setter
    def tw(self, value):
        """
        tw : section web_thickness
        """
        self.web_thickness = value
    #
    @property
    def bf(self):
        """
        bft : section top flange width
        """
        #if self.top_flange_width.value < self.bottom_flange_width.value:
        #    return self.top_flange_width
        #return self.bottom_flange_width
        return min(self.top_flange_width, self.bottom_flange_width)
    @bf.setter
    def bf(self, value):
        """
        Assuming double symmetrical section
        bft : section top flange width
        """
        self.top_flange_width = value
        self.bottom_flange_width = value
    #
    @property
    def tf(self):
        """
        tft : section top flange thickness
        """
        #if self.top_flange_thickness.value < self.bottom_flange_thickness.value:
        #    return self.top_flange_thickness
        #return self.bottom_flange_thickness
        return min(self.top_flange_thickness, self.bottom_flange_thickness)
    @tf.setter
    def tf(self, value):
        """
        Assuming double symmetrical section
        tf : section top flange thickness
        """
        self.top_flange_thickness = value
        self.bottom_flange_thickness = value
    #
    # Single symetrical section
    #
    @property
    def bft(self):
        """
        bft : section top flange width
        """
        return self.top_flange_width
    @bft.setter
    def bft(self, value):
        """
        bft : section top flange width
        """
        self.top_flange_width = value
    #
    @property
    def tft(self):
        """
        tft : section top flange thickness
        """
        return self.top_flange_thickness
    @tft.setter
    def tft(self, value):
        """
        tft : section top flange thickness
        """
        self.top_flange_thickness = value
    #
    #
    @property
    def bfb(self):
        """
        bfb : section bottom flange width
        """
        return self.bottom_flange_width
    @bfb.setter
    def bfb(self, value):
        """
        bfb : section bottom flange width
        """
        self.bottom_flange_width = value
    #
    @property
    def tfb(self):
        """
        tfb : section bottom flange thickness
        """
        return self.bottom_flange_thickness
    @tfb.setter
    def tfb(self, value):
        """
        tfb : section bottom flange thickness
        """
        self.bottom_flange_thickness = value
    #
    @property
    def r(self):
        """
        fillet radious
        """
        return self.root_radius
    @r.setter
    def r(self, value):
        """
        fillet radious
        """
        self.root_radius = value
    #
    #
    def _shape(self):
        """
        """
        _section = []
        if 'symmetrical' in self.type.lower():
            _section.append("{:2s}+   bf    +{:33s}{:1.3E} {:1.3E}  {:1.3E} {:1.3E}\n"
                            .format("", "", 
                                    self.d.convert('millimetre').value, 
                                    self.tw.convert('millimetre').value, 
                                    self.bf.convert('millimetre').value, 
                                    self.tf.convert('millimetre').value))
            _section.append("+ +====+====+ tf\n")
            _section.append("{:}| tw\n".format(7*" ")) 
            _section.append("d{:6s}|{:6s}Z\n".format("", ""))
            _section.append("{:7s}|{:6s}^\n".format("", ""))
            _section.append("+ +====+====+ + > Y\n")
            _section.append("{:2s}+   bf    +\n".format(""))            
        else:
            _section.append("{:4s}+ bft +{:34s}{:1.3E} {:1.3E}  {:1.3E} {:1.3E}\n"
                            .format("","", 
                                    self.d.convert('millimetre').value, 
                                    self.tw.convert('millimetre').value, 
                                    self.bft.convert('millimetre').value, 
                                    self.tft.convert('millimetre').value))
            _section.append("+   +==+==+{:55s}{:1.3E} {:1.3E}\n"
                            .format("", 
                                    self.bfb.convert('millimetre').value, 
                                    self.tfb.convert('millimetre').value))
            _section.append("{:7s}|\n".format("")) 
            _section.append("d{:6s}|{:6s}Z\n".format("", ""))
            _section.append("{:7s}|{:6s}^\n".format("", ""))
            _section.append("+ +====+====+ + > Y\n")
            _section.append("{:2s}+   bfb   +\n".format(""))
        
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
        stress points = [1 2 3 4 5 6 7 8 9]
        """
        # get section's coordinates
        coord_y = self.section_coordinates.y # lateral
        coord_z = self.section_coordinates.z # vertical
        # get shear stress
        tau_y, tau_z = self.shear_stress(actions.Fy, actions.Fz, 
                                           stress_type=self.shear_stress_type)
        # FIXME: don't know what to do here
        tau_x = [tau_y[x] * 0 for x in range(len(tau_y))]
        # get bending stress
        sigma_x = [(actions.Fx / self.area) for _ in coord_y]
        sigma_y = [(actions.My * _coord / self.Iy) for _coord in coord_z]
        sigma_z = [(actions.Mz * _coord / self.Iz) for _coord in coord_y]
        #
        if stress:
            if isinstance(stress.tau_x, list):
                # assuming section stress already calculated
                #print('---> list')
                stress.tau_x = self._combine_stress(tau_x, stress.tau_x)
                stress.tau_y = self._combine_stress(tau_y, stress.tau_y)
                stress.tau_z = self._combine_stress(tau_z, stress.tau_z)
                #
                stress.sigma_x = self._combine_stress(sigma_x, stress.sigma_x)
                stress.sigma_y = self._combine_stress(sigma_y, stress.sigma_y)
                stress.sigma_z = self._combine_stress(sigma_z, stress.sigma_z)
            else:
                # Assuming global stress
                stress_tau_x = [stress.tau_x for x in range(9)]
                stress.tau_x = self._combine_stress(tau_x, stress_tau_x)
                #
                stress_tau_y = [stress.tau_y for x in range(9)]
                _index = [4] # 3, 4, 5
                for x in _index:
                    stress_tau_y[x] *= 0
                stress.tau_y = self._combine_stress(tau_y, stress_tau_y)
                #
                stress_tau_z = [stress.tau_z for x in range(9)]
                _index = [0, 2, 6, 8]
                for x in _index:
                    stress_tau_z[x] *= 0
                stress.tau_z = self._combine_stress(tau_z, stress_tau_z)
                #
                stress_sigma_x = [stress.sigma_x for _ in coord_y]
                stress.sigma_x = self._combine_stress(sigma_x, stress_sigma_x)
                #
                _factor_z = [1, 1, 1, 1, 0, -1, -1, -1, -1]
                stress_sigma_y = [stress.sigma_y * _coord for _coord in _factor_z]
                stress.sigma_y = self._combine_stress(sigma_y, stress_sigma_y)
                #
                _factor_y = [1, 0, -1, 0, 0, 0, 1, 0, -1]
                stress_sigma_z = [stress.sigma_z * _coord  for _coord in _factor_y]
                stress.sigma_z = self._combine_stress(sigma_z, stress_sigma_z)
                
                # print(stress.sigma_x[x].convert('megapascal').value) #aldh6850
                # print(stress.sigma_y[x].convert('megapascal').value) #aldh6850
                # print(stress.sigma_z[x].convert('megapascal').value) #aldh6850
                # print(stress.tau_x[x].convert('megapascal').value) #aldh6850
                # print(stress.tau_y[x].convert('megapascal').value) #aldh6850
                # print(stress.tau_z[x].convert('megapascal').value) #aldh6850
        else:
            stress = BeamStress(sigma_x, sigma_y, sigma_z, 
                                tau_x, tau_y, tau_z)
        #
        return stress
    #
    def _add_global_stress(self, stress_local, stress_global):
        """
        """  
        # _new_stress = [ _item + math.copysign(1, _item.value) * stress_global  
        #                 if _item.value != 0  else _item for _item in stress_local] #aldh6850
        
        #aldh6850 - update to ensure the "global" stress has the same sign as the "local" stress to be conservative
        #aldh6850 - update to ensure when the "local" stress is zero the "global" stress is used
        
        _new_stress = [ _item + math.copysign(1, _item.value) * abs(stress_global)  
                        if _item.value != 0  else stress_global for _item in stress_local] #aldh6850
        
        
        return _new_stress
    #
    def _combine_stress(self, stress_1, stress_2):
        """
        """
        # change * by +
        _new_stress = [stress_1[x] + math.copysign(1, stress_1[x].value) * abs(stress_2[x]) 
                       for x in range(9)]
        return _new_stress
    #
    def _get_section_coordinates(self):
        """
        1    2     3
        +----+-----+
        |____+_____| 4    ^ z
             |            |
             + 5          +--> y
             |
         ____+_____  6
        |          |
        +----+-----+      
        7    8     9
        """
        # horizontal
        coord_y = [-1 * self.bft/2.0, 0 * self.bft, self.bft/2.0, 
                   0 * self.bft, 0 * self.bft, 0 * self.bfb, 
                   -1 * self.bfb/2.0, 0 * self.bfb, self.bfb/2.0]
        # vertical
        _Zcb = self.Zc - self.d
        coord_z = [self.Zc, self.Zc, self.Zc, 
                   self.Zc - self.tft, 0 * _Zcb , _Zcb + self.tfb,
                   _Zcb, _Zcb, _Zcb]
        
        self.section_coordinates = points(coord_y, coord_z)
        #print('ok')
    #
    #
    def set_default(self):
        """ """
        self.cls._default = self.name    
#
#
#