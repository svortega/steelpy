# Copyright (c) 2015-2020 steelpy
#

# Python stdlib imports

# package imports
from steelpy.f2uModel.material.main import Materials
import steelpy.roarks.formulas.circular_ring as ring
from steelpy.process.units.main import Units
from steelpy.f2uModel.sections.main import Sections
#from steelpy.sections.tee import Tee
from steelpy.f2uModel.sections.shapes.solid import Rectangle
from steelpy.codes.api.design import API_design


#-------------------------------------------------
#                Supporting Section
#-------------------------------------------------
#
class RoarkFormulas:
    
    __slots__ = ('case','ring_load', 'phase', 'theta', 'phi',
                 'istep', 'xloc', 'xmom', 'xhop', 'xshr')
    #
    def __init__(self, case, load, phase = 0., theta = 0., phi = 0., step = 32):
        """
        """
        self.case = case
        self.ring_load = load	
        self.phase  = phase
        self.theta  = theta
        self.phi  = phi
        #
        self.istep = step
        self.xloc = []
        self.xmom = []
        self.xhop = []
        self.xshr = []
        #
        for i in range(self.istep):
            self.xloc.append(0)
            self.xmom.append(0)
            self.xhop.append(0)
            self.xshr.append(0)
#
#
#
#
class  CircularRing:
    
    def __init__(self, istep:int = 24):
        """
        """
        self._units = Units()
        _material = Materials()
        _material[1] = 'plastic'
        self._material = _material[1]
        # Number of points around ring
        self.istep  = istep 
        #
        self.RingCase = {}
        #
        self.iload  = 0
        #
        #
        self.xloc = []
        self.xm = []
        self.xh = []
        self.xs = []
        self.bstri = []
        self.bstro = []
        self.hstr = []
        self.sstr = []
        self.cstri = []
        self.cstro = []
        self.cstmax = 0
        self.shmax = 0
        #
        # Outside radius of ring
        self.orad = 0 * self._units.mm
        # inside flange
        self.fiw = 0 * self._units.mm
        self.fit = 0 * self._units.mm
        # outside flange
        self.fow = 0 * self._units.mm 
        self.fot = 0 * self._units.mm
        # web
        self.wlh = 0 * self._units.mm
        self.wth = 0 * self._units.mm
        #
        self.LocalAxis = 'x y z'
        #self.MaterialType = 'steel'
        #
        self._chord = API_design()
        self._section = Sections()
        self._stiffness_flag = False
    #
    @property
    def units(self):
        """
        """
        return self._units
        
    #     
    @property
    def material(self):
        """
        """
        return self._material
    
    #@material.setter
    #def material(self, value):
    #    """
    #    """
    #    self._material = value
    #
    @property
    def geometry(self):
        """
        """
        return self._chord.section
    
    @property
    def chord(self):
        """
        """
        return self._chord
    #
    def case(self, case, load, phase, theta, phi):
        """
        """
        self.iload += 1
        #print (self.iload)
        #
        self.RingCase[self.iload] = RoarkFormulas(case, load, phase, 
                                                  theta, phi, self.istep)
    #
    #
    def input_radius(self, orad):
        """
        """
        #
        #self.SectionType = SectionType
        self.orad = orad
    #
    #    
    @property
    def stiffener(self):
        """
        """
        self._stiffness_flag = True
        return self._section
    #
    def print_results(self, PrintOption = 'MONITOR'):
        """
        """
        #
        #---------------------
        for i in range(self.istep):
            self.xm.append(0)
            self.xh.append(0)
            self.xs.append(0)
            self.bstri.append(0)
            self.bstro.append(0)
            self.hstr.append(0)
            self.sstr.append(0)
            self.cstri.append(0)
            self.cstro.append(0)
        #
        # Open report.txt For Append As #1
        #
        print (" ")
        print ("        STRESSES IN CIRCULAR RINGS 5th edition")
        print ("                R.J.ROARK & W.C.YOUNG")
        print ("******************************************************")
        print (" ")
        print ('RINGS DATA ECHO')
        print ('number of load cases: {:}'.format(self.iload))
        print ('number of points around ring: {:}'.format(self.istep))
        print ("Cross Section")
        print ('outside radius of ring: {:}'.format(self.orad))
        print ('outside flange: {:} x {:} thk'.format(self.fow,self.fot))
        print ('inside flange:  {:} x {:} thk'.format(self.fiw,self.fit))
        print ('web             {:} x {:} thk'.format(self.wlh,self.wth))
        print ("")
        print ('CALCULATED SECTION PROPERTIES AND FACTORS')
        # sectional properties of curved beam are calculated
        # and output to file
        if self._stiffness_flag :
            # _section.height = 1.1 * (self.geometry.D * self.geometry.t)**0.50
            #_D = _D + _Tft + _Tfb
            self._section.curved(self.orad)
            #(self.section.area, self.section.R, self.section.Iy, self.section.F, 
            # self.section.e, _c, _c1, _ki, _ko, _shearFactor) = ring.cibeam2(self.SectionType, self.orad,
            #                                                     _D, _Tw, _Bft, _Tft, _Bfb, _Tfb)
            #
            print ('section area :',self._section.area,'mm**2')
            print ('second moment of area @ centroid :',_section.Iy,'mm**4')
            print ('centroidal radius ',self._section.R)
            #print ('extreme fibre distances   ci: ',_c,',  co: ', _c1)
            print ('stress factors            ki: ', self._section.ki,',  ko: ',self._section.ko)
            print(' ')
            #
            # constant factors calculated and written to file
            
            _k1, _k2 = ring.factors2(self.E, self.G, self.nu, 
                                     _section.area, _section.R, 
                                     _section.Iy, _section.F, _section.e,
                                     _section.d)
        else:
            R = self.geometry.D
            _k1, _k2 = 1, 1
            _section = Rectangle()
            _section.d = 1.0 * self.units.m
            _section.t = self.geometry.t
            _section._get_properties()
            _section.curved(R)
        #
        # start of main subroutineme loop - each load case is called 
        # in turn control is passed to roark routines 
        #  insert extra roark cases here
        #
        #print('-->')
        for _ring in self.RingCase.values():
            _xloc, _xm, _xh, _xs = ring.circular_ring(_ring.case, _ring.ring_load.value,
                                                      _ring.theta.value, _ring.phi.value, 
                                                      _ring.phase.value, self.istep,
                                                      R.value, _k1, _k2,)
            #
            print ('LOCATION    MOMENT      THRUST      SHEAR')
            print ('             N.mm         N           N')
            for i in range (len(_xloc)):
                print("{:>5.1f}    {: 1.4e} {: 1.4e} {: 1.4e}"
                      .format(_xloc[i], _xm[i], _xh[i], _xs[i]))
                #
                self.xm[i] += _xm[i]
                self.xh[i] += _xh[i]
                self.xs[i] += _xs[i]
            print ('')
        #
        _tau_y = max([abs(_item.value) for _item in _section.tau_y])
        # stresses are calculated for superposition of 
        (self.cstmax, self.shmax, self.bstri, self.bstro,
         self.hstr, self.sstr, self.cstri, self.cstro) = ring.stress2(_section.Iy.value, _section.area.value, 
                                                                      _section.c.value, _section.c1.value, 
                                                                      _section.ki.value, _section.ko.value, 
                                                                      _tau_y, self.istep, 
                                                                      _xloc, self.xm, self.xh, self.xs)
        # SUPERPOSITION OF ALL LOAD CASES
        print ("")
        print ('                                                              SUPERPOSITION OF ALL LOAD CASES')
        print ('          -------------Forces-------------    -----------------------Stresses (n/mm^2)---------------------------')
        print ('Location  Moment     Thrust        Shear     Bending     Bending      Hoop        Shear       Bend.in.   Bend.out',)
        print ('Degrees    N.mm        N             N        Inner       Outer                               +Hoop       +Hoop',)
        for i in range(len(_xloc)):
            print ("{:>5.1f}  {: 1.4e} {: 1.4e} {: 1.4e} {: 1.4e} {: 1.4e} {: 1.4e} {: 1.4e} {: 1.4e} {: 1.4e}"
                   .format(_xloc[i], self.xm[i], self.xh[i], self.xs[i], self.bstri[i],
                           self.bstro[i], self.hstr[i],
                           self.sstr[i], self.cstri[i], self.cstro[i]))
        
        print ("")
        print ('max. combined stress    max. shear stress')
        print ('      N/mm**2                N/mm**2')
        print ("   {: 1.3f}            {: 1.3f}".format(self.cstmax, self.shmax))
        
        lo = min(self.xs)
        hi = max(self.xs)
        Vy = max(abs(lo), abs(hi))
        print ('')
        print  ("results written to specified filename")
#
#
#
if __name__=="__main__": 
    #
    # Input data etc
    #
    #
    print (" ")
    print ("STRESSES IN CIRCULAR RINGS-AFTER R.J.ROARK & W.C.YOUNG")
    print ("******************************************************")
    #
    #
    #filenm = input(" Enter output file name: ")
    #
    #open (7, file=filenm, status="new")
    #
    #ititle = input(" Enter job title: ")
    #
    #
    print (" ")
    print('DATA ENTRY')
    #
    iload = int(input('Number of load cases: '))
    #
    istep = int(input('Number of points around ring: '))
    orad = float(input('Outside radius of ring: '))
    #
    fow = float(input('Outside flange wdth: '))
    fot = float(input('Outside flange thk: '))
    #
    fiw = float(input('inside flange wdth: '))
    fit = float(input('inside flange thk: '))
    #
    wlh = float(input('web dpth: '))
    wth = float(input('web thk: '))
    #
    # Input Data
    ring1 = RoakRing()
    ring1.RingData(orad, istep)
    ring1.outside_flange_section(fow, fot)
    ring1.inside_flange_section(fiw, fit)
    ring1.web_section(wlh, wth)
    #
    print (' ')
    for i in range(iload) :
        #
        print('Load Case :', i+1)
        roak_case = int(input('Roark case: '))
        apld = float(input('Applied load: '))
        phase = float(input('Phase: '))
        theta = float(input('Theta: '))
        thi = float(input('Thi: '))
        #
        ring1.case(roak_case, apld, phase, theta, thi)
        print(' ')
        #
    #
    # Printing Results
    ring1.print_results()
    #
    #