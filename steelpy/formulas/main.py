# Copyright (c) 2009-2023 steelpy
#

# Python stdlib imports

# package imports
#from steelpy.material.main import Materials
from steelpy.process.units.main import Units
#from steelpy.sections.main import Sections
#from steelpy.sections.tee import Tee
#from steelpy.sections.shapes.solid import Rectangle
#from steelpy.codes.api.main import API_design
from .roarks.chapter9.ring_main import Ring
#import steelpy.formulas.roarks.circular_ring as ring
#
#from .pilkey.chapter11.table2 import SimpleBeam
from .pilkey.main import BeamBasic
#
#-------------------------------------------------
#                Supporting Section
#-------------------------------------------------
#
class Formulas:
    __slots__ = ['_units']
    
    def __init__(self):
        """
        """
        self._units = Units()
        #self._ring = Ring()
    #
    def ring(self, istep:int = 360):
        """ """
        return Ring(istep)
#
#if __name__=="__main__": 
#    #
#    # Input data etc
#    #
#    #
#    print (" ")
#    print ("STRESSES IN CIRCULAR RINGS-AFTER R.J.ROARK & W.C.YOUNG")
#    print ("******************************************************")
#    #
#    #
#    #filenm = input(" Enter output file name: ")
#    #
#    #open (7, file=filenm, status="new")
#    #
#    #ititle = input(" Enter job title: ")
#    #
#    #
#    print (" ")
#    print('DATA ENTRY')
#    #
#    iload = int(input('Number of load cases: '))
#    #
#    istep = int(input('Number of points around ring: '))
#    orad = float(input('Outside radius of ring: '))
#    #
#    fow = float(input('Outside flange wdth: '))
#    fot = float(input('Outside flange thk: '))
#    #
#    fiw = float(input('inside flange wdth: '))
#    fit = float(input('inside flange thk: '))
#    #
#    wlh = float(input('web dpth: '))
#    wth = float(input('web thk: '))
#    #
#    # Input Data
#    ring1 = RoakRing()
#    ring1.RingData(orad, istep)
#    ring1.outside_flange_section(fow, fot)
#    ring1.inside_flange_section(fiw, fit)
#    ring1.web_section(wlh, wth)
#    #
#    print (' ')
#    for i in range(iload) :
#        #
#        print('Load Case :', i+1)
#        roak_case = int(input('Roark case: '))
#        apld = float(input('Applied load: '))
#        phase = float(input('Phase: '))
#        theta = float(input('Theta: '))
#        thi = float(input('Thi: '))
#        #
#        ring1.case(roak_case, apld, phase, theta, thi)
#        print(' ')
#        #
#    #
#    # Printing Results
#    ring1.print_results()
#    #
#    #
#
#
#