# 
# Copyright (c) 2019-2020 steelpy
#
# Python stdlib imports
#import math
#import datetime
from typing import Union, ClassVar

# package imports
from steelpy.codes.api.wsd_22ed import APIwsd22ed
from steelpy.codes.process.process import BeamDesignParameters, CodeResults
from steelpy.f2uModel.load.actions import Actions
from steelpy.f2uModel.sections.process.stress import BeamStress
import steelpy.codes.process.print_report as print_report

#
#
class API_design(BeamDesignParameters):
    """
    """
    def __init__(self):
        """
        """
        BeamDesignParameters.__init__(self)
        self._actions:ClassVar = Actions()
        #
        # Default output
        self.member_name:Union[str:int] = 'N/A'
        #self.file_out = 'APIwsdCodeCheck.out'
        self.member_number:int = 0
        #
        #      Factors
        #self.flag_axial = 'COMPRESSION'
        self.design_condition:str = 'OPERATING'
        self.SFxt:float =  1.67
        self.SFb:float = 1.333
        self.SFxc:float = 2.00
        self.SFh:float = 2.00
        #
        self.hydro_check:bool = False
        #
        # internal stress 
        self.stress = BeamStress(0 * self.units.MPa, 0 * self.units.MPa, 0 * self.units.MPa, 
                                 0 * self.units.MPa, 0 * self.units.MPa, 0* self.units.MPa)        
    #
    @property
    def actions(self):
        """
        """
        return self._actions
    
    @actions.setter
    def actions(self, value):
        """
        """
        self._actions = value
    #
    #
    @property
    def wave(self):
        """
        """
        self.hydro_check = True
        return self._wave
    
    @wave.setter
    def wave(self, value):
        """
        """
        self.hydro_check = True
        self._wave = value
    #
    #
    def applied_loads(self, Px, Vy, Vz, BMx, BMy, BMz):
        """
        """
        # Axial        
        self.actions.Fx = Px
        # Shear        
        self.actions.Fz = Vy # ip
        self.actions.Fy = Vz # op
        # Bending
        self.actions.Mx = BMx
        self.actions.My = BMy # ip
        self.actions.Mz = BMz # op
    #
    #
    def reduction_factor(self):
        """
        6.3.2.5
        """
        pass
    #
    def safety_factors(self, SFxt, SFb, SFxc, SFh):
        """
        """
        self.SFxt = SFxt
        self.SFb = SFb
        self.SFxc = SFxc
        self.SFh = SFh
        self.design_condition = 'USER'
    #
    #
    def WSD(self):
        """
        """
        APIwsd = APIwsd22ed()
        APIwsd.get_load_stress(self.section, self.actions)
        #
        #print ("+++++++++++++++++++++++++++++")
        #print ("              Calculation: ", self.Header)
        #
        APIwsd.axial_tension(self.section, self.material)
        APIwsd.axial_compression(Klr=max(self.Klrz, self.Klry), 
                                          section=self.section, 
                                          material=self.material)
        APIwsd.bending(self.section, self.material)
        APIwsd.shear(self.section, self.material)
        APIwsd.torsional_shear(self.section, self.material)
        #
        # Hydro Check
        if self.hydro_check :
            #print ('Wave Length = {:1.3f} m'.format(self._wave.wave_length.value))
            # Calc Hydro Pressure
            APIwsd.safety_factors(self.design_condition, self.material,
                                  self.SFxt, self.SFb, self.SFxc, self.SFh)            
            APIwsd.hydrostatic_pressure(self.z, self.wave)
            APIwsd.hoop_buckling(self.Lr, self.section, self.material)
            APIwsd.combination_hydro(self.material)
            # Print out Results
            #APIwsd.PrintSectionProperties(self)
            #APIwsd.PrintShear(self)
            #APIwsd.PrintHydroPressureCalc(self)
            #APIwsd.PrintAxialBendingAndHP(self) 
            #
            #self.Header += 1
            comb = APIwsd.print_hydro_pressure(self)
            comb.extend(APIwsd.print_axial_bending_andHP(self))
        else :
            APIwsd.combination(self.Klrz, self.Klry, 
                                        self.Cmz, self.Cmy,
                                        self.material)
            #
            #APIwsd22ed.PrintSectionProperties(self)
            #APIwsd22ed.PrintShear(self)
            comb = APIwsd.print_axial_bending_noHP(self)
        #
        #
        #self.Header += 1
        code_detail = ["API RP2A-WSD-Ed22 [Nov 2014]", "Strength Of Tubular Members"]
        output = print_report.print_header(code_detail)
        #
        # Calc Section Properties
        output.extend(self.section.print_properties())
        output.extend(self.material.print_properties())
        output.extend(APIwsd.print_shear())
        output.extend(comb)
        #
        self.results = CodeResults(axial= APIwsd.axial_results,
                                   shear= APIwsd.shear_results,
                                   bending= APIwsd.bending_results,
                                   combined= APIwsd.combined_results,
                                   report=output)
    #
    def LRFD(self):
        """
        """
        pass
    #
    def print_results(self):
        """
        """
        #
        self.section.properties
        #
        print("Ly : {:1.3f} m".format(self.Ly.value))
        print("Lz : {:1.3f} m".format(self.Lz.value))
        #
        self.Klrz = self.Kz * self.Lz.value / self.section.rz.value
        self.Klry = self.Ky * self.Ly.value / self.section.ry.value
        #
        try:
            self.Lr
        except AttributeError:
            self.Lr = max(self.Ly, self.Lz)
        #
        self.WSD()
        #
        for line in self.results.report:
            print(line.rstrip())
    #
    #