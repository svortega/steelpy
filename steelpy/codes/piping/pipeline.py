# Copyright (c) 2015-2016 steelpy

# Python stdlib imports
import sys
import math
import re

# package imports
#from steelpy.sectionproperty import section
#from steelpy.material.material import Material
#from steelpy.wave.theory import WaveStokes5
#import steelpy.codes.piping.process as process
#from steelpy.codes.piping.ASME_B313 import ASME
#from steelpy.codes.piping.DNV_F101 import DNV
#from steelpy.codes.piping.PD8010_2015_1 import PD8010_2015_1
#from steelpy.codes.piping.PD8010_2015_2 import PD8010_2015_2
#import steelpy.codes.piping.printout as print_results
#
#
#
#-------------------------------------------------
#                Supporting Section
#-------------------------------------------------
#
#
def search_line(lineIn, key, keyWord=None, count=1):
    """
    search key word anywere in the the string
    """
    lineOut = lineIn
    _match = False
    for _key, _item in key.items():
        rgx = re.compile(_item, re.IGNORECASE)
        keys = rgx.search(lineIn)

        if keys:
            keyWord = _key
            lineOut = re.sub(keys.group(), " ", lineIn, count)
            _match = True

    lineOut = lineOut.strip()
    return keyWord, lineOut, _match
#
def get_code(lineIn):
    """
    """
    _key = {"PD8010_1": r"\b((PD|BS)\s*((\_|\-)\s*)?8010\s*((\_|\-)\s*)?(part)?\s*((\_|\-)\s*)?1)\b",
            "PD8010_2": r"\b((PD|BS)\s*((\_|\-)\s*)?8010\s*((\_|\-)\s*)?(part)?\s*((\_|\-)\s*)?2)\b",
            "DNV_F101": r"\b(DNV\s*((\_|\-)\s*)?(OS)?\s*((\_|\-)\s*)?F101)\b",
            "ASME_B31": r"\b(ASME\s*((\_|\-)\s*)?B\s*((\_|\-)\s*)?31\s*((\_|\-|\.)\s*)?3)\b"}

    keyWord, lineOut, _match = search_line(lineIn, _key)

    return keyWord  
#
def get_pipe_type(lineIn):
    """
    """
    _key = {"Riser": r"\b((riser|landfall))\b",
            "Seabed": r"\b((sea\s*bed|mud\s*line)(including\s*tie\s*((\_|\-)\s*)?in)?)\b",
            "AboveGround": r"\b(above\s*((\_|\-)\s*)?ground)\b",
            "Buried": r"\b(below\s*((\_|\-)\s*)?ground|buried)\b"}

    keyWord, lineOut, _match = search_line(lineIn, _key)

    return keyWord
#
#
def get_design_condition(lineIn):
    """
    """
    _key = {"hydrotest": r"\b(hydro(test)?)\b",
            "operating": r"\b(op(erating)?)\b",
            "ultimate": r"\b((ultimate|user))\b"}

    keyWord, lineOut, _match = search_line(lineIn, _key)

    return keyWord
#
def get_design_method(lineIn):
    """
    """
    _key = {"stress": r"\b(stress)\b",
            "strain": r"\b(strain)\b",
            "upheaval" : r"\b(upheaval\s*(buckling)?)\b"}

    keyWord, lineOut, _match = search_line(lineIn, _key)

    return keyWord
#
#
#-------------------------------------------------
#               Main Driver Section
#-------------------------------------------------
#
class Pipeline_Assessment:
    """
    """
    #
    def __init__(self, RootSearch="FAST"):
        """
        """
        #
        self.root_search = RootSearch.upper()
        #
        self.header = 1
        #
        # Default General Data
        #
        self.pipe_name = 'Pipe'
        self.pipe_type = None
        self.design_condition = None
        #self.file_out = 'Pipe_PD8010.out'
        #self.design_factors_type = 'AUTO'
        #
        #
        # Default Pipe Section
        # Corrosion Allowance
        self.tcorr = 0
        # 
        self.tfab = 0
        # Maximum Ovality
        self.fo = 0.025
        # Manufacturing Tolerance (-ve)
        # (include percentage sign)
        self.tol = 0.125
        #
        #
        # Default Material
        #
        #       [N/mm2]
        self.SMYS = 248.0
        #        [N/mm2]
        self.E = 205000.0
        #        [N/mm2]
        self.G = 77200.0
        #        [kg/m3]
        self.rho_s = 7850
        #        (MN/m3)
        self.rho_w = 0.00001005
        #         (ksi)[N/mm2]
        self.SMTS = self.SMYS/0.75
        # poisson
        self.Poisson = 0.30
        # gravity [m/s^2]
        self.g = 9.810
        # Coefficent of termal exp
        # (C^-1)
        self.alpha_T = 1.170E-05
        #
        #
        self.material_derating = False
        self.material_type = "CMN"
        self.derate_method = "DNV"
        self.Tmax = 0.0
        #
        # Default Length
        self.Lc = 0
        #
        # Default Pressure
        #
        # Internal pressure
        self.Pi = None
        # External Pressure
        self.Po = None
        #
        self.pressure = False
        #
        # Default Temperature
        self.T1 = None
        self.T2 = None
        
        #self.Ta = 20.0
        self.alpha = 1.10E-05
        self.temperature = False
        self.Stemp = 0.0
        #
        # Default Loading
        self.load_type = None
        #
        #
        # Default Flexibility Factors 
        self.FS_code = None
        self.SIFs_type = 'straight'
        self.pipe_description = 'pipe'
        #
        #
        #
        # Default Hidro Data
        self.wave_data = False
        #
        self.flanges = None
        self.T_ = None
        self.r2 = None
        self.S = None
        self.theta = None
        self.rx = 0
        self.Tc = 0
        self.R1 = 1
        self.Tr = 0
        # ------
        #
        self.load_type = None
        self.Fx = None
        self.F = 0
        self.Fs = 0
        self.Mb = 0
        self.T = 0
        #
        self.T_inlet = False
        self.ovality = False
        self.delta_T = None
        #
        self.sigma_e = 0
        self.sigma_L = 0
        self.tau = 0
        self.sigma_h = 0
        self.Sb = 0
        #
        self.Px = 0
    #
    def design_code(self, code):
        """
        code : PD8010 Part 1 (2015) Steel pipelines on land
               PD8010 Part 2 (2015) Subsea pipelines
               DNV-OS-F101   (2007) Submarine Pipeline Systems
               ASME B31.3    (2006) Process Piping
        """
        self.code = get_code(code)
        if not self.code:
            print('   *** error design code {:} not recognised'.format(code))
            sys.exit()
        #
    #
    def design_data(self, design_condition, design_method):
        """
        desing_condition : hydrotest/operating
        method : stress/strain
        """
        #
        self.design_condition = get_design_condition(design_condition)
        if not self.design_condition:
            print('   *** error design condition {:} not recognised'
                  .format(desing_condition))
            sys.exit()
        
        self.design_method = get_design_method(design_method)
        if not self.design_method:
            print('   *** error design method {:} not recognised'
                  .format(method))
            sys.exit()
        #
    #
    def pipe_data(self, pipe_type, pipe_name=None,
                  pipe_restraint=False, pipe_history=False):
        """
        pipe_type : riser/landfall/seabed (subsea)
                    above ground/buried    (land)
        
        name      : pipe name/identification
        
        restrain  : True/False (A pipeline is deemed to be totally restrained 
                                when axial movement and bending resulting from
                                temperature or pressure change is totally prevented)
        
        pipe_history : unknown = False
                         known = True --> weld joint factor = 1.0
        
        """
        self.pipe_type = get_pipe_type(pipe_type)
        if not self.pipe_type:
            print('   *** error pipe type {:} not recognised'.format(pipe_type))
            sys.exit()
        
        self.pipe_restrained = pipe_restraint
        self.pipe_history = pipe_history
        self.pipe_name = pipe_name
    #
    def pipe_section(self, Do, tnom, tcorr,
                     tol=0.125, fo=0.025, 
                     Lc=0, coating=0):
        """
        Do    : Outside diameter of a pipe
        tnom  : Nominal wall thickness
        tcorr : Corrosion allowance
        tol   : Manufacturing Tolerance (-ve) (include percentage sign)
        fo    : Maximum Ovality
        Lc    : Pipe section length
        coating : pipe over coating
        """
        #
        self.Do = float(Do)
        self.tnom = float(tnom)
        self.tcorr = float(tcorr)
        self.tol = float(tol)
        self.fo = float(fo)
        self.Lc = Lc
        self.coating =  abs(coating)
        #        
    #
    def material(self, SMYS, SMTS=0, E=200000, alpha=1.10E-05,
                 Poisson=0.30, G=77200.0):
        """
        SMYS : Specified minimum yield strength
        E    : Young's modulus of Elasticity
        alpha : Linear coefficient of thermal expansion
        Poisson : Poisson's ratio
        """
        # 
        self.SMYS = float(SMYS)
        self.E = float(E)
        self.alpha_T = float(alpha)
        self.Poisson = float (Poisson)
        self.G = float(G)
        self.SMTS = float(SMTS)
        #
        if self.SMTS == 0 :
            self.SMTS= self.SMYS/0.90
            #
    #
    def material_derate(self, material_type="CMN", derate_code="DNV", Tmax=0.0):
        """
        material_type
        derate_code
        Tmax
        """
        #
        self.material_type = material_type.upper()
        self.derate_method = derate_code.upper()
        self.Tmax = Tmax
        #
        self.material_derating = True
        #
    #
    #
    def design_pressure(self, Pi, Po=0):
        """
        Pi : Internal pressure
        Po : External pressure
        """
        self.Pi = float(Pi)
        self.Po = float(Po)
        self.pressure = True
    #
    def design_temperature(self, T1=0, T2=0, delta_T=None,
                           T_sw=0, T_inlet=None):
        """
        T1 : Installation temperature
        T2 : Maximum or minimum metal temperature
        delta_T : 
        Tw : Seawater temperature
        """
        #
        self.T1 = float(T1)
        self.T2 = float(T2)
        self.T_sw = float(T_sw)
        
        if delta_T:
            self.delta_T = delta_T
        #
        if T_inlet:
            self.T_inlet = []
            for j in range(T_inlet):
                self.T_inlet.append(j)
        #
        self.temperature = True
        #
    #
    def functional_load(self, F=0, Fs=0, Mb=0, T=0):
        """
        F  : Axial force
        Fs : Shear Force applied to a pipeline
        Mb : bending moment applied to a pipeline
        T  : Torque applied to the pipeline
        
        Loading conditions that should be classified as functional loads include:
        a) weight of the pipeline system and its contents
        b) thermal effects
        c) pressure effects
        d) transient operational effects
        e) hydrostatic pressure of the environment
        f) residual installation load remaining after hydrotest
        """      
        self.F = F
        self.Fs = abs(Fs)
        self.Mb = abs(Mb)
        self.T = abs(T)
        #
    #
    def enviromental_load(self, Px=0, Vip=0, Vop=0, BMip=0, BMop=0, BMt=0):
        """
        Px   : Axial force
        Vip  : Shear Force In Plane
        Vop  : Shear Force Out of Plane
        BMip : Bending moment In Plane
        BMop : Bending moment Out of Plane
        BMvt : Torsional bending moment
        
        Environmental loads can be due to wind, waves, currents, earthquakes and
        other environmental phenomena.
        """
        #
        # Axial        
        self.Px = Px
        # Shear        
        self.Vip = Vip
        self.Vop = Vop
        # Bending
        self.BMip = BMip
        self.BMop = BMop
        self.BMt = BMt
        #
        #
        self.load_type = "load"
        #
    #
    #
    def stress_input(self, sigma_h, Tau, sigma_L=0):
        """
        sigma_h : Hoop stress
        tau     : Shear stress
        sigma_L : Longitudinal stress
        """
        #
        self.sigma_h = float(sigma_h)
        #
        self.sigma_L = float(SigmaL)
        self.tau = float(Tau)
        #
        self.load_type = "fea_stress"
        #
    #
    def strain_input(self, epsilon_b, epsilon_ph=None, 
                     epsilon_pL=None, epsilon_pr=None):
        """
        epsilon_b : Maximum bending strain
        epsilon_ph : Principal circunferential (hoop) strain
        epsilon_pL : Longitudinal plastic strain
        epsilon_pr : Radial plastic strain
        """
        #
        self.epsilon_b = epsilon_b
        self.epsilon_pL = epsilon_pL
        self.epsilon_ph = epsilon_ph
        self.epsilon_pr = epsilon_pr
        #
        self.load_type = "fea_strain"
        #
    #
    def hydrostatic_pressure(self, Hw, Tw, d, z, wave_length=None,
                             Rhow = 0.0000010250, g = 9.810):
        """
        Hw : Wave height
        Tw : Wave period
        d  : Water depth
        z  : Depth of the member relative to still water level
             Measure positive upwards (mm)
        wave_length = Wave Length  (mm)
        """
        #
        self.Hw = Hw
        self.Tw = Tw
        self.d = d
        self.z = z
        #
        self.wave_length = wave_length
        self.rho_w = Rhow
        self.g = g
        #
        self.wave_data = True
        #
    #
    #
    def bend_data(self, FS_code="ASME_B31.3", description='PIPE-BEND',
                 flanges=0, T_= 0, r2=0, R1S=0, theta=0, rx=0):
        """
        FS_code :
        description :
        flanges :
        T   : Nominal wall thickness of the fitting (mm)
        r2  : Mean radius of matching piping (mm)
        R1S : Bend radius of piping (mm)
        theta : 
        rx  :
        """
        #
        self.FS_code = FS_code
        self.pipe_description = description
        self.flanges = flanges
        self.T_ = T_
        self.r2 = r2
        self.rx = rx
        #
        if self.pipe_description == 'pipe-bend':
            self.R1 = R1S
        #
        else:
            self.S = R1S
            self.theta = theta
            #
        #
        self.SIFs_type = 'bend'
        #
    #
    def tee_data (self, FS_code="ASME_B31.3", description='WELDING-TEE',
                  flanges=0, T_=0, r2=0, Tc=0, rx=0):
        """
        FS_code :
        description :
        flanges :
        T :
        r2 :
        R1S
        theta
        """
        #
        self.FS_code = FS_code.upper()
        self.pipe_description = description.upper()
        self.flanges = flanges
        self.T_ = T_
        self.r2 = r2
        #
        if self.pipe_description == 'REINFORCED-TEE':
            self.Tr = Tc
        #
        else:
            self.Tc = Tc
            self.rx = rx
        #
        self.SIFs_type = 'TEE'
        #
    #
    def SIFs(self, li=1.0, lo=1.0, C1=1.0, k=1.0):
        """
        li : Stress Intensification In-Plane
        lo : Stress Intensification Out-of-Plane
        C1 :
        k  : Flexibility Factor
        """
        #
        # 
        self.li = self(li)
        self.lo = float(lo)
        self.C1 = float(C1)
        self.k = float(k)
        #
        self.SIFs_type = "user"
        #
    #
    #
    def submerged_weight(self, w_operation, w_installation):
        """
        w_opearion  (N/m)
        w_installation (N/m)
        """
        self.w_operation = w_operation
        self.w_installation =  w_installation
    #
    def soil_conditions(self, rho_sub, f):
        """
        rho_sub = Cover submerged unit weight (N/m^3)
        f : uplift coefficient
        """
        self.rho_sub = rho_sub
        self.f = f
    #
    def design_factors(self, fa, fa_hs, e=1.0):
        """
        e    : weld_joint_factor
        fd   : Design factor
        fd_h : Hoop stress desigh factor
        """
        #
        self.fa = fa
        self.fa_hs = fa_hs
        self.design_factors_type = 'USER'
        self.e = e
        #
    #
    def set_ovality(self, j):
        """
        """
        self.ovality = j
        
    #
    def substance_category(self, category):
        """
        fluids_category = A/B/C/D/E
        
        Categorization of fluids according 
        to hazard potential
        """
        #
        self.substance_category = category
        #        
    #
    #
    def get_results(self, output=True):
        """
        """
        #
        Fx = 0
        # Flag Tension or Compression
        try:
            if self.F/abs(self.F) == -1.0:
                Fx = self.F
        except ZeroDivisionError:
            pass
        #
        try:
            if self.Px/abs(self.Px) == -1.0:
                Fx += self.Px
        except ZeroDivisionError:
            pass 
        #
        if Fx != 0:
            self.Fx = Fx
        #
        # Temperature
        if not self.delta_T:
            self.delta_T = [(self.T2 - self.T1) - self.T_sw]
        #
        if self.Tmax == 0.0:
            self.Tmax = max(self.T1, self.T2)
        #
        # Material Derate
        #
        # Derate Material
        if self.material_derating:
            Fy, Fu = process.material_derating(self)
            self.sigma_y = Fy
            self.sigma_u = Fu
        # No Materail Derating
        else:
            self.sigma_y = self.SMYS
            self.sigma_u = self.SMTS
        #
        #
        #if self.header == 1:
        #    print('header')

        #
        # Po is calculated based on hydro pressure if
        # card was activated by user.
        if self.wave_data :
            self.Po = process.get_external_pressure(self, output)
        #
        #
        if "PD8010_1" in self.code :
            process.get_PD80101(self)
    
        elif "PD8010_2" in self.code :
            process.get_PD80102(self, output)
        #
        elif self.code == "DNV":
            get_DVN(self)
        #
        elif self.code == "ASME":
            get_ASME(self)
        #
        else:
            print("CODE NOT RECOGNIZED")
            exit()
    #
    def print_results(self, file_out=None):
        """
        """
        # User Name output file
        if not file_out:
            self.file_out = str(self.pipe_name) +'_pipe.out'
        
        else:
            self.file_out = file_out        
        #
        # ==========================
        #     Write Out Report
        # ==========================
        #        
        output = print_results.header(self)
        output.extend(print_results.pipe_geometry(self))
        #output.extend(print_results.section_properties(self))
        output.extend(print_results.material_properties(self))
        #
        if self.load_type:
            output.extend(print_results.design_data(self))
        
        if self.FS_code:
            output.extend(print_results.flexibility_stress_factors(self))
        
        if self.load_type:
            output.extend(print_results.stress_calculation(self))
        
        #
        if "PD8010" in self.code :
            # Part 1 : Steel pipelines on land
            if '_1' in self.code :
                output.extend(print_results.allowable_stress(self))
                self.header = self.header + 1
            # Part 2
            else:
                if self.wave_data:
                    output.extend(print_results.hydrostatic_pressure(self))
                
                output.extend(print_results.allowable_stress(self))
                
                try:
                    output.extend(print_results.buckling(self))
                    
                    output.extend(print_results.propagation_upheaval(self))
                
                    #output.extend(print_results.upheaval_buckling(self))  
                
                except AttributeError:
                    pass
                

                
                self.header = self.header + 1
                #                
        #
        elif "DNV" in self.code :
            pass
        #
        elif "ASME" in self.code :
            pass
        #
        else:
            print("CODE NOT RECOGNIZED")
            exit()
        #
        output_file = open(file_out,'w+')
        output_file.write("".join(output))
        output_file.close()
        print('--->')
#
#
#