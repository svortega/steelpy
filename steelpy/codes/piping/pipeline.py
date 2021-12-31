# Copyright (c) 2015-2016 steelpy

# Python stdlib imports
from __future__ import annotations
import sys
import math
import re

# package imports
from steelpy.process.units.main import Units
from steelpy.codes.piping.PD8010_2015_1 import PD8010_2015_1
from steelpy.codes.piping.PD8010_2015_2 import PD8010_2015_2
from steelpy.codes.piping.process import expansion_flexibility, PD8010_buckling
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
    match = False
    for _key, _item in key.items():
        rgx = re.compile(_item, re.IGNORECASE)
        keys = rgx.search(lineIn)
        if keys:
            keyWord = _key
            lineOut = re.sub(keys.group(), " ", lineIn, count)
            match = True
    lineOut = lineOut.strip()
    return keyWord, lineOut, match
#
def get_code(lineIn):
    """
    """
    key = {"PD8010_1": r"\b((PD|BS)\s*((\_|\-)\s*)?8010\s*((\_|\-)\s*)?(part)?\s*((\_|\-)\s*)?1)\b",
            "PD8010_2": r"\b((PD|BS)\s*((\_|\-)\s*)?8010\s*((\_|\-)\s*)?(part)?\s*((\_|\-)\s*)?2)\b",
            "DNV_F101": r"\b(DNV\s*((\_|\-)\s*)?(OS)?\s*((\_|\-)\s*)?F101)\b",
            "ASME_B31": r"\b(ASME\s*((\_|\-)\s*)?B\s*((\_|\-)\s*)?31\s*((\_|\-|\.)\s*)?3)\b"}

    keyWord, lineOut, _match = search_line(lineIn, key)

    return keyWord  
#
def get_pipe_type(lineIn):
    """
    """
    key = {"Riser": r"\b((riser|landfall))\b",
            "Seabed": r"\b((sea\s*bed|mud\s*line)(including\s*tie\s*((\_|\-)\s*)?in)?)\b",
            "AboveGround": r"\b(above\s*((\_|\-)\s*)?ground)\b",
            "Buried": r"\b(below\s*((\_|\-)\s*)?ground|buried)\b"}
    keyWord, lineOut, _match = search_line(lineIn, key)
    return keyWord
#
#
def get_design_condition(lineIn):
    """
    """
    key = {"hydrotest": r"\b(hydro(test)?)\b",
            "operating": r"\b(op(erating)?)\b",
            "ultimate": r"\b((ultimate|user))\b"}

    keyWord, lineOut, _match = search_line(lineIn, key)
    return keyWord
#
def get_design_method(lineIn):
    """
    """
    key = {"stress": r"\b(stress)\b",
            "strain": r"\b(strain)\b",
            "upheaval" : r"\b(upheaval\s*(buckling)?)\b"}

    keyWord, lineOut, _match = search_line(lineIn, key)
    return keyWord
#
#
def functional_stress(self, t_check, PD8010, output):
    """
    """
    # External overpressure
    P = abs(self.Pi - self.Po)
    # hoop stress
    sigma_h = PD8010.hoop_stress(P, tmin=t_check)
    # longitudinal stress
    deltaT = self.delta_T[0]
    sigma_L, Sb = PD8010.longitudinal_stress(t_check, 
                                             sigma_h, deltaT, 
                                             P, self.Fa, self.Mb,
                                             self.pipe_restrained,
                                             output=output)
    # shear stress
    tau = PD8010.shear_stress(t_check, self.T, self.Fs, output=output)
    # 
    sigma_e = PD8010.equivalent_stress(sigma_h, sigma_L, tau, output=output)   
    return sigma_e
#
#
def get_PD80101(self, output):
    """
    """
    PD8010 = PD8010_2015_1()

    PD8010.section_input(self.Do, self.tnom, 
                         self.tcorr, self.tol,
                         output=output)

    PD8010.material_input(self.sigma_y,
                          self.sigma_u,
                          self.material.E,
                          self.material.poisson,
                          self.material.alpha,
                          output=output)
    #
    self.code_name = "PD8010-1:2015 Part 1 : Steel Pipelines on Land"
    #
    # Start Calculation
    #
    # 6.4 Strength Check
    t_check = self.tnom
    #calculate_stress(self, t_check, PD8010)
    self.sigma_e += functional_stress(self, t_check, PD8010, output=output)
    if self.load_type :
        # Modify stress from detail FEA
        # to include thck variation
        if 'fea_stress' in self.load_type:
            self.sigma_e += PD8010.stressFEA(self)
        else:
            self.sigma_e += PD8010.enviromental_stress(self.Px, self.Vip, self.Vop,
                                                       self.BMip, self.BMop, self.BMt,
                                                       output=output)
            #self.M = math.sqrt(self.BMip**2 + self.BMop**2)
            #self.Mt = math.sqrt(self.Vip**2 + self.Vop**2)         
    #
    #self.sigma_e += sigmae_1 #+ sigmae_2
    # code check stresses
    self.Sae, self.fd = PD8010.allowable_equivalent_stress()
    self.a = PD8010.substances_categorization(self.substance_category)
    self.Sah, self.fd_hs = PD8010.allowable_hoop_stress(self.a, self.pipe_history)
    # External overpressure
    P = abs(self.Pi - self.Po)    
    self.sigma_h = PD8010.hoop_stress(P, tmin=t_check)
    self.UR_h, self.UR_eq = PD8010.limits_calculated_stress(self.sigma_h, self.Sah, 
                                                            self.sigma_e, self.Sae,
                                                            self.design_method)
    #
    # 6.4.4 Buckling
    PD8010_buckling(self, PD8010, output=output)    
    #
    #print('-->')
    
#
def get_PD80102(self, output):
    """
    """
    # Code
    self.code_name = "PD8010-2:2015 Part 2 : Subsea Pipelines"
    #
    PD8010 = PD8010_2015_2()
    # Section
    self.tmin, self.Di = PD8010.section_input(self.Do, self.tnom,
                                              self.tcorr, self.tol,
                                              output=output)
    # Material
    PD8010.material_input(self.sigma_y,
                          self.sigma_u,
                          self.material.E,
                          self.material.poisson,
                          self.material.alpha,
                          output=output)
    # 6.4.1 Design Factors
    self.fd, self.fd_hs = PD8010.design_factors(self.design_condition,
                                                self.pipe_type, 
                                                output=output)
    #6.4.2.3 Expansion and flexibility
    if self.SIFs_type != "user":
        (self.h, self.k, self.lo, 
         self.li, self.C1) = expansion_flexibility(self, output=output)
    #
    # 6.4.3 Strain Design Method
    if 'strain' in self.design_method:
        get_PD8010_strain(self, output)
    # 6.4.2 Stress Design Method
    else:
        # 6.4 Strength Check
        if self.load_type:
            # Modify stress from detail FEA
            # to include thck variation
            if 'fea_stress' in self.load_type:
                self.sigma_e += PD8010.stressFEA(self)
            else:
                self.sigma_e += PD8010.enviromental_stress(self.Px, self.Vip, self.Vop,
                                                           self.BMip, self.BMop, self.BMt,
                                                           output=output)
                #self.M = math.sqrt(self.BMip**2 + self.BMop**2)
                #self.Mt = math.sqrt(self.Vip**2 + self.Vop**2)                
        # if Load input in forces
        #else:
        t_check = self.tnom
        self.sigma_e += functional_stress(self, t_check, PD8010, output=output)
        #
        #
        #self.sigma_e += sigmae_1 + sigmae_2
        if output:
            print('')
            print('Total equivalent stress : {:2.3f}'
                  .format(self.sigma_e.convert('megapascal').value))
        #
        # Overpressure
        P = abs(self.Pi - self.Po)
        #if P.value < 0:
        #    P = self.Po
        self.sigma_h = PD8010.hoop_stress(P, tmin=self.tmin, output=output)
        # code check stresses
        self.UR_h, self.UR_eq = PD8010.allowable_stress(self.sigma_h,
                                                        self.sigma_e,
                                                        self.fd,
                                                        self.fd_hs,
                                                        self.design_method,
                                                        output=output)
    # 6.4.4 Buckling
    PD8010_buckling(self, PD8010, output=output)
#
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
        self.units = Units()
        #
        self.root_search = RootSearch.upper()
        #
        self.header = 1
        #
        ## Default General Data
        ##
        ##
        ##
        self.material_derating = False
        #self.material_type = "CMN"
        #self.derate_method = "DNV"
        #self.Tmax = 0.0
        ##
        ##
        ## Default Loading
        self.load_type = False
        ##
        ##
        ## Default Flexibility Factors 
        self.FS_code = None
        self.SIFs_type = 'straight'
        self.pipe_description = 'pipe'
        ##
        ##
        ##
        ## Default Hidro Data
        self.wave_data = False
        ##
        #self.flanges = None
        #self.T_ = None
        #self.r2 = None
        #self.S = None
        #self.theta = None
        #self.rx = 0
        #self.Tc = 0
        #self.R1 = 1
        #self.Tr = 0
        ## ------
        ##
        #self.load_type = None
        self.Fa = 0 * self.units.N
        self.Fs = 0 * self.units.N
        self.Mb = 0 * self.units.N * self.units.m
        self.T  = 0 * self.units.N * self.units.m
        ##
        #self.T_inlet = False
        #self.ovality = False
        #self.delta_T = None
        ##
        # Equivalent stress (8010-2)
        self.sigma_e = 0 * self.units.MPa
        #self.sigma_L = 0
        #self.tau = 0
        #self.sigma_h = 0
        #self.Sb = 0
        ##
        ## Env Forces
        ##
        self.Px = 0 * self.units.N
        # Shear        
        self.Vip = 0 * self.units.N
        self.Vop = 0 * self.units.N
        # Bending
        self.BMip = 0 * self.units.N * self.units.m
        self.BMop = 0 * self.units.N * self.units.m        
        self.BMt = 0 * self.units.N * self.units.m        
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
    def pipe_section(self, tcorr:Units, tol:float=0.125, fo=0.025, 
                     Lc:Units|None=None, coating:Units|None=None):
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
        #self.Do = Do
        #self.tnom = tnom
        self.tcorr = tcorr
        self.tol = tol
        self.fo = fo
        #
        self.Lc = Lc
        if not self.Lc:
            self.Lc = 0 * self.units.m
        #
        self.coating = coating
        if not self.coating:
            self.coating = 0 * self.units.m
    #
    #
    @property
    def geometry(self):
        """
        """
        return self._geometry
    
    @geometry.setter
    def geometry(self, value:Sections):
        """
        """
        self._geometry = value  
    #    
    #
    @property
    def material(self):
        """ """
        return self._material
    
    @material.setter
    def material(self, value):
        """ """
        self._material = value
    #
    #
    #def material(self, SMYS:Units, SMTS:Units|None=None, 
    #             E:Units|None=None, G:Units|None=None, 
    #             alpha:float=1.10E-05, Poisson:float=0.30):
    #    """
    #    SMYS : Specified minimum yield strength
    #    E    : Young's modulus of Elasticity
    #    alpha : Linear coefficient of thermal expansion
    #    Poisson : Poisson's ratio
    #    """
    #    # 
    #    self.SMYS = SMYS
    #    #
    #    self.E = E
    #    if not self.E:
    #        self.E = 200000 * self.units.MPa
    #    #
    #    self.G = G
    #    if not self.G:
    #        self.G = 77200.0 * self.units.MPa   
    #    #
    #    self.alpha_T = alpha
    #    self.Poisson = float (Poisson)
    #    #
    #    self.SMTS = SMTS
    #    if not self.SMTS:
    #        self.SMTS= self.SMYS / 0.90
    #
    #
    def material_derate(self, material_type:str="CMN", 
                        derate_code:str="DNV", Tmax:Units|None=None):
        """
        material_type
        derate_code
        Tmax
        """
        #
        self.material_type = material_type.upper()
        self.derate_method = derate_code.upper()
        self.Tmax = Tmax
        if not self.Tmax:
            self.Tmax = 0*self.units.m
        self.material_derating = True
    #
    #
    def design_pressure(self, Pi, Po:Units|None=None):
        """
        Pi : Internal pressure
        Po : External pressure
        """
        self.Pi = Pi
        self.Po = Po
        if not self.Po:
            self.Po = 0 * self.units.bar
        self.pressure = True
    #
    def design_temperature(self, T1:Units|None=None, 
                           T2:Units|None=None,
                           T_sw:Units|None=None, 
                           delta_T:Units|None=None,
                           T_inlet:list|None=None):
        """
        T1 : Installation temperature
        T2 : Maximum or minimum metal temperature
        delta_T : 
        Tw : Seawater temperature
        """
        #
        self.T1 = T1
        if not self.T1:
            self.T1 = 0 * self.units.K 
        
        self.T2 = T2
        if not self.T2:
            self.T2 = 0 * self.units.K 
        #
        self.T_sw = T_sw
        if not self.T_sw:
            self.T_sw = 0 * self.units.K      
        
        self.delta_T = []
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
    def functional_load(self, Fa:Units, Fs:Units, Mb:Units, T:Units):
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
        self.Fa = Fa
        self.Fs = Fs
        self.Mb = Mb
        self.T = T
    #
    #
    def enviromental_load(self, Px:Units|None=None, 
                          Vip:Units|None=None, Vop:Units|None=None, 
                          BMip:Units|None=None, BMop:Units|None=None,
                          BMt:Units|None=None):
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
        if not self.Px:
            self.Px = 0 * self.units.N
        # Shear        
        self.Vip = Vip
        if not self.Vip:
            self.Vip = 0 * self.units.N
        self.Vop = Vop
        if not self.Vop:
            self.Vop = 0 * self.units.N
        # Bending
        self.BMip = BMip
        if not self.BMip:
            self.BMip = 0 * self.units.N * self.units.m
        self.BMop = BMop
        if not self.BMop:
            self.BMop = 0 * self.units.N * self.units.m        
        self.BMt = BMt
        if not self.BMt:
            self.BMt = 0 * self.units.N * self.units.m
        
        self.load_type = "load"
    #
    #
    def stress_input(self, sigma_h:Units, Tau:Units, 
                     sigma_L:Units|None=None):
        """
        sigma_h : Hoop stress
        tau     : Shear stress
        sigma_L : Longitudinal stress
        """
        self.sigma_h = sigma_h
        self.tau = Tau
        #
        self.sigma_L = sigma_L
        if not self.sigma_L:
            self.sigma_L = 0 * self.units.MPa
        #
        self.load_type = "fea_stress"
    #
    #
    def strain_input(self, epsilon_b:float, epsilon_ph=None, 
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
    def bend_data(self, FS_code:str="ASME_B31.3", description:str='PIPE-BEND',
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
        else:
            self.S = R1S
            self.theta = theta
        #
        self.SIFs_type = 'bend'
    #
    #
    def tee_data (self, FS_code:str="ASME_B31.3", description:str='WELDING-TEE',
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
        else:
            self.Tc = Tc
            self.rx = rx
        #
        self.SIFs_type = 'TEE'
    #
    #
    def SIFs(self, li:float=1.0, lo:float=1.0, 
             C1:float=1.0, k:float=1.0):
        """
        li : Stress Intensification In-Plane
        lo : Stress Intensification Out-of-Plane
        C1 :
        k  : Flexibility Factor
        """
        self.li = self(li)
        self.lo = lo
        self.C1 = C1
        self.k = k
        self.SIFs_type = "user"
    #
    #
    def submerged_weight(self, w_operation:Units, w_installation:Units):
        """
        w_opearion  (N/m)
        w_installation (N/m)
        """
        self.w_operation = w_operation
        self.w_installation =  w_installation
    #
    def soil_conditions(self, rho_sub:Units, f:float):
        """
        rho_sub = Cover submerged unit weight (N/m^3)
        f : uplift coefficient
        """
        self.rho_sub = rho_sub
        self.f = f
    #
    def design_factors(self, fa:float, fa_hs:float, e:float=1.0):
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
    def set_ovality(self, j:float):
        """
        """
        self.ovality = j
    #
    def set_substance_category(self, category:str):
        """
        fluids_category = 
        A : Typically non-flammable water-based fluids / Water, brine, dilute effluents
        B : Flammable and/or toxic fluids that are liquids / Oil and petroleum products, Methanol
        C : Non-flammable fluids that are non-toxic gases / Nitrogen, oxygen, argon and air
        D : Non-toxic, single-phase natural gas.
        E : Flammable and/or toxic fluids that are gases / Hydrogen, carbon dioxide, natural gas, 
                                                           ethane, ethylene, liquefied petroleum,
                                                           gas liquids, ammonia and chlorine, 
                                                           Spiked or live crude oil
        
        Categorization of fluids according 
        to hazard potential
        """
        self.substance_category = category       
    #
    #
    def get_results(self, output:bool=True):
        """
        """
        # section data
        self.Do = self.section.D * self.units.m
        self.tnom = self.section.tw * self.units.m
        #
        # Temperature
        if not self.delta_T:
            self.delta_T = [abs(self.T2 - self.T1) - self.T_sw]
        #
        #if self.Tmax == 0.0:
        self.Tmax = max(self.T1, self.T2)
        #
        # Material Derate
        #
        # Derate Material
        if self.material_derating:
            Fy, Fu = process.material_derating(self)
            self.sigma_y = Fy
            self.sigma_u = Fu
        else:
            # No Materail Derating
            self.sigma_y = self.material.Fy
            self.sigma_u = self.material.Fu
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
            get_PD80101(self, output)
    
        elif "PD8010_2" in self.code :
            get_PD80102(self, output)
        #
        #elif self.code == "DNV":
        #    get_DVN(self)
        ##
        #elif self.code == "ASME":
        #    get_ASME(self)
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
        elif "ASME" in self.code :
            pass
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