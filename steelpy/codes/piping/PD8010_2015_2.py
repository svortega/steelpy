# Copyright (c) 2015-2016 steelpy

# Python stdlib imports
import math
import datetime

# package imports
#from steelpy.codes.piping.process import pipe_section
from steelpy.math.rootsearch import GoalSeeker
#from steelpy.codes.piping.ASME_B313 import ASME
#from steelpy.codes.piping.DNV_F101 import DNV
from steelpy.wave.theory import WaveStokes5
from steelpy.sectionproperty.shapes.elliptical import hollow_ellipse

#
def pipe_section(Do, tnom, tcheck):
    """
    """
    #
    _Dinom = (Do - 2 * tnom)
    
    # Internal Area
    Ai = (math.pi / 4.0 *_Dinom**2)    
    
    # External Area
    Ae = (math.pi / 4.0 * Do**2)    
    
    # Pipe section area (nominal thickness)
    CSA = (math.pi / 4.0 *
           (Do**2 -(Do - 2*tnom)**2))
    
    # Pipe section area (Tnom - Tcorr thickness)
    Anomfab = (math.pi / 4.0 *
               (Do**2 -(Do - 2 * tcheck)**2))
    
    #
    I = (math.pi / 64.0 * 
          (Do**4 - (Do - 2 * tcheck)**4))
    #
    Ze = ( I / (Do / 2.0))
    #
    Ip = ((math.pi/32.0)*
          (Do**4 - (Do - 2*tcheck)**4))
    #
    return Ai, Ae, CSA, Anomfab, I, Ze, Ip
#
def ellipse_section(D, t, fo):
    """
    D 
    t
    fo
    
    Calculation of elastic moment of inertia for an ellipse
    """
    Dmin = (1.0 - 0.50 * fo) * D
    Dmax = (1.0 + 0.50 * fo) * D
    
    a = (Dmin - t)/ 2.0
    b = (Dmax - t) / 2.0
    
    Area, Zc, Yc, Iy, Zey, Iz, Zez = hollow_ellipse(a, b, t)
    
    return min(Iy, Iy)
#
#
# ===================================================
#      PD8010-2:2004 Part 2 : Subsea Pipelines
# ===================================================
class PD8010_2015_2:
    """
    PD 8010-1 & 2: 2015
    """
    #
    def __init__(self):
        """
        """
        self.g = 9.810 # m/s^2
        self.rho_w = 0.0000010250 # kg/mm^2
        
        # Prop height and natural half-length
        h = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50]
        self.h = [hi*1000.0 for hi in h]
    #
    # ===================================================
    # Section Properties
    # ===================================================
    #     
    def section_input(self, Do, tnom, tcorr, tol, output=False): 
        """
        Do
        tnom
        tol
        """
        #
        self.Do = Do
        self.tnom = tnom
        self.tcorr = tcorr 
        self.tol = tol
        
        self.tfab = self.tnom * self.tol
        self.tmin = self.tnom - self.tfab - self.tcorr
        self.Di = (self.Do - 2 * self.tmin)        
        #
        #-------------------------------------------------
        #   Cross-Sectional Area
        #
        self.A = ((math.pi / 4.)*
                  (self.Do**2 -(self.Do - 2*self.tmin)**2))
        #
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        #   Second Moment of Area about Mayor Axis
        #   --------------------------------------
        #
        self.I = ((math.pi / 64.0)*
                  (self.Do**4 - (self.Do - 2*self.tmin)**4))
        # print ('circular Second Moment',self.Iip)
        # 
        #
        self.Ze = ( self.I / (Do / 2.0))
        #-------------------------------------------------
        #   Plastic Modulus about Mayor Axis
        #    def Zx(self):
        #
        self.Zp = ((self.Do**3 - 
                    (self.Do - 2 * self.tmin)**3) / 6.0 )
        # print ('circular Plastic Modulus',self.Zpip)
        #
        #-------------------------------------------------
        #   Radius of gyration about Mayor Axis
        #
        self.r = (math.sqrt(self.I / self.A))
        # print ('circular Radius of gyration about Mayor Axis',self.rip)
        # 
        #self.rop = self.rip
        #
        #-------------------------------------------------
        #   Torsional Constant
        #
        self.J = (2 * self.I)
        # print ('Torsional Constant',self.J)
        # 
        #-------------------------------------------------
        #   Polar Moment of Inertia
        self.Ip = ((math.pi / 32.0) *
                   (self.Do**4 - 
                    (self.Do - 2*self.tmin)**4))
        # print ('Polar Moment of Inertia Ip:',self.Ip)
        #
        #
        #
        #-------------------------------------------------
        #  Mass
        #self.mass = (self.A * self.Rhos)/1000**2
        #
        if output:
            print ('')
            print ('Section Properties')
            print ('Do = {:} mm'.format(self.Do))
            print ('t nom = {:} mm'.format(self.tnom))
            print ('t corr = {:} mm'.format(self.tmin))
            print ('t min = {:} mm'.format(self.tmin))
            #print (' ')        
        #
        return self.tmin, self.Di
    #
    def material_input(self, sigma_y, sigma_u, E, Poisson, alpha_T, output=False):
        """
        """
        self.sigma_y = sigma_y
        self.E = E
        self.Poisson = Poisson
        self.sigma_u = sigma_u
        self.alpha_T = alpha_T
        
        if output:
            print ('')
            print ('Material Properties')
            print ('Fy = {:}'.format(self.sigma_y))
            print ('Fu = {:}'.format(self.sigma_u))
            print ('E = {:}'.format(self.E))
            print ('Coefficient of Thermal Expansion = {:}'.format(self.alpha_T))
    #    
    # 6.4 Strength
    def design_factors(self, design_condition, pipe_type, output=False):
        """
        fd   : Design factor
        fd_h : Hoop stress desigh factor
        
        
        6.4.1 Design factors
        
        The design factors, fd, should be used for the stress-based design
        method in 6.4.2. The design factors for strain-based design method 
        are given in 6.4.3.
        
        *NOTE BS EN 14161 allows a design factor of up to 0.83
              to be used. The design factors recommended for UK
              use are shown in Table 2. Higher design factors may
              be used for all or part of a pipeline system, provided
              that an equivalent level of safety is achieved throughout 
              the system under consideration, and across all relevant
              limit states. A full risk assessment (see Annex D) is 
              recommended if higher design factors are used, and might 
              be subject to regulatory review.
        
        Where the application of a reliability-based limit state design 
        method leads to a reduction in the wall thickness that is 
        necessary to meet the recommendations for pressure containment, 
        particular attention should be paid to installation and operability 
        considerations, i.e. it should be demonstrated that the selected 
        wall thickness is appropriate for all the load conditions and 
        combinations that can reasonably be expected throughout the life 
        of the pipeline.
        
        Account should be taken of the variation of strength with 
        temperature on the basis of verifiable test data appropriate to 
        the material under consideration.
        
        """
        #
        # Table 2 - Design factor, fd
        #
        #
        # Equivalent stresses arising from
        # construction or hydrotest loads
        if "hydrotest" in design_condition:
            #
            # Seabed including tie-in
            # & Riser/landfall
            self.fd = 1.0
            self.fd_hs = 1.0 # Hoop stress
        #
        # Equivalent stresses resulting from
        # functional and environmental or
        # accidental loads
        else:
            # Riser/landfall
            if 'riser' in pipe_type.lower():
                self.fd = 0.72 
                self.fd_hs = 0.60 # Hoop stress
            # Seabed including tie-in
            else: 
                self.fd = 0.96
                self.fd_hs = 0.72 # Hoop stress
        #
        if output :
            print ("")
            print ("Design Factors")
            print ("Test Data  = {:}".format(design_condition))
            print ("Pipe Type  = {:}".format(pipe_type))
            print ("fd         = {:2.3f}".format(self.fd))
            print ("fd hs      = {:2.3f}".format(self.fd_hs))
        #
        return self.fd, self.fd_hs
        #
    #
    #
    # 6.4.2 Stress-based design
    # -------------------------
    def allowable_stress(self, sigma_h, sigma_e, fd, fd_hs,
                         design_method, output=False):
        """
        
        
        6.4.2.1 Allowable stress
        
        Stress in the pipeline system should meet the 
        inequality shown in equation (2).
        where SigmaA is determined in accordance with 
        6.4.2.2 or 6.4.2.4 as appropriate.
        The effect of temperature re-rating on the SMYS should be included.
        
        NOTE In the absence of more directly applicable data, the derating data
        in DNV-OS-F101 may be used.
        """
        #
        URhoop = sigma_h / (fd_hs * self.sigma_y)
        # 
        if 'stress' in design_method:
            UReq = sigma_e / (fd * self.sigma_y)
            checkout = ("UR eq: {:2.3f}".format(UReq))
        
        else:
            UReq = self.epsilon_p / 0.0050
            
            if self.epsilon_p > 0.0050:
                checkout = ("Fail  Epsilonp > 0.0010: {:2.3f}"
                            .format(self.epsilon_p))
                
            else:
                checkout = ("Pass  Epsilonp < 0.0010: {:2.3f}"
                            .format(self.epsilon_p))
        #
        if output:
            print("")
            print("Allowable Stress")
            print("UR hoop : {:2.3f}".format(URhoop))
            print("{:}".format(checkout))
        #
        return URhoop, UReq
    #
    def hoop_stress(self, P, Do=None, tmin=None, 
                    output=False):
        """
        Do   :
        tmin : minimimum wall thickness (m)
        P : External Overpressure
        
        6.4.2.2 Hoop stress
        
        The wall thickness used for hoop stress calculation should be the minimum value
        allowing for permitted wall thickness variations, such as fabrication tolerances,
        and subtracting any corrosion allowance, i.e. as shown in equation (4).
        tnom = tmin + tfab + tcorr    (4)
        """
        #
        _tmin = self.tmin
        if tmin:
            _tmin = float(tmin)
        
        _Do = self.Do
        if Do:
            _Do = float(Do)
        #
        _Dih = (_Do - 2*_tmin)
        #
        # NOTE For clad or lined pipelines, the strength contribution
        # of the cladding or lining is usually not taken into account,
        # unless it is necessary to contribute to the structural integrity.
        #
        # For all other stress checks in this section, tnom should
        # be used in the calculation of component stresses.
        #
        if _Do/_tmin > 20:
            # Hoop stress should be calculated using equation (3)
            # (thin wall) when the ratio of Do/tmin is greater than 20.
            sigma_h = P * (_Do/(2 * _tmin))
        #
        else:
            # Equation (5) (thick wall) should be used when the ratio
            # of Do/tmin is less than or equal to 20.
            sigma_h = (P* ((_Do**2 + _Dih**2) / 
                           (_Do**2 - _Dih**2)))
            #
        #
        if output:
            print("")
            print("External Overpressure P : {:2.3f}".format(P))
            print("Hoop Stress, tnom : {:2.3f}".format(_tmin))
            print("Sigma h  = {:2.3f}".format(sigma_h))
        #
        # self.sigma_h = sigma_h
        #
        return sigma_h
        #
    #
    #
    def equivalent_stress(self, sigma_h, sigma_L, tau=0, output=False):
        """
        sigma_h
        sigma_L
        tau
        
        6.4.2.4 Equivalent stress
        
        Unless a strain-based design approach is adopted (see 6.4.3), equivalent stresses
        should be evaluated using the von Mises stress criterion shown in equation (6).
        
        NOTE 1 - Nominal wall thickness may be used in the evaluation.
        The total component longitudinal stress should be the sum of the longitudinal
        stresses arising from pressure, bending, temperature, weight (force), other
        sustained loadings and occasional loadings. Accidental loads should be taken
        into account as indicated in F.7. Account should be taken of the variation in
        axial restraint throughout the pipeline.
        
        NOTE 2 - A pipeline is deemed to be totally restrained when axial movement and
        bending resulting from temperature or pressure change is totally prevented.
        
        """
        # 
        #
        #sigma_h =  hoop_stress(self, tnom)
        #
        # Unless a strain-based design approach is adopted 
        # (see 6.4.3), equivalent stresses should be evaluated
        # using the von Mises stress criterion shown in equation(6)
        # (6)
        sigma_e = math.sqrt(sigma_h**2 + sigma_L**2
                            - sigma_h * sigma_L
                            + 3*tau**2)
        #
        if output:
            print (" ")
            print ("Functional Stress Calculation")        
            print ("Von Mises Equivalent Stress = {:2.3f}".format(sigma_e))
        #
        #self.SFeq = (self.fd * self.sigma_y / self.sigma_e)
        #print (("Equiv Stress Usage Factor: " +"%2.3f")%(self.SFeq))
        #
        return sigma_e
    #
    # 
    def stressFEA(self, tnom=None, output=False):
        """
        tnom : Nominal wall thickness (m)
        
        Stress based on FE analysis inlcuding Temp & Int Pressure
        """
        #
        if tnom:
            _tnomcorr = float(tnom)
        else:
            _tnomcorr = self.tnom - self.tcorr
        #
        # NOTE Nominal wall thickness may be used in the evaluation.
        #
        _Sigmah = self.sigma_h * (self.tnom/_tnomcorr)
        #
        self.sigma_L = self.sigma_L * (self.tnom/_tnomcorr)
        #
        self.tau = self.tau * (self.tnom/_tnomcorr)
        #
        #
        # Unless a strain-based design approach is adopted 
        # (see 6.4.3), equivalent stresses should be evaluated
        # using the von Mises stress criterion shown in equation(6)
        # (6)
        self.sigma_e = math.sqrt(_Sigmah**2 + self.sigma_L**2 -
                                (_Sigmah * self.sigma_L) + 
                                (3*self.tau**2))
        #
        # NOTE Minimum wall thickness may be used in the evaluation
        #
        self.sigma_h = self.sigma_h * (self.tnom/self.tmin)
        #
        if output:
            print (" ")
            print ("Equivalent Stress")
            print ("T fact: {:}".format(self.tnom/_tnomcorr))
            print ("Hoop Stress  = {:2.3f}".format(_Sigmah))
            print ("Longitudinal Stress = {:2.3f}".format(self.sigma_L))
            print ("Tau  = {:2.3f}".format(self.tau))
            print ("Sigma e  = {:2.3f}".format(self.sigma_e))
            print ("Hoop Stress  = {:2.3f}".format(self.sigma_h))
        #
        return self.sigma_e
        #
    #
    #
    def equivalent_strain(self, epsilon_pL, epsilon_ph, epsilon_pr, output=False):
        """
        
        epsilon_pL : Principal longitudinal plastic strain
        epsilon_ph : Principal circunferential (hoop) strain
        epsilon_pr : Radial plastic strain
        
        6.4.3 Strain-based design
        
        The limit on equivalent stress recommended in 6.4.2.4 
        may be replaced by a limit on allowable strain,
        provided that all the following conditions are met:
        
        a) The allowable hoop stress criterion (see 6.4.2.1 and
           6.4.2.2) is met.
        
        b) Under the maximum operating temperature and pressure,
           the plastic component of the equivalent strain does 
           not exceed 0.005 (0.5 %).
           
        c) The reference state for zero strain is the as-built 
          state (after pressure test). The plastic component 
          of the equivalent uniaxial tensile strain should be
          calculated using equation (8).
          This analysis can be performed conservatively by assuming
          a linearly elastic - perfectly plastic stress/strain curve.
          Other, more realistic stress/strain curves may be used. 
          However, it is essential that the assumed curve is validated 
          as being conservative by material stress/strain curves from
          the manufactured pipe.
        
        d) Any plastic deformation occurs only when the pipeline is
           first raised to its maximum operating pressure and temperature,
           but not during subsequent cycles of depressurization, reduction
           in temperature to the minimum operating temperature, or return 
           to the maximum operating pressure and temperature. This should 
           be determined via analytical methods or an appropriate finite 
           element analysis. The analysis should include an estimate of 
           the operational cycles that the pipeline is likely to 
           experience during the operational lifetime.
        
        e) The Do/tnom ratio does not exceed 60.
        
        f) Axial or angular misalignment at welds is maintained within defined
           tolerances.
        
        g) A fracture analysis is carried out in accordance with 6.4.5.
        
        h) A fatigue analysis is carried in accordance with 6.4.6
        
        i) The weld metal yield stress matches or overmatches the longitudinal yield
           stress of the pipe.
           
        j) For welds where allowable defect sizes are based on an ECA,
           UT supplements radiographic testing, unless automated ultrasonic
           testing (AUT) is performed.
        
        k) Additional limit states are analysed as follows:
           1) bending failure resulting from application of a moment in
              excess of the moment capacity of the pipe;
           2) ovalization - distortion of the pipe wall associated with
              bending to high strain levels (see 6.4.4.2 and Annex G);
           3) local buckling (see 6.4.4.1 and Annex G);
           4) global buckling - lateral or upheaval buckling due to 
              overall axial compression (see 6.4.4.1 and Annex G).
        
        Plastic deformation reduces pipeline flexural rigidity; this effect can reduce
        resistance to upheaval buckling and should be checked if upheaval buckling
        might occur. The effects of strain localization should be taken into account in
        the strain-based design.
        
        NOTE Strain localization is associated with discontinuities in stiffness of the pipeline
            (bending or axial) and can therefore develop in the following locations:
          - changes in wall thickness;
          - buckle arrestor locations;
          - locally thinned regions, e.g. due to corrosion;
          - field joints and coatings;
          - welds, due to undermatching of the strength of the weld.
        
        """
        #
        if self.Do/self.tnom > 60:
            print ("Do/tnom > 60  --> Strain Design Not Applicable")
            sys.exit()
        
        else:
            epsilon_p = math.sqrt((2.0/3.0) * 
                                  (epsilon_pL**2 +
                                   epsilon_ph**2 +
                                   epsilon_pr**2))
        #
        if output:
            print (" ")
            print ("Strain Design")
            print ("Equivalent Uniaxial Tensile Strain : {:2.5f}".format(epsilon_p))
        #
        return epsilon_P
    #
    #
    def longitudinal_stress(self, tcheck, 
                            sigma_h, delta_T, P,
                            Mb=None, restrained=False,
                            output=False):
        """
        restrained : False/True
        temperature_pressure : False/True

        6.4.2.3 Longitudinal stress

        The total longitudinal stress should be the sum of the longitudinal stress arising
        from pressure, bending, temperature, mass, other sustained loadings and
        occasional loadings (see Annex G).

        NOTE A - pipeline is deemed to be totally restrained when axial movement and
        bending resulting from temperature or pressure change is totally prevented.
        """
        #
        checkout = []
        checkout.append("Hoop Stress = {:2.3f}".format(sigma_h))
        #   --------------------------------------
        #
        Ai, Ae, CSA, Anomfab, I, Ze, Ip =  pipe_section(self.Do, self.tnom, tcheck)
        #
        sigma_b = Mb / Ze
        
        # Restrained
        if restrained:
            # Temperature induced part of Force
            # - As*E*AlphaT*DeltaT
            Stemp = 0
            if delta_T:
                Stemp = self.E * self.alpha_T * delta_T 
                checkout.append("Temperature Long Stress = {: 2.3f}".format(Stemp))
            #
            # Pressure induced part of Force
            # -DeltaPi*Ai*(1-2*poisson)
            _Sp = 0
            if P:
                #v = self.Poisson #Ai * (1 - 2 * self.Poisson)
                if self.Do/tcheck > 20:
                    _Sp = sigma_h * self.Poisson 
                else:
                    #p = (Pi - Po)
                    _Sp = (sigma_h - P) * self.Poisson 
                
                #temperature_pressure = True
                checkout.append("Pressure Long Stress = {:2.3f}".format(_Sp))
            #
            #
            _SL = _Sp - Stemp
        
        # Unrestrained
        else:
            # Bending Stress due to temperature, weight
            # of pipe contents, insulation, snow and ice,
            # wind or earthquake is calculated by the 
            # following equation:
            #
            if self.Do/tcheck > 20:
                k = self.Do / (self.Do - 2 * tcheck)
                
            else:
                k = 1.0
            #
            _SL = sigma_h/(k**2 + 1) + sigma_b
            #
            checkout.append('k = {:2.3f}'.format(k))
            checkout.append("SL = {:2.3f}".format(_SL))
        #
        sigma_L =   _SL 
        #
        if output:
            print("")
            print("Longitudinal Stress")      
            
            for line in checkout:
                print("{:}".format(line))
                
            print("Longitudinal Stress = {:2.3f}".format(sigma_L))
        #
        return sigma_L, sigma_b
    #
    def enviromental_stress(self, Px, Vip, Vop, BMip, BMop, BMt,
                            tcheck=None, output=False):
        """
        """
        if not tcheck:
            tcheck = self.tnom
        
        Ai, Ae, CSA, Anomfab, I, Ze, Ip =  pipe_section(self.Do, self.tnom, tcheck)
        
        M = math.sqrt(BMip**2 + BMop**2)
        
        sigma_I = abs(Px / Anomfab) + (M / Ze)
        
        sigma_II = 0
        
        V = math.sqrt(Vip**2 + Vop**2)
        
        sigma_tau =  BMt / (2*Ip) + V / CSA
        
        sigma_e = math.sqrt(sigma_I**2 + - sigma_I * sigma_II 
                            + sigma_II**2 + sigma_tau**2)
        
        if output:
            print("")
            print("Enviromental Stress Calculation")
            print("Equivalent Stress = {:2.3f}".format(sigma_e))
        
        return sigma_e    
    #
    def shear_stress(self, tcheck, T, Fs, output=False):
        """
        T  : torque
        Fs : Shear Force
        
        
        6.4.2.4 Shear stree
        
        The shear stress, tau, should be calculated from the torque and shear force applied
        to the pipeline using equation (9)(7).
        """
        #
        Ai, Ae, CSA, Anomfab, I, Ze, Ip =  pipe_section(self.Do, self.tnom, tcheck)
        #
        tau = ((T / (2.0 * Ze)) + (2 * Fs / Anomfab))
        if output:
            print("")
            print("Shear stress")
            print("Tau  = {:2.3f}".format(tau))
        #
        return tau
    #    
    #    
    #
    # ===================================================
    # 6.4.4 Buckling
    #
    # 6.4.4.1 General
    # 
    # The following buckling modes should be taken into account:
    #
    # a) local buckling of the pipe wall due to external pressure,
    #    axial tension or compression, bending and torsion or a 
    #    combination of these loads (see G.1);
    #
    # NOTE 1 For fabrication processes which introduce cold 
    #        deformations giving different strength in tension
    #        and compression, a fabrication factor, Mufab, should 
    #        be determined. Guidance on the selection of a 
    #        suitable fabrication factor is given in 
    #        DNV-OS-F101:2000, Section 5.
    # 
    # b) propagation buckling due to external pressure, following
    #    the formation of local buckles or localized damage (see G.2);
    # 
    # c) restrained pipe buckling due to axial compressive forces, 
    #    induced by high operating temperatures and pressures. This
    #    can take the form of horizontal snaking of pipelines, or 
    #    vertical upheaval of trenched or buried pipelines (see G.3).
    # 
    # NOTE 2 The formulae given in Annex G define one approach to 
    #        analysis. Alternative approaches are available and may
    #        be used where justified.
    # 
    # In all buckling analyses, the nominal wall thickness should be used.
    #
    #
    # 6.4.4.2 Ovality
    #
    # Ovality, or out-of-roundness, of pipes or a section of pipeline
    # that could cause buckling or interference with pigging operations 
    # should be avoided.
    #
    # NOTE 1 In some situations, where loading is dominated by bending,
    #        buckling might not occur but unacceptable levels of 
    #        ovalization can result.
    # 
    # NOTE 2 Ovalization may be calculated in accordance with G.4 in 
    #        the absence of a more rigorous evaluation.
    #
    #
    # ===================================================
    #
    # Annex G (informative)
    # Buckling
    #
    # G.1 Local buckling
    #
    # -------------------------
    # G.1.1 General
    # NOTE 1 Local buckling of the pipe wall can be avoided if the various loads to
    # which the pipe is subjected are less than the characteristic values in G.1.2 to G.1.7.
    #
    # NOTE 2 Guidance on buckling is given in DNV-OS-F101.
    # Where the concrete cladding is thick enough and reinforced to provide a
    # structural member conforming to BS 6349-1-4 and BS EN 1992-1-1, it may be
    # used to provide support against buckling provided that appropriate justification
    # is given.
    #
     # -------------------------
    #
    # G.1 Local Buckling
    def external_pressure(self, fo, sigma_u,
                          root_search='FAST', factor=1.0,
                          output=False):
        """
        factor = 2.0 (Legacy from BS8010-3-1993)
        
        
        G.1.2 External Pressure
        
        The characteristic value, Pc, that causes collapse
        when the external pressure is acting alone, can be
        calculated using equations (G.1) to (G.4).
        
        """
        #
        #
        # (G.4) Maximun Ovality
        # self.fo = (self.Dmax - Dmin)/self.Do
        #
        # (G.3) Yield Pressure
        self.Py = ((2*self.sigma_y)*(self.tnom/self.Do))
        #
        # (G.2) Elastic Critical Pressure
        self.Pe = (((2*self.E)/(1.0 - self.Poisson**2)) * 
                   (self.tnom/self.Do)**3)
        #
        #self.Pe = (((2*self.E)/(1.0 - self.Poisson**2)) * 
        #           (self.tnom/self.Do)**3)        
        #
        # (G.1)
        # Define Funtion 
        # TODO : Check why 2?
        def f(x): return (((x/self.Pe) - 1.0)*
                          ((x/self.Py)**2 - 1.0) - 
                          factor * (x/self.Py)*(fo*self.Do/self.tnom))
        #
        # Find first root (Note that it may not be the minimum
        # use 'full' to find all roots, but process is slow)
        #self.Search = 'FAST'
        # Find Pc
        Pc = GoalSeeker(f, sigma_u*10, root_search)
        #
        if output:
            print("")
            print("G.1.2 External Pressure")            
            print("Maximun Ovality (fo) = {:2.3f}".format(fo))
            print("Yield Pressure (Py) = {:2.3f}".format(self.Py))
            print("Elastic Critical Pressure (Pe) = {:2.3f}".format(self.Pe))
            #print("Po max  = {:2.3f}".format(Po))
            print("Pc  = {:2.3f}".format(Pc))
            #print("URhp  = {:2.3f}".format(abs(Po/Pc)))
        #
        return Pc
    # 
    def axial_compression(self, Fx, D, t, output=False):
        """
        G.1.3 Axial Compression
        
        If D/tnom is less than 60, local buckling under
        axial compression does not occur until the mean 
        axial compression load, Fxc, reaches the yield
        load, Fy, i.e. as shown in equation (G.5).
        """
        #
        _Fy = math.pi*(D - t) * t * self.sigma_y
        #
        self.Fxc = _Fy
        #
        if output:
            print ("")
            print ("G.1.3 Axial Compresion")
            print ("Fxc  = {:2.3f}".format(self.Fxc))
            print ("URaxial  = {:2.3f}".format(Fx/self.Fxc))
        #
        return self.Fxc
    #
    def bending(self, Mb, D, t, output=False):
        """
        D
        t
        sigma_y :
        
        
        G.1.4 Bending
        
        The characteristic bending moment value, Mc, 
        required to cause buckling when bending moments
        are acting alone, can be obtained using equations
        (G.6) and (G.7).
        """
        #
        # (G.7)
        _Mp = ((D - t)**2 * (t * self.sigma_y))
        # (G.6)
        Mc = ((1.0 - 0.0024 * (D/t)) * _Mp)
        
        # The characteristic bending strain, Epsilon_bc, at  
        # which buckling due to bending moments acting alone 
        # occurs, can be obtained using equation (G.8).
        
        # (G.8)
        epsilon_bc = 15.0 * (t/D)**2
        #
        if output:
            print ("")
            print ("G.1.4 Bending")
            print ("Mc  = {:2.3f}".format(Mc))
            print ("URbm  = {:2.3f}".format(Mb/Mc))
            #print (("URBMop  = "+"%2.3f")%(self.Mz/self.Mc))            
            print ("Epsilon bc  = {:2.3f}".format(epsilon_bc))
        #
        return Mc, epsilon_bc
        #
    #
    def torsion(self, T,  D, t, output=False):
        """
        G.1.5 Torsion
        
        The characteristic value, Tauc, that causes buckling
        when torsion is acting alone, can be obtained using
        equations (G.9) to (G.13).
        
        """
        #
        # (G.12)
        _tau_y = self.sigma_y / math.sqrt(3.0)
        # (G.13)
        _alpha_t = (self.E/ _tau_y)*(t/D)**(3.0/2.0)
        # (G.9)
        if _alpha_t < 1.5:
            self.tau_c = (0.542 * _alpha_t)*_tau_y
        # (G.11)
        elif _alpha_t > 9.0:
            self.tau_c = _tau_y
        # (G.10)
        else:
            self.tau_c = (0.813 + 
                          0.068 * math.sqrt(_alpha_t - 1.50))
        #
        _fvt = abs((T * D)/(2 * self.Ip))
        #
        if output:
            print("")
            print("G.1.5 Torsion")
            print ("Tau y  = {:2.3f}".format(_tau_y))
            print ("Alpha t = {:2.3f}".format(_alpha_t))
            print ("Tau c = {:2.3f}".format(self.tau_c))
            print ("Torsion  = {:2.3f}".format(_fvt))
            print ("URt  = {:2.3f}".format(_fvt/self.tau_c))
        #
        return self.tau_c
    #
    def load_combination(self, Mb, Mc, Fx, Fxc, fo, Po, Pc, output=False):
        """
        G.1.6 Load Combinations
        
        The maximum external overpressure, P, in the 
        presence of compressive axial force, Fx, and/or
        bending moment, M, when fo is less than 0.05 (5%),
        can be calculated using equation (G.14), where:
        - Gamma is calculated using equation (G.15);
        - Sigmahb is calculated using equation (G.16);
        - Sigmahcr is calculated using equation (G.17) 
          or equation (G.18) as appropriate.
        """
        #
        #
        if fo <= 0.05:
            #
            _sigma_hE = self.E * (self.tnom/(self.Do - self.tnom))**2
            
            # (G.16)
            _sigma_hb = Po * self.Do / (2*self.tnom)
            
            # (G.17)
            if _sigma_hE <= (2.0/3.0) * self.sigma_y :
                _sigma_hcr = _sigma_hE
            
            # (G.18)
            else:
                _sigma_hcr = (self.sigma_y * 
                             (1 - ((1.0/3.0)*(2*self.sigma_y / 
                                              (3*_sigma_hE)))))
            
            # (G.15)
            _gamma = (1 + 300.0 * 
                      (self.tnom/self.Do)*(_sigma_hb/_sigma_hcr))
            
            #
            try:
                pressure = abs(Po/Pc)
            except ZeroDivisionError:
                pressure = 0
            
            if Fx:
                axial = abs(Fx/Fxc)
            else:
                axial = 0
            
            try:
                bending = abs(Mb/Mc)
            
            except ZeroDivisionError:
                bending = 0
            
            self.UR_lc = (((bending + axial)**_gamma) + pressure)
            #print("")
            #print('bending + axial = {:}'.format((bending + axial)**_gamma))
            #print('pressure        = {:}'.format(pressure))
        
        else:
            print ("fo > 0.05 (5%)")
            print ("Section not applicable")
            sys.exit()
        #
        if output:
            print("")
            print("G.1.6 Load Combination")
            print("Sigma hE  = {:2.3f}".format(_sigma_hE))
            print("Sigma hb  = {:2.3f}".format(_sigma_hb))
            print("Gamma  = {:2.3f}".format(_gamma))
            print("UR load comb  = {:1.6f}".format(self.UR_lc))
        #
        return self.UR_lc
    # 
    def strain_criteria(self, P, Pc, epsilon_bc, output=False):
        """
        G.1.7 Strain criteria
        
        The bending strain, Epsilonb, required to cause
        buckling, in the presence of external overpressure,
        P, can be calculated using equation (G.19).
        """
        #
        # Values for Epsilonbc and Pc can be obtained from 
        # equations (G.8) and (G.1) respectively.        
        # (G.19)
        epsilon_b = (1 - (P/Pc))*epsilon_bc
        #
        if output:
            print("")
            print("G.1.7 Strain criteria")
            print("Epsilon b = {:2.3f}".format(epsilon_b))
        #
        return epsilon_b
    #         
    # -------------------------
    # G.2  Propagation buckling
    def propagation_buckling(self, Pi, Po, t, output=False):
        """
        G.2  Propagation Buckling
        
        The potential for a pipeline to propagate local 
        buckles is dependent on the external overpressure, 
        P, and its relationship with the propagation pressure Pp.
        """
        #
        # The propagation pressure, Pp, can be calculated using 
        # equation (G.21).
        # (G.21)
        self.Pp = (10.70 * self.sigma_y * (t/self.Do)**2.25)
        #
        # The external overpressure, P, can be calculated using
        # equation (G.20).
        self.P = Po - Pi
        # Check if external pressure is not significant
        if self.P < 0:
            self.P = Po
        #
        # If P is less than Pp, then, even though it is possible
        # for the pipe to develop a local buckle, the buckle will
        # not propagate.
        #
        # If P is greater than or equal to Pp and a local buckle
        # or local damage has occurred, then the pipeline is likely
        # to undergo propagation buckling. It can be advisable to 
        # provide buckle arresters
        #
        if output:
            print("")
            print("G.2 Propagation Buckling")
            print("Pp = {:2.3f}".format(self.Pp))
            print("P  = {:2.3f}".format(self.P))
            print("URpp  = {:2.3f}".format(abs(self.P/self.Pp)))
        #
        return self.P, self.Pp
    #
    # -------------------------
    #  G.3 Upheaval buckling
    def upheaval_buckling(self, output=False):
        """
        
        G.3 Upheaval Buckling
        
        Two major factors contribute towards the upheaval of subsea pipelines: the
        effective axial driving force, arising from the internal pressure and temperature,
        and the presence of vertical out-of-straightness (OOS) in the seabed profile.
        
        """
        
        # No yet Included
        print ("No yet Included")
        #
    #
    def buckling_force(self, tcheck, Pi, d, delta_T, SF=1.2, output=False):
        """
        """
        #
        Ai, Ae, CSA, Anomfab, I, Ze =  pipe_section(self.Do, self.tnom, tcheck)
        
        # pressure difference
        self.delta_P = Pi - self.rho_w * d * self.g
        
        # fully build-in in Force
        self.F = []
        for _deltaT in delta_T:
            self.F.append((self.E * self.alpha_T * CSA * _deltaT 
                           + Ai * (1.0 - 2.0 * self.Poisson) * self.delta_P) * SF)
        print(" F full = {:}".format(self.F))
        #
        
        # factored temperature difference
        self.delta_Tfactor = []
        for F in self.F:
            self.delta_Tfactor.append((F - Ai * (1.0 - 2.0 * self.Poisson) * self.delta_P) 
                                      / (self.E * self.alpha_T * CSA))
        
        if output:
            print(" Factored temperature difference = {:}".format(self.delta_Tfactor))
        #
        return self.F
    #
    def imperfection_shape(self, tcheck, w_ins):
        """
        """
        Ai, Ae, CSA, Anomfab, I, Ze =  pipe_section(self.Do, self.tnom, tcheck)
        
        L = [(hi * 72.0 * self.E * I / w_ins)**0.25 for hi in self.h]
        
        # Define height and half length between inflection points
        self.hr = [hi/2.4545 for hi in self.h]
        
        self.Lr = [Li/3.0 for Li in L]
        #print(' Lr : ', self.Lr)
        
        return self.Lr, self.hr
    #
    def ovality_bending_effect(self, D, t, j, output=False):
        """
        D 
        t
        j
        """
        #
        f = [(i+1)/100.0 for i in range(j)]
        # 
        I = [ellipse_section(D, t, fo) for fo in f]
        
        EI = [self.E * Ii for Ii in I]
        
        if output:
            print('EIo = {:1.4e} N mm^2'.format(EI[0]))
        
        return EI
    #
    def download_response(self, F, Lr, EI, hr, w_op, rho_sub, f, Dc):
        """
        F : Buckling design force
        Lr : Inflection points height
        EIj : 
        """
        Phi_L = []
        for i in range(len(Lr)):
            Phi_L.append([])
            for j in range(len(EI)):
                Phi_L[i].append(Lr[i] * math.sqrt(F / EI[j]))
        
        Phi_W = []
        for i in range(len(Lr)):
            Phi_W.append([])
            for j in range(len(EI)):
                
                if Phi_L[i][j] < 2.75:
                    Phi_W[i].append(0.180)
                
                else:
                    if Phi_L[i][j] < 10.0:
                        calc = (2.402 / (Phi_L[i][j])**2) - (7.874 / (Phi_L[i][j])**4)
                        Phi_W[i].append(calc)
                    
                    else:
                        Phi_W[i].append(0.023)
        # Required total download
        Wreq = []
        for i in range(len(Lr)):
            Wreq.append([])
            for j in range(len(EI)):
                Wreq[i].append(Phi_W[i][j] * F**2 * hr[i] / EI[j])
        # require cover download
        Wcover = []
        for i in range(len(Lr)):
            Wcover.append([])
            for j in range(len(EI)):
                Wcover[i].append(Wreq[i][j] - w_op)
        # require cover height
        Hreq = []
        for i in range(len(Lr)):
            Hreq.append([])
            for j in range(len(EI)):
                calc = ((math.sqrt(rho_sub**2 * Dc**2 
                                   + 4 * rho_sub * f * Wcover[i][j] ) 
                         - rho_sub * Dc) / (2.0 * rho_sub * f))
                Hreq[i].append(calc)
        #
        return Hreq
    #
    # -------------------------
    # G.4 Ovalization
    def ovalization(self, P, fo, epsilon_b, tnom=None, output=False):
        """
        G.4  Ovalization
        
        The total ovalization, f, of a pipe due to
        the combined effects of unidirectional bending
        and external pressure can be calculated 
        using equations (G.25) to (G.27).
        
        Values for Pe and fo can be obtained from 
        equations (G.2) and (G.4) respectively.
        
        NOTE If cyclic or reversed bending is applied,
        the resulting ovalization can be considerably 
        greater than that predicted by the equation.
        """
        #
        t = self.tnom
        if tnom:
            t = tnom
        # (G.27)
        self.Cf = 0.120 * (1 + self.Do/(120 * t))
        
        # (G.26)
        self.Cp = 1.0 / (1.0 - P/self.Pe)
        
        # (G.25)
        self.f = (self.Cp 
                  * (self.Cf * (epsilon_b * self.Do /t)**2 + fo))
        
        if output:
            print("")
            print("G.4 Ovalizacion")
            print("Cp  = {:2.3f}".format(self.Cp))
            print("Cf  = {:2.3f}".format(self.Cf))            
            print("f   = {:2.3f}".format(self.f))
        #
        return self.Cf, self.Cp, self.f
    #
    #
    # -------------------------
    #    
#
