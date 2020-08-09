# Copyright (c) 2015-2016 steelpy

# Python stdlib imports
import math
import datetime

# package imports
from steelpy.math.rootsearch import GoalSeeker
from steelpy.codes.piping.ASME_B313 import ASME
from steelpy.codes.piping.DNV_F101 import DNV
from steelpy.wave.theory import WaveStokes5
#

#
#
def pipe_section(Do, tnom, tcheck):
    """
    """
    #
    _Dinom = (Do - 2 * tnom)
    
    # Internal Area
    Ai = (math.pi*_Dinom**2/4.0)    
    
    # External Area
    Ae = (math.pi*Do**2/4.0)    
    
    # Pipe section area (nominal thickness)
    CSA = ((math.pi / 4.)*
           (Do**2 -(Do - 2*tnom)**2))
    
    # Pipe section area (Tnom - Tcorr thickness)
    Anomfab = ((math.pi / 4.0)*
               (Do**2 -(Do - 2 * tcheck)**2))
    
    #
    Ze = ((math.pi / 64.0) * 
          (Do**4 - (Do - 2 * tcheck)**4) /
          (Do / 2.0))
    #
    return Ai, Ae, CSA, Anomfab, Ze
#
#-------------------------------------------------
#                  PD8010 Section
#-------------------------------------------------
#
class PD8010_2015_1:
    """
    PD 8010-1 & 2: 2015
    """
    #
    def __init__(self):
        """
        """
        pass
    #
    # 
    # ===================================================
    # Section Properties
    # ===================================================
    #     
    def section_input(self, Do, tnom, tmin): 
        """
        Do
        tnom
        tmin
        """
        #
        self.Do = float(Do)
        self.tnom = float(tnom)
        self.tmin = (tmin)
        #
        print (' ')
        print ('Section Properties ')
        print ('Do =', self.Do)
        print ('Tnom =',self.tnom)
        print ('Tmin =',self.tmin)
        print (' ')
        #
        #     
        #-------------------------------------------------
        #   Cross-Sectional Area
        #
        self.A = ((math.pi / 4.)*
                  (self.Do**2 -(self.Do - 2*self.tmin)**2))
        #   
        #
        #-------------------------------------------------
        #               Section Properties
        #-------------------------------------------------
        #
        #
        #   Second Moment of Area about Mayor Axis
        #   --------------------------------------
        #
        self.I = ((math.pi / 64.0)*
                  (self.Do**4 - (self.Do - 2*self.tmin)**4))
        # print ('circular Second Moment',self.Iip)
        # 
        #   
        #
        #
        #-------------------------------------------------
        #   Plastic Modulus about Mayor Axis
        #    def Zx(self):
        #
        self.Zp = ((self.Do**3 - 
                    (self.Do - 2 * self.tmin)**3) / 6.0 )
        # print ('circular Plastic Modulus',self.Zpip)
        #
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
        #
        #-------------------------------------------------
        #   Polar Moment of Inertia
        self.Ip = ((math.pi / 32.0) *
                   (self.Do**4 - 
                    (self.Do - 2*self.tmin)**4))
        # print ('Polar Moment of Inertia Ip:',self.Ip)
        #
        #
        #-------------------------------------------------
        #  Mass
        #self.mass = (self.A * self.Rhos)/1000**2
        #
    #
    def material(self, sigma_y, E, Poisson, alpha_T):
        """
        """
        self.sigma_y = sigma_y
        self.E = E
        self.Poisson = Poisson
        print ('Material Properties ')
        print ('Fy =',self.sigma_y)
        #print ('Fu =',self.Sigmau)
        print ('E =',self.E)
        #print ('G =',self.G)
        self.alpha_T = alpha_T
    #
    #
    # ===================================================
    #    PD8010-1:2004 Part 1 : Steel Pipelines on Land
    # ===================================================
    # 
    def substances_categorization(self, category):
        """
        category : A - Typically non-flammable water-based fluids
                   B - Flammable and/or toxic fluids that are liquids
                       at ambient temperature and at atmospheric
                       pressure conditions
                   C - Non-flammable fluids that are non-toxic
                       gases at ambient temperature and
                       atmospheric pressure conditions
                   D - Non-toxic, single-phase natural gas
                   C - Flammable and/or toxic fluids that are gases
                       at ambient temperature and atmospheric
                       pressure conditions and are conveyed as
                       gases and/or liquids
                       Mixtures of petroleum or chemical
                       substances, having a Reid vapour pressure
                       greater than 31 kPa absolute
        
        
        5.2 Categorization of fluids
        
        The substances to be transported should be categorized, with respect to
        hazard potential in respect of public safety, into one of the five categories given
        in Table 1.
        
        NOTE 1 - Attention is drawn to the Pipelines Safety Regulations 1996 [15] for the
        definition and classification of dangerous fluids (hazardous substances).
        Gases or liquids not specifically included by name should be classified in the
        category containing substances most closely similar in hazard potential to those
        quoted. If the category is not clear, the more hazardous category should be
        assumed.
        
        NOTE 2 - Guidance is given in a number of publications including
        HSE publication L82 [22] and ICE publication Nomenclature for hazard and risk
        assessment [23].
        
        NOTE 3 - Additional guidance on CO2 pipelines is given in DNV-RP-J202.
        
        """
        #
        a = 0.72
        #
        #if category == "B":
        #    a = 0.72
        #
        if category == "C":
            a = 0.30
        #
        return a
    #
    def substance_factor(self, substance, Q=None):
        """
        substance : ammonia
                    carbon dioxide dense/gas phase
                    ethylene
                    hydrogen
                    liquid petroleum gas
                    natural gas liquid
                    user
        
        Q         : Substance Factor
        
        
        Table 3 Substance factors
        """
        
        if 'ammonia' in substance:
            self.Q = 2.50
        
        elif 'carbon' in substance:
            self.Q = 2.0
            if 'gas' in substance:
                self.Q = 1.0
        
        elif 'ethylene' in substance:
            self.Q = 0.80
        
        elif 'hydrogen' in substance:
            self.Q = 0.45
        
        elif 'petroleum' in substance:
            self.Q = 1.0
    
        elif 'natural' in substance:
            self.Q = 1.25
        
        else:
            self.Q = Q
    #
    def minimum_routeing_distance(self):
        """
        c) For pipelines having a design factor not exceeding 0.72, the minimum
           distance for routeing purposes between the pipeline and occupied
           buildings, Y, should be determined using equation (1).
        """
        self.Y = (self.Q * 
                  ((self.Do**2 / 32000.0) + (self.Do / 160.0) + 11.0) 
                  * (self.p / 3.20 + 1.40))
    #
    def bends(self, n):
        """
        n :
        
        
        6.2.2.3 Bends
        
        Changes in direction may be made by bending pipe or installing factory-made
        bends or elbows. All bends should be free from buckling, cracks or other
        evidence of mechanical damage. The nominal internal diameter of a bend
        should not be reduced in ovality by more than 2.5% at any point around the
        bend. Sufficient tangent lengths should be left at each end of a bend to ensure
        good alignment and to facilitate welding. Pipes bent cold should not contain a
        girth weld within the bent section.
        """
        
        self.thin = 50.0 / (n + 1)
    #
    # 
    def limits_calculated_stress(self, sigma_h, Sah, sigma_e, Sae, design_method):
        """
        sigma_h
        Sah
        sigma_e
        Sae
        
        6.4.3 Limits of calculated stress
        

        """
        #
        print (" ")
        print ("Allowable Stress")
        #
        URhoop = (sigma_h / Sah )
        print (("URhoop : {:2.3f}").format(URhoop))
        if 'stress' in design_method:
            UReq = (sigma_e / Sae)
        else:
            UReq = self.epsilon_p / 0.0050
            
            if self.epsilon_p > 0.0050:
                print (("Fail  Epsilonp > 0.0010:" +"%2.3f")%(self.epsilon_p))
                
            else:
                print (("Pass  Epsilonp < 0.0010:" +"%2.3f")%(self.epsilon_p))            
        #
        print (("UReq   : {:2.3f}").format(UReq))
        #
        return URhoop, UReq
    #
    def allowable_hoop_stress(self, a, pipe_history=False):
        """
        a : 
        pipe_history : unknown = False
                         known = True
        
        
        6.4.3.1 Allowable hoop stress
        
        The allowable hoop stress (Sah) should be 
        calculated using equation (11).
        
        The weld joint factor, e, should be 1.0 for pipe conforming to
        BS EN ISO 3183:2012 and/or API 5L:2012 when supplied as seamless,
        longitudinally welded or spirally welded pipe. If the pipe history is unknown,
        the weld joint factor e should not exceed 0.60 for pipe of 0.114 m outside
        diameter or smaller, or 0.80 for pipe larger than 0.114 m outside diameter.
        
        NOTE - The effect of temperature de-rating on the SMYS of carbon steel is included
        in the design factors for temperatures up to 120C.
        """
        e = 1.0
        if not pipe_history :
            #
            if self.Do > 114.0:
                e = 0.80
            
            else:
                e = 0.60
        #
        #
        Sah = (a * e * self.sigma_y)
        fd_hs = (a * e)
        print ("Allowable hoop stress Sah = ", Sah)
        #
        return Sah, fd_hs
    #
    def allowable_equivalent_stress(self):
        """
        sigma_y : 
        
        6.4.3.2 Allowable Equivalent Stress
        
        NOTE Further guidance is given in BS EN 13480, ASME B31.3 and ASME B31.8.
        """
        #
        Sae = 0.90 * self.sigma_y
        fd =  0.90
        print ("Allowable equivalent stress Sah = {:}".format(Sae))
        
        return Sae, fd
    #
    def anchor_blocks(self):
        """
        
        6.17 Anchor blocks
        
        The design of anchor blocks to prevent axial movement of a pipeline should
        take into account the pipeline expansion force and any pipe-to-soil
        friction-preventing movement. The axial compressive force necessary to restrain
        a pipeline should be calculated using equation (14) for thin wall or
        equation (15) for thick wall.
        """
        # Thin wall
        self.F = (A * (self.E * self.alpha * (self.T2 - self.T1) 
                       + 0.50*self.sigma_hl 
                       - self.Poisson * self.sigma_hl))
        # Thick wall
        self.F = (A * (self.E * self.alpha * (self.T2 - self.T1) 
                       + (self.sigma_hl / (K**2 + 1))
                       - self.Poisson * (self.sigma_hl - P)))
        #
    #
    def hydrostatic_test(self, D, t, Tf):
        """
        D  : pipe diametre
        t  : wall thickness
        Tf : Temperature factor change
        
        
        11.7.2 Method of assessment for hydrostatic test
        
        The relationship between pressure and temperature should be calculated in
        accordance with equation (20).
        """
        self.delta_p = 0.10 * (264.70 * Tf / (D/t + 100))
    #
    def hoop_stress(self, Pi, Po,
                    Do=None, tmin=None):
        """
        Do   :
        tmin : minimimum wall thickness (m)
        Pi : Internal pressure
        Po : External pressure

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

        _Dih = (_Do - 2*_tmin)

        Pi = float(Pi)
        Po = float(Po)
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
            sigma_h = (Pi - Po)*(_Do/(2 * _tmin))
        #
        else:
            # Equation (5) (thick wall) should be used when the ratio
            # of Do/tmin is less than or equal to 20.
            sigma_h = ((Pi - Po)*
                       ((_Do**2 + _Dih**2) / 
                        (_Do**2 - _Dih**2)))
            #
        #
        print (" ")
        print ("Hoop Stress", _tmin)
        print (("Sigmah  = "+"%2.3f")%(sigma_h))
        #
        # self.sigma_h = sigma_h
        #
        return sigma_h
        #
    #
    def shear_stress(self, tcheck, T, Fs):
        """
        T  : torque
        Fs : Shear Force
        
        
        6.4.2.4 Shear stree
        
        The shear stress, tau, should be calculated from the torque and shear force applied
        to the pipeline using equation (9)(7).
        """
        #
        Ai, Ae, CSA, Anomfab, Ze =  pipe_section(self.Do, self.tnom, tcheck)
        #
        tau = ((T / (2.0 * Ze)) + (2 * Fs / Anomfab))
        #
        print (("Tau  = {:2.3f}").format(tau))
        #
        return tau
    #    
    def longitudinal_stress(self, tcheck, 
                            Fx, Mb, sigma_h, 
                            T1, T2, Pi, Po,
                            restrained=False,
                            temperature_pressure=False):
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
        #   --------------------------------------
        #
        Ai, Ae, CSA, Anomfab, Ze =  pipe_section(self.Do, self.tnom, tcheck)
        #
        print (" ")
        print ("Longitudinal Stress Calculation")
        #
        # ------------------------------------------
        # Hoop Stress Component of von Mises Stress
        # S = H -DeltaPi*Ai*(1-2*poisson) - As*E*AlphaT*DeltaT
        # ------------------------------------------
        #
        # Restrained
        if restrained:
            # Temperature induced part of Force
            # - As*E*AlphaT*DeltaT
            Stemp = 0
            if T2:
                delta_T = (T2 - T1)
                Stemp = self.E * self.alpha_T * delta_T # CSA * 
                temperature_pressure = True
                print (("Temperature Long Stress = {: 2.3f}").format(Stemp))
            #
            # Pressure induced part of Force
            # -DeltaPi*Ai*(1-2*poisson)
            _Sp = 0
            if Pi:
                v = self.Poisson #Ai * (1 - 2 * self.Poisson)

                if self.Do/self.tnom > 20:
                    _Sp = sigma_h * v 
                else:
                    p = (Pi - Po)
                    _Sp = (sigma_h - p) * v 

                temperature_pressure = True
                print (("Pressure Long Stress = {: 2.3f}").format(_Sp))
            #
            #            
            _SL = _Sp - Stemp
            print("SL = {:2.3f}".format(_SL))

        # Unrestrained
        else:
            # Bending Stress due to temperature, weight
            # of pipe contents, insulation, snow and ice,
            # wind or earthquake is calculated by the 
            # following equation:
            #
            if self.Do/self.tnom > 20:
                k = Do / (Do - 2 * tcheck)

            else:
                k = 1.0
            #
            SL2 = sigma_h/(k**2 + 1) + (Mb / Ze)
            #
            print('k = {:2.3f}'.format(k))
            print ("SL = {:2.3f}".format(SL2))            
            _SL = (_Sp  + SL2)
        #
        # ------------------------------------------
        # True Wall Axial Force
        # N = S - PiAi + PeAe
        _PiAi = 0
        _PeAe = 0
        #if Pi:
            # Internal Pressure (N)
            #_PiAi = Pi * Ai
            # External Pressure (N)
            #_PeAe = Po * Ae
        #
        #print ("Pi, Pe",_PiAi, _PeAe ) 
        #
        sigma_b = Mb / Ze
        # ------------------------------------------        
        # Summing all components of longitudinal
        # normal stress:
        #self.sigma_L = self.Sb + self.Sdl + self.Sp - self.Stemp
        #
        _Fx = ( _SL + _PiAi - _PeAe + Fx )
        # Include Temperature + hoop stress + given axial load
        sigma_L =   _SL + Fx / Anomfab
        #Fx = abs(min(( _SL + _PiAi - _PeAe + Fx ),0))

        #
        #print(' Fx = {:}'.format(_Fx))
        print (("Longitudinal Stress = {:2.3f}").format(sigma_L))
        return sigma_L, sigma_b
    # 
    def equivalent_stress(self, sigma_h, sigma_L, tau):
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
        print (" ")
        print ("Equivalent Stress")
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
        print (("sigmae  = {:2.3f}").format(sigma_e))
        #
        #self.SFeq = (self.fd * self.sigma_y / self.sigma_e)
        #print (("Equiv Stress Usage Factor: " +"%2.3f")%(self.SFeq))
        #
        return sigma_e
    #    
    def expansion_flexibility(self, FS_code, SIFs_type, pipe_description,
                              flanges=0, T_= 0, r2=0, R1S=0, theta=0):
        """ 
        FS_code : flexibility & Stress Code
        SIFs_type :
        pipe_description : 

        Return :
        h  : Flexibility Characteristic
        k  : Flexibility Factor
        lo : Stress Intensification Out-Plane
        li : Stress Intensification In-Plane
        C1 : Flanges

        6.4.2.3 Expansion and flexibility

        Pipelines and piping should be designed with sufficient
        flexibility to prevent expansion or contraction causing
        excessive forces or stresses in pipe material, joints, 
        equipment, anchors or supports.

        Expansion calculations should be carried out on pipelines
        where flexibility is in doubt, and where temperature changes
        are expected. Thermal and pressure expansion or contraction
        can cause movement at termination points, changes in direction 
        or changes in size. The necessary flexibility should be provided
        if such movements are unrestrained. Account should be taken of 
        buckling forces that can be imposed on pipelines (see 6.4.4).

        The effect of restraints, such as support friction, branch 
        connections and lateral interferences should be taken into account.
        Calculations should take into account stress intensification 
        factors found to be present in components other than plain 
        straight pipe. Account should be taken of any extra flexibility
        of such components.

        NOTE In the absence of more directly applicable data, the 
        flexibility factors and stress intensification factors given in
        BS EN 13480 or ASME B31.3 may be used.

        Pipelines can be restrained so that the longitudinal movement 
        owing to thermal and pressure changes is absorbed by direct 
        axial compression or tension of the pipe. In such cases expansion
        calculations should be carried out taking into account all the 
        forces acting on the pipeline. Consideration should be given to 
        elastic instability due to longitudinal compressive forces.

        Where movement is restrained, flexibility should be provided by
        means of loops, offsets or special fittings. The total operating 
        temperature range should be taken as the difference between the
        maximum and minimum metal temperatures for the operating cycle 
        under consideration and should be used in calculating stresses in 
        loops, bends and offsets.

        The temperature range used in the calculation of reactions on 
        anchors and equipment should be taken as the difference between
        the maximum or minimum metal temperatures and the installation 
        temperature, whichever gives the greater reaction.

        Where there is a likelihood of repeated stress changes (including
        thermal stress) giving rise to fatigue conditions, the stress 
        range and allowable number of cycles should be calculated in 
        accordance with 6.4.6. Nominal pipe wall thickness (including any
        corrosion allowance) and nominal outside diameter should be
        used for expansion and flexibility calculations.
        """
        #
        if FS_code == "ASME_B31.3":
            asme = ASME()
            asme.bend_data(pipe_description,
                           flanges, 
                           T_, 
                           r2, 
                           R1S, 
                           theta)

            h, k, lo, li, C1 = asme.AppendixD()
        #
        # BS EN 13480
        elif FS_code == "BS_EN_13480":
            #
            # Bend Flexibility
            if SIFs_type == 'BEND':
                #
                print ("No implement yet")
            #
            # Tee Flexibility
            elif SIFs_type == 'TEE':
                #
                print ("No implement yet")
            #
            # Unknow
            else:
                print ("No implement yet")
            #
            # 
            # Flanges
            if self.Flanges == 0 or self.Flanges > 2:
                self.C1 = 1.0
            #
            # Rest
            else:
                print ("No implement yet")
            #
            #
        #        
        #
        # No flexibility & Stress Factors
        else:
            # Flexibility Characteristic
            h = 1.0
            # Flexibility Factor
            k = 1.00
            # Stress Intensification
            # Out-of-Plane
            lo = 1.0
            # In-Plane
            li = 1.0
            #
            C1 = 1.0
        #
        #
        print (" ")
        print ("Expansion Flexibility ")
        print ("Based on ",FS_code)
        print ("Member Description: ", pipe_description)
        print ("Member Type: ",SIFs_type)
        print (("h  = "+"%2.3f")%(h))
        print (("k  = "+"%2.3f")%(k))
        print (("lo = "+"%2.3f")%(lo))
        print (("li = "+"%2.3f")%(li))
        print (("C1 = "+"%2.3f")%(C1))
        #
        return h, k, lo, li, C1
    #    
#     

