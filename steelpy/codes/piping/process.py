# Copyright (c) 2015-2016 steelpy

# Python stdlib imports
import math
#import datetime

# package imports
from steelpy.wave.theory import WaveStokes5
from steelpy.codes.piping.ASME_B313 import ASME
from steelpy.codes.piping.DNV_F101 import DNV
from steelpy.codes.piping.PD8010_2015_1 import PD8010_2015_1
from steelpy.codes.piping.PD8010_2015_2 import PD8010_2015_2


#
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
def functional_stress(self, t_check, PD8010, output):
    """
    """
    # External overpressure
    P = (self.Pi - self.Po)
    # hoop stress
    sigma_h = PD8010.hoop_stress(P, tmin=t_check)

    # longitudinal stress
    deltaT = self.delta_T[0]
    sigma_L, Sb = PD8010.longitudinal_stress(t_check, 
                                             sigma_h,
                                             deltaT, P,
                                             self.Mb,
                                             self.pipe_restrained,
                                             output=output)

    # shear stress
    tau = PD8010.shear_stress(t_check, self.T, self.Fs, output=output)
    
    sigma_e = PD8010.equivalent_stress(sigma_h, sigma_L, tau, output=output)   

    return sigma_e
#
#
def get_PD8010_upheaval(self):
    """
    """
    t_check = self.tnom
    self.F = PD8010.buckling_force(t_check,
                                   self.Pi, self.d,
                                   self.delta_T)

    self.Lr, self.hr = PD8010.imperfection_shape(t_check,
                                                 self.w_installation)

    if self.ovality:
        self.EI = PD8010.ovality_bending_effect(self.Do, self.tnom,
                                                self.ovality)
    else:
        self.EI = [self.E * PD8010.I for i in self.F]

    Dc = self.Do - self.coating
    self.Hreq = {}
    for x in range(len(self.delta_T)):
        _name = round(self.delta_T[x], 2)
        self.Hreq[_name] = PD8010.download_response(self.F[x], self.Lr, 
                                                    self.EI, self.hr,
                                                    self.w_operation, 
                                                    self.rho_sub, 
                                                    self.f, Dc)    
#
def get_PD8010_strain(self):
    """
    """
    # 6.4 Strength Check
    # 6.4.3.a
    #
    if "load" in self.load_type:
        t_check = self.tnom
        self.sigma_h = PD8010.hoop_stress(self, t_check)
    #
    elif "fea_stress" in self.load_type:
        PD8010.allowable_stress(self.sigma_h,
                                self.sigma_e,
                                self.design_method)
    #
    # 6.4.3.b & 6.4.3.d
    PD8010.equivalent_strain(self)
    #
    # 6.4.3.c (?)
    #
    # # 6.4.3.i.1.- Bending Failure (?)
    #
    #    
#
#
def get_PD80101(self):
    """
    """
    PD8010 = PD8010_2015_1()

    PD8010.section_input(self.Do, self.tnom)

    PD8010.material_input(sigma_y = self.sigma_y,
                          E=self.E,
                          Poisson=self.Poisson,
                          alpha_T= self.alpha_T)
    #
    self.code_name = "PD8010-1:2015 Part 1 : Steel Pipelines on Land"
    #
    self.Po = 0.0
    #
    # Start Calculation
    #
    # 6.4 Strength Check
    if "load" in self.load_type :
        t_check = self.tnom
        calculate_stress(self, t_check, PD8010)

    #PD8010.SectionProperties(self)
    self.sigma_h = PD8010.hoop_stress(self.Pi, self.Po)
    # code check stresses
    self.Sae, self.fd = PD8010.allowable_equivalent_stress()
    self.a = PD8010.substances_categorization(self.substance_category)
    self.Sah, self.fd_hs = PD8010.allowable_hoop_stress(self.a, self.pipe_history)

    self.UR_h, self.UR_eq = PD8010.limits_calculated_stress(self.sigma_h, self.Sah, 
                                                            self.sigma_e, self.Sae,
                                                            self.design_method)
    
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
                          self.E,
                          self.Poisson,
                          self.alpha_T,
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
                sigma_e += PD8010.stressFEA(self)
        # if Load input in forces
        #else:
        t_check = self.tnom
        sigmae_1 = functional_stress(self, t_check, PD8010, output=output)
        #
        try:
            sigmae_2 = PD8010.enviromental_stress(self.Px, self.Vip, self.Vop,
                                                  self.BMip, self.BMop, self.BMt,
                                                  output=output)
            self.M = math.sqrt(self.BMip**2 + self.BMop**2)
            self.Mt = math.sqrt(self.Vip**2 + self.Vop**2)
        
        except AttributeError:
            sigmae_2 = 0
        #
        self.sigma_e = sigmae_1 + sigmae_2
        
        if output:
            print('')
            print('Total equivalent stress : {:2.3f}'.format(self.sigma_e))
        #
        # External overpressure
        P = (self.Pi - self.Po)        
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
def PD8010_buckling(self, PD8010, output):
    """
    """
    #
    if output:
        print("")
        print("6.4.4 Buckling Check (AnnexG)")
    
    # 6.4.4.1.a
    # G.1.2 External pressure
    if self.Po > 0:
        self.Pc = PD8010.external_pressure(self.fo,
                                           self.sigma_u,
                                           self.root_search,
                                           output=output)
        
        # 6.4.4.1.b
        # G.2 Propagation buckling
        self.P, self.Pp = PD8010.propagation_buckling(self.Pi, self.Po,
                                                      self.tnom, output=output)
    else:
        self.Pc = 0
    
    if self.load_type:
        if "load" in self.load_type :
            # G.1.3 Axial compression
            if self.Fx:
                self.Fxc = PD8010.axial_compression(self.Fx, 
                                                    self.Do, self.tnom,
                                                    output=output)
            else:
                self.Fxc = 0
            # G.1.4 Bending
            self.Mc, self.epsilon_bc = PD8010.bending(self.Mb + self.M, 
                                                      self.Do, self.tnom,
                                                      output=output)
            # G.1.5 Torsion
            self.tau_c = PD8010.torsion(self.T + self.Mt, 
                                        self.Do, self.tnom,
                                        output=output)
            # G.1.6 Load Combination
            
            self.UR_lc = PD8010.load_combination(self.Mb + self.M, self.Mc, 
                                                 self.Fx, self.Fxc,
                                                 self.fo, 
                                                 self.Po, self.Pc,
                                                 output=output)

    #
    try:
        # 6.4.4.1.c
        #PD8010.upheaval_buckling(self)
        
        # G.1.7 Strain criteria
        self.epsilon_b = PD8010.strain_criteria(self.P, self.Pc, self.epsilon_bc)
        # 6.4.4.2 Ovality
        # G.4 Ovalization
        self.Cf, self.Cp, self.f = PD8010.ovalization(self.P, self.fo, 
                                                      self.epsilon_b)                
    except AttributeError:
        pass
    #
    # upheaval buckling
    if 'upheaval' in self.design_method:
        get_PD8010_upheaval(self)
    #       
    
#
def get_DNV(self):
    """
    """
    # DNV Section
    #
    #
    self.code_name = "DNV-OS-F101 (2007) Submarine Pipeline Systems"
    print("NO YET IMPLEMENT")
    exit()
#
def get_ASME(self):
    """
    """
    # ASME Section
    #
    #
    self.code_name = "ASME B31.3-2006 Process Piping"
    print("NO YET IMPLEMENT")
    exit()
#
def expansion_flexibility(self, output):
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
    if self.FS_code == "ASME_B31.3":
        asme = ASME()
        asme.bend_data(self.pipe_description,
                       self.flanges, 
                       self.T_, 
                       self.r2, 
                       self.R1, 
                       self.theta)

        h, k, lo, li, C1 = asme.AppendixD()
    #
    # BS EN 13480
    elif self.FS_code == "BS_EN_13480":
        #
        # Bend Flexibility
        if self.SIFs_type == 'BEND':
            #
            print ("No implement yet")
        #
        # Tee Flexibility
        elif self.SIFs_type == 'TEE':
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
    if output:
        print ("")
        print ("Expansion Flexibility ")
        print ("Based on : {:}".format(self.FS_code))
        print ("Member Description: {:}".format(self.pipe_description))
        print ("Member Type: {:}".format(self.SIFs_type))
        print ("h  = {:2.3f}".format(h))
        print ("k  = {:2.3f}".format(k))
        print ("lo = {:2.3f}".format(lo))
        print ("li = {:2.3f}".format(li))
        print ("C1 = {:2.3f}".format(C1))
    #
    return h, k, lo, li, C1
#
#
def material_derating(self):
    """
    """
    if self.derate_method == "DNV":
        derate = DNV()
        derate.temperature(self.Tmax)
        derate.material(self.SMYS, self.SMTS)
        Fy, Fu = derate.characteristic_material_properties(self.design_condition, 
                                                           self.material_type)
        #self.sigma_y = Fy
        #self.sigma_u = Fu
    
    return Fy, Fu 
#
#
def get_external_pressure(self, output):
    """
    """
    if not self.wave_length:
        wave = WaveStokes5.Stoke5()
        wave.Data(self.Hw/1000., self.Tw, self.d/1000.)
        #print ('Wave Length =',Wave.WaveLength)
        self.wave_length = wave.WaveLength * 1000.0
    #
    self.k = (2*math.pi)/self.wave_length
    
    # Maximum External Pressure
    self.Hz = (-self.z + (self.Hw/2)*
               ((math.cosh(self.k * (self.d + self.z)))/
                (math.cosh(self.k * self.d))))
    
    # (13.2-20)
    Po = (self.rho_w * self.g * self.Hz)
    #
    if output:
        print("")
        print("External Pressure")
        print('Wave Length               = {:2.3f}'.format(wave.WaveLength))
        print("Maximum External Pressure = {:2.3f}".format(self.Hz))
        print("                      Po  = {:2.3f}".format(Po))
    #
    return Po
#
#