# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
from typing import NamedTuple
#import sys
import math


# package imports
from steelpy.utils.dataframe.main import DBframework
import numpy as np

#
#-------------------------------------------------
#
#
#
#-------------------------------------------------
#
#
#
#
@dataclass(kw_only=True)
class ShapeStressBasic:
    name: str|int
    build: str = 'welded'
    #
    # -------------------------------------
    # Operations
    # -------------------------------------
    #
    #
    def stress(self, E: float, G: float, poisson: float, 
               actions=None, stress=None, df=None):
        """return cross section stress"""
        #print('-->')
        #try:
        df.columns
        dfres = DBframework()
        # -------------------------------------
        try:
            stress_df = df[['node_end', 'tau_x', 'tau_y', 'tau_z',
                            'sigma_x', 'sigma_y', 'sigma_z']]
            stress = self._stress(stress=stress_df, G=G, E=E)
        except KeyError:
            actions_df = df[['node_end', # 'load_title', 
                             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                             'Psi', 'B', 'Tw']]
                             #'theta1', 'theta2', 'theta3']]
            stress = self._stress(actions=actions_df,
                                  G=G, E=E, poisson=poisson)
            # -------------------------------------
            header = ['load_name', 'mesh_name', 
                      'load_level', 'load_system',
                      'element_name', 
                      'node_end','stress_points', 'y', 'z']
            # -------------------------------------
            coord = stress.stress_point
            items = [[row.load_name, row.mesh_name, 
                      row.load_level, row.load_system, 
                      row.element_name, row.node_end,
                      x+1, coord.y[x], coord.z[x]]
                     for x in range(len(coord.y))
                     for row in df.itertuples()]
            df_stress = dfres.DataFrame(data=items, columns=header, index=None)
            # -------------------------------------
            # axial stress
            df_stress['tau_x'] = np.array(stress.tau_x).flatten()
            # In Plane shear stress
            df_stress['tau_y'] = np.array(stress.tau_y).flatten()
            # Out Plane shear stress
            df_stress['tau_z'] = np.array(stress.tau_z).flatten()
            # torsion stress
            df_stress['sigma_x'] = np.array(stress.sigma_x).flatten()
            # In plane bending stress
            df_stress['sigma_y'] = np.array(stress.sigma_y).flatten()
            # Out plane bending stress
            df_stress['sigma_z'] = np.array(stress.sigma_z).flatten()
            # return dataframe
            return df_stress

        #except AttributeError:
        #    if stress:
        #        print("stress")
        #    elif actions:
        #        print("actions")
        #    else:
        #        print('--> ??')
        #        1 / 0
    #
    #
    #
    def shear_stress(self, Vy, Vz,
                     #stress_type:str ='average',
                     alpha: float = 2.0):
        """
        Vy :
        Vz : 
        alpha: Shear correction factor
        Qy,z : First moment of area
        """
        # -------------------------------------------------
        #            Shear Stress Calculation
        I = self.I
        #coord =  self.section_coordinates()
        #
        # Q area shear stress
        Qby, Qbz = self.Qb()
        #
        # InPlane
        tau_z = [Vy * item / I.y
                 for item in Qbz]
        #
        # Out Plane
        tau_y = [Vz * item / I.z
                 for item in Qby]
        #
        #1 / 0
        return tau_y, tau_z
    #
    def torsional_stress(self, Mt, B=None):
        """
        """
        coord =  self.section_coordinates()
        tau_max = self.taux_max(Mt=Mt)
        tau_x = [tau_max # * abs(item)
                 for item in coord.z]
        return tau_x
    #
    def _stress(self, actions: list, G: float, E: float, 
                poisson: float, stress=None):
        """
        z
        ^
        + > Y
        
        sigma_x : Normal stress due to axial force alone.
        sigma_y : Normal stress due to bending moment about y axis.
        sigma_z : Normal stress due to bending moment about z axis.
        tau_y : Shear stress due to shear force in y direction.
        tau_z : Shear stress due to shear force in z direction.
        tau_x : Shear stress due to torsional moment.
        """
        #
        # get section's coordinates
        coord =  self.section_coordinates()
        prop = self.properties(poisson=poisson)
        #
        # ----------------------------------------------
        # get shear stress
        tau_y, tau_z = self.shear_stress(Vy=actions.Fy,
                                         Vz=actions.Fz)
        #
        #
        tau_x = [ actions.Mx * 0
                  for item in coord.z]
        #
        # ----------------------------------------------
        # get bending stress
        sigma_x = [actions.Fx / prop.area for item in coord.y]
        sigma_y = [actions.Mz * item / prop.Iy for item in coord.z]
        sigma_z = [actions.My * item / prop.Iz for item in coord.y]
        #
        # ----------------------------------------------
        # Thin walled cross sections
        try: # shape available : I & C sections
            # 
            tau_t = self.tau_t(psi=actions.Psi, G=G)
            #
            sigma_w, tau_w = self.warping_stress(B=actions.B,
                                                 Tw=actions.Tw)
            #
            # shear + torsion
            tau_y = [tau_t[x] + item
                     for x, item in enumerate(tau_y)]
            # shear + torsion + warping torsion
            tau_z = [tau_t[x] + tau_w[x] + item
                     for x, item in enumerate(tau_z)]
            #
            # bending + warping torsion
            sigma_y = [sigma_w[x] + item
                       for x, item in enumerate(sigma_y)]
            #
        except AttributeError:
            # ----------------------------------------------
            # FIXME: Torsion - solid cross section 
            tau_x = self.torsional_stress(Mt=actions.Mx)
        #
        # ----------------------------------------------
        # stress process 
        stress_out = BeamStress(sigma_x, sigma_y, sigma_z, 
                                tau_x, tau_y, tau_z, coord)
        #
        if stress:
            stress_out = self.add_stress(stress=stress, other=stress_out)
        #
        return stress_out
    #
    def add_stress(self, stress, other):
        """ """
        if isinstance(stress.tau_x, list):
            stress.tau_x = self._combine_stress(other.tau_x, stress.tau_x)
            stress.tau_y = self._combine_stress(other.tau_y, stress.tau_y)
            stress.tau_z = self._combine_stress(other.tau_z, stress.tau_z)
            #
            stress.sigma_x = self._combine_stress(other.sigma_x, stress.sigma_x)
            stress.sigma_y = self._combine_stress(other.sigma_y, stress.sigma_y)
            stress.sigma_z = self._combine_stress(other.sigma_z, stress.sigma_z)
        else:
            # Assuming global stress
            stress.tau_x = self._add_global_stress(other.tau_x, stress.tau_x)
            stress.tau_y = self._add_global_stress(other.tau_y, stress.tau_y)
            stress.tau_z = self._add_global_stress(other.tau_z, stress.tau_z)
            #
            stress.sigma_x = self._add_global_stress(other.sigma_x, stress.sigma_x)
            stress.sigma_y = self._add_global_stress(other.sigma_y, stress.sigma_y)
            stress.sigma_z = self._add_global_stress(other.sigma_z, stress.sigma_z)
        #
        return stress
    #
    def _add_global_stress(self, stress_local, stress_global):
        """
        """  
        # _new_stress = [ _item + math.copysign(1, _item) * stress_global  
        #                 if _item.value != 0  else _item for _item in stress_local] #aldh6850
        #
        #aldh6850 - update to ensure the "global" stress has the same sign as the "local" stress to be conservative
        #aldh6850 - update to ensure when the "local" stress is zero the "global" stress is used
        #
        _new_stress = [ item + math.copysign(1, item) * abs(stress_global)  
                        if item != 0  else stress_global
                        for item in stress_local] #aldh6850
        #
        return _new_stress
    #
    def _combine_stress(self, stress_1, stress_2):
        """
        stress_1: 
        stress_2:
        """
        # change * by +
        _new_stress = [stress_1[x] + math.copysign(1, stress_1[x].value) * abs(stress_2[x]) 
                       for x in range(9)]
        return _new_stress   
    #
    #  
    #    
    # -------------------------------------
    #
    #@property
    def properties(self, poisson: float):
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
        return self._properties(poisson=poisson)
    #
    #
    #def push_property(self):
    #    """ """
    #    self.properties
    #    

#
#
#
#
class Sigma(NamedTuple):
    """
    Cauchy stress tensor
    https://en.wikipedia.org/wiki/Stress_(mechanics)
    """
    sigma_11:float
    sigma_21:float
    sigma_31:float
    #
    sigma_12:float
    sigma_22:float
    sigma_32:float
    #
    sigma_13:float
    sigma_23:float 
    sigma_33:float
    

# TODO: class duplicated
#@dataclass
class Stress:
    """
    General stress class
    """
    __slots__ = ('sigma_x', 'sigma_y', 'sigma_z', 
                 'tau_x', 'tau_y', 'tau_z','_units')    
    
    def __init__(self, sigma_x, sigma_y, sigma_z,
                 tau_x, tau_y, tau_z):
        """
        """ 
        #print('-->')
        self.sigma_x: list = sigma_x
        self.sigma_y: list = sigma_y
        self.sigma_z: list = sigma_z
        self.tau_x:   list = tau_x
        self.tau_y:   list = tau_y
        self.tau_z:   list = tau_z
    
    def von_mises(self) -> list:
        """
        Returns Von-Mises stress
        """
        items = []
        for x in range(len(self.sigma_x)):
            items.append((0.50 * ((self.sigma_x[x] - self.sigma_y[x])**2 +
                                   (self.sigma_y[x] - self.sigma_z[x])**2 +
                                   (self.sigma_z[x] - self.sigma_x[x])**2)
                           + 3 * (self.tau_x[x]**2 + self.tau_y[x]**2
                                  + self.tau_z[x]**2))**0.50)
        return items
    #
    def sigma(self) -> list:
        """
        Cauchy stress tensor
        """
        sigma_items = []
        for x in range(len(self.sigma_x)):
            sigma_items.append(Sigma(self.sigma_x[x], self.tau_x[x], self.tau_y[x],
                                      self.tau_x[x], self.sigma_y[x], self.tau_z[x],
                                      self.tau_y[x], self.tau_z[x], self.sigma_z[x]))
        return sigma_items
    #
    def analyze_stress_state(self, index) -> None:
        """
        """
        import numpy as np

        _sigma_items = self.sigma()
        _items = _sigma_items[index]
        # load the stresses into our matrix and compute the 
        # deviatoric and isotropic stress matricies
        sigma = np.array([[_items.sigma_11.value, _items.sigma_21.value, _items.sigma_31.value],
                          [_items.sigma_12.value, _items.sigma_22.value, _items.sigma_32.value],
                          [_items.sigma_13.value, _items.sigma_23.value, _items.sigma_33.value]])
    
        self.isotropic_stress = 1.0/3.0*np.trace(sigma)*np.eye(3)
        self.deviatoric_stress = sigma - self.isotropic_stress
        # compute principal stresses
        eigvals = list(np.linalg.eigvalsh(sigma))
        eigvals.sort()
        eigvals.reverse()
        # compute max shear stress
        self.maximun_shear = (max(eigvals) - min(eigvals)) / 2.0
        # compute the stress invariants
        I1 = np.trace(sigma)
        J2 = (1.0 / 2.0 * np.trace(np.dot(self.deviatoric_stress,
                                          self.deviatoric_stress)))
        J3 = (1.0 / 3.0 * np.trace(np.dot(self.deviatoric_stress, 
                                          np.dot(self.deviatoric_stress,
                                                 self.deviatoric_stress))))
        # compute other common stress measures
        self.mean_stress = 1.0/3.0*I1
        self.equivalent_stress  = np.sqrt(3.0*J2)
        # compute lode coordinates
        self.lode_r = np.sqrt(2.0*J2)
        self.lode_z = I1/np.sqrt(3.0)

        dum = (3.0 * np.sqrt(6.0)
               * np.linalg.det(self.deviatoric_stress / self.lode_r))
        self.lode_theta = 1.0 / 3.0 * np.arcsin(dum)
        # compute the stress triaxiality
        self.triaxiality = self.mean_stress / self.equivalent_stress
        # Print out what we've found
        headerprint(" Stress State Analysis ")
        matprint("Input Stress",sigma)
        headerprint(" Component Matricies ")
        matprint("Isotropic Stress",self.isotropic_stress)
        matprint("Deviatoric Stress",self.deviatoric_stress)
        headerprint(" Scalar Values ")
        valprint("P1",eigvals[0])
        valprint("P2",eigvals[1])
        valprint("P3",eigvals[2])
        valprint("Max Shear", self.maximun_shear)
        valprint("Mean Stress",self.mean_stress)
        valprint("Equivalent Stress", self.equivalent_stress)
        valprint("I1",I1)
        valprint("J2",J2)
        valprint("J3",J3)
        valprint("Lode z",self.lode_z)
        valprint("Lode r",self.lode_r)
        valprint("Lode theta (rad)",self.lode_theta)
        valprint("Lode theta (deg)",np.degrees(self.lode_theta))
        valprint("Triaxiality",self.triaxiality)
        headerprint(" End Output ")
    #
    #
#
#
def headerprint(string):
    """ Prints a centered string to divide output sections. """
    mywidth = 64
    mychar = "="
    numspaces = mywidth - len(string)
    before = int(math.ceil(float(mywidth-len(string))/2))
    after  = int(math.floor(float(mywidth-len(string))/2))
    print("\n"+before*mychar+string+after*mychar+"\n")

def valprint(string, value):
    """ Ensure uniform formatting of scalar value outputs. """
    print("{0:>30}: {1: .10e}".format(string, value))

def matprint(string, value):
    """ Ensure uniform formatting of matrix value outputs. """
    print("{0}:".format(string))
    print(value)
  
#
#
@dataclass()
class BeamStress:
    """
    beam element stress
    """
    #__slots__ = ['sigma_x', 'sigma_y', 'sigma_z',
    #             'tau_x', 'tau_y', 'tau_z','x']
    
    #def __init__(self, sigma_x, sigma_y, sigma_z,
    #             tau_x, tau_y, tau_z, x):
    """
    """
    sigma_x : list
    sigma_y : list
    sigma_z : list
    tau_x   : list
    tau_y   : list
    tau_z   : list
    stress_point : list|tuple
    
    #def von_mises(self) -> list:
    #    """
    #    Returns Von-Mises stress
    #    """
    #    _items = []
    #    for x in range(len(self.sigma_x)):
    #        _items.append(((self.sigma_x[x] + self.sigma_y[x] + self.sigma_z[x])**2 +
    #                       3 * (self.tau_x[x]**2 + self.tau_y[x]**2 + self.tau_z[x]**2))**0.50)
    #    return _items
#
#
class PlateStress:
    """
    plate element stress
    """
    __slots__ = ['sigma_x', 'sigma_y', 'sigma_z', 
                 'tau_x', 'tau_y', 'tau_z','_units']    
    
    def __init__(self, sigma_x, sigma_y, sigma_z,
                 tau_x, tau_y, tau_z):
        """
        """ 
        #print('-->')
        self.sigma_x: list = sigma_x
        self.sigma_y: list = sigma_y
        self.sigma_z: list = sigma_z
        self.tau_x:   list = tau_x
        self.tau_y:   list = tau_y
        self.tau_z:   list = tau_z
    
    def von_mises(self) -> list:
        """
        Returns Von-Mises stress
        """
        _items = []
        for x in range(len(self.sigma_x)):
            _items.append(((self.sigma_y[x]**2 + self.sigma_z[x]**2 - self.sigma_y[x] * self.sigma_z[x]) +
                           3 * (self.tau_y[x]**2 + self.tau_z[x]**2))**0.50)
        return _items
    #
    def factor_by(self, other:float)->None:
        """ """
        for x in range(len(self.sigma_x)):
            self.sigma_x[x] *= other
            self.sigma_y[x] *= other
            self.sigma_z[x] *= other
            self.tau_x[x] *= other
            self.tau_y[x] *= other
            self.tau_z[x] *= other
    #
    def principal_stresses(self, index:int):
        """ """
        import numpy as np
        #
        _sigma_items = self.sigma()
        _items = _sigma_items[index]
        # load the stresses into our matrix and compute the
        # deviatoric and isotropic stress matricies
        sigma = np.array([[_items.sigma_11.convert('pascal').value,
                           _items.sigma_21.convert('pascal').value,
                           _items.sigma_31.convert('pascal').value],
                          [_items.sigma_12.convert('pascal').value,
                           _items.sigma_22.convert('pascal').value,
                           _items.sigma_32.convert('pascal').value],
                          [_items.sigma_13.convert('pascal').value,
                           _items.sigma_23.convert('pascal').value,
                           _items.sigma_33.convert('pascal').value]])
        # compute principal stresses
        eigvals = list(np.linalg.eigvalsh(sigma))
        eigvals.sort()
        eigvals.reverse()
        # Get principal stresses
        P1 = eigvals[0]
        P2 = eigvals[1]
        P3 = eigvals[2]
        # compute max shear stress
        tau_max = (max(eigvals) - min(eigvals)) / 2.0
        #
        return [P1, P2, P3, tau_max]
    #
    def sigma(self) -> list:
        """
        Cauchy stress tensor
        """
        _sigma_items = []
        for x in range(len(self.sigma_x)):
            _sigma_items.append(Sigma(self.sigma_x[x], self.tau_x[x], self.tau_y[x],
                                      self.tau_x[x], self.sigma_y[x], self.tau_z[x],
                                      self.tau_y[x], self.tau_z[x], self.sigma_z[x]))
        return _sigma_items