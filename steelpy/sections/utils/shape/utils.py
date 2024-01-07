# 
# Copyright (c) 2009 steelpy


# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Mapping
from typing import NamedTuple
#import re
import math
#

# package imports
#
#from .operations import get_sect_properties,  get_Isection
from steelpy.utils.dataframe.main import DBframework
import numpy as np
#
#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------
#
#
#
#
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
class ShapeBasic:
    name: str|int
    build: str = 'welded'
    #
    # -------------------------------------
    # Operations
    # -------------------------------------
    #
    def stress(self, actions=None, stress=None, df=None):
        """return cross section stress"""
        #print('-->')
        try:
            df.columns
            dfres = DBframework()
            # -------------------------------------
            try:
                stress_df = df[['node_end', 'tau_x', 'tau_y', 'tau_z',
                                'sigma_x', 'sigma_y', 'sigma_z']]
                stress = self._stress(stress=stress_df)
            except KeyError:
                actions_df = df[['node_end', # 'load_title', 
                                 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                                 'B', 'Tw']]
                stress = self._stress(actions=actions_df)
                # -------------------------------------
                header = ['load_name', 'component_name', 
                          'load_level', 'load_system',
                          'element_name','node_end','stress_points', 'y', 'z']
                # -------------------------------------
                coord = stress.stress_point
                items = [[row.load_name, row.component_name, 
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

        except AttributeError:
            if stress:
                print("stress")
            elif actions:
                print("actions")
            else:
                print('--> ??')
                1 / 0
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
    #
    # -------------------------------------
    #
    #def _print_section_properties(self):
    #    """
    #    """
    #    file = shape_io.print_header()
    #    file.extend(self._shape())
    #    file.extend(shape_io.print_properties(self))
    #    return file    
    #    
    # -------------------------------------
    #
    #@property
    def properties(self):
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
        return self._properties()
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
#
class ShapeProperty(NamedTuple):
    """ """
    area:float
    Zc:float
    Yc:float
    Iy:float
    Zey:float
    Zpy:float
    ry:float
    Iz:float
    Zez:float
    Zpz:float
    rz:float
    J:float
    Cw:float
    #
    def __str__(self) -> str:
        """ """
        #print('--->')
        output = "{:1.4e} {:1.4e} {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n"\
            .format(self.area, self.Iy, self.Iz, self.Zc, self.ry, 1)
        output += "{:25s} {:1.4e} {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n"\
            .format("", self.Zey, self.Zez, 1, self.rz, self.Cw)
        output += "{:25s} {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n"\
            .format("", self.Zpy, self.Zpz, 1, self.area*0)
        return output
#
#
class SectionPropertyXX(NamedTuple):
    """
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
    area: float
    Zc  : float
    Yc  : float
    Iy  : float
    Zy : float
    Sy : float
    SFy : float
    ry  : float
    Iz  : float
    Zz : float
    Sz : float
    SFz : float
    rz  : float
    SCy  : float
    SCz  : float
    Cw  : float
    #
    def __str__(self) -> str:
        """ """
        #print('--->')
        output = "{:<14s} {:1.4e} {:1.4e} {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n"\
            .format(self.area, self.Iy, self.Iz, self.Zc, self.ry, 1)
        output += "{:<14s} {:1.4e} {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n"\
            .format(" "*14, self.Zy, self.Zz, self.SCy, self.rz, self.Cw)
        output = "{:<14s} {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n"\
            .format(self.area*0, self.Sy, self.Sz, self.SCz)        
        return output
#
#
#
# ---------------------------------
#
def get_sect_properties(properties:list[float], steps:int=10):
    """ """
    #sect_prop = [None for _ in range(steps)]
    sect_prop = []
    for x, item in enumerate(properties):
        try:
            #sect_prop[x] = item.value
            sect_prop.append(item.value)
        except AttributeError:
            #sect_prop[x] = item
            sect_prop.append(item)
            raise IOError('units required')
    return sect_prop
#
#
def get_sect_prop_df(df):
    """ """
    for cname in df:
        try:
            df[cname] = df[cname].apply(lambda x: x.convert('metre').value)
        except AttributeError:
            pass
    #print('---')
    return df   
#
#
#
def get_Isection(parameters: list):
    """ [d, tw, bf, tf, bfb, tfb, r, title] """
    # basic information
    section = parameters[:4] # d, tw, bf, tf
    # check if str title at the end
    if isinstance(parameters[-1], str):
        title = parameters.pop()
    else:
        title = None
    #
    # check if root radius
    if len(parameters) == 4:
        section.append(parameters[2]) # bfb
        section.append(parameters[3]) # tfb
        section.append(0)             # root radius
        
    elif len(parameters) == 5:
        r = parameters.pop()
        section.append(parameters[2]) # bfb
        section.append(parameters[3]) # tfb
        section.append(r)
        
    elif len(parameters) == 6:
        section.append(0)
        
    else:
        if len(parameters) != 7:
            raise IOError('Error Ibeam input data ')
    #
    section.append(title)
    #  
    return section

