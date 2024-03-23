# 
# Copyright (c) 2009 steelpy


# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
#from collections.abc import Mapping
from typing import NamedTuple
#import re
#import math
#

# package imports
#
#from steelpy.sections.utils.shape.stress import BeamStress
#from steelpy.utils.dataframe.main import DBframework
#import numpy as np
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
#
#
class ShapeProperty(NamedTuple):
    """
    area : Area [length^2]
    Zc,Yc : Centroid [length]
    Iy,Iz : Moment of Inertia [length^4]
    Sy,Sz : Elastic section modulus [length^3]
    Zy,Zz : Plastic section modulus [length^3]
    rx,rz : Radious of gyration [length]
    J : Torsional Constant
    Cw : Warping Constant
    alpha_sy, alpha_sz : Shear correction coefficient 
    """
    area:float #[length^2]
    # Centroid [length]
    Zc:float 
    Yc:float
    # y-y axis
    Iy:float # Moment of Inertia [length^4]
    Sy:float # Elastic section modulus [length^3]
    Zy:float # Plastic section modulus [length^3]
    ry:float # Radious of gyration [length]
    # z-z axis
    Iz:float
    Sz:float
    Zz:float
    rz:float
    # Torsion
    J:float
    Cw:float
    # Shear correction coefficient 
    alpha_sy: float = 0
    alpha_sz: float = 0
    #
    #
    @property
    def Asy(self):
        """ """
        if self.alpha_sy == 0:
            return 0
        return self.area / self.alpha_sy
    #
    @property
    def Asz(self):
        """ """
        if self.alpha_sz == 0:
            return 0
        return self.area / self.alpha_sz
    #
    def __str__(self) -> str:
        """ """
        #print('--->')
        output = "{:1.4e} {:1.4e} {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n"\
            .format(self.area, self.Iy, self.Iz, self.Zc, self.ry, 1)
        output += "{:25s} {:1.4e} {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n"\
            .format("", self.Sy, self.Sz, 1, self.rz, self.Cw)
        output += "{:25s} {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n"\
            .format("", self.Zy, self.Zz, 1, self.area*0)
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


