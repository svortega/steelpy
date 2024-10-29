# 
# Copyright (c) 2009 steelpy


# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
#from collections.abc import Mapping
from typing import NamedTuple
import re
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
    desity: float = 7850 # kg/m^3
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
            .format(self.area, self.Iy, self.Iz, self.Zc, self.ry, self.J)
        output += "{:25s} {:1.4e} {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n"\
            .format("", self.Sy, self.Sz, 1, self.rz, self.Cw)
        output += "{:25s} {:1.4e} {:1.4e} {:1.4e} {:1.4e}\n"\
            .format("", self.Zy, self.Zz, 1, self.area*self.desity)
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
def get_sect_list(parameters: list|tuple,
                  number:int = 9, step:int=3)->list:
    """Return : [diameter/height, base, a,
                 FAvy, FAvz, shear_stress,
                 build, compactness, title]"""
    # basic information
    section = [None] * number
    section[step+0] = 1.0       # FAvy
    section[step+1] = 1.0       # FAvz
    section[step+2] = 'maximum' # shear_stress
    section[step+3] = 'welded'  # build
    #
    for x in range(step):
        try:
            section[x] = parameters[x].value
        except AttributeError:
            raise IOError("units required")
        except IndexError:
            break
    return section
#
#
def get_prop_dict(parameters: dict) -> list:
    """Return : [FAvy, FAvz, shear_stress,
                 build, compactness, title]"""
    #
    section = [1.0, 1.0, 'maximum', 'welded', None, None]
    for key, item in parameters.items():
        if re.match(r"\b(SA(_|-|\s*)?inplane)\b", key, re.IGNORECASE):
            section[0] = item
        elif re.match(r"\b(SA(_|-|\s*)?outplane)\b", key, re.IGNORECASE):
            section[1] = item
        elif re.match(r"\b(shear(_|-|\s*)?stress)\b", key, re.IGNORECASE):
            section[2] = item
        elif re.match(r"\b(build)\b", key, re.IGNORECASE):
            section[3] = item
        elif re.match(r"\b(compactness)\b", key, re.IGNORECASE):
            section[4] = item
        elif re.match(r"\b(title)\b", key, re.IGNORECASE):
            section[5] = item
    #
    return section
#
#
def get_sect_dict(parameters: dict,
                  number:int = 11, step:int = 4)->list:
    """Return : [d/h, tw, b, tb, r,
                 FAvy, FAvz, shear_stress,
                 build, compactness, title]"""
    # basic information
    section = [None] * step
    #
    for key, item in parameters.items():
        # h
        if re.match(r"\b(h(eight)?|d(epth)?)\b", key, re.IGNORECASE):
            try:
                section[0] = item.value
            except AttributeError:
                raise IOError("units required")
        # tw
        elif re.match(r"\b(t(hickness|hk)?(_|-|\s*)?w(eb)?)\b", key, re.IGNORECASE):
            try:
                section[1] = item.value
            except AttributeError:
                raise IOError("units required")
        # b
        if re.match(r"\b((b(ase)?|w(idth)?)(_|-|\s*)?(f(lange)?)?)\b", key, re.IGNORECASE):
            try:
                section[2] = item.value
            except AttributeError:
                raise IOError("units required")
        # bt
        elif re.match(r"\b(t(hickness|hk)?(_|-|\s*)?(b(ase)?|f(lange)?))\b", key, re.IGNORECASE):
            try:
                section[3] = item.value
            except AttributeError:
                raise IOError("units required")
        # r
        elif re.match(r"\b(r(oot)?(_|-|\s*)?(r(adius|atio)?)?)\b", key, re.IGNORECASE):
            try:
                section[4] = item.value
            except AttributeError:
                raise IOError("units required")
    # bf
    if not section[3]:
        section[3] = section[1]
    # r
    if not section[4]:
        section[4] = 0.0
    #
    section.extend(get_prop_dict(parameters))
    return section
#
#
class ShapeDim(NamedTuple):
    """ """
    shape:str
    d:float
    tw:float
    b:float
    tb:float
    r: float
    #
    FAvy:float
    FAvz:float
    shear_stress:str
    build:str|None
    compactness:str|None
    title:str|None
    #
    @property
    def h(self) -> float:
        return self.d
