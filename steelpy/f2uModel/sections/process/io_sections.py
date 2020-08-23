# 
# Copyright (c) 2018-2020 steelpy


# Python stdlib imports
#import re
from typing import NamedTuple
#

# package imports
from steelpy.process.io_module.text import search_line

#
#
class PropertyOut(NamedTuple):
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
#
class SectionProperty(NamedTuple):
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
def find_section_dimensions(line_in: str) -> str:
    """
    """
    _key = {"height":r"\b((h(eight)?|d(epth)?)(z)?)\b",
            "web_thickness":r"\b(web(\_)?t(hickness)?(\_)?(y|w(all)?)?)\b",
            
            "top_flange_width":r"\b((t(op)?|b(t)?)(\_)?(f(lange)?)?((\_)?width)?)\b",
            "top_flange_thickness":r"\b(t(op|hickness)?(\_)?(f(lange)?)?(\_)?t(hickness)?)\b",
            
            "bottom_flange_width":r"\b(b(ottom)?(\_)?(f(lange)?)?(\_)?(width|b))\b",
            "bottom_flange_thickness" : r"\b((t(hickness)?)?(\_)?b(ottom)?(\_)?f(lange)?((\_)?thickness)?)\b",
            
            "theta" : r"\b(theta|angle)\b",
            "r" : r"\b(r)\b"}

    keyWord, line_out, _match = search_line(line_in, _key)
    
    if not _match:
        raise TypeError('    ** error parameter {:} not identified'.format(line_in))

    return keyWord
#
#
def get_dimension(self, _dim: str, value: float):
    """
    """
    #  Beam Section Definition
    if 'height' in _dim.lower():
        self.height = float(value)

    elif 'web_thickness' in _dim.lower():
        self.web_thickness = float(value)

    elif 'top_flange_width' in _dim.lower():
        self.top_flange_width = float(value)

    elif 'top_flange_thickness' in _dim.lower():
        self.top_flange_thickness = float(value)

    elif 'bottom_flange_width' in _dim.lower():
        self.bottom_flange_width = float(value)

    elif 'bottom_flange_thickness' in _dim.lower():
        self.bottom_flange_thickness = float(value)
    
    elif 'angle' in _dim.lower():
        self.angle = float(value)
    # root radius
    elif 'r' in _dim.lower():
        self.r = float(value) 
    
    #return self
#
