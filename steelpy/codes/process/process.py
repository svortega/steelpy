# 
# Copyright (c) 2019-2020 iLift
#

# Python stdlib imports
from typing import NamedTuple, Tuple, ClassVar, Union, List

# package imports
#from steelpy.material.material import Materials
#from steelpy.sections.shape import Sections
from steelpy.process.units.main import Units
from steelpy.process.io_module.text import match_line


#
#
def find_stress_item(word_in):
    """
    """
    _key = {"sigma_x": r"\b((sigma|stress)(\_)?(x))\b",
            "sigma_y": r"\b((sigma|stress)(\_)?(y))\b",
            "sigma_z": r"\b((sigma|stress)(\_)?(z))\b",
            #
            "tau_x": r"\b((tau)(\_)?(x))\b",
            "tau_y": r"\b((tau)(\_)?(y))\b",
            "tau_z": r"\b((tau)(\_)?(z))\b"}
    
    _match = match_line(word_in, _key)
    if not _match:
        raise IOError('  ** item {:} not recognized'.format(word_in))
    return _match
#
def assign_stress_item(self, mat_items):
    """
    Assign stress 
    
    """
    #if not self.stress:
    #    self.stress = Stress(0, 0, 0, 0, 0, 0)
    
    for key, value in mat_items.items():
        _item = find_stress_item(key)
        #
        if _item == 'sigma_x':
            self.stress.sigma_x += value
        elif _item == 'sigma_y':
            self.stress.sigma_y += value
        elif _item == 'sigma_z':
            self.stress.sigma_z += value
        elif _item == 'tau_x':
            self.stress.tau_x += value
        elif _item == 'tau_y':
            self.stress.tau_y += value
        elif _item == 'tau_z':
            self.stress.tau_z += value
        else:
            raise IOError('error stress item : {:} not recognized'
                          .format(_item))
    #
    #self.stress = Stress(_sigma_x, _sigma_y, _sigma_z, _tau_x, _tau_y, _tau_z)
    #print('ok')
#
#
#
class CodeResults(NamedTuple):
    """
    """
    axial:Tuple
    shear:Tuple
    bending:Tuple
    combined:Tuple
    report:List[str]
    #
    @property
    def total(self):
        """"""
        total_UR = self.combined.UR
        UR_flag = self.combined.UR_flag
        shear_res = max(self.shear[:2])
        if shear_res > total_UR:
            total_UR = shear_res
            UR_flag = self.shear.UR_flag
        return SummaryResults(total_UR, UR_flag)
#
#
class ChapterResults(NamedTuple):
    """
    """
    URy : float
    URz : float
    UR_flag : str
    allowable_y:float
    allowable_z:float
    #
    @property
    def status(self):
        if max(self.URy, self.URz) > 1.0:
            return "fail"
        return "pass"
#
#class AxialResults(NamedTuple):
#    """
#    """
#    UR : float
#    UR_flag : str
#    allowable:float
#    warning:str
#
class SummaryResults(NamedTuple):
    """
    """
    UR : float
    UR_flag : str
    
    @property
    def status(self):
        if self.UR > 1.0:
            return "fail"
        return "pass"
#
#
#
#
#
class BeamDesignParameters:
    """Beam Desing Patameters """
    __slots__ = ["units", "_material", "_section",
                 "g", "C", "theta",
                 "Ly", "Lz", "Ky", "Kz",
                 "Cmy", "Cmz"]
    
    def __init__(self):
        """
        """
        self.units = Units()
        #
        self.g = 9.810 * self.units.m / self.units.second**2
        self.C = 0.30
        #
        self.Ky = 1.0 
        self.Kz = 1.0
        #
        # Moment Modifiers
        self.Cmy = 0.85
        self.Cmz = 0.85
        # id this the section rotation?
        self.theta = 0 * self.units.degrees
    #
    @property
    def material(self):
        """
        """
        return self._material
    
    @material.setter
    def material(self, value):
        """
        """
        self._material = value
    
    @property
    def section(self):
        """
        """
        return self._section
    
    @section.setter
    def section(self, value):
        """
        """
        self._section = value    
    #
    