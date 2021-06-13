# 
# Copyright (c) 2019-2021 steelpy
#

# Python stdlib imports
from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, ClassVar

# package imports
from steelpy.process.units.main import Units
from steelpy.process.io_module.text import match_line


#
#
@dataclass
class SelfWeight:
    """ """
    __slots__ = [ 'x', 'y', 'z' ]

    def __init__(self):
        """ """
        self.x: float = 0
        self.y: float = 0
        self.z: float = 0
#
#
class NodeForce(NamedTuple):
    """
    """
    fx: float
    fy: float
    fz: float
    mx: float
    my: float
    mz: float
    number: int
    load_name: str
    system:str
    load_complex:int
    #
    def __str__(self,units:str="si") -> str:
        """ """
        output  = (f"{str(self.number):12s} {10*' '} "
                   f"{self.fx: 1.3e} {self.fy: 1.3e} {self.fy: 1.3e}"
                   f"{0: 1.3e} {0: 1.3e} {0: 1.3e}\n")
        #step = 12*" "
        output += (f"{self.coordinate_system.upper():12s} {10*' '} "
                   f"{self.mx: 1.3e} {self.my: 1.3e} {self.mz: 1.3e}\n")
        return output
    #
    @property
    def coordinate_system(self):
        if self.system != 0:
            return "local"
        return "global" 
#
#
#
def find_force_item(word_in):
    """
    """
    _key = {"Fx": r"\b((f|axial|p)(\_)?(x|1)?)\b",
            "Fy": r"\b((f|v|shear)(\_)?(y|i(n)?|2(\_)?(p(lane)?)?))\b",
            "Fz": r"\b((f|v|shear)(\_)?(z|o(ut)?|3(\_)?(p(lane)?)?))\b",
            #
            "Mx": r"\b((bending(\_)?)?m(oment)?(\_)?(x|t(orsion(al)?)?|1))\b",
            "My": r"\b((bending(\_)?)?m(oment)?(\_)?(y|i(n)?(\_)?p(lane)?|2))\b",
            "Mz": r"\b((bending(\_)?)?m(oment)?(\_)?(z|o(ut)?(\_)?p(lane)?|3))\b"}
    
    _match = match_line(word_in, _key)
    if not _match:
        raise IOError('  ** item {:} not recognized'.format(word_in))
    return _match
#
# FIXME: badly done to assign forces
def assign_force_item(self, mat_items):
    """
    Assign material's input data to a dictionary by name/number
    """
    #
    #if not self.actions:
    #    self.actions = Actions(Fx = 0 * self._units.N, 
    #                           Fy = 0 * self._units.N, 
    #                           Fz = 0 * self._units.N, 
    #                           Mx = 0 * self._units.N * self._units.m, 
    #                           My = 0 * self._units.N * self._units.m, 
    #                           Mz = 0 * self._units.N * self._units.m)
    #
    for key, value in mat_items.items():
        _item = find_force_item(key)
        #
        if _item == 'Fx':
            self.Fx = value
        elif _item == 'Fy':
            self.Fy = value        
        elif _item == 'Fz':
            self.Fz = value
        #
        elif _item == 'Mx':
            self.Mx = value        
        elif _item == 'My':
            self.My = value
        elif _item == 'Mz':
            self.Mz = value
        else:
            raise IOError('error Force item : {:} not recognized'
                          .format(key))
    #
    #self.actions = Actions(_P, _Vy, _Vx, _Mt, _Mx, _My)
    #print('ok')
#
#
class Actions:
    """
    Force & bending moments
    """
    __slots__ = ['_Fx', '_Fy', '_Fz', '_Mx', '_My', '_Mz']
    
    def __init__(self):
        """
        """
        #self._units = Units()
        self._Fx = 0 #* self._units.N
        self._Fy = 0 #* self._units.N
        self._Fz = 0 #* self._units.N
        self._Mx = 0 #* self._units.N * self._units.m
        self._My = 0 #* self._units.N * self._units.m
        self._Mz = 0 #* self._units.N * self._units.m
        #
        # ----- Load offsets -----
        #self.Xe = 0 * units.m # load offset from support 1
    
    #@property
    #def units(self):
    #    """
    #    units [length, mass, time, temperature, force, pressure/stress]/n
    #    """
    #    return self._units    
    #
    def __setattr__(self, name, value):
        """ """
        if name=="device":
            print("device test")
        else:
            super().__setattr__(name, value)
    #
    #def __setitem__(self, **kwargs):
    #    """
    #    """
    #    print('here')
    def member_actions(self, **kwargs):
        """
        """
        assign_force_item(self, kwargs)
    #
    #
    #def __str__(self, units:str="si") -> str:
    #    """ """
    #    output  = (f"{str(self.number):12s} {self.fx: 1.3e} "
    #               f"{self.fy: 1.3e} {self.fy: 1.3e}\n")
    #    step = 12*" "
    #    output += (f"{step} {10*' '} {self.mx: 1.3e} "
    #               f"{self.my: 1.3e} {self.mz: 1.3e}\n")
    #    return output
#