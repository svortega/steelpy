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
class PointNode(NamedTuple):
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
    def __str__(self):
        print('-->')
        return "{: 14.5f} {: 14.5f} {: 14.5f} {: 14.5f} {: 14.5f} {: 14.5f}".format(self.fx, self.fy, self.fz,
                                                                                    self.mx, self.my, self.mz)
    #
    #@property
    #def name(self):
      #  """
      #  """
      #  return self.load_name
    #@name.setter
    #def name(self, load_name:str):
    #    """
    #    """
    #    self.load_name = load_name
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
    __slots__ = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                 '_units']
    
    def __init__(self):
        """
        """
        self._units = Units()
        self.Fx: ClassVar = 0 * self._units.N
        self.Fy: ClassVar = 0 * self._units.N
        self.Fz: ClassVar = 0 * self._units.N
        self.Mx: ClassVar = 0 * self._units.N * self._units.m
        self.My: ClassVar = 0 * self._units.N * self._units.m
        self.Mz: ClassVar = 0 * self._units.N * self._units.m
        #
        # ----- Load offsets -----
        #self.Xe = 0 * units.m # load offset from support 1
    
    @property
    def units(self):
        """
        units [length, mass, time, temperature, force, pressure/stress]/n
        """
        return self._units    
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