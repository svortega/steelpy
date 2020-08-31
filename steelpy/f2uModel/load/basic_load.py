#
# Copyright (c) 2009-2020 fem2ufo
# 

# Python stdlib imports
from array import array
from collections.abc import Mapping
from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, Iterator, Dict, Iterable, ClassVar, Union


# package imports
from steelpy.f2uModel.load.element import BeamDistributed, BeamPoint
from steelpy.f2uModel.load.node import NodeLoad
#
#
# ---------------------------------
#
@dataclass
class SelfWeight:
    """ """
    __slots__ = ['x', 'y', 'z']
    
    def __init__(self):
        """ """
        self.x:float = 0
        self.y:float = 0
        self.z:float = 0
#
#
class LoadTypes:
    """
    """
    __slots__ = ['_nodal_load', '_nodal_mass', '_nodal_displacement',
                 '_beam_line', '_beam_point', '_selfweight',
                 'name', 'number', 'title']
    
    def __init__(self):
        """
        """
        self._nodal_load = NodeLoad()
        self._nodal_displacement = NodeLoad()
        self._nodal_mass = BeamDistributed()
        self._beam_line = BeamDistributed()
        self._beam_point = BeamPoint()
        self._selfweight = SelfWeight()
    #
    @property
    def selfweight(self):
        """
        The self weight form allows you to specify multipliers to 
        acceleration due to gravity (g) in the X, Y, and Z axes. 
        If switched on, the default self weight acts in the Y axis 
        with a magnitude and sign of -1."""
        return self._selfweight
    #
    @property
    def point_node(self):
        """
        """
        return self._nodal_load
    
    @point_node.setter
    def point_node(self, values:List):
        """
        Point Load
        """
        for value in values:
            self._nodal_load[value[0]] = value[1:]
    #
    @property
    def mass_node(self):
        """
        """
        return self._nodal_mass
    
    @mass_node.setter
    def mass_node(self, values:List):
        """
        """
        for value in values:
            self._nodal_mass[value[0]] = value[1:]
    #
    @property
    def displacement_node(self):
        """
        """
        return self._nodal_displacement
    
    @displacement_node.setter
    def displacement_node(self, values:List):
        """
        """
        for value in values:
            self._nodal_displacement[value[0]] = value[1:]
    #
    #
    @property
    def udl_beam(self):
        """
        Uniformly Distributed Load (udl)
        """
        return self._beam_line
    
    @udl_beam.setter
    def udl_beam(self, values:List):
        """
        Uniformly Distributed Load (udl)
        value : [qx1, qy1, qz1, qx2, qy2, qz2, L1, L2]
        
                        |
             q1         | q2
        o------|        |----------o
        |                          |
        +  L1  +        +    L2    +
        """
        for value in values:
            self._beam_line[value[0]] = value[1:]
    #
    @property
    def line_beam(self):
        """
        Linear Varying Load (lvl) - Non Uniformly Distributed Load
        """
        return self._beam_line
    
    @line_beam.setter
    def line_beam(self, values:List):
        """
        Linear Varying Load (lvl) - Non Uniformly Distributed Load
                value : [qx1, qy1, qz1, qx2, qy2, qz2, L1, L2]
        
                        |
             q1         | q2
        o------|        |----------o
        |                          |
        +  L1  +        +    L2    +
        """
        for value in values:
            self._beam_line[value[0]] = value[1:]
    #
    #
    @property
    def point_beam(self):
        """ Concentrated force """
        return self._beam_point
    
    @point_beam.setter
    def point_beam(self, values:List):
        """
        Concentrated force
        """
        for value in values:
            self._beam_point[value[0]] = value[1:]    
    
#
class BasicLoad(Mapping):
    """
    FE Load Cases
    
    LoadType
        |_ name
        |_ number
        |_ basic
        |_ combination_level
        |_ time_series
        |_ 
        |_ temperature
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    #
    __slots__ = ['_load', 'gravity']
    #
    def __init__(self):
        """
        """
        self._load: Dict = {}
        self.gravity = 9.80665 # m/s^2
    #
    def __setitem__(self, load_name:Union[str,int],
                    load_title:str) -> None:
        """
        """
        self._load[load_name] = LoadTypes()
        self._load[load_name].name = load_name
        self._load[load_name].title = load_title
        self._load[load_name].number = len(self._load)
    
    def __getitem__(self, load_name:Union[str,int]):
        """
        """
        try:
            return self._load[load_name]
        except KeyError:
            raise IOError("load case not defined")
    #
    def __delitem__(self, load_name:Union[str,int]):
        """
        """
        del self._load[load_name]
    #
    def __len__(self) -> int:
        return len(self._load)
    #
    def __iter__(self):
        """
        """
        return iter(self._load)
    #
    #
    @property
    def g(self, load_name):
        """"""
        return self.gravity
#
  
    
    
