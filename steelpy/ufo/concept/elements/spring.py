# 
# Copyright (c) 2009-2018 fem2ufo
# 


# Python stdlib imports
from typing import Tuple, Dict, List, ClassVar, Iterable, MutableMapping
from collections.abc import Mapping


class BasicSpring:
    """
    """
    __slots__ = ('_material', '_section', '_nodes',
                 '_elements', '_type',
                 '_points') 
    
    def __init__(self, spring_type:str, points: ClassVar) -> None:
        """
        """
        self._type: str = spring_type
        self._points: ClassVar = points
        #self._mesh: ClassVar = mesh
        #
        # self._properties = BeamProperties()
    
    #
    @property
    def material(self):
        """
        """
        return self._material[0]

    @material.setter
    def material(self, material):
        """
        """
        self._material = [material]

    #    
#
class Springs:
    """
    """
    __slots__ = ('_springs', '_points', '_spring_type',
                 '_material')
    
    def __init__(self, spring_type:str, points:ClassVar, 
                 materials:MutableMapping) -> None:
        """
        """
        self._spring_type = spring_type
        self._points: ClassVar = points
        self._springs: Dict = {}
        self._material:MutableMapping = materials
    #
    def __getitem__(self, spring_name: str):
        """
        """
        return self._springs[spring_name]
    #
    def __setitem__(self, spring_name: str, corner_nodes: List[float]) -> None:
        """
        """
        self._springs[spring_name] = BasicSpring(spring_type= self._spring_type,
                                                 points = self._points)
        if corner_nodes:
            self._springs[spring_name].corner_points = corner_nodes
    #
    #
    def __len__(self) -> int:
        return len(self._springs)

    
    def __iter__(self)-> Iterable:
        """
        """
        return iter(self._springs)
        #

    def __contains__(self, value) -> bool:
        return value in self._springs
    #     
#
class SpringElement:
    """
    FE concept spring 
    
    Spring
        |_ name
        |_ number
        |_ force [f0, f1, f2,..., fn]
        |_ displacement [d0, d1, d2,..., dn]
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
    """
    #
    __slots__ = ('number', 'name', 'force', 'displacement',
                 'Pdelta', 'material', 'node', 'type',
                 'hydrodynamics', 'stress', 'section',  # need to remove this line
                 'releases', 'guidepoint', 'offsets')  # need to remove this line

    def __init__(self, Name, Number=None):
        self.number = Number
        self.name = Name
        self.force = []
        self.displacement = []
        self.Pdelta = []
        #
        self.material = []
        self.node = []
        self.type = 'spring'
        #
        self.stress = []
        self.section = []
        self.releases = []
        self.offsets = []
        self.guidepoint = []

    @property
    def length_node2node(self):
        """
        """
        # TODO: if not beam return false
        # if self.type.element ==
        #
        _dx = self.node[0].x - self.node[1].x
        _dy = self.node[0].y - self.node[1].y
        _dz = self.node[0].z - self.node[1].z
        _x1 = (_dx * _dx + _dy * _dy + _dz * _dz) ** 0.50

        return _x1

    #


#