# 
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
from collections.abc import Mapping
from typing import NamedTuple, Tuple, Union, List, Dict, Iterator


# package imports
from steelpy.f2uModel.mesh.boundary import BoundaryNodes

#
#
class BoundaryItem(NamedTuple):
    """
    """
    x: float
    y: float
    z: float
    rx: float
    ry: float
    rz: float
    number: int
    name: Union[str, int]


#
#
#
class BoundaryType:
    
    __slots__ = ["_cls", "_boundary_name"]
    
    def __init__(self, cls, boundary_name):
        """
        """
        self._cls = cls
        self._boundary_name = boundary_name
    #
    @property
    def support(self):
        """ """
        self._cls._supports.node[self._boundary_name]
    
    @support.setter
    def support(self, conditions):
        """ """
        self._cls._supports.node[self._boundary_name] = conditions
        #print('--')
    #
    #
#
#
class Boundaries(Mapping):
    
    __slots__ = ["_labels", "_supports",
                 "f2u_points", "_points"]
    
    def __init__(self, points, boundary_points):
        """
        """
        global f2u_points
        f2u_points = points
        #
        self._supports = boundary_points  # BoundaryNodes()
        #
        self._labels: List[Union[str,int]] = []
        self._points: List[Union[str,int]] = []
    
    def __setitem__(self, support_name: Union[int, str],
                    coordinates: Union[List[float], Dict[str, float]]) -> None:
        """
        """
        try:
            self._labels.index(support_name)
            raise Exception('boundary name {:} already exist'.format( support_name))
        except ValueError:
            self._labels.append(support_name)
            node = f2u_points.get_point_name(coordinates)
            self._points.append(node)
    
    def __getitem__(self, support_name: int) -> Tuple:
        """
        node
        """
        return BoundaryType(cls=self, boundary_name=support_name)
    
    #
    #
    def __len__(self) -> float:
        return len(self._labels)

    def __iter__(self) -> Iterator:
        """
        """
        return iter(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    #