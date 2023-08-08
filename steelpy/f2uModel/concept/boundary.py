# 
# Copyright (c) 2009-2023 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
from collections.abc import Mapping
from typing import NamedTuple, Iterator


# package imports
from .elements.boundary import BoundaryNodes
#
#
class Boundary:

    def __init__(self) -> None:
        """
        """
        self._nodes = BoundaryNodes()

    #
    @property
    def node(self):
        """"""
        return self._nodes

    @node.setter
    def node(self, values):
        """"""
        for value in values:
            self._nodes[value[0]] = value[1:]


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
    name: str|int
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
        return self._cls._supports.node[self._boundary_name]
    
    @support.setter
    def support(self, conditions):
        """ """
        self._cls._supports.node[self._boundary_name] = conditions
        #print('--')
    #
    @property
    def point(self):
        """ """
        index = self._cls._labels.index(self._boundary_name)
        return self._cls._points[index]
    #
#
#
class ConceptBoundaries(Mapping):
    
    __slots__ = ["_labels", "_supports",
                 "f2u_points", "_points"]
    
    def __init__(self):
        """
        """
        self._supports = Boundary()
        self._labels: list[str|int] = []
        self._points: list[tuple[float]] = []
    
    def __setitem__(self, support_name: int|str,
                    coordinates: list[float]|dict[str, float]) -> None:
        """
        """
        try:
            self._labels.index(support_name)
            raise Exception('boundary name {:} already exist'.format( support_name))
        except ValueError:
            self._labels.append(support_name)
            try:
                self._points.append((coordinates[0],
                                     coordinates[1],
                                     coordinates[2]))
            except IndexError:
                self._points.append((coordinates[0],
                                     coordinates[1], 0))
    
    def __getitem__(self, support_name:int) -> tuple:
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
    def df(self, df, columns:dict|None=None):
        """ """
        self._labels = df.name.tolist()
        self._points = df[["x", "y", "z"]].values.tolist()
        points = df[["name", "support"]].values.tolist()
        self._supports.node = points
        #print('---')
    #
    #