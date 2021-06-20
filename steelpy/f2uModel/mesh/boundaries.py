# 
# Copyright (c) 2009-2021 fem2ufo
# 


# Python stdlib imports
from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Iterable, Union
#

# package imports
from steelpy.f2uModel.mesh.sqlite.boundary import BoundaryNodeSQL
from steelpy.f2uModel.mesh.inmemory.boundary import BoundaryNodes
#
#
#
#
class Boundaries:
    """
    """
    __slots__ = ['_nodes']
    
    def __init__(self, mesh_type:str,
                 db_file:Union[str,None]):
        """
        """
        if mesh_type != "inmemory":
            self._nodes = BoundaryNodeSQL(db_file=db_file)
        else:
            self._nodes = BoundaryNodes()
    #
    #def __setitem__(self, node_name: int,
    #                coordinates: Union[List[float], Dict[str, float]]) -> None:
    #    """
    #    """
    #    self._boundaries[node_name] = coordinates
    #
    #def __getitem__(self, node_name:str):
    #    """
    #    """
    #    #print('<-- mat' )
    #    return self._boundaries[node_name]
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
    def __str__(self, units:str="si") -> str:
        """ """
        lenght = ' m'
        space = " "
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += f"{33*space}BOUNDARIES\n"
        output += "\n"
        output += (f"Node {14*space} x {6*space} y {6*space} z {5*space} mx {5*space} my {5*space} mz {5*space} title")
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for key, node in self._nodes.items():
            output += node.__str__()
        return output
#
#