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