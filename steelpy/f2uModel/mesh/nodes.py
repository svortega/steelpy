# 
# Copyright (c) 2009-2021 fem2ufo
# 


# Python stdlib imports
from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Iterable, Union
#

# package imports
from steelpy.f2uModel.mesh.sqlite.nodes import NodeSQL
from steelpy.f2uModel.mesh.inmemory.nodes import NodesInMemory
from steelpy.f2uModel.mesh.operations.nodes import node_renumbering
#
#
#
#
class Nodes(Mapping):
    """
    """
    __slots__ = ['_nodes']
    
    def __init__(self, mesh_type:str,
                 db_file:Union[str,None]):
        """
        """
        if mesh_type != "inmemory":
            self._nodes = NodeSQL(db_file=db_file,
                                  db_system=mesh_type)
        else:
            self._nodes = NodesInMemory()
    #
    def __setitem__(self, node_name: int,
                    coordinates: Union[List[float], Dict[str, float]]) -> None:
        """
        """
        self._nodes[node_name] = coordinates
    #
    def __getitem__(self, node_name:int):
        """
        """
        #print('<-- mat' )
        return self._nodes[node_name]
    #
    def __len__(self) -> float:
        return len(self._nodes)

    def __iter__(self):
        """
        """
        _nodes = sorted(self._nodes)
        return iter(_nodes)

    def __contains__(self, value) -> bool:
        return value in self._nodes
    #
    def __str__(self, units:str="si") -> str:
        """ """
        lenght = ' m'
        space = " "
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += f"{33*space}NODES\n"
        output += "\n"
        output += (f"Node {12*space} x  [{lenght}] {4*space} y  [{lenght}] {4*space} z  [{lenght}] ")        
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"        
        for key, node in self._nodes.items():
            output += node.__str__()
        return output
    #
    def get_node_name(self, coordinates, tol: float = 0.01):
        """ """
        return self._nodes.get_point_name(coordinates, tol)
    #
    #
    def get_new_node(self, coordinates):
        """ """
        return self._nodes.get_new_point(coordinates)
    #
    #
    def update_item(self, node_name:int, item:str, value:Union[float,int]):
        """ """
        self._nodes.update_item(node_name, item, value)
        #print('---')
    #
    #
    #def update_number(self, node_name:int, number:int):
    #    """ """
    #    self._nodes.update_number(node_name, number)
    #
    #
    def renumbering(self, elements):
        """ """
        new_node_list = node_renumbering(self._nodes, elements)
        #new_node_list = list(reversed(new_node_list))
        self._nodes.renumbering(new_node_list)
        #print('end')
    #
    #
#
#