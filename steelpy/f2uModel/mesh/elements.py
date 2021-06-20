# 
# Copyright (c) 2009-2021 fem2ufo
# 


# Python stdlib imports
from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Iterable, Union
#

# package imports
from steelpy.f2uModel.mesh.sqlite.element import ElementSQL
from steelpy.f2uModel.mesh.inmemory.element import ElementInMemory
#
#
#
#
class Elements(Mapping):
    """
    """
    __slots__ = ['_elements']
    
    def __init__(self, nodes, materials, sections,
                 mesh_type:str="inmemory",
                 db_file:Union[str,None]=None):
        """
        """
        if mesh_type != "inmemory":
            self._elements = ElementSQL(db_file=db_file,
                                        db_system=mesh_type)
        else:
            self._elements = ElementInMemory(nodes=nodes,
                                             materials=materials,
                                             sections=sections)
    #
    def __setitem__(self, element_number: int,
                    parameters: Union[List[float], Dict[str, float]]) -> None:
        """
        parameters = ['beam', node1, node2, material, section, roll_angle]
        """
        if isinstance(parameters, (list, tuple)):
            element_type = parameters[0]
            node_1 = parameters[1]
            node_2 = parameters[2]
            material = parameters[3]
            section = parameters[4]
            #
            if isinstance(parameters[-1], str):
                title = parameters.pop()
            else:
                title = "NULL"
            #
            try:
                roll_angle = parameters[5]
            except IndexError:
                roll_angle = 0.0
        elif isinstance(parameters, dict):
            1/0
        else:
            raise Exception('   *** Element input format not recognized')
        #
        self._elements[element_number] = [element_type, node_1, node_2, 
                                          material, section, 
                                          roll_angle, title]
    #
    def __getitem__(self, element_number:str):
        """
        """
        return self._elements[element_number]
    #
    def __len__(self) -> float:
        return len(self._elements)

    def __iter__(self):
        """
        """
        return iter(self._elements)

    def __contains__(self, value) -> bool:
        return value in self._elements
    #
    def __str__(self, units:str="si") -> str:
        """ """
        lenght = ' m'
        space = " "
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += f"{33*space}ELEMENTS\n"
        output += "\n"
        output += (f"Element     node1    node2 {4*space} material {4*space} section {4*space} beta {4*space} title")        
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for key, element in self._elements.items():
            output += element.__str__()
        return output
    #
    @property
    def get_free_nodes(self):
        """
        find nodes not sharing elements
        """
        return self._elements.get_free_nodes
    #
    #def iter_elements(self, arraysize=1000):
    #    """
    #    """
    #    return self._elements.iter_elements(arraysize)
    #
    #
    def get_number(self, start:int=0)-> Iterable[int]:
        """
        """
        return self._elements.get_number(start)
    #
    @property
    def get_connectivities(self):
        """ """
        return self._elements.get_connectivities
    #
    #
    #
    #def update_item(self, element_number:int, item:str, value:Union[float,int]):
    #    """ """
    #    if "number" in item:
    #        self._elements[element_number].number = value
    #    else:
    #        raise IOError("item not recognized")
    #    #self._elements.update_item(element_number, item, value)
    #    #print('here')
    ##
    ##
    #def update_number(self, element_number:int, number:int):
    #    """ """
    #    self._elements[element_number].number = number
#   #
#