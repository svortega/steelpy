# 
# Copyright (c) 2009-2021 fem2ufo
#

# Python stdlib imports
import math 
from itertools import chain
from array import array
#from collections import OrderedDict
from collections import Counter
from collections.abc import Mapping
import functools
#import logging
from typing import Dict, List, ClassVar, Tuple, Iterable, Union


# package imports
#from steelpy.interface.properties.beam import HydroBeamProperties
#from steelpy.f2uModel.mesh.inmemory.geometry import DirectionCosines, Eccentricities
#from steelpy.f2uModel.mesh.inmemory.geometry import Releases
#from steelpy.f2uModel.mesh.operations.elements  import (beam_Klocal, trans_3d_beam, Rmatrix, trans3Dbeam, beam_stiffness)
#from steelpy.beam.main import Beam
from steelpy.f2uModel.mesh.operations.elements import BeamElement
#
@functools.lru_cache(maxsize=2048)
def get_node_end(_number_nodes: int, line: int) -> Tuple[int, int]:
    """
    """
    if _number_nodes == 3:
        # print('triangle')
        if line == 1:
            return 1, 2
        elif line == 2:
            return 0, 2
        else:
            return 0, 1
    else:
        # print('quad')
        if line == 1:
            return 0, 1
        elif line == 2:
            return 1, 2
        elif line == 3:
            return 2, 3
        else:
            return 3, 0
#
#
#
class ElementInMemory(Mapping):
    __slots__ = ['_labels', '_number','_type', '_connectivity', '_element',
                 '_sections', '_materials', '_mesh', '_releases',  '_roll_angle',
                 '_direction_cosines', '_eccentricities', '_offset_index',
                 '_f2u_nodes', '_f2u_materials', '_f2u_sections', '_title']

    
    def __init__(self, nodes, materials, sections) -> None:
        """
        Manages f2u elements
        """
        #
        # f2u classes
        self._f2u_nodes = nodes
        self._f2u_materials: ClassVar = materials
        self._f2u_sections: ClassVar = sections
        #self._f2u_releases: ClassVar = releases
        # TODO: element's code check parameters missing
        #self._design_parameters: ClassVar = CodeCheck()
        #self._f2u_direction_cosines: ClassVar = DirectionCosines()
        #self._f2u_eccentricities: ClassVar = Eccentricities()
        #self._f2u_releases: ClassVar = Releases()
        #
        self._labels: array = array('i', [])
        self._number: array = array('i', [])
        self._title: List = []
        #self._materials = array('i', [])
        self._sections:List[Union[str,int]] = []
        self._materials:List[Union[str,int]] = []
        #
        self._roll_angle: array = array('f', [])
        #self._direction_cosines = array('i', [])
        #self._offset_index = array('i', [])
        #
        self._type: List = []
        self._connectivity: List = []
        #self._eccentricities: List = []
        self._releases: List = []

    #
    def __setitem__(self, element_number: int, parameters: List) -> None:
        """
        parameters = ['beam', node1, node2, material, section, roll_angle]
        """
        try:
            self._labels.index(element_number)
            raise Exception('element {:} already exist'.format(element_number))
        except ValueError:
            # default
            self._labels.append(element_number)
            self._number.append(self._labels.index(element_number))
            #
            #if isinstance(parameters, (list, tuple)):
            #    element_type = parameters[0]
            #    node_1 = parameters[1]
            #    node_2 = parameters[2]
            #    material = parameters[3]
            #    section = parameters[4]
            #    try:
            #        self._roll_angle.append(parameters[5])
            #    except IndexError:
            #        self._roll_angle.append(0.0)
            #elif isinstance(parameters, dict):
            #    pass
            #else:
            #    raise Exception('   *** Element input format not recognized')
            #
            #_hinge = self.spcl_bc(element_type)
            #self._releases.append([_hinge, _hinge])
            self._releases.append([])
            #
            self._type.append(parameters[0])
            self._connectivity.append(parameters[1:3])
            self._materials.append(parameters[3])
            self._sections.append(parameters[4])
            self._roll_angle.append(parameters[5])
            self._title.append(parameters[6])
            #self._offset_index.append(-1)
            #self._direction_cosines.append(-1)
            # should be in demand
            #self._eccentricities.append([])
    
    def __getitem__(self, element_number: int) -> ClassVar:
        """
        """
        try:
            return BeamElement(self, element_number)
        except ValueError:
            raise IndexError(' ** element {:} does not exist'.format(element_number))

    #
    def __iter__(self) -> Iterable:
        """
        """
        return iter(self._labels)

    def __delitem__(self, element_number: Union[str,int]) -> None:
        """
        """
        try:
            _nodes_empty = []
            _index = self._labels.index(element_number)
            # remove element form node's set
            for _item in self._connectivity[_index]:
                _node = self._f2u_nodes[_item]
                try:
                    _node.sets.elements.remove(element_number)
                except ValueError:
                    pass
                # capture nodes with no more elements
                if not _node.sets.elements:
                    _nodes_empty.append(_node.number)
            # remove element form list
            i = self._labels.index(element_number)
            self._labels.pop(i)
            self._sections.pop(i)
            self._materials.pop(i)
            self._roll_angle.pop(i)
            #self._direction_cosines.pop(i)
            self._connectivity.pop(i)
            #
            #offset_index = self._offset_index[i]
            #self._offset_index.pop(i)
            ## FIXME
            #if offset_index != -1:
            #    self._eccentricities.pop(offset_index)
            #    1/0
            # delete empty nodes
            for _node_number in _nodes_empty:
                del self._f2u_nodes[_node_number]
            #
            # FIXME: number should be updated according new index
            self._number.pop( i )
        except ValueError:
            raise KeyError('    *** warning element {:} does not exist'
                           .format(element_number))

    def __contains__(self, value) -> bool:
        return value in self._labels

    def __len__(self) -> float:
        return len(self._labels)
    #
    #
    def get_number(self, start:int=0)-> Iterable[int]:
        """
        """
        try:
            n = max(self._labels)
        except ValueError:
            n = start
        #
        while True:
            n += 1
            yield n
    #
    #def iter_elements(self, arraysize=1000):
    #    """
    #    """
    #    for element_name in self._labels:
    #        yield BeamElement(self, element_name)
    #
    def spcl_bc(self, element_type:str):
        """
        Impose condition for special global shapes
        1 fix (0 free)
        """
        if element_type == 'truss':
            return [1,1,1,1,0,0]
        else: # beam
            return [1,1,1,1,1,1]
    #
    @property
    def get_free_nodes(self):
        """
        find nodes not sharing elements
        """       
        #columns = list(zip(*self._connectivity))
        #column = columns[0]
        #column
        flat = list(chain.from_iterable(self._connectivity))
        return [k for k, v in Counter(flat).items() if v == 1]
    #
    #
    #
    #def update_item(self, element_number:int, item:str, value:Union[float,int]):
    #    """ """
    #    _index = self._labels.append(element_number)
    #    print('here')
    #
    @property
    def get_connectivities(self):
        """ """
        return self._connectivity
#
#
#
