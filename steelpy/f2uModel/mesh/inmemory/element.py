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
from dataclasses import dataclass
import functools
#import logging
from typing import Dict, List, ClassVar, Tuple, Iterable, Union


# package imports
#from steelpy.interface.properties.beam import HydroBeamProperties
#from steelpy.f2uModel.mesh.inmemory.geometry import DirectionCosines, Eccentricities
#from steelpy.f2uModel.mesh.inmemory.geometry import Releases
from steelpy.trave3D.preprocessor.assemble import beam_Ks, trans_3d_beam, Rmatrix
#from steelpy.beam.main import Beam
from steelpy.f2uModel.mesh.operations.element import BeamBasic
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
@dataclass
class BeamElement:
    """
    """
    __slots__ = ['name', 'index', '_elements']

    def __init__(self, cls, element_name: Union[int, str]) -> None:
        """
        """
        self.index: int = cls._labels.index(element_name)
        self._elements: ClassVar = cls
        self.name: Union[int, str] = element_name
    #
    @property
    def number(self) -> int:
        return self._elements._number[self.index]

    @number.setter
    def number(self, value:int) -> None:
        """"""
        self._elements._number[ self.index ] = value
    #
    @property
    def type(self)-> str:
        """
        """
        return self._elements._type[self.index]

    #
    @property
    def connectivity(self) -> List:
        """
        """
        return self._elements._connectivity[self.index]

    @connectivity.setter
    def connectivity(self, nodes:List[int]) -> List:
        """
        """
        self._elements._connectivity[self.index] = nodes
    #
    @property
    def nodes(self) -> List:
        """
        """
        _nodes = []
        for _node in self._elements._connectivity[self.index]:
            _nodes.append(self._elements._f2u_nodes[_node])
        return _nodes

    @nodes.setter
    def nodes(self, connectivity: List) -> None:
        """
        """
        for _node in connectivity:
            self._elements._f2u_nodes[_node].sets.elements.append(self.number)
        self._elements._connectivity[self.index] = list(connectivity)

    #
    # TODO: offset should be set directly in fem file
    @property
    def offsets(self):
        """
        return eccentricities
        """
        return self.eccentricities

    @offsets.setter
    def offsets(self, eccentricities: List) -> None:
        """
        input
        eccentricities : list [eccentricities number per node]
        """
        _offsets = []
        for _item in eccentricities:
            try:
                self._elements._f2u_eccentricities[_item]
                _offsets.append(_item)
            except KeyError:
                _offsets.append(None)
        if any(_offsets):
            self._elements._eccentricities.append(_offsets)
            _index = len(self._elements._eccentricities) - 1
            self._elements._offset_index[self.index] = _index
        else:
            raise ValueError(' no valid eccentricities were given')

    @property
    def eccentricities(self):
        """
        return eccentricities
        """
        _index = self.offset_index
        _list = []
        for _item in self._elements._eccentricities[_index]:
            _list.append(self._elements._f2u_direction_cosines[_item])
        return _list

    @property
    def offset_index(self):
        """
        """
        # _index = self._elements._labels.index(self.number)
        _index_ecc = self._elements._offset_index[self.index]
        if _index_ecc == -1:
            raise ValueError(' no eccentricity defined')
        else:
            return _index_ecc

    #
    @property
    def releases(self) -> Tuple:
        """
        """
        _list = []
        for _item in self._elements._releases[self.index]:
            _list.append(self._elements._f2u_releases[_item])
        return _list

    #
    @property
    def hinges(self) -> Tuple:
        """
        """
        return self.releases

    #
    @property
    def global_offset(self) -> List:
        """
        return eccentricities in the global system
        """
        1 / 0
        return self._mesh.eccentricities[self.eccentricities]

    @global_offset.setter
    def global_offset(self, eccentricities: List) -> None:
        """
        set offset in the global system
        input
        eccentricities : list [x, y, z]
        """
        1 / 0
        self._elements._eccentricities[self.index] = eccentricities
    #
    #
    @property
    def material(self) -> List:
        """
        """
        material_name = self._elements._materials[self.index]
        #return self._elements._f2u_materials.get_item_by_number(_material_number)
        return material_name

    @material.setter
    def material(self, material_name: str) -> None:
        """
        """
        #self._elements._f2u_materials[material_name].sets.elements.append(self.number)
        #material_number = self._elements._f2u_materials[material_name].number
        #self._elements._materials[self.index] = material_number
        self._elements._materials[self.index] = material_name

    #
    #
    @property
    def section(self) -> List:
        """
        """
        section_name = self._elements._sections[self.index]
        #return self._elements._f2u_sections.get_item_by_number(section_number)
        #return self._elements._f2u_sections[section_name]
        return section_name

    @section.setter
    def section(self, section_name: str) -> None:
        """
        """
        # self._mesh.sections[section_name].sets.elements.append(self.number)
        #section_number = self._elements._f2u_sections[section_name].number
        #self._elements._sections[self.index] = section_number
        self._elements._sections[self.index] = section_name

    #
    #
    @property
    def unit_vector(self) -> List[float]:
        """
        """
        if self.type in [ 'beam', 'truss' ]:
            _node1 = self._elements._f2u_nodes[ self._elements._connectivity[ self.index ][ 0 ] ]
            _node2 = self._elements._f2u_nodes[ self._elements._connectivity[ self.index ][ -1 ] ]
            dx = _node2.x - _node1.x
            dy = _node2.y - _node1.y
            dz = _node2.z - _node1.z
            # direction cosines
            L = math.dist(_node1[:3], _node2[:3])
            l = dx / L
            m = dy / L
            n = dz / L
            return [l, m, n]
        else:
            raise IOError("no yet included")
    #
    @property
    def direction_cosines(self) -> List[float]:
        """
        """
        return self.unit_vector
    #
    #
    @property
    def length(self) -> float:
        """
        """
        node1, node2 = self.nodes
        return math.dist(node1[:3], node2[:3])
    #
    @property
    def beta(self):
        """beta angle roll"""
        return self._elements._roll_angle[self.index]
    
    @beta.setter
    def beta(self, value):
        """beta angle roll"""
        self._elements._roll_angle[self.index] = value
    #
    #
    def __str__(self) -> str:
        if (title := self._elements._title[self.index]) == "NULL":
            title = ""
        return "{:8d} {:8d} {:8d} {:>12s} {:>12s} {: 6.4f} {:>6.3f} {:>12s}\n"\
               .format(self.name, *self.connectivity,
                       self.material, self.section, self.beta,
                       self.length, title)
    #
    #
    #@property
    def beam(self) -> BeamBasic:
        """
        """
        section = self.section
        section = self._elements._f2u_sections[section]
        section = section.properties
        material = self.material
        material = self._elements._f2u_materials[material]
        beam = BeamBasic(L=self.length, 
                         E=material.E.convert('pascal').value, 
                         Iy=section.Iy, Iz=section.Iz)  
        return beam
    #
    # @properties.setter
    # def properties(self, property_number:int) -> None:
    #    """
    #    """
    #    section_number = self._sections[section_name].number
    #    self._properties[self.index] = property_number
    #
    #
    #def geometric_matrix(self, s):
    #    """
    #    """
    #    if self.type in ['beam', 'truss']:
    #        eg = beam_geom(self.length_node2node,
    #                       s)
    #        return trans_3d_beam(eg, self.r)
    #    else:
    #        raise IOError("no yet included")        
    #
    #@property
    #def mass_matrix(self):
    #    """
    #    """
    #    if self.type in ['beam', 'truss']:
    #        em = beam_mass(self.length_node2node,
    #                       self.section, self.material,
    #                       ilump=2)
    #        return trans_3d_beam(em, self.r)
    #    else:
    #        raise IOError("no yet included")     
    #
    @property
    def unit_vector(self) -> List[ float ]:
        """
        """
        node1 = self._elements._f2u_nodes[self._elements._connectivity[self.index][0]]
        node2 = self._elements._f2u_nodes[self._elements._connectivity[self.index][-1]]
        dx = node2[0] - node1[0]
        dy = node2[1] - node1[1]
        dz = node2[2] - node1[2]
        # direction cosines
        L = math.dist(node1[:3], node2[:3])
        l = dx / L
        m = dy / L
        n = dz / L
        return [l, m, n]    
    #
    @property
    def Kmatrix(self):
        """
        K stiffness matrix
        """
        #FIXME: material 
        section = self._elements._f2u_sections[self.section].properties
        material = self._elements._f2u_materials[self.material]
        # solve K matrix
        K = beam_Ks(self.length,
                    section.area, section.J,
                    section.Iy, section.Iz,
                    material.E.convert("pascal").value, 
                    material.G.convert("pascal").value,
                    section.area, section.area)
        return trans_3d_beam(K, self.R)
        #return self.trans3d(K)
    #
    @property
    def R(self):
        """
        Rotation matrix
        """
        if self.type in ['beam', 'truss']:
            return Rmatrix(*self.unit_vector, self.beta)
        else:
            raise IOError("no yet included")
    #
    #def trans3d(self,M):
    #    """ 3-d coordinate transformations"""
    #    return trans_3d_beam(M, self.R)
    #
    @property
    def DoF(self):
        """
        """
        dof = []
        for _node in self._elements._connectivity[self.index]:
            dof.append(self._elements._f2u_nodes[_node].index)
        return dof
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
