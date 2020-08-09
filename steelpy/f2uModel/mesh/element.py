# 
# Copyright (c) 2009-2020 fem2ufo
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
from steelpy.f2uModel.mesh.geometry import DirectionCosines, Eccentricities
from steelpy.f2uModel.mesh.geometry import Releases
from steelpy.trave3D.preprocessor.assemble import (beam_stiffness, beam_Ks,
                                                   trans_3d_beam, Rmatrix)
# from steelpy.properties.codecheck.codecheck import CodeCheck

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
class Element:
    """
    """
    __slots__ = ['name', 'index', '_elements']

    def __init__(self, cls, element_name) -> None:
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
    #
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
        #
        # try:
        _node1 = self._elements._nodes[self._elements._connectivity[self.index][0]]
        _node2 = self._elements._nodes[self._elements._connectivity[self.index][-1]]
        # except Exception as error:
        #    print(error)
        #    return None
        #
        # TODO: if not beam return false
        # if self.offsets:
        try:
            _offsets = self.offsets
            _case = 2
            try:
                node1 = _offsets[0]
                _case = _offsets[0].system
            except AttributeError:
                node1 = [0, 0, 0]
            #
            try:
                node2 = _offsets[1]
                _case = _offsets[1].system
            except AttributeError:
                node2 = [0, 0, 0]

            # global case
            if _case == 'global':
                _dx = ((_node1.x + node1[0])
                       - (_node2.x + node2[0]))
                _dy = ((_node1.y + node1[1])
                       - (_node2.y + node2[1]))
                _dz = ((_node1.z + node1[2])
                       - (_node2.z + node2[2]))
                #
                _x1 = (_dx * _dx + _dy * _dy + _dz * _dz) ** 0.50
            # local case
            else:
                # beam
                #_dx = _node1.x - _node2.x
                #_dy = _node1.y - _node2.y
                #_dz = _node1.z - _node2.z
                #_x1 = (_dx * _dx + _dy * _dy + _dz * _dz) ** 0.50
                #_x1 = self.length_node2node
                # considering shortening only (ASAS)
                # TODO : this should work for any format
                _x1 = self.length_node2node - node1[0] + node2[0]
            #
            return _x1
        except ValueError:
            return self.length_node2node
        #

    @property
    def length_node2node(self) -> float:
        """
        """
        _node1 = self._elements._f2u_nodes[self._elements._connectivity[self.index][0]]
        _node2 = self._elements._f2u_nodes[self._elements._connectivity[self.index][-1]]
        #_dx = _node1.x - _node2.x
        #_dy = _node1.y - _node2.y
        #_dz = _node1.z - _node2.z
        #dist2 = (_dx * _dx + _dy * _dy + _dz * _dz)**0.50
        return math.dist(_node1[:3], _node2[:3])

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
    #@property
    #def properties(self) -> List:
    #    """
    #    """
    #    # property_number = self._elements._properties[self.index]
    #    return self._elements._f2u_properties[self.index]
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
    #@property
    #def Kmatrix(self):
    #    """
    #    K stiffness matrix
    #    """
    #    if self.type in ['beam', 'truss']:
    #        section = self._elements._f2u_sections[self.section]
    #        material = self._elements._f2u_materials[self.material]._get_data()
    #        #beam_length: float, area, J, Iy, Iz,
    #        #                   emod:float, gmod:float
    #        prop = section._get_properties()
    #        #ek = beam_stiffness(self.length_node2node,
    #        #                    self.section.area.value, 
    #        #                    self.section.J.value,
    #        #                    self.section.Iy.value,
    #        #                    self.section.Iz.value,
    #        #                    self.material.E.convert("pascal").value,
    #        #                    self.material.G.convert("pascal").value)
    #        ek = beam_Ks(self.length_node2node,
    #                     prop.area, prop.J, prop.Iy, prop.Iz,
    #                     material.E, material.G,
    #                     prop.area, prop.area)
    #        return trans_3d_beam(ek, self.R)
    #    else:
    #        raise IOError("no yet included")
    ##
    #@property
    #def R(self):
    #    """
    #    Rotation matrix
    #    """
    #    if self.type in ['beam', 'truss']:
    #        beta = self._elements._roll_angle[self.index]
    #        return Rmatrix(*self.unit_vector, beta)
    #    else:
    #        raise IOError("no yet included")        
    #
    @property
    def DoF(self):
        """
        """
        dof = []
        for _node in self._elements._connectivity[self.index]:
            dof.append((self._elements._f2u_nodes[_node].number) * 6)
        return dof
#
#
class Elements(Mapping):
    """
    element[name] = [name, connectivity, material, section, type, group]
    connectivity[number] = [name, node1, node2,..., nodei]
    """
    __slots__ = ('_labels', '_number','_type', '_connectivity', '_element',
                 '_sections', '_materials', '_mesh', '_releases',  '_roll_angle',
                 '_direction_cosines', '_eccentricities', '_offset_index',
                 '_f2u_nodes', '_f2u_materials', '_f2u_sections')
                 #'_f2u_direction_cosines', '_f2u_eccentricities', '_f2u_releases')

    
    def __init__(self) -> None: # properties
        """
        Manages f2u elements
        """
        #
        # f2u classes
        #self._f2u_nodes = nodes
        #self._f2u_materials: ClassVar = materials
        #self._f2u_sections: ClassVar = sections
        #self._f2u_releases: ClassVar = releases
        # TODO: element's code check parameters missing
        #self._design_parameters: ClassVar = CodeCheck()
        #self._f2u_direction_cosines: ClassVar = DirectionCosines()
        #self._f2u_eccentricities: ClassVar = Eccentricities()
        #self._f2u_releases: ClassVar = Releases()
        #
        self._labels = array('i', [])
        self._number = array('i', [])
        #self._sections = array('i', [])
        #self._materials = array('i', [])
        self._sections:List[Union[str,int]] = []
        self._materials:List[Union[str,int]] = []
        #
        self._roll_angle = array('f', [])
        #self._direction_cosines = array('i', [])
        #self._offset_index = array('i', [])
        #
        self._type: List = []
        self._connectivity: List = []
        #self._eccentricities: List = []
        self._releases: List = []

    #
    def __setitem__(self, element_name: Union[int, int], 
                    parameters: Union[List[float], Dict[str, float]]) -> None:
        """
        farg = [name, connectivity, material, section, type, group]
        """
        try:
            self._labels.index(element_name)
            raise Exception('element {:} already exist'.format(element_name))
        except ValueError:
            # default
            self._labels.append(element_name)
            self._roll_angle.append(0.0)
            self._number.append(self._labels.index(element_name))
            #
            if isinstance(parameters, (list, tuple)):
                element_type = parameters[0]
                node_1 = parameters[1]
                node_2 = parameters[2]
                material = parameters[3]
                section = parameters[4]
            elif isinstance(parameters, dict):
                pass
            else:
                raise Exception('   *** Element input format not recognized')
            #
            #_hinge = self.spcl_bc(element_type)
            #self._releases.append([_hinge, _hinge])
            self._releases.append([])
            #
            self._type.append(element_type)
            self._connectivity.append([node_1, node_2])
            self._sections.append(section)
            self._materials.append(material)
            #self._properties.append(-1)
            #self._offset_index.append(-1)
            #self._direction_cosines.append(-1)
            # should be in demand
            #self._eccentricities.append([])
    
    def __getitem__(self, element_name: Union[str,int]) -> ClassVar:
        """
        """
        try:
            #_index = self._labels.index(element_name)
            return Element(self, element_name)
        except ValueError:
            raise IndexError(' ** element {:} does not exist'.format(element_name))

    #
    def __iter__(self) -> Iterable:
        """
        """
        return iter(self._labels)

    def __delitem__(self, element_name: Union[str,int]) -> None:
        """
        """
        try:
            _nodes_empty = []
            _index = self._labels.index(element_name)
            # remove element form node's set
            for _item in self._connectivity[_index]:
                _node = self._f2u_nodes[_item]
                try:
                    _node.sets.elements.remove(element_name)
                except ValueError:
                    pass
                # capture nodes with no more elements
                if not _node.sets.elements:
                    _nodes_empty.append(_node.number)
            # remove element form list
            i = self._labels.index(element_name)
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
                           .format(element_name))

    def __contains__(self, value) -> bool:
        return value in self._labels

    def __len__(self) -> float:
        return len(self._labels)
    #
    #
    #@property
    #def direction_cosines(self) -> ClassVar:
    #    """
    #    """
    #    return self._f2u_direction_cosines
    #
    #@property
    #def unit_vectors(self) -> ClassVar:
    #    """
    #    """
    #    return self._f2u_direction_cosines
    #
    #@property
    #def eccentricities(self) -> ClassVar:
    #    """
    #    """
    #    return self._f2u_eccentricities
    #
    #@property
    #def offsets(self) -> ClassVar:
    #    """
    #    """
    #    return self._f2u_eccentricities
    #
    #@property
    #def releases(self) -> ClassVar:
    #    """
    #    """
    #    return self._f2u_releases
    #
    #
    def get_number(self, start:int=1)-> Iterable[int]:
        """
        """
        try:
            n = max(self._number)
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
    #
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
