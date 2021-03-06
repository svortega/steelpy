# 
# Copyright (c) 2009-2021 fem2ufo
#

# Python stdlib imports
import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Dict, List, ClassVar, Tuple, Iterable, Union


# package imports
from steelpy.process.units.main import Units
#import steelpy.process.io_module.text as common
#
#
#
@dataclass
class Element:
    """
    """
    __slots__ = ['name', 'index', '_cls']

    def __init__(self, cls, index) -> None:
        """
        """
        self.index: int = index
        self._cls: ClassVar = cls
        self.name: Union[int, str] = cls._labels[index]
    #
    @property
    def number(self) -> int:
        return self._cls._number[self.index]

    @number.setter
    def number(self, value:int) -> None:
        """"""
        self._cls._number[self.index] = value
    #
    @property
    def type(self)-> str:
        """
        """
        return self._cls._type[self.index]

    #
    @property
    def connectivity(self) -> List:
        """
        """
        #conn = self._cls._connectivity[self.index]
        jnt = []
        for conn in self._cls._connectivity[self.index]:
            jnt.append(self._cls.f2u_points[conn])
        return jnt
    #
    @property
    def material(self) -> List:
        """
        """
        material_name = self._cls._materials[self.index]
        return self._cls.f2u_materials[material_name]

    @material.setter
    def material(self, material) -> None:
        """
        """
        try:
            self._cls._materials[self.index] = material.name
        except AttributeError:
            self._cls.f2u_materials[material]
            self._cls._materials[self.index] = material
    #
    @property
    def section(self) -> List:
        """
        """
        section_name = self._cls._sections[self.index]
        return self._cls.f2u_sections[section_name]

    @section.setter
    def section(self, section) -> None:
        """
        """
        try:
            self._cls._sections[self.index] = section.name
        except AttributeError:
            self._cls.f2u_sections[section]
            self._cls._sections[self.index] = section
    #
    #@property
    #def segment_length(self):
    #    """
    #    segment length from start node
    #    """
    #    return self._cls._segment[self.index] * f2u_units.m
#
#
@dataclass
class SegmentedBeam:
    """ """
    __slots__ = ['indices', '_cls', '_step_name']
    
    def __init__(self, cls, step_name, indices):
        """
        """
        self._cls = cls
        self.indices = indices
        self._step_name = step_name
    #
    #
    @property
    def name(self):
        """ """
        index = self._get_index()
        return self._cls._step_label[index[0]]
    #
    @property
    def length(self):
        """ """
        index = self._get_index()
        return self._cls._segment[index[0]] * self._cls.f2u_units.m

    @length.setter
    def length(self, item):
        """ """
        index = self.indices[-1]
        length = item.value
        self._cls._segment[index] = length
        new_index = self._cls._duplicate_element(index)
        self._cls._step_label[new_index] = self._step_name
    #
    #
    @property
    def material(self):
        """ """
        index  = self._get_index()
        name = self._cls._materials[index[0]]
        return self._cls.f2u_materials[name]

    @material.setter
    def material(self, value):
        """ """
        index = self.indices[-1]
        try:
            self._cls._materials[index] = value.name
        except AttributeError:
            self._cls.f2u_materials[value]
            self._cls._materials[index] = value
    #
    @property
    def section(self):
        """ """
        index = self._get_index()
        return self._cls.f2u_sections[self._cls._sections[index[0]]]

    @section.setter
    def section(self, value):
        """ """
        index = self.indices[-1]
        try:
            self._cls._sections[index] = value.name
        except AttributeError:
            self._cls.f2u_sections[value]
            self._cls._sections[index] = value
    #
    #
    @property
    def _mesh(self):
        """
        """
        index = self._get_index()
        return self._cls._mesh[index[0]]
    
    @_mesh.setter
    def _mesh(self, element:Union[str,int]):
        """
        """
        index = self._get_index()
        self._cls._mesh[index[0]] = element
    #     
    #
    def _get_index(self):
        """ """
        index = [_index for _index in self.indices
                 if self._step_name == self._cls._step_label[_index]]
        if not index:
            raise IndexError
        return index
    #
    #def __iter__(self):
    #    """ """
    #    for index in self.indices[1:]:
    #        yield Element(self._cls, index)
    #   
#
#
class Steps:
    __slots__ = ['indices', '_cls']
    
    def __init__(self, cls):
        """
        """
        self._cls = cls._cls
        self.indices = [index for index, name in enumerate(self._cls._labels) 
                        if cls.name == name]
    
    def __setitem__(self, step_name, coord):
        """ """
        index = self.indices[-1]
        nodes = self._cls._connectivity[index]
        node1 = self._cls.f2u_points[nodes[0]]
        node1 = self._cls.f2u_points._get_coordinates(node1)
        node2 = self._cls.f2u_points._get_coordinates(coord)
        length = math.dist(node1, node2)
        self._cls._segment[index] = length
        new_index = self._cls._duplicate_element(index)
        self._cls._step_label[new_index] = step_name
    #    
    # 
    def __getitem__(self, step_name: Union[int,str]) -> Tuple:
        """
        step_name : node number
        """
        return SegmentedBeam(self._cls, step_name, self.indices)
    #
    def __iter__(self):
        """ """
        for index in self.indices:
            step_name = self._cls._step_label[index]
            yield SegmentedBeam(self._cls, step_name, self.indices)
    #
    def __len__(self):
        return len(self.indices)
#
#
@dataclass
class Beam(Element):
    __slots__ = ['name', 'index', '_cls', '_steps']

    def __init__(self, cls, element_index) -> None:
        """
        """
        Element.__init__(self, cls, element_index)
        self._steps = Steps(self)
    #
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
                self._cls._f2u_eccentricities[_item]
                _offsets.append(_item)
            except KeyError:
                _offsets.append(None)
        if any(_offsets):
            self._cls._eccentricities.append(_offsets)
            _index = len(self._cls._eccentricities) - 1
            self._cls._offset_index[self.index] = _index
        else:
            raise ValueError(' no valid eccentricities were given')

    @property
    def eccentricities(self):
        """
        return eccentricities
        """
        _index = self.offset_index
        _list = []
        for _item in self._cls._eccentricities[_index]:
            _list.append(self._cls._f2u_direction_cosines[_item])
        return _list

    @property
    def offset_index(self):
        """
        """
        # _index = self._cls._labels.index(self.number)
        _index_ecc = self._cls._offset_index[self.index]
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
        for _item in self._cls._releases[self.index]:
            _list.append(self._cls._f2u_releases[_item])
        return _list

    #
    @property
    def hinges(self) -> Tuple:
        """
        """
        return self.releases

    #
    @property
    def type(self) -> str:
        """
        """
        return self._cls._type[self.index]
    
    @type.setter
    def type(self, beam_type:str):
        """
        """
        self._cls._type[self.index] = beam_type
    #
    @property
    def beta(self):
        """beta angle roll"""
        return self._cls._roll_angle[self.index]
    
    @beta.setter
    def beta(self, value):
        """beta angle roll"""
        self._cls._roll_angle[self.index] = value
    #
    #
    @property
    def step(self):
        """
        """
        return self._steps
    #
    @property
    def length(self) -> Units:
        """
        """
        _nodes = self.connectivity
        length = math.dist([_nodes[0].x, _nodes[0].y, _nodes[0].z], 
                           [_nodes[1].x, _nodes[1].y, _nodes[1].z])
        return length * self._cls.f2u_units.m
    #
    def find_coordinate(self, node_distance:float, node_end:int=0) -> Tuple:
        """
        """
        _node = self.connectivity
        _nodeNo3 = [0, 0, 0]
        #
        if math.isclose(node_distance.value, 0, rel_tol=0.01):
        #if distance <= 0.0001:
            _nodeNo3[0] = _node[node_end].x
            _nodeNo3[1] = _node[node_end].y
            _nodeNo3[2] = _node[node_end].z
        else:
            if node_end == 1:
                _v1 = (_node[0].x - _node[1].x)
                _v2 = (_node[0].y - _node[1].y)
                _v3 = (_node[0].z - _node[1].z)
            else:
                _v1 = (_node[1].x - _node[0].x)
                _v2 = (_node[1].y - _node[0].y)
                _v3 = (_node[1].z - _node[0].z)
            #
            _norm = (_v1**2 + _v2**2 + _v3**2)**0.50
            _v1 /= _norm
            _v2 /= _norm
            _v3 /= _norm
            _nodeNo3[0] = (_node[node_end].x + _v1 * node_distance)
            _nodeNo3[1] = (_node[node_end].y + _v2 * node_distance)
            _nodeNo3[2] = (_node[node_end].z + _v3 * node_distance)
        #
        #Coordinates = get_coordinate_system(_node[0].system)
        #return Coordinates(*_nodeNo3)
        return _nodeNo3
    #
    #
    @property
    def unit_vector(self) -> List[float]:
        """
        """
        _node1, _node2 = self.connectivity
        dx = _node2.x.value - _node1.x.value
        dy = _node2.y.value - _node1.y.value
        dz = _node2.z.value - _node1.z.value
        # direction cosines
        L = math.dist([_node1.x.value, _node1.y.value, _node1.z.value], 
                      [_node2.x.value, _node2.y.value, _node2.z.value])
        l = dx / L
        m = dy / L
        n = dz / L
        return [l, m, n]  
#
#
class ConceptElements(Mapping):
    """
    element[name] = [name, connectivity, material, section, type, group]
    connectivity[number] = [name, node1, node2,..., nodei]
    """
    __slots__ = ['_labels', '_number','_type', '_connectivity', '_element', '_element_type',
                 '_sections', '_materials', '_mesh', '_releases',  '_roll_angle', 
                 '_direction_cosines', '_eccentricities', '_offset_index', 
                 '_segment', '_step_label', 
                 'f2u_points', 'f2u_materials', 'f2u_sections', 'f2u_units']

    
    def __init__(self, element_type:str, points, 
                 materials, sections) -> None: # properties
        """
        Manages f2u elements
        """
        #global f2u_materials, f2u_sections, f2u_units, f2u_points
        self.f2u_materials = materials
        self.f2u_sections =  sections
        self.f2u_units = Units()
        self.f2u_points = points
        #
        self._element_type = element_type
        #
        self._labels:List[Union[str,int]] = []
        self._sections:List[Union[str,int]] = []
        self._materials:List[Union[str,int]] = []
        self._roll_angle:List[float] = []
        self._number:List[int] = []
        #
        self._type: List[str] = []
        self._connectivity: List = []
        self._segment:List[float] = []
        self._step_label:List[Union[str,int]] = []
        self._mesh:List[Union[str,int]] = []
        #
        #self._direction_cosines = array('i', [])
        #self._offset_index = array('i', [])        
        #self._eccentricities: List = []
        #self._releases: List = []        
    #
    def __setitem__(self, element_name: Union[int, int], 
                    parameters: Union[List[float], Dict[str, float]]) -> None:
        """
        farg = [name, connectivity, material, section, type, group]
        """
        try:
            index = self._labels.index(element_name)
            raise Exception('{:} {:} already exist'.format(self._element_type, element_name))
        except ValueError:
            # default
            self._labels.append(element_name)
            self._roll_angle.append(0.0)
            index = self._labels.index(element_name)
            self._number.append(index)
            self._type.append(self._element_type)
            # set connectivity
            try:
                node_1 = self.f2u_points.get_point_name(parameters[0])
            except IOError:
                node_1 = self.f2u_points.get_new_point(parameters[0])

            try:
                node_2 = self.f2u_points.get_point_name(parameters[1])
            except IOError:
                node_2 = self.f2u_points.get_new_point(parameters[1])
            #
            self._connectivity.append([node_1, node_2])
            # set blank data
            self._sections.append(-1)
            self._materials.append(-1)
            #
            # set deafult material, section
            if self.f2u_materials.default:
                self._materials[index] = self.f2u_materials.default
            #
            if self.f2u_sections.default:
                self._sections[index] = self.f2u_sections.default
            #
            self._segment.append(0)
            self._step_label.append(element_name)
            self._mesh.append(-1)
            # to be defined
            #self._properties.append(-1)
            #self._offset_index.append(-1)
            #self._direction_cosines.append(-1)
            # should be in demand
            #self._eccentricities.append([])
            #self._releases.append([])
    
    def __getitem__(self, element_name: Union[str,int]) -> ClassVar:
        """
        """
        try:
            _index = self._labels.index(element_name)
            _indices = [index for index, name in enumerate(self._labels) 
                        if element_name == name]
            _index = _indices[0]
            if self._type[_index] in ["beam", "truss"]:
                return Beam(self, _index)
            else:
                return Element(self, _index)
        except ValueError:
            raise IndexError(' ** element {:} does not exist'.format(element_name))

    #
    def __iter__(self) -> Iterable:
        """
        """
        labels = list(dict.fromkeys(self._labels))
        return iter(labels)

    def __delitem__(self, element_name: Union[str,int]) -> None:
        """
        """
        try:
            _nodes_empty = []
            _index = self._labels.index(element_name)
            # remove element form node's set
            for _item in self._connectivity[_index]:
                _node = f2u_points[_item]
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
            self._segment.pop(i)
            self._step_label.pop(i)
            self._mesh.pop(i)            
            #
            #offset_index = self._offset_index[i]
            #self._offset_index.pop(i)
            ## FIXME
            #if offset_index != -1:
            #    self._eccentricities.pop(offset_index)
            #    1/0
            # delete empty nodes
            for _node_number in _nodes_empty:
                del f2u_points[_node_number]
            #
            # FIXME: number should be updated according new index
            self._number.pop(i)
        except ValueError:
            raise KeyError('    *** warning element {:} does not exist'
                           .format(element_name))

    def __contains__(self, value) -> bool:
        return value in self._labels

    def __len__(self) -> float:
        labels = list(dict.fromkeys(self._labels))
        return len(labels)
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
    def _duplicate_element(self, index) -> int:
        """
        """
        element_name = self._labels[index]
        step = index + 1
        self._labels.insert(step, element_name)
        self._roll_angle.insert(step, self._roll_angle[index])
        self._type.insert(step, self._type[index])
        # set connectivity 
        self._connectivity.insert(step, self._connectivity[index])
        # set blank data
        self._sections.insert(step, self._sections[index])
        self._materials.insert(step, self._materials[index])
        #
        self._segment.insert(step, 0)
        self._step_label.insert(step, -1)
        self._mesh.insert(step, -1)
        # renumber
        self._number.insert(step, self._labels.index(element_name))
        self._number= [index for index, _ in enumerate(self._labels)]
        return step
    #
    #
#
#
#
