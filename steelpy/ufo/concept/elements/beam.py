# 
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
import math
#from collections.abc import Mapping
from dataclasses import dataclass
import re
from operator import sub #, add

# package imports
from steelpy.ufo.utils.beam import BeamBasic
from steelpy.ufo.utils.node import NodePoint
from steelpy.ufo.utils.element import get_beam_df
from steelpy.ufo.mesh.process.brotation import Rmatrix2
from steelpy.utils.units.main import Units

#
from steelpy.sections.utils.shape.main import is_section
from steelpy.material.utils.operations import is_material
from steelpy.utils.dataframe.main import DBframework
#
#
@dataclass
class ElementConcept:
    """
    """
    __slots__ = ['name', 'index', '_cls']

    def __init__(self, cls, index) -> None:
        """
        """
        self.index: int = index
        self._cls = cls
        self.name: int|str = cls._labels[index]
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
    def nodes(self) -> list:
        """
        """
        jnt = []
        for conn in self._cls._connectivity[self.index]:
            nname =self._cls.f2u_points.get_name(conn)
            jnt.append(self._cls.f2u_points[nname])
        return jnt
    #
    @property
    def material(self) -> list:
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
    def section(self) -> list:
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
    #@property
    #def connectivity(self):
    #    """ """
    #    connodes = self._cls._connectivity[self.index]
    #    connodes
    #
    #
    def __str__(self) -> str:
        """ """
        beam_name = self.name
        node1, node2 = self._cls._connectivity[self.index]
        title = ""
        return "{:>8s} {:>8s} {:>8s} {:>12s} {:>12s} {: 1.2e} {:>1.3e} {:}\n"\
               .format(str(beam_name), str(node1), str(node2),
                       str(self.material.name), str(self.section.name),
                       self.beta, self.L, title)    
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
        return self._cls._segment[index[0]] #* self._cls.f2u_units.m

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
        materials = self._cls.f2u_materials
        index  = self._get_index()
        name = self._cls._materials[index[0]]
        return materials[name]

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
        sections = self._cls.f2u_sections
        index = self._get_index()
        return sections[self._cls._sections[index[0]]]

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
    def _mesh(self, element: str|int):
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
        # [x, y, z, boundary]
        node1 = self._cls.f2u_points.get_coordinates(node1)
        node2 = self._cls.f2u_points.get_coordinates(coord)
        length = math.dist(node1[:3], node2[:3])
        self._cls._segment[index] = length
        new_index = self._cls._duplicate_element(index)
        self._cls._step_label[new_index] = step_name
    #    
    # 
    def __getitem__(self, step_name: int|str) -> tuple:
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
#
class ConceptBeam(BeamBasic):
    """
    element[name] = [name, connectivity, material, section, type, group]
    connectivity[number] = [name, node1, node2,..., nodei]
    """
    __slots__ = ['_labels', '_number','_type', '_connectivity', '_beam_type',
                 '_sections', '_materials', '_mesh', '_releases',  '_roll_angle', 
                 '_direction_cosines', '_eccentricities', '_offset_index', 
                 '_segment', '_step_label', 
                 'f2u_points', 'f2u_materials', 'f2u_sections']

    
    def __init__(self, beam_type:str,
                 points, materials, sections,
                 labels, element_type,
                 component: int) -> None: # properties
        """
        Manages f2u elements
        """
        super().__init__(component)
        #
        self.f2u_materials = materials
        self.f2u_sections = sections
        self.f2u_points = points
        #
        self._beam_type = beam_type
        self._labels = labels
        self._type = element_type
        #
        self._sections:list[str|int] = []
        self._materials:list[str|int] = []
        self._roll_angle:list[float] = []
        self._number:list[int] = []
        #
        self._connectivity: list = []
        self._segment:list[float] = []
        self._step_label:list[str|int] = []
        self._mesh:list[str|int] = []
        #
        #self._direction_cosines = array('i', [])
        #self._offset_index = array('i', [])        
        #self._eccentricities: List = []
        #self._releases: List = []        
    #   
    # -------------------------------------------
    #
    def _get_data(self, parameters:list|tuple|dict)-> list:
        """
        parameters: node1, node2, material, section
        output = [connectivity, material, section ]
        """
        #
        #output = [None] * 3
        #
        # set connectivity
        # end 1
        try:
            node_1 = self.f2u_points.get_point_name(parameters[0])
        except IOError:
            node_1 = self.f2u_points.get_new_point(parameters[0])
        # end 2
        try:
            node_2 = self.f2u_points.get_point_name(parameters[1])
        except IOError:
            node_2 = self.f2u_points.get_new_point(parameters[1])
        #
        #
        material = self.f2u_materials.default
        section = self.f2u_sections.default
        for item in parameters[2:]:
            if is_material(item):
                material = item.name
            elif is_section(item):
                section = item.name
            elif isinstance(item, str):
                try:
                    material = self.f2u_materials[item].name
                except KeyError:
                    try:
                        section = self.f2u_sections[item].name
                    except KeyError:
                        raise IOError(f' beam input data {item} not valid')
            else:
                raise IOError(f' beam input data {item} not valid')
        #
        # Checks
        if not material:
            raise IOError('material missing')
        
        if not section:
            raise IOError('section missing')
        #
        conn = [node_1, node_2]
        return [conn, material, section]
    #
    def __setitem__(self, name: int|str,
                    parameters: list|tuple|dict) -> None:
        """
        nmae : beam name
        parameters = connectivity, material, section
        """
        try:
            index = self._labels.index(name)
            raise Exception(f'beam {name} already exist')
        except ValueError:
            # default
            self._labels.append(name)
            index = self._labels.index(name)
            self._type.append(self._beam_type)
            #
            self._roll_angle.append(0.0)
            self._number.append(index)
            #
            values = self._get_data(parameters)
            self._connectivity.append(values[0])
            self._materials.append(values[1])
            self._sections.append(values[2])
            #
            self._segment.append(0)
            self._step_label.append(name)
            self._mesh.append(-1)
            #
            # TODO: to be defined
            #self._properties.append(-1)
            #self._offset_index.append(-1)
            #self._direction_cosines.append(-1)
            # should be in demand
            #self._eccentricities.append([])
            #self._releases.append([])
    
    def __getitem__(self, element_name: str|int):
        """
        """
        try:
            _index = self._labels.index(element_name)
            _indices = [index for index, name in enumerate(self._labels) 
                        if element_name == name]
            _index = _indices[0]
            if self._type[_index] in ["beam", "truss"]:
                return BeamItemConcept(self, _index)
            else:
                return ElementConcept(self, _index)
        except ValueError:
            raise IndexError(' ** element {:} does not exist'.format(element_name))
    #
    # -------------------------------------------
    #
    def __delitem__(self, element_name: str|int) -> None:
        """
        """
        try:
            _nodes_empty = []
            _index = self._labels.index(element_name)
            # remove element form node's set
            for _item in self._connectivity[_index]:
                _node = self.f2u_points[_item]
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
            for _node_id in _nodes_empty:
                del self.f2u_points[_node_id]
            #
            # FIXME: number should be updated according new index
            self._number.pop(i)
        except ValueError:
            raise KeyError('    *** warning element {:} does not exist'
                           .format(element_name))
    #
    # -------------------------------------------
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
    # -------------------------------------------
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
    # -------------------------------------------
    #
    @property
    def df(self):
        """ """
        nodes = list(zip(*self._connectivity))
        data = {'name': self._labels,
                'type': self._beam_type,
                'material': self._materials,
                'section': self._sections,
                'node_1': nodes[0],
                'node_2': nodes[1],
                'node_3' : None, 
                'node_4' : None,
                 'roll_angle': self._roll_angle,
                 'title': self._step_label,}
        db = DBframework()
        membdf = db.DataFrame(data=data)     
        return membdf

    @df.setter
    def df(self, df):
        """ """
        # clean df
        df = get_beam_df(df)
        columns = list(df.columns)
        coord = list(filter(lambda v: re.match('(x|y|z)(1|2)',
                                               v, re.IGNORECASE), columns))

        if coord:
            coord1 = ['x1', 'y1', 'z1']
            coord2 = ['x2', 'y2', 'z2']
            #
            for index, row in df.iterrows():
                member_name = row['name']
                coord1 = [row['x1'], row['y1'], row['z1']]
                coord2 = [row['x2'], row['y2'], row['z2']]
                parameters = [coord1, coord2, row['material'], row['section']]
                              #row['roll_angle'], row['title']]
                self.__setitem__(member_name, parameters)
        else: # node, point
            1/0
        #
        #print('-->')
    #
    # -------------------------------------------
    #
#
#
def alphaNumOrder(string):
    """ Returns all numbers on 5 digits to let sort the string with numeric order.
    Ex: alphaNumOrder("a6b12.125")  ==> "a00006b00012.00125"
    """
    return ''.join([format(int(x), '05d') if x.isdigit()
                   else x for x in re.split(r'(\d+)', string)])
#
#
@dataclass
class BeamItemConcept(ElementConcept):
    __slots__ = ['name', 'index', '_cls', '_steps']

    def __init__(self, cls, element_index) -> None:
        """
        """
        super().__init__(cls, element_index)
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
    def offsets(self, eccentricities: list) -> None:
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
    def releases(self) -> tuple:
        """
        """
        _list = []
        for _item in self._cls._releases[self.index]:
            _list.append(self._cls._f2u_releases[_item])
        return _list

    #
    @property
    def hinges(self) -> tuple:
        """
        """
        return self.releases

    #
    #@property
    #def type(self) -> str:
    #    """
    #    """
    #    return self._cls._type[self.index]
    #
    #@type.setter
    #def type(self, beam_type:str):
    #    """
    #    """
    #    self._cls._type[self.index] = beam_type
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
    def L(self) -> Units:
        """
        """
        node1, node2 = self.nodes
        return math.dist(node1[:3], node2[:3])
        
    #
    def find_coordinate(self, node_distance:float, node_end:int=0) -> tuple:
        """
        """
        node = self.nodes
        nodeNo3 = [0, 0, 0]
        #
        if math.isclose(node_distance, 0, rel_tol=0.01):
            nodeNo3[0] = node[node_end].x
            nodeNo3[1] = node[node_end].y
            nodeNo3[2] = node[node_end].z
        else:
            if node_end == 1:
                v1 = (node[0].x - node[1].x)
                v2 = (node[0].y - node[1].y)
                v3 = (node[0].z - node[1].z)
            else:
                v1 = (node[1].x - node[0].x)
                v2 = (node[1].y - node[0].y)
                v3 = (node[1].z - node[0].z)
            #
            norm = (v1**2 + v2**2 + v3**2)**0.50
            v1 /= norm
            v2 /= norm
            v3 /= norm
            #
            nodeNo3[0] = (node[node_end].x + v1 * node_distance)
            nodeNo3[1] = (node[node_end].y + v2 * node_distance)
            nodeNo3[2] = (node[node_end].z + v3 * node_distance)
        #
        nodeNo3 = NodePoint(name=0,
                            component=0, 
                            number=0, 
                            coord_system='cartesian', 
                            x=nodeNo3[0], 
                            y=nodeNo3[1], 
                            z=nodeNo3[2],
                            r=None, theta=None, phi=None,
                            title=None, 
                            index=0,
                            boundary=None)
        return nodeNo3.system()
    #
    #
    @property
    def unit_vector(self) -> list[float]:
        """
        """
        node1, node2 = self.nodes
        L = math.dist(node1[:3], node2[:3])
        uv = list(map(sub, node2[:3], node1[:3]))
        return [item / L for item in uv]
    #
    #
    def T3D(self):
        """
        Returns the transformation matrix for the member
        """
        #if m2D:
            # Element length and orientation
            #node1,  node2 = self.nodes
            #Tlg = Rmatrix2D(node1, node2)
            #return Tlg
            #raise NotImplementedError()
        #else:        
        #if self.type in ['beam', 'truss']:
        #return Rmatrix(*self.unit_vector, self.beta)
        nodei, nodej = self.nodes
        #L = self.L #.value
        #r3 = Rmatrix2(nodei, nodej, L=L)
        #return Rmatrix(*self.unit_vector, self.beta)
        return Rmatrix2(nodei, nodej)
#
#
