#
# Copyright (c) 2009-2022 fem2ufo
#

# Python stdlib imports
from dataclasses import dataclass
import math
from typing import NamedTuple, Dict, Union, Tuple, List

# package imports
from steelpy.f2uModel.mesh.operations.beam.stiffness import (beam_Klocal, trans_3d_beam, 
                                                             Rmatrix, Rmatrix_new,
                                                             trans3Dbeam, beam_stiffness)
from steelpy.f2uModel.mesh.operations.beam.beam_sol import BeamBasic
#from steelpy.process.math.vector import Vector

#
#
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
    def Kglobal(self):
        """
        K stiffness matrix
        """
        #FIXME: material 
        section = self._elements._f2u_sections[self.section].properties
        material = self._elements._f2u_materials[self.material]
        # solve K matrix
        #K = beam_Klocal(self.length,
        #            section.area, section.J,
        #            section.Iy, section.Iz,
        #            material.E.convert("pascal").value, 
        #            material.G.convert("pascal").value,
        #            section.area, section.area)
        #
        K = beam_stiffness(self.length,
                           section.area, section.J,
                           section.Iy, section.Iz,
                           material.E.convert("pascal").value, 
                           material.G.convert("pascal").value)
        #
        #self.beta = 30
        #disb = self.R
        node1 = self._elements._f2u_nodes[self._elements._connectivity[self.index][0]]
        node2 = self._elements._f2u_nodes[self._elements._connectivity[self.index][-1]]
        #dirc = Rmatrix_new(node1[:3], node2[:3], self.beta)
        #xxx = trans3Dbeam(K, node1[:3], node2[:3])
        #return trans3Dbeam(K, node1[:3], node2[:3])
        #return trans_3d_beam(K, dirc)
        return trans_3d_beam(K, self.R)
        #return self.trans3d(K)
    #
    #
    @property
    def R(self):
        """
        Rotation matrix
        """
        #if self.type in ['beam', 'truss']:
        return Rmatrix(*self.unit_vector, self.beta)
        #else:
        #    raise IOError("no yet included")
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
@dataclass
class BeamElementSQL:
    """ """
    __slots__ = ['name', 'db_file', 'type']
    
    def __init__(self, element_name:int, db_file:str) -> None:
        """
        """    
        self.name = element_name
        self.db_file = db_file
        self.type: str = "beam"
    #
    @property
    def number(self) -> int:
        """ """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name)
        return data[1]
    
    @number.setter
    def number(self, value:int) -> None:
        """"""
        1/0
        conn = create_connection(self.db_file)
        item = "number"
        with conn:
            update_element_item(conn, self.name, item, value)
    #
    @property
    def connectivity(self) -> List:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            connodes = get_connectivity(conn, self.name)
        return connodes

    @connectivity.setter
    def connectivity(self, nodes:List[int]) -> List:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            #push_connectivity(conn, self.name, nodes)
            update_connectivity(conn, self.name, nodes)
        #self._connectivity[self.index] = nodes
    #
    @property
    def material(self) -> List:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name)
        return data[4]

    @material.setter
    def material(self, material_name: str) -> None:
        """
        """
        conn = create_connection(self.db_file)
        item = "material"
        with conn:
            update_element_item(conn, self.name, item, material_name)

    #
    @property
    def section(self) -> List:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name)
        return data[5]

    @section.setter
    def section(self, section_name: str) -> None:
        """
        """
        conn = create_connection(self.db_file)
        item = "section"
        with conn:
            update_element_item(conn, self.name, item, self.name)

    #
    @property
    def beta(self):
        """beta angle roll"""
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name)
        return data[3]
    
    @beta.setter
    def beta(self, value):
        """beta angle roll"""
        conn = create_connection(self.db_file)
        item = "roll_angle"
        with conn:
            update_element_item(conn, self.name, item, self.name)
    #
    #
    def __str__(self) -> str:
        """ """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name)
        #title =  data[-1]
        if (title := data[-1]) == "NULL":
            title = ""        
        #
        return "{:8d} {:8d} {:8d} {:>12s} {:>12s} {: 6.4f} {:>6.3f} {:>12s}\n"\
               .format(self.name, *self.connectivity,
                       self.material, self.section, self.beta,
                       self.length, title)
    #
    #
    @property
    def DoF(self) -> List[ int ]:
        """
        """
        conn = create_connection(self.db_file)
        dof = [ ]
        for node_name in self.connectivity:
            node = get_node(conn, node_name=node_name)
            number = node[0] - 1
            dof.append(number) #  * 6
        return dof
    #
    @property
    def length(self) -> float:
        """
        """
        conn = create_connection(self.db_file)
        nodes = self.connectivity
        node1 = get_node(conn, node_name=nodes[0])
        node2 = get_node(conn, node_name=nodes[1])
        # _dx = _node1.x - _node2.x
        # _dy = _node1.y - _node2.y
        # _dz = _node1.z - _node2.z
        # dist2 = (_dx * _dx + _dy * _dy + _dz * _dz)**0.50
        return dist(node1[3:6], node2[3:6])
    #
    @property
    def unit_vector(self) -> List[ float ]:
        """
        """
        # TODO: get_node should be aligned with inmemmory
        conn = create_connection(self.db_file)
        node1 = get_node(conn, node_name=self.connectivity[0])
        node2 = get_node(conn, node_name=self.connectivity[1])
        dx = node2[3] - node1[3]
        dy = node2[4] - node1[4]
        dz = node2[5] - node1[5]
        # direction cosines
        L = dist(node1[3:6], node2[3:6])
        l = dx / L
        m = dy / L
        n = dz / L
        return [l, m, n]
    #
    @property
    def Kglobal(self):
        """ """
        #conn = create_connection(self.db_file)
        material, section, beta = self._K_data()
        #section = get_sectionSQL(conn, self.section)
        #material = get_materialSQL(conn, self.material)
        # solve K matrix
        R = Rmatrix(*self.unit_vector, beta)
        # R = Rmatrix(*self.direction_cosines, self.beta)
        # K = beam_stiffness(self.length,
        #                   section.area,
        #                   section.J,
        #                   section.Iy,
        #                   section.Iz,
        #                   material.E,
        #                   material.G)
        K = beam_Klocal(self.length,
                    section.area, section.J,
                    section.Iy, section.Iz,
                    material.E, material.G,
                    section.area, section.area)
        return trans_3d_beam(K, R)
    #
    def _K_data(self):
        """ """
        conn = create_connection(self.db_file)
        cur = conn.cursor()
        cur.execute ("SELECT * FROM tb_Elements\
                    WHERE tb_Elements.name = {:};".format(self.name))
        row = cur.fetchone()
        #
        #connodes = get_connectivity(conn, self.name)
        #data = [*row[4:], connodes]
        material = get_materialSQL(conn, row[4]) 
        section = get_sectionSQL(conn, row[5])
        beta = row[6]
        conn.close()
        return material, section, beta
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
#
#