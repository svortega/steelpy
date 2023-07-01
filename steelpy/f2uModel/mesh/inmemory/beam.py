#
# Copyright (c) 2009-2023 fem2ufo
#

# Python stdlib imports
from __future__ import annotations
from array import array
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from itertools import chain
import math
#from typing import NamedTuple
from operator import sub, add

# package imports
# steelpy.f2uModel.mesh.

from steelpy.f2uModel.mesh.process.elements.bstiffness import (beam_Klocal, trans_3d_beam,
                                                               Rmatrix, Rmatrix_new,
                                                               trans3Dbeam, beam_stiffness)
#
#from ...process.beam.beam_sol import BeamBasic
#from steelpy.formulas.pilkey.main import BasicBeam
#from steelpy.beam.main import BasicCalcs
#
#from steelpy.process.math.vector import Vector
#
from ..process.elements.beam import BeamBasic, BeamItemBasic
#
#
#
#
class BeamIM(BeamBasic):
    __slots__ = ['_labels', '_number',  '_title', '_connectivity',
                 '_sections', '_materials', '_releases',  '_roll_angle',
                 '_direction_cosines', '_eccentricities', '_offset_index',
                 '_f2u_nodes', '_f2u_sections', '_f2u_materials']

    
    def __init__(self, nodes, materials, sections) -> None:
        """
        Beam element 
        """
        super().__init__()
        #
        self._f2u_nodes = nodes
        self._f2u_sections = sections
        self._f2u_materials = materials
        #
        self._sections:list[str|int] = []
        self._materials:list[str|int] = []
        #
        self._roll_angle: array = array('f', [])
        #self._direction_cosines = array('i', [])
        #self._offset_index = array('i', [])
        #self._eccentricities: List = []
        #
        self._connectivity: list = []
        self._releases: list = []
    #
    def __setitem__(self, element_name: int, parameters: list) -> None:
        """
        parameters = [node1, node2, material, section, roll_angle]
        """
        try:
            self._labels.index(element_name)
            raise Exception('element {:} already exist'.format(element_name))
        except ValueError:
            # check 
            mat = self._f2u_materials[parameters[2]]
            sect = self._f2u_sections[parameters[3]]
            #
            self._labels.append(element_name)
            mnumber = next(self.get_number())
            self._number.append(mnumber)
            #
            self._connectivity.append(parameters[:2])
            self._materials.append(parameters[2])
            self._sections.append(parameters[3])
            self._roll_angle.append(parameters[4])
            try:
                self._title.append(parameters[5])
            except IndexError:
                self._title.append("NULL")
            #
            self._releases.append([])
    #
    def __getitem__(self, element_name: int):
        """
        """
        try:
            return BeamItemIM(self, element_name)
        except ValueError:
            raise IndexError(' ** element {:} does not exist'.format(element_name))    
    #
    def __delitem__(self, element_name: str|int) -> None:
        """
        """
        try:
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
            #
            # FIXME: number should be updated according new index
            self._number.pop( i )
        except ValueError:
            raise KeyError('    *** warning element {:} does not exist'
                           .format(element_name))    
    #   
    #
    # f2u process
    #
    def get_number(self, start:int=1):
        """
        """
        try:
            n = max(self._number) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
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
    @property
    def get_connectivities(self):
        """ """
        return self._connectivity
    #
    # Beam process
    #
#
#
@dataclass
class BeamItemIM(BeamItemBasic):
    """
    """
    __slots__ = ['name', 'index', '_cls']

    def __init__(self, cls, element_name: int|str) -> None:
        """
        """
        super().__init__(element_name=element_name)
        self._cls = cls
        self.index: int = cls._labels.index(element_name)
        #self.name: int|str = element_name
    #
    @property
    def number(self) -> int:
        return self._cls._number[self.index]

    @number.setter
    def number(self, value:int) -> None:
        """"""
        self._cls._number[ self.index ] = value
    #
    @property
    def type(self)-> str:
        """
        """
        return self._cls._type[self.index]

    #
    @property
    def connectivity(self) -> list:
        """
        """
        return self._cls._connectivity[self.index]

    @connectivity.setter
    def connectivity(self, nodes:list[int]) -> list:
        """
        """
        self._cls._connectivity[self.index] = nodes
    #
    @property
    def nodes(self) -> list:
        """
        """
        _nodes = []
        for _node in self._cls._connectivity[self.index]:
            _nodes.append(self._cls._f2u_nodes[_node])
        return _nodes

    @nodes.setter
    def nodes(self, connectivity: list) -> None:
        """
        """
        for _node in connectivity:
            self._cls._f2u_nodes[_node].sets.elements.append(self.number)
        self._cls._connectivity[self.index] = list(connectivity)

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
    @property
    def global_offset(self) -> list:
        """
        return eccentricities in the global system
        """
        1 / 0
        return self._cls.eccentricities[self.eccentricities]

    @global_offset.setter
    def global_offset(self, eccentricities: list) -> None:
        """
        set offset in the global system
        input
        eccentricities : list [x, y, z]
        """
        1 / 0
        self._cls._eccentricities[self.index] = eccentricities
    #
    #
    @property
    def material(self) -> list:
        """
        """
        material_name = self._cls._materials[self.index]
        #return self._cls._f2u_materials.get_item_by_number(_material_number)
        #return material_name
        return self._cls._f2u_materials[material_name]

    @material.setter
    def material(self, material_name: str) -> None:
        """
        """
        #self._cls._f2u_materials[material_name].sets.elements.append(self.number)
        #material_number = self._cls._f2u_materials[material_name].number
        #self._cls._materials[self.index] = material_number
        self._cls._materials[self.index] = material_name

    #
    #
    @property
    def section(self) -> list:
        """
        """
        section_name = self._cls._sections[self.index]
        #return self._cls._f2u_sections.get_item_by_number(section_number)
        #
        #return section_name
        return self._cls._f2u_sections[section_name]

    @section.setter
    def section(self, section_name: str) -> None:
        """
        """
        # self._mesh.sections[section_name].sets.elements.append(self.number)
        #section_number = self._cls._f2u_sections[section_name].number
        #self._cls._sections[self.index] = section_number
        self._cls._sections[self.index] = section_name

    #
    #
    #
    #
    @property
    def L(self) -> float:
        """
        """
        node1, node2 = self.nodes
        return math.dist(node1[:3], node2[:3])
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
    def __str__(self) -> str:
        """ """
        if (title := self._cls._title[self.index]) == "NULL":
            title = ""
        return "{:8d} {:8d} {:8d} {:>12s} {:>12s} {: 1.2e} {:>1.3e} {:>12s}"\
               .format(self.name, *self.connectivity,
                       str(self.material.name), str(self.section.name),
                       self.beta, self.L, title)
        #return (f"{self.name} {self.connectivity[0]} {self.connectivity[1]} \
        #        {self.material.name} {self.section.name} {self.beta}")
    #
    @property
    def unit_vector(self) -> list[ float ]:
        """
		Direction cosines for the local x-axis
        """
        node1 = self._cls._f2u_nodes[self._cls._connectivity[self.index][0]]
        node2 = self._cls._f2u_nodes[self._cls._connectivity[self.index][-1]]
        # direction cosines
        L = math.dist(node1[:3], node2[:3])
        #
        uv = list(map(sub, node2[:3], node1[:3]))
        return [item / L for item in uv]
        #return [l, m, n]
    #
    #
    @property
    def direction_cosines(self) -> list[float]:
        """
        """
        return self.unit_vector    
    #
    #
    # @properties.setter
    # def properties(self, property_number:int) -> None:
    #    """
    #    """
    #    section_number = self._sections[section_name].number
    #    self._properties[self.index] = property_number
    #
    #
    @property
    def DoF(self):
        """
        """
        dof = []
        for _node in self._cls._connectivity[self.index]:
            dof.append(self._cls._f2u_nodes[_node].index)
        return dof
    #
    #

    
#
#