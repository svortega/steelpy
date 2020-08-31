# 
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#from collections import Counter
from typing import List, ClassVar, Dict
#from itertools import chain



# package imports
from steelpy.f2uModel.mesh.element import Elements
#from steelpy.f2uModel.mesh.geometry import Releases
from steelpy.f2uModel.mesh.sets import Groups
from steelpy.f2uModel.mesh.boundary import Boundaries
from steelpy.f2uModel.mesh.node import Nodes
#from steelpy.f2uModel.mesh.operations import dump_f2u_mesh
#
#
from steelpy.f2uModel.mesh.operations import get_nodes_connected
#from steelpy.f2uModel.sql.sqlmodule import f2uDB
#

#
class Mesh:
    """
    mesh[beam_name] = [number, element1, element2, elementn]
    """
    __slots__ = ['_nodes', '_elements', '_materials', '_sections',
                 '_eccentricities', '_boundaries', '_groups']

    def __init__(self):
        """
        """
        self._boundaries: ClassVar = Boundaries()
        self._nodes = Nodes()
        self._elements: ClassVar = Elements()
        # groups
        self._groups = Groups()
    #
    @property
    def nodes(self) -> ClassVar:
        """
        """
        return self._nodes
    
    @nodes.setter
    def nodes(self, values:List) -> ClassVar:
        """
        """
        for value in values:
            self._nodes[value[0]] = value[1:4]
            try:
                if isinstance(value[4], str):
                    self._boundaries.node[value[0]] = value[4]
                else:
                    self._boundaries.node[value[0]] = value[4:]
            except IndexError:
                pass

    @property
    def elements(self) -> ClassVar:
        """
        """
        return self._elements
    
    @elements.setter
    def elements(self, values) -> ClassVar:
        """
        """
        for value in values:
            self._elements[value[0]] = value[1:]

    @property
    def boundaries(self) -> ClassVar:
        """
        """
        return self._boundaries
    #
    @property
    def sections(self) -> ClassVar:
        """
        """
        return self._sections
    
    @sections.setter
    def sections(self, value) -> None:
        """
        """
        self._sections = value
    #
    @property
    def materials(self) -> ClassVar:
        """
        """
        return self._materials
    
    @materials.setter
    def materials(self, value) -> None:
        """
        """
        self._materials = value
    #
    @property
    def groups(self) -> ClassVar:
        """
        """
        return self._groups
    #
    #@property
    def renumbering(self):
        """
        """
        print("** Renumbering Nodes")
        get_nodes_connected(self._elements, self._nodes)
        self._nodes._renumber()
        print("** End Renumbering Nodes")
    #
    #
    #
  
