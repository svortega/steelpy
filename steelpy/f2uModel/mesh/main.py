# 
# Copyright (c) 2009-2021 fem2ufo
#

# Python stdlib imports
from typing import List, ClassVar, Dict, Union
#

# package imports
#from steelpy.f2uModel.mesh.inmemory.element import Elements
#from steelpy.f2uModel.mesh.geometry import Releases
from steelpy.f2uModel.mesh.inmemory.sets import Groups
#from steelpy.f2uModel.mesh.boundary import Boundaries
#
from steelpy.f2uModel.mesh.nodes import Nodes
from steelpy.f2uModel.mesh.boundaries import Boundaries
from steelpy.f2uModel.mesh.elements import Elements
#
#from steelpy.f2uModel.mesh.operations.nodes import node_renumbering
#from steelpy.f2uModel.sql.sqlmodule import f2uDB
#

#
class Mesh:
    """
    mesh[beam_name] = [number, element1, element2, elementn]
    """
    __slots__ = ['_nodes', '_elements', # '_materials', '_sections',
                 '_eccentricities', '_boundaries', '_groups', 'db_file']

    def __init__(self, materials, sections,
                 mesh_type:str="inmemory",
                 db_file:Union[str,None]=None):
        """
        """
        self.db_file = db_file
        self._nodes = Nodes(mesh_type=mesh_type,
                            db_file = self.db_file)
        self._boundaries = Boundaries(mesh_type=mesh_type,
                                      db_file = self.db_file)
        #
        #self._elements: ClassVar = Elements()
        self._elements: ClassVar = Elements(nodes= self._nodes,
                                            materials=materials,
                                            sections=sections,
                                            mesh_type=mesh_type,
                                            db_file = self.db_file)
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
        self._nodes.renumbering(self._elements)
        #for node in single_nodes:
        #    boundary = self._boundaries.node[node]
        #    if not boundary:
        #        self._boundaries.node[ node ] = 'free'
        print("** End Renumbering Nodes")
    #
    #
    #
  
