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
                 #'_load_case', '_load_combination',  '_sql']

    def __init__(self):
        """
        """
        #self._units = Units()
        #self._materials: ClassVar = Materials()
        #self._sections: ClassVar = Sections()
        #self._properties: ClassVar = properties
        self._boundaries: ClassVar = Boundaries()
        self._nodes = Nodes()
        self._elements: ClassVar = Elements()
        # groups
        self._groups = Groups()
        # loads
        #self._load_case = LoadCase()
        #self._load_combination = LoadCombination()
        #self._sql = f2uDB(model=self)
    
    #@property
    #def units(self):
    #    """
    #    """
    #    return self._units
    
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
        get_nodes_connected(self._elements)
        self._nodes._renumber()
        print("** End Renumbering Nodes")
    #
    #def dump_sql_model(self):
    #    """
    #    """
    #    print("** Writing SQL to Disk")
    #    populate_sql(component_name="f2u_model", 
    #                 component=self)
    #    populate_loading_sql(component_name="f2u_model", 
    #                         component=self)
    #    print("** End Writing SQL to Disk")
    ##
    #def load_sql_model(self):
    #    """
    #    """
    #    component_name = "f2u_model.db"
    #    load_sql_model(component_name)
    #    print("-->")  
    ##
    #def dump_mesh(self):
    #    """
    #    """
    #    print("** Writing Mesh to Disk")
    #    dump_f2u_mesh(self)
    #    print("** End Writing Mesh to Disk")
    #
    #
  
