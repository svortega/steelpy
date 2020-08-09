# 
# Copyright (c) 2009-2020 fem2ufo
# 


# Python stdlib imports
import sys
from typing import NamedTuple, Tuple, List, Dict, Iterable, TypeVar, ClassVar

# package imports
#from steelpy.f2uModel.femodel.mesh.mesh import Mesh
from steelpy.f2uModel.concept.concept import Concepts
#from steelpy.f2uModel.material.material import Materials
#from steelpy.f2uModel.sections.section import Sections
#from steelpy.f2uModel.sets.sets import Groups
#from steelpy.f2uModel.properties.properties import Properties
#from steelpy.f2uModel.model.assemble import InputData
#from steelpy.f2uModel.femodel.boundary.boundary import Boundaries
#
#from fem2ufo.f2u_model.femodel.mesh.node import Nodes


#
def get_number(items: dict) -> int:
    """
    return maximum consecutive number of items
    """
    _number = [item.number for item in items.values()]
    try:
        return max(_number)
    except ValueError:
        return 0


# GEOMETRY Section
class Component:
    """
    FE Geometry model class
    
    Parameters:
      :number: integer
      :name: string
      :nodes: list of node's class
      :elements: list of element's class
      :materials: list of material's class
      :sets: list of groups (elements, nodes, etc)
      :sections: list of section's class
      :vectors: list of guide points
      :eccentricities: list of eccentricities
      :joints: list of joint's class
      :hinges: list of hinges definitios
      :loads: list of load's class
      :data: FE model data
      :units: FE model units
      :soil: list of soil's class
      :hydrodynamics: hydrodynamic's data
      :boundaries: list of FE model boundary
    
    Parameters:  
      :number:  integer internal number 
      :name:  string node external name
    """
    __slots__ = ('number', 'name', 'units', 'type', 'data',
                 '_materials', '_sections', '_properties',
                 'sets', 'load', 'mesh', '_concept',
                 '_material_default', '_section_default', 
                 '_set_dafault', '_boundaries', '_nodes')

    def __init__(self, name: str, number: int = None) -> None:
        """
        """
        self.number: int = number
        self.name: str = name
        self.type: str = 'substructure'
        # set main components
        #self._materials: ClassVar = Materials()
        #self._sections: ClassVar = Sections()
        #self._boundaries: ClassVar = Boundaries()
        #self._properties: ClassVar = Properties()
        # set mesh components
        #self._nodes = Nodes()
        # set mesh
        self.mesh: ClassVar = Mesh(materials= self._materials,
                                   sections= self._sections,
                                   #properties= self._properties,
                                   boundaries=self._boundaries)
        #
        # set concepts
        self._concept: ClassVar = Concepts(mesh= self.mesh,
                                           properties= self._properties)
                                           #materials= self._materials, 
                                           #sections= self._sections,
                                           #
                                           #boundaries= self._boundaries)
        
        #self.sets: ClassVar = Groups(concepts=self._concept,
        #                             mesh= self.mesh)
        #self.data: Tuple = InputData()
        # start defaults
        self._material_default: bool = None
        self._section_default: bool = None
        self._set_dafault: bool = None

    #
    @property
    def materials(self) -> ClassVar:
        """
        """
        return self._materials

    #
    @property
    def sections(self) -> ClassVar:
        """
        """
        return self._sections

    #
    @property
    def properties(self) -> ClassVar:
        """
        """
        return self._properties

    @property
    def concept(self) -> ClassVar:
        """
        """
        return self._concept

    #
    @property
    def groups(self):
        """ """
        return self.sets
    #
    @property
    def default_material(self, material_name: str) -> None:
        """
        """
        _material_default = self.material[material_name]

    #
    @property
    def default_section(self, section_name: str) -> None:
        """
        """
        _section_dafault = self.sets[section_name]

    #
    @property
    def default_group(self, group_name: str) -> None:
        """
        """
        _set_dafault = self.sets[group_name]

    #
    def add_joint(self, joint_name: str, node_number: int) -> None:
        """
        """
        self.joints.map_joint_with_node(joint_name, node_number)
        #
        # self.joints[joint_name] = geometry.Joint(joint_name)
        # self.joints[joint_name].node = node_number

    #
    def add_support(self, name: str, node_name: int = False,
                    element_name: str = False, boundary: str = False) -> None:
        """
        """
        self.boundaries[name] = geometry.Support(name)

        if node_name:
            try:
                # print('-->')
                self._mesh.nodes[node_name]
            except IndexError:
                print('   *** error : Node {:} not found'.print(node_name))
                sys.exit
            # self._mesh.nodes[node_name].set_boundary(self.boundaries[name])
            self.boundaries[name].nodes.append(node_name)
        elif element_name:
            print('fix boundary element here')
            sys.exit()
        #
        if boundary:
            print(' add here boundary')

    #
    def add_material(self, name: str, number: int = None) -> None:
        """
        """
        if not number:
            number = get_number(self.material) + 1

        self.material[name] = Materials(name, number)
    #
