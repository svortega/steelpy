# 
# Copyright (c) 2009-2021 fem2ufo
#

# Python stdlib imports
from dataclasses import dataclass
import os
from typing import Tuple, Dict, List, ClassVar, Union

# package imports
from steelpy.f2uModel.mesh.main import Mesh
from steelpy.f2uModel.concept.main import Concepts
from steelpy.f2uModel.sections.main import Sections
from steelpy.f2uModel.material.main import Materials
from steelpy.f2uModel.properties.main import Properties
from steelpy.f2uModel.load.main import Load
from steelpy.f2uModel.results.main import Results
from steelpy.f2uModel.concept.meshing import Meshing

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


#
class f2uModel:
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
    __slots__ = ['component', 'type', 'data', '_results',
                 'db_file',
                 '_materials', '_sections', '_properties',
                 'sets', 'mesh', '_concept', '_concept_flag',
                 '_boundaries', '_nodes', '_load', '_meshing']

    def __init__(self, component: Union[str, int],
                 mesh_type:str='sqlite') -> None:
        """
        mesh_type : sqlite/inmemory
        """
        print("-- module : fem2ufo Version 5.00dev")
        print('{:}'.format(52*'-'))
        #
        self.component: Union[str, int] = component
        self.type: str = 'substructure'
        self.db_file:str = component + "_f2u.db"
        if mesh_type != "inmemory":
            try: # remove file if exist
                os.remove(self.db_file)
            except FileNotFoundError:
                pass
        #
        # set main components
        self._materials: ClassVar = Materials(mesh_type=mesh_type,
                                              db_file=self.db_file)
        self._sections: ClassVar = Sections(mesh_type=mesh_type,
                                            db_file=self.db_file)
        self._properties: ClassVar = Properties()
        # set mesh components
        self.mesh: ClassVar = Mesh(materials=self._materials,
                                   sections=self._sections,
                                   mesh_type=mesh_type,
                                   db_file=self.db_file)
        # loads component
        self._load = Load(mesh_type=mesh_type,
                          db_file=self.db_file)
        self._results = Results(basic_loads=self._load,
                                mesh_type=mesh_type,
                                db_file=self.db_file)
        #
        #self.sets: ClassVar = Groups(concepts=self._concept,
        #                             mesh= self.mesh)
        #self.data: Tuple = InputData()
        #
        # set concepts
        self._concept: ClassVar = Concepts(mesh=self.mesh, 
                                           materials=self.materials,
                                           sections=self.sections,
                                           properties= self._properties)
        self._concept_flag = False
        #
        self._meshing = Meshing(mesh=self.mesh, 
                                concept=self._concept,
                                load=self._load)
        #
        # start defaults
        #self._material_default: bool = None
        #self._section_default: bool = None
        #self._set_dafault: bool = None

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
    #
    @property
    def concept(self) -> ClassVar:
        """
        """
        self._concept_flag = True
        return self._concept

    #
    @property
    def groups(self):
        """ """
        return self.sets
    #
    #
    #
    @property
    def load(self) -> ClassVar:
        """
        """
        return self._load
    #
    #
    def build(self) -> None:
        """
        """
        #_sql = f2uDB(component_name= self.component)
        #_sql.materials = self._materials
        self.sections.get_properties()
        if self._concept_flag:
            self._meshing.get_mesh()
            self.mesh.renumbering()
            #_sql.write_concept(self._concept)
        #_sql.write_geometry(self.mesh)
        #_sql.write_load(self._load)
        #self._sql = _sql
        print('end meshing')
    #

