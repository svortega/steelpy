# 
# Copyright (c) 2009-2023 fem2ufo
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
import os
#
# package imports
# steelpy.f2uModel
from .mesh.main import Mesh
from .concept.main import Concepts
from .properties.main import Properties
from .process.meshing import Meshing
#from .plot.main import PlotModel
# 
from steelpy.sections.main import Sections
from steelpy.material.main import Materials
#
#from steelpy.material.process.operations import get_isomat_prop_df
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
    __slots__ = ['component', 'type', 'data',
                 'db_file', '_plot', 'mesh_type',
                 '_materials', '_sections', '_properties',
                 '_mesh', '_concept', '_concept_flag']
                 # '_nodes', '_meshing', '_boundaries',  'sets', '_load', '_results',

    def __init__(self, component:str|int) -> None:
        """
        mesh_type : sqlite/inmemory
        """
        print("-- module : fem2ufo Version 5.00dev")
        print('{:}'.format(52*'-'))
        #
        self.component: str|int = component
        self.type: str = 'substructure'
        self.mesh_type:str = 'sqlite'
        #
        #BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        filename = component + "_f2u.db"
        path = os.path.abspath(filename)
        self.db_file = path
        #directory = os.path.dirname(path)
        #
        #self.db_file:str = component + "_f2u.db"
        #if mesh_type != "ttt": #"inmemory":
        try: # remove file if exist
            os.remove(self.db_file)
        except FileNotFoundError:
            pass
        #
        # set main components
        #mesh_type2 = 'sqlite'
        self._materials = Materials(mesh_type=self.mesh_type,
                                    db_file=self.db_file)

        self._sections = Sections(mesh_type=self.mesh_type,
                                  db_file=self.db_file)

        self._mesh = Mesh(materials=self._materials,
                          sections=self._sections,
                          mesh_type=self.mesh_type,
                          db_file=self.db_file)
        #
        self._properties = Properties()
        #
        # set concepts
        #self._concept = Concepts(materials=self._materials,
        #                         sections=self._sections,
        #                         properties= self._properties)
        self._concept_flag = False
        #
        #self._meshing = Meshing(concept=self._concept,
        #                        mesh_type=mesh_type,
        #                        db_file=self.db_file)
        #
        # start defaults
        #self._material_default: bool = None
        #self._section_default: bool = None
        #self._set_dafault: bool = None
        #self._plot = PlotModel(mesh=self._mesh)
        #self._plot._concept = self._concept
    #
    def materials(self, values:None|list=None,
                  df=None):
        """
        """
        if isinstance(values, list):
            if isinstance(values[0], list):
                for item in values:
                    self._materials[item[0]] = item[1:]
            else:
                self._materials[values[0]] = values[1:]
        #
        try:
            df.columns
            self._materials.df = df
            #group = df.groupby("type")
            # Elastic type
            #try:
            #    elastic = group.get_group("elastic")
            #    #elastic = get_isomat_prop_df(elastic)
            #    #elastic = elastic.drop_duplicates(['name'])
            #    self._materials.elastic(df=df)
            #except KeyError:
            #    # nonlin = group.get_group("plastic")
            #    raise IOError('Material type not valid')
        except AttributeError:
            pass
        #
        return self._materials   

    #
    def sections(self, values:None|list=None,
                 df=None):
        """
        """
        #if values:
        if isinstance(values, list):
            if isinstance(values[0], list):
                for item in values:
                    self._sections[item[0]] = item[1:]
            else:
                self._sections[values[0]] = values[1:]
        #
        # dataframe input
        try:
            df.columns   
            self._sections.df = df
        except AttributeError:
            pass            
        #
        return self._sections

    #
    def properties(self):
        """
        """
        return self._properties
    #
    def concept(self):
        """
        """
        self._concept = Concepts(materials=self._materials,
                                 sections=self._sections,
                                 properties= self._properties)
        #self._plot._concept = self._concept
        self._concept_flag = True
        return self._concept

    #
    def groups(self):
        """ """
        return self.sets
    #
    def mesh(self):
        """ """
        #self._mesh = Mesh(materials=self._materials,
        #                  sections=self._sections,
        #                  mesh_type=self.mesh_type,
        #                  db_file=self.db_file)
        return self._mesh
    #
    #
    #
    def build(self) -> None:
        """
        """
        #
        #
        #
        self._sections.get_properties()
        #
        if self._concept_flag:
            meshing = Meshing(concept=self._concept,
                              mesh=self._mesh)
            meshing.get_mesh()
            self._mesh.renumbering()
            #mesh._load._basic.FER()
            #return mesh
            #_sql.write_concept(self._concept)
        #
        # check wave load case
        #
        self._mesh._load._basic.wave_process()
        #
        # TODO : remove second _load for simplification
        self._mesh._load._basic.FER(elements= self._mesh._elements)        
        #
        #
        print('end meshing')
        return self._mesh
    #
    #@property
    #def plot(self):
    #    """ """
    #    return self._plot
#
#
#
#