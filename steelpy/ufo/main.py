# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
#
# package imports
from steelpy.ufo.mesh.main import Mesh
from steelpy.ufo.concept.main import Concept
from steelpy.ufo.properties.main import Properties
#
# 
#
#
#
class UFOmodel:
    """
    UFO Model class
    
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
    __slots__ = ['_name', # 'type', 'data',
                 #'db_file', '_plot', 'mesh_type',
                 '_properties', # '_materials', '_sections', 
                 '_mesh', '_concept']

    def __init__(self, name:str|int) -> None:
        """
        mesh_type : sqlite/inmemory
        """
        print("-- module : ufo Version 6.00dev")
        print('{:}'.format(52*'-'))
        #
        #if not component:
        #    component = "f2u_model"
        self._name = name
        #self.type: str = 'substructure'
        #self.mesh_type:str = 'sqlite'
        #
        # set main components
        #
        #self._mesh:dict = {}
        #self._mesh = Mesh(db_name=self._name,
        #                  sql_file=sql_file)
        #
        self._properties = Properties()
        #
        self._concept = Concept(name=self._name,
                                properties=self._properties)
        #self._concept_flag = False
    #
    # -------------------
    #
    #def materials(self, values:None|list=None,
    #              df=None):
    #    """
    #    """
    #    if isinstance(values, list):
    #        if isinstance(values[0], list):
    #            for item in values:
    #                self._materials[item[0]] = item[1:]
    #        else:
    #            self._materials[values[0]] = values[1:]
    #    #
    #    try:
    #        df.columns
    #        self._materials.df = df
    #    except AttributeError:
    #        pass
    #    #
    #    return self._materials   
    #
    #def sections(self, values:None|list=None,
    #             df=None):
    #    """
    #    """
    #    #if values:
    #    if isinstance(values, list):
    #        if isinstance(values[0], list):
    #            for item in values:
    #                self._sections[item[0]] = item[1:]
    #        else:
    #            self._sections[values[0]] = values[1:]
    #    #
    #    # dataframe input
    #    try:
    #        df.columns   
    #        self._sections.df = df
    #    except AttributeError:
    #        pass            
    #    #
    #    return self._sections
    #
    def properties(self):
        """
        """
        return self._properties
    #
    # -------------------
    #
    #@property
    def concept(self):
        """
        """
        #if not name:
        #    name = self.component
        #
        #self._concept = Concepts(name=name, 
        #                         #materials=self._materials,
        #                         #sections=self._sections,
        #                         properties= self._properties)
        return self._concept
    #
    # -------------------
    #
    def mesh(self, #name:str|None = None,
             sql_file:str|None = None):
        """ """
        #if not name:
        name = self._name
        if sql_file:
            name = None
        
        return Mesh(name=name,
                    sql_file=sql_file)

    #
    # -------------------
    #
    def build(self, name:str|None = None) -> None:
        """
        """
        #
        if not name:
            name = self._name        
        #
        #self._sections.get_properties()
        #
        if self._concept:
            self._mesh[name] = self._concept.mesh()
        #    meshing = Meshing(concept=self._concept,
        #                      component=name, 
        #                      mesh_type=self.mesh_type)
        #    self._mesh[name] = meshing.get_mesh()
        #    self._mesh[name].renumbering()
        #    #mesh._load._basic.FER()
        #    #return mesh
        #    #_sql.write_concept(self._concept)
        #
        # check wave load case
        #
        for key, item in self._mesh.items():
            item.build()
            #item._load._basic.wave_process()
            # TODO : remove second _load for simplification
            #item._load._basic.FER(elements= item._elements)        
        #
        #
        print('end meshing')
        return self._mesh[name]
    #
    #@property
    #def plot(self):
    #    """ """
    #    return self._plot
#
#
#
#