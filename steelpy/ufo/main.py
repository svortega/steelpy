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
from steelpy.ufo.mesh.sqlite.main import UFOmain
# 
#
#
#
class UFOmodel(UFOmain):
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
    __slots__ = ['_name', '_properties', 
                 '_mesh', '_concept']

    def __init__(self, name:str|int|None = None,
                 sql_file:str|None = None) -> None:
        """
        mesh_type : sqlite/inmemory
        """
        print("-- module : ufo Version 6.50dev")
        print('{:}'.format(52*'-'))
        #
        super().__init__(name=name,
                         sql_file=sql_file)        
        #
        # mesh basic
        self._mesh = Mesh(component=self._component, 
                          sql_file=self.db_file)
        #self._name = name
        self._properties = Properties()
        #
        self._concept = Concept(component=self._component,
                                mesh=self._mesh, 
                                properties=self._properties)

    #
    #
    # -------------------
    #
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
        return self._concept
    #
    # -------------------
    #
    def mesh(self): #name:str|None = None,
             #sql_file:str|None = None):
        """ """
        return self._mesh

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