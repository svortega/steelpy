# 
# Copyright (c) 2009 steelpy
# 
# Python stdlib imports
from __future__ import annotations


# package imports
from steelpy.ufo.mesh.main import ConceptMesh
from steelpy.ufo.process.main import ufoBasicModel, ModelClassBasic
#
from steelpy.ufo.concept.process.joint import Connection
from steelpy.ufo.concept.process.boundary import BoundaryConcept
from steelpy.ufo.concept.process.geometry import Releases
from steelpy.ufo.concept.process.meshing import MeshingConcept
#
from steelpy.ufo.concept.elements.sets import Groups
from steelpy.ufo.concept.elements.points import NodesIM
from steelpy.ufo.concept.elements.main import ConceptElements
#
from steelpy.ufo.load.main import ConceptLoad
#
from steelpy.ufo.plot.main import PlotConcept
#
#from steelpy.sections.main import Section
#from steelpy.material.main import Material
#
#
#
class Concept(ModelClassBasic):
    """ f2u Concept model Class """
    __slots__ = ['_labels', '_item', '_name',
                 '_properties', '_mesh']
    #
    def __init__(self, name:str, properties):
        """
        """
        super().__init__()
        self._name = name
        #self._mesh = mesh
        self._properties = properties
        self._item : dict = {}
        #
        self._mesh = ConceptMesh(name=self._name)
    #
    def __setitem__(self, name: int|str, parameters: list|str) -> None:
        """
        """
        try:
            self._labels.index(name)
            raise Exception('Concept {:} already exist'.format(name))
        except ValueError:
            self._labels.append(name)
            #self._item.append(name)
            #
            self._mesh[name] = parameters
            #
            self._item[name] = ConceptItem(component=name,
                                           title=parameters,
                                           mesh=self._mesh[name], 
                                           #materials=self._mesh[name]._materials,
                                           #sections=self._mesh[name]._sections,
                                           properties= self._properties)
    
    #
    def __getitem__(self, name: int|str):
        """ """
        try:
            self._labels.index(name)
            return self._item[name]
            #
            #
            #return Concepts(name=name,
            #                title=parameters, 
            #                materials=self._mesh[name]._materials,
            #                sections=self._mesh[name]._sections,
            #                properties= self._properties)            
        except ValueError:
            raise IndexError(f' Concept {name} not valid')
    #
    # --------------------
    # Mesh
    # --------------------
    #    
    def mesh(self):
        """ Meshing"""
        #1 / 0
        for name in self._labels:
            self._item[name].mesh()
#
#
class ConceptItem(ufoBasicModel):
    """
    Structural elements
    beam
    truss
    shell
    solid
    cable
    pile
    
    Utilities
    connection
    mass
    intersection
    spring
    
    Contact
    """
    __slots__ = ['_name', '_title', '_mesh',
                 'joints', '_beams', 'pile',
                 '_nodes','_materials', '_sections', '_load',
                 #'shells', 'membranes', 'solids', 'springs',
                 '_hinges', '_boundaries', '_properties',
                 '_groups', '_elements', '_point'] 
    
    def __init__(self, component:str,
                 title: str, mesh, 
                 #materials, sections, 
                 properties) -> None: 
        """
        """
        #
        #super().__init__(name, title)
        #
        self._component = component
        self._title = title
        self._mesh = mesh
        #
        self._materials = mesh._materials
        self._sections = mesh._sections
        self._nodes = mesh._nodes
        #
        #self._materials = Material(mesh_type="inmemory")
        #
        #self._sections = Section(mesh_type="inmemory")        
        #
        # 
        self.joints = Connection(points=self._nodes)
        #
        self._hinges = Releases()
        #
        #
        self._boundaries = BoundaryConcept(points=self._nodes,
                                           component=self._component)
        #
        # Points
        self._point = NodesIM(component=self._component,
                              boundary=self._boundaries)        
        #
        self._elements = ConceptElements(points=self._nodes,
                                         materials=self._materials,
                                         sections=self._sections)
        #
        # groups
        self._groups = Groups()
        #
        #self._load = load
        self._load = ConceptLoad(points=self._nodes,
                                 elements=self._elements,
                                 component=self._component, 
                                 boundaries=self._boundaries)
    #
    #
    #
    # --------------------
    # Element
    # --------------------
    #
    def point(self, values: None|list|dict=None,
              df=None):
        return self._point
    #
    #
    def boundary(self, values: None|list|dict = None,
                 df = None):
        """
        """
        if values:
            if isinstance(values, list):
                1/0
                for item in values:
                    self._boundaries[item[0]] = item[1:]
            else:
                raise IOError('boundary input not valid')
        #
        try:
            df.columns
            self._boundaries.df(df)
        except AttributeError:
            pass
            
        return self._boundaries
    #
    #
    #
    # --------------------
    # Loading
    # --------------------
    #    
    #
    # --------------------
    # Operations
    # --------------------
    #
    #
    # --------------------
    # Mesh
    # --------------------
    #    
    #
    def mesh(self):
        """ Meshing"""
        print('{:}'.format(52*'-'))
        meshing = MeshingConcept(concept=self)
        mesh = meshing.get_mesh()
        #mesh.renumbering()
        mesh.build()        
        return mesh
    #
    # --------------------
    # Plotting
    # --------------------
    #
    def plot(self, figsize:tuple = (10, 10)):
        """ """
        #print('--')
        return PlotConcept(cls=self, figsize=figsize)
    #
    #    
