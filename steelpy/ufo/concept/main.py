# 
# Copyright (c) 2009 steelpy
# 
# Python stdlib imports
from __future__ import annotations


# package imports
from steelpy.ufo.utils.main import ufoBasicModel, ModelClassBasic
#
from steelpy.ufo.concept.elements.joint import Connection
from steelpy.ufo.concept.elements.boundary import BoundaryConcept
from steelpy.ufo.concept.elements.geometry import Releases
#
from steelpy.ufo.concept.meshing.main import MeshingConcept
#
from steelpy.ufo.concept.elements.sets import Groups
from steelpy.ufo.concept.elements.points import NodesIM
from steelpy.ufo.concept.elements.main import ConceptElements
#
from steelpy.ufo.load.main import ConceptLoad
#
from steelpy.ufo.plot.main import PlotConcept
#
from steelpy.utils.dataframe.main import DBframework
#
#
#
class Concept(ModelClassBasic):
    """ f2u Concept model Class """
    __slots__ = ['_labels', '_item', '_component',
                 '_properties', '_mesh']
    #
    def __init__(self, component:str|int,
                 mesh, properties):
        """
        """
        super().__init__(component)
        #
        self._mesh = mesh
        self._properties = properties
        self._item : dict = {}
        self._labels:list = []
    #
    def __setitem__(self, name: int|str, title: int|str) -> None:
        """
        """
        try:
            self._labels.index(name)
            raise Exception('Concept {:} already exist'.format(name))
        except ValueError:
            self._labels.append(name)
            self._mesh[name] = title
            self._item[name] = ConceptItem(component=name,
                                           title=title,
                                           mesh=self._mesh[name], 
                                           properties= self._properties)
    
    #
    def __getitem__(self, name: int|str):
        """ """
        try:
            self._labels.index(name)
            return self._item[name]
        except ValueError:
            raise IndexError(f' Concept {name} not valid')
    #
    # --------------------
    # Mesh
    # --------------------
    #    
    def mesh(self):
        """ Meshing"""
        for name in self._labels:
            self._item[name].mesh()
    #
    #
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
                 properties) -> None: 
        """
        """
        #
        super().__init__(component)
        #
        self._title = title
        self._mesh = mesh
        #
        self._materials = mesh._materials
        self._sections = mesh._sections
        self._nodes = mesh._nodes
        #
        # 
        self.joints = Connection(points=self._nodes)
        #
        self._hinges = Releases()
        #
        # Points
        self._point = NodesIM(component=self._component)
        #                      boundary=self._boundaries)
        #
        #
        self._boundaries = BoundaryConcept(points=self._point,
                                           component=self._component)
        #        
        #
        self._elements = ConceptElements(points=self._nodes,
                                         materials=self._materials,
                                         sections=self._sections,
                                         component=self._component, )
        #
        # groups
        self._groups = Groups()
        #
        #self._load = load
        self._load = ConceptLoad(points=self._point,
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
        """ """
        if values:
            if isinstance(values, dict):
                mname = values['name']
                if isinstance(mname, (list | tuple)):
                    db = DBframework()
                    dfnew = db.DataFrame(values)
                    self._point.df = dfnew
                else:
                    mname = values.pop('name')
                    self._point[mname] = values
            
            elif isinstance(values, (list,tuple)):
                if isinstance(values[0], (list,tuple)):
                    for value in values:
                        self._point[value[0]] = value[1:]
                elif isinstance(values[0], dict):
                    for value in values:
                        self._point[value['name']] = value
                else:
                    self._point[values[0]] = values[1:]
        #
        # dataframe input
        try:
            df.columns   
            self._point.df = df
        except AttributeError:
            pass        
        return self._point
    #
    #
    def boundary(self, values: None|list|dict = None,
                 df = None):
        """
        """
        if values:
            1 / 0
            if isinstance(values, dict):
                bname = values['name']
                if isinstance(bname, (list | tuple)):
                    db = DBframework()
                    self._boundaries.df = db.DataFrame(values)
                else:
                    sname = values.pop('name')
                    self._boundaries[sname] = values   

            elif isinstance(values, (list, tuple)):
                if isinstance(values[0], (list, tuple)):
                    for item in values:
                        self._boundaries[item[0]] = item[1:]

                elif isinstance(values[0], dict):
                    for item in values:
                        self._boundaries[item['name']] = item
                else:
                    self._boundaries[values[0]] = values[1:]

            else:
                raise IOError('Boundary input data not valid')
        #
        try:
            df.columns
            self._boundaries.df = df
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
