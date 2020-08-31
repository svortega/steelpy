# 
# Copyright (c) 2009-2020 fem2ufo
# 


# Python stdlib imports
from typing import Tuple, Dict, List

# package imports
#
from steelpy.f2uModel.concept.joint import Connection
from steelpy.f2uModel.concept.element import ConceptElements
from steelpy.f2uModel.concept.boundary import ConceptBoundaries
from steelpy.f2uModel.concept.load import ConceptLoad
#from steelpy.f2uModel.concept.beam import Beams
#from steelpy.f2uModel.concept.shell import Shells
#from steelpy.f2uModel.concept.spring import Springs
from steelpy.f2uModel.mesh.geometry import Releases




class Concepts:
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
    __slots__ = ['joints', 'beam', 'pile', 'points', #'_load',
                 'shells', 'membranes', 'solids', 'springs', 'load',
                 '_hinges', 'boundary', '_properties', '_name', '_mesh']
    
    def __init__(self, mesh, properties): # load, 
        """
        """
        self._mesh = mesh
        self.points = mesh.nodes
        self.joints = Connection(points=self.points)
        self._hinges = Releases()
        #
        self.boundary = ConceptBoundaries()
        #
        self.beam = ConceptElements(element_type="beam", points = self.points,
                                    materials=mesh.materials, sections=mesh.sections)

        #self.truss = Elements(element_type='truss', points = self.points,
        #                      materials=mesh.materials, sections=mesh.sections)
        #
        #self.piles = Beams(beam_type='pile', 
        #                   points=self.points, elements=mesh.elements,
        #                   materials=mesh.materials, sections=mesh.sections, 
        #                   properties=properties,hinges=self._hinges)
        #
        #self.shells = Shells(shell_type='shell',
        #                     points=self.points, materials=mesh.materials)
        #
        #self.springs = Springs(spring_type='spring', 
        #                       points=self.points, materials=mesh.materials)
        #
        #self._load = load
        self.load = ConceptLoad(self.points, self.beam) # self._load,
    #
    #
    #
    #
    #
    #def __iter__(self):
        #"""
        #"""
        #for _name in self.beams.keys():
        #    yield self.beams[_name]
        #
        #for _name in self.trusses.keys():
        #    yield self.trusses[_name]
        #
        #for _name in self.piles.keys():
        #    yield self.piles[_name]
    #
    #
    #@property
    #def hinges(self):
    #    """
    #    """
    #    return self._hinges
    #
    #@property
    #def get_name(self):
    #    """
    #    """
    #    return self._name
    #
    #
    #def __delattr__(self, name:str) -> None:
    #    """
    #    """
    #    _joints = []
    #    if 'piles' in name:
    #        _tobe_deleted = [key for key in self.piles.keys()]
    #        for _item in _tobe_deleted:
    #            for _element in self.piles[_item].elements:
    #                _joints.extend(_element.connectivity)
    #            del self.piles[_item]
    #    #
    #    # deleting redundant link joints
    #    _joints = set(_joints)
    #    for _number in _joints:
    #        try:
    #            del self.joints[_number]
    #        except KeyError:
    #            print('-->', _number)
    #
