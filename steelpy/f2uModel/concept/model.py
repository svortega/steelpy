# 
# Copyright (c) 2009-2019 fem2ufo
# 


# Python stdlib imports
from copy import copy
import logging
from typing import Dict, ClassVar

# package imports
from fem2ufo.f2u_model.model.component import Component
from fem2ufo.f2u_model.femodel.metocean.metocean import Metocean
from fem2ufo.f2u_model.femodel.foundation.soil import Foundation
from fem2ufo.f2u_model.femodel.load.load import Load
from fem2ufo.f2u_model.femodel.analysis.analysis import Analysis
from fem2ufo.f2u_model.model.assemble import Assemble
from fem2ufo.f2u_model.model.assemble import InputData
from fem2ufo.process.control import euclid


#
#
class F2U_model:
    """
    fem2ufo fe concept model
    
    Model
        |_ name
        |_ number
        |_ data
        |_ component
        |_ load
        |_ foundation
        |_ metocean
        |_ analysis
        |_ dataset [time_history, etc.]
    
    """
    __slots__ = ('_components', '_component_name',
                 '_loading', '_load_name',
                 '_analysis', '_assemble', 'data',
                 '_metocean', '_metocean_name',
                 '_foundation', '_foundation_name',
                 '_dataset', '_dataset_name')

    def __init__(self) -> None:
        """
        """
        self._components: Dict = {}
        self._loading: Dict = {}
        #
        self._foundation: Dict = {}
        self._metocean: Dict = {}
        self._analysis: Dict = {}
        #
        self._component_name: str = None
        self._load_name: str = None
        self._metocean_name: str = None
        self._foundation_name: str = None
        #
        self.data: ClassVar = InputData()
        self._assemble: ClassVar = Assemble()

    #
    @property
    def assemble(self) -> ClassVar:
        """
        """
        return self._assemble

    #
    @property
    def components(self) -> ClassVar:
        """
        """
        return self._components

    #
    @property
    def component(self) -> ClassVar:
        """
        """
        return self._components[self._component_name]

    @component.setter
    def component(self, name: str) -> None:
        """
        component
            |_ number: integer
            |_ name: string
            |_ nodes: list of node's class
            |_ elements: list of element's class
            |_ materials: list of material's class
            |_ sets: list of groups (elements, nodes, etc)
            |_ sections: list of section's class
            |_ vectors: list of guide points
            |_ eccentricities: list of eccentricities
            |_ joints: list of joint's class
            |_ hinges: list of hinges definitios
            |_ loads: list of load's class
            |_ data: FE model data
            |_ units: FE model units
            |_ boundaries: list of FE model boundary
        """
        self._component_name = name
        try:
            return self._components[name]
        except KeyError:
            number = len(self._components) + 1
            self._components[name] = Component(name, number)
            # self._loading[name] = Load(component=self._component[name])
            self.assemble[self._component_name] = number

    #
    @property
    def load(self) -> ClassVar:
        """
        """
        return self._loading[self._load_name]

    @load.setter
    def load(self, load_name: str) -> None:
        """
        """
        self._load_name = load_name
        try:
            return self._loading[load_name]
        except KeyError:
            self._loading[load_name] = Load(component=self.component)
            self.assemble[self._component_name].load = load_name

    #
    @property
    def loads(self) -> ClassVar:
        """
        """
        return self._loading

    #
    @property
    def metocean(self) -> ClassVar:
        """
        """
        return self._metocean[self._metocean_name]

    @metocean.setter
    def metocean(self, name: str) -> None:
        """
        metocean
            |_ name
            |_ number
            |_ data
            |_ type
            |_ wave
            |_ current
            |_ wind
            |_ cdcm
            |_ air_drag
            |_ non_hydro
            |_ elevation
            |_ hydro_diameter
            |_ buoyancy
            |_ flooded
            |_ seastate
        """
        self._metocean_name = name
        try:
            return self._metocean[name]
        except KeyError:
            number = len(self._metocean) + 1
            self._metocean[name] = Metocean(name, number)
            self.assemble[self._component_name].metocean = name

    #
    #
    @property
    def analysis(self) -> ClassVar:
        """
        """
        return self._analysis

    @analysis.setter
    def analysis(self, name: str = None) -> None:
        """
        """
        if not name:
            name = self.name
        number = len(self._analysis) + 1
        self._analysis[name] = Analysis(name, 1)

    #
    @property
    def foundation(self) -> ClassVar:
        """
        """
        return self._foundation[self._foundation_name]

    @foundation.setter
    def foundation(self, name: str) -> None:
        """
        """
        self._foundation_name = name
        try:
            return self._foundation[name]
        except KeyError:
            number = len(self._foundation) + 1
            self._foundation[name] = Foundation(name, number)
            # FIXME: 
            try:
                self.assemble[self._component_name].foundation = name
            except KeyError:
                pass
    #
    #
    def _get_link_nodes(self, component) -> Dict:
        """
        """
        _bound_node: Dict = {}
        for key, _joint in component.concept.joints.items():
            try:
                if 'link' in _joint.type:
                    _bound_node[key] = _joint.node
            except AttributeError:
                continue
        #
        # return {key:item if item.type == 'link' else None
        #        for key, item in _bound_node.items()}
        return _bound_node

    #
    # TODO: need to be checked
    def find_super_node(self, _tol:float) -> None:
        """
        Find supernodes 
        """
        logging.info('    * Finding link nodes')
        # remove redundant boundaries 
        _components = {}
        for _component in self._components:
            _components[_component] = self._get_link_nodes(self._components[_component])
        _subcomponent = copy(_components)
        # loop all components
        for key1, item1 in _components.items():
            _boundary1 = self._components[key1]._boundaries
            del _subcomponent[key1]  # delete duplicate component in subcomponent
            for key2, item2 in _subcomponent.items(): # loop rest components
                _boundary2 = self._components[key2]._boundaries
                _nodes = self._components[key2].mesh.nodes # get nodes of subcomponent
                for _node1 in item1.values(): # loop nodes of components
                    _point = euclid.Point3(_node1.x, _node1.y, _node1.z)
                    _flag = False
                    for _node2 in item2.values(): # loop nodes of subcomponents
                        # set tolerance of link node
                        _sphere = euclid.Sphere(euclid.Point3(_node2.x, _node2.y, _node2.z), _tol)
                        # check if node within the boundary
                        if _sphere.intersect(_point):
                            _member2_name = _node2.sets.elements[0]
                            _element2 = self._components[key2].mesh.elements[_member2_name]
                            _nodes_element2 = _element2.nodes
                            for x, _node_link in enumerate(_nodes_element2):
                                if _node2.number == _node_link.number:
                                    logging.info('    * Link node found {:6.0f} --> {:}'
                                                 .format(_node1.number, _node2.number))
                                    print('    * Link node found {:6.0f} --> {:}'
                                          .format(_node1.number, _node2.number))
                                    _element2.connectivity[x] = _node1.number
                                    # include node componenent 1 (to be confirmed)
                                    _nodes[_node1.number] = _node1[:3]
                                    _boundary2[_node1.number] = _boundary1[_node1.number][:6]
                                    # delete node from component 2
                                    del self._components[key1].concept.joints[_node1.number]
                                    _flag = True
                                    break
                            #
                            if _flag:
                                break
                            logging.warning('   ** --> Link node {:} not found'
                                            .format(_node1.number))
                            print('   ** --> Link node {:} not found'
                                  .format(_node1.number))
        #
        #print('--->')
    #
    #
    #    
#
# class F2U_foundation:
#    __slots__ = ('data', '_foundation', '_foundation_name')
#    #
#    def __init__(self) -> None:
#        """
#        """
#        self._foundation: Dict={}
#        self._foundation_name:str = None
#        self.data = InputData()   
#
#
