#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#
#
# package imports
from array import array
# steelpy.f2uModel.load
from steelpy.ufo.load.process.actions import SelfWeight
from steelpy.ufo.load.process.basic_load import BasicLoadMain
from steelpy.ufo.load.concept.beam import BeamLoadItemIM
from steelpy.ufo.load.concept.node import NodeLoadItemIM
from steelpy.ufo.load.concept.wave_load import WaveLoadItemIM

#
#
# ---------------------------------
#
#
class BasicLoadConcept(BasicLoadMain):
    __slots__ = ['_load', '_labels', '_title', '_number',
                 'f2u_points', 'f2u_elements', '_component']

    def __init__(self, points, elements, component: int):
        """
        """
        super().__init__()
        self._labels: list = []
        self._title: list[str] = []
        self._number: array = array("I", [])        
        #
        self._load: dict = {}
        #
        # FIXME: reduce dependency
        self.f2u_points = points
        self.f2u_elements= elements
        self._component = component
    #
    def __setitem__(self, load_name: int, load_title: str) -> None:
        """
        load_name :
        load_title :
        """
        try:
            self._labels.index(load_name)
            self._title.index(load_title)
            raise Warning("Basic Load title {:} already defined".format(load_title))
        except ValueError:
            self._labels.append(load_name)
            self._title.append(load_title)
            # TODO: fix numbering
            load_id = len(self._load) + 1
            self._load[load_name] = LoadTypesConcept(name=load_name,
                                                     number=load_id,
                                                     title=load_title,
                                                     component=self._component, 
                                                     points=self.f2u_points,
                                                     beams=self.f2u_elements._beams)
            self._number.append(load_id)

    def __getitem__(self, load_name: str|int):
        """
        """
        try:
            return self._load[load_name]
        except KeyError:
            raise IOError("load case not defined")
    #
    #def __contains__(self, value) -> bool:
    #    return value in self._labels

    #def __len__(self) -> float:
    #    return len(self._labels)

    #def __iter__(self):
    #    """
    #    """
    #    items = list(set(self._labels))
    #    return iter(items)
    #
    def __delitem__(self, load_name: str|int):
        """
        """
        del self._load[load_name]
#
#
class LoadTypesConcept:
    """
    """
    __slots__ = ['_node', '_node_id', '_beam', '_beam_id',
                 '_selfweight', '_line', '_line_id', '_wave',
                  'name', 'number', 'title', 'f2u_points', 'f2u_beams']

    def __init__(self, name: str|int, number: int, title: str,
                 component: int, 
                 points, beams):
        """
        """
        self.f2u_points = points
        self.f2u_beams = beams
        #
        self.name = name
        self.number = number
        self.title = title
        #
        self._selfweight = SelfWeight()
        
        self._node = NodeLoadItemIM(load_name=self.name,
                                    load_title=self.title,
                                    nodes=self.f2u_points)
        
        self._beam = BeamLoadItemIM(load_name=self.name,
                                    load_title=self.title,
                                    component=component)
                                    #beams=self.f2u_beams)
        #
        self._wave = WaveLoadItemIM(load_name=self.name)

    #
    @property
    def gravity(self):
        """
        The self weight form allows you to specify multipliers to
        acceleration due to gravity (g) in the X, Y, and Z axes.
        If switched on, the default self weight acts in the Y axis
        with a magnitude and sign of -1."""
        return self._selfweight
    
    @gravity.setter
    def gravity(self, values):
        """ """
        self._selfweight[self.name] = [*values, self.title]
    #
    #
    @property
    def points(self):
        """ return all points"""
        return self._node
    #
    @property
    def point(self):
        """ return current point"""
        return self._node[self._node_id]

    @point.setter
    def point(self, values):
        """ set point"""
        # set connectivity
        try:
            node_id = self.f2u_points.get_point_name(values)
        except IOError:
            node_id = self.f2u_points.get_new_point(values)
        self._node_id = node_id
    #
    @property
    def line(self):
        """ """
        return self._line[self._line_id]
    #
    @line.setter
    def line(self, values):
        """ """
        1 / 0
    #
    @property
    def beams(self):
        """ return all beam"""
        return self._beam
    #
    @property
    def beam(self):
        """ return current beam"""
        return self._beam[self._beam_id]

    @beam.setter
    def beam(self, values):
        """ """
        #if isinstance(values[0], list):
        #    for value in values:
        #        self._beam[value[ 0 ] ] = value[ 1: ]
        #else:
        #    self._beam[ values[ 0 ] ] = values[ 1: ]
        if isinstance(values, str):
            try:
                self.f2u_beams[values]
                self._beam_id = values
            except IndexError:
                1 / 0
                raise IOError(f'Beam {values} not found')
            
        else:
            self._beam_id = values.name
    #
    #
    @property
    def wave(self):
        """ """
        return self._wave
    
    @wave.setter
    def wave(self, wave_data):
        """
        design_load : max_BSOTM
        wave_load=None, 
        """
        design_load: str = 'max_BS'
        criterion: str = 'local'
        
        if isinstance(wave_data, (str, int)):
            wave_data = [wave_data, design_load, criterion, self.title]
        
        elif isinstance(wave_data, (list, tuple)):
            try:
                wave_data.extend([criterion, self.title])
            except AttributeError:
                wave_data = wave_data + (criterion, self.title, )
        
        elif isinstance(wave_data, dict):
            raise NotImplementedError
        
        else:
            raise IOError('wave data not valid')
        #
        self._wave[self.name] = wave_data
#
#
#class TimeHistoryConcept(Mapping):
#    __slots__ = ['_load', '_labels', '_title', '_number', 'f2u_points']
#
#
# ---------------------------------
#
#
#class BasicLoadIM(BasicLoadBasic):
#    """
#    FE Load Cases
#    
#    LoadType
#        |_ name
#        |_ number
#        |_ basic
#        |_ combination_level
#        |_ time_series
#        |_ 
#        |_ temperature
#    
#    **Parameters**:  
#      :number:  integer internal number 
#      :name:  string node external name
#    """
#    __slots__ = ['_labels', '_title','_number', 'gravity',
#                 '_basic', '_f2u_elements', '_f2u_nodes']
#
#    def __init__(self, nodes, elements):
#        """
#        """
#        super().__init__()
#        self._labels: list = []
#        self._title: list[str] = []
#        self._number: array = array("I", [])        
#        #
#        self._basic: dict = {}
#        self._f2u_elements = elements
#        self._f2u_nodes = nodes
#    #
#    def __setitem__(self, load_name: str|int, load_title: str) -> None:
#        """
#        load_name :
#        load_title :
#        """
#        1 / 0
#        try:
#            self._labels.index(load_name)
#            self._title.index(load_title)
#            raise Warning("Basic Load {:} already defined".format(load_title))
#        except ValueError:
#            self._labels.append(load_name)
#            self._title.append(load_title)
#            load_id = next(self.get_number())
#            self._basic[load_name] = LoadTypeInMemory(name=load_name,
#                                                      number=load_id,
#                                                      title=load_title,
#                                                      nodes=self._f2u_nodes,
#                                                      elements=self._f2u_elements)
#            self._number.append(load_id)
#
#    def __getitem__(self, load_name: str|int):
#        """
#        """
#        try:
#            return self._basic[load_name]
#        except KeyError:
#            raise IOError("Basic load case {:} not defined".format(load_name))
#
#    #
#    def __delitem__(self, load_name: str|int):
#        """
#        """
#        del self._basic[load_name]
#    #
#    #
#    # def get_basic_load(self, elements, nodes, materials,
#    #                   sections):
#    #    """
#    #    """
#    #    1/0
#    #    #return get_basic_load(self, elements, nodes, 
#    #    #                      materials, sections)
#    #
##
##
##
#class LoadTypeInMemory(LoadTypeBasic):
#    """
#    """
#    __slots__ = ['_node', '_beam', '_selfweight',
#                 'name', 'number', 'title', '_wave']
#                 #'_f2u_beams', '_f2u_nodes']
#
#    def __init__(self, nodes, elements):
#        """
#        """
#        #super().__init__(name, number, title)
#        #self.name = name
#        #self.number = number
#        #self.title = title
#        self._selfweight = SelfWeight()
#        self._node = NodeLoadItemIM(load_name=name,
#                                    load_title = title, 
#                                    nodes=nodes)
#        beams = elements.beams()
#        self._beam = BeamLoadItemIM(load_name=name,
#                                    load_title = title,
#                                    beams=beams)
#        #
#        self._wave = WaveLoadItemIM(load_name=name)
#    #
#
