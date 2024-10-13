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
                 component: int, points, beams):
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
                                    points=self.f2u_points)
        
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
    #
    @property
    def beam(self):
        """ return current beam"""
        return self._beam[self._beam_id]

    @beam.setter
    def beam(self, values):
        """ """
        if isinstance(values, (str, int)):
            try:
                self.f2u_beams[values]
                self._beam_id = values
            except IndexError:
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

