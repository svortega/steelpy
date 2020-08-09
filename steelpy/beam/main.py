# 
# Copyright (c) 2019-2020 steelpy
# 

# Python stdlib imports
#from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Tuple, Union


# package imports
from steelpy.beam.static.operations import BeamOps
from steelpy.beam.static.load import Load
from steelpy.beam.static.support import Support
from steelpy.f2uModel.sections.main import Sections
from steelpy.f2uModel.material.main import Materials



#
class Beam:
    __slots__ = ['_labels', '_fixity', '_supcoord',
                 '_load', '_section', '_material', 'steps',
                 '_beam_length', '_reactions', '_response']

    def __init__(self, beam_length):
        """
        beam_length : total length of the beam 
        section_type:str
        """
        print("-- module : Beam Version 0.10dev")
        print('{:}'.format(52*'-'))
        #
        self._beam_length:float = beam_length.value
        self.steps: int = 10
        #
        self._material = Materials()
        self._material['beam'] = 'elastic'
        self._section = Sections()
        #self._section['beam'] = section_type
        #
        self._labels:List = []
        self._fixity:List = []
        self._supcoord:List = []
        self._reactions:List = []
        #
        self._load = Load(cls=self)
        self._response = BeamOps(cls=self)
    #
    @property
    def section(self):
        """
        """
        return self._section['beam']

    @section.setter
    def section(self, section_type:str):
        """
        section_type : select beam's cross section shape (tubular/rectagulat,I_section)
        """
        #self._section = Sections ()
        self._section["beam"] = section_type
    #
    @property
    def material(self):
        """
        """
        return self._material['beam']

    #@material.setter
    #def material(self, value):
    #    """
    #    """
    #    self._material = value
    #
    @property
    def support(self):
        """
        """
        return Support(self)
    #
    @property
    def load(self):
        """
        """
        return self._load
    #
    #
    def join(self, other, join_type):
        """
        b1 = Beam()
        b2 = Beam()
        b = b1.join(b2, 'fixed')
        """
        pass
    #
    @property
    def shear(self):
        """ """
        return self._response.shear(self.steps)
    #
    @property
    def bending_moment(self):
        """ """
        return self._response.bending(self.steps)
    #
    @property
    def slope(self):
        """ """
        return self._response.slope(self.steps)
    #
    @property
    def deflection(self):
        """ """
        return self._response.deflection(self.steps)
    #
    @property
    def response(self):
        """ """
        return  self._response
#
#
