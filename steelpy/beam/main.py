# 
# Copyright (c) 2019-2022 steelpy
# 

# Python stdlib imports
#from collections.abc import Mapping
from typing import NamedTuple, Dict, List, Tuple, Union


# package imports
from steelpy.beam.static.operations import BeamResponse
from steelpy.beam.static.load import Load, Combination
from steelpy.beam.static.support import Support
from steelpy.sections.main import Sections
from steelpy.material.main import Materials
from steelpy.design.beam.main import BeamDesign


#
class Beam:
    __slots__ = ['_support', '_load', '_section', '_material', 
                 'steps', 'beam_length', '_response', 
                 '_load_combination', '_design', 'beam_name'] 

    def __init__(self, name:Union[int,str]='beam_design',
                 mesh_type:str='inmemory'):
        """
        beam_length : total length of the beam 
        section_type:str
        """
        print("-- module : Beam Version 0.15dev")
        print('{:}'.format(52*'-'))
        #
        self.beam_name = name
        self.steps: int = 10
        #
        db_file:str = "beam_f2u.db"
        self._material = Materials(mesh_type=mesh_type,
                                   db_file=db_file)
        self._material['beam'] = 'elastic'
        self._section = Sections(mesh_type=mesh_type,
                                 db_file=db_file)
        #
        self._load = Load(cls=self)
        self._response = BeamResponse(cls=self)
        self._support = Support(cls=self)
        self._load_combination = Combination(cls=self)
        #
        component = "-"*12
        self._design = BeamDesign(name=name, component=component)
    #
    @property
    def length(self):
        """ """
        return self.beam_length
    
    @length.setter
    def length(self, beam_length):
        """ """
        try:
            self.beam_length = beam_length.value
        except AttributeError: # unit must be in metres
            self.beam_length = beam_length
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
    #    self._material['beam'] = value
    #
    @property
    def support(self):
        """
        """
        return self._support #Support(self)
    
    #@support.setter
    #def support(self, value:list):
    #    """
    #    """
    #    self._support = value
    #    #Support(self)
    #
    @property
    def load(self):
        """
        """
        return self._load
    #
    @property
    def load_combination(self):
        """
        """
        return self._load_combination    
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
    @property
    def design(self):
        """
        """
        self._design.material = self._material['beam']
        self._design.section = self._section["beam"]
        self._design.L = self.length
        self._design._stress = self._response.stress()
        return self._design
#
#
