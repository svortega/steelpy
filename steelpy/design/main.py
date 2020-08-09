# Copyright (c) 2019-2020 steelpy

# Python stdlib imports

# package imports
from steelpy.design.clamp.clamp import ClampDesign
from steelpy.design.beam.design import BeamDesign

class Design:
    """
    """
    __slots__ = ['clamp_design', 'beam_design']
    
    def __init__(self):
        """
        """
        self.clamp_design = ClampDesign()
        self.beam_design = BeamDesign()
    
    @property
    def clamp(self):
        """
        """
        return self.clamp_design
    
    @property
    def beam(self):
        """
        """
        return self.beam_design