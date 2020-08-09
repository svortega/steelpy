# Copyright (c) 2019-2020 steelpy

# Python stdlib imports

# package imports
from steelpy.process.units.main import Units
from steelpy.beam.main import Beam
#

class BeamDesign:
    
    __slots__ = ['_units', '_mesh', '_elements']
    
    def __init__(self):
        """
        """
        self._units = Units()
    #
    @property
    def mesh(self):
        return self._mesh
    
    @mesh.setter
    def mesh(self, mesh):
        self._mesh = mesh
    
    #
    @property
    def section(self):
        """
        """
        pass