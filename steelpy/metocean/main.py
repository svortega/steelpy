# 
# Copyright (c) 2019-2021 steelpy
#
# 
# Python stdlib imports
from typing import Union, Dict

# package imports
#from steelpy.metocean.irregular.main import WaveIrregular
from steelpy.metocean.regular.main import RegularWaves
#from steelpy.metocean.irregular.spectrum import Sprectrum
#from steelpy.process.units.units import Units
#from steelpy.metocean.interface.wave import RegularWaves, IregularWaves
#from steelpy.metocean.interface.wind import Winds
#from steelpy.metocean.interface.current import Currents
#from steelpy.f2uModel.load.combination import MetoceanCombination
#
#
class Metocean:
    """
    FE Metocean Class
    
    Metocean
        |_ name
        |_ number
        |_ data
        |_ type
        |_ wave
        |_ current
        |_ wind
        |_ cdcm
        |_ non_hydro
        |_ elevation
        |_ hydro_diameter
        |_ buoyancy
        |_ flooded
        |_ seastate
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
      
      : seastate : metocean combination
    """
    #
    __slots__ = ['_regular_wave', '_iregular_wave',
                 '_current', '_wind', '_units',
                 '_combination', '_spectrum']
    #
    def __init__(self):
        """
        """
        self._regular_wave = RegularWaves()
        #self._iregular_wave = WaveIrregular()
        #self._wind = Winds()
        #self._current = Currents()
        #self._units = Units()
        #self._combination = MetoceanCombination()
        #self._spectrum = Sprectrum()
    #
    #@property
    #def spectrum(self):
    #    """
    #    """
    #    return self._spectrum
    #
    #@property
    #def irregular_wave(self):
    #    """
    #    """
    #    return self._iregular_wave
    #
    #@property
    def regular_wave(self):
        """
        """
        return self._regular_wave
    #
    #@property
    #def linear_wave(self):
    #    """
    #    """
    #    return self.regular_waves
    #
    #@property
    #def deterministic_wave(self):
    #    """
    #    """
    #    return self.regular_waves     
    ##
    #@property
    #def wind(self):
    #    """
    #    """
    #    return self._wind  
    ##
    #@property
    #def current(self):
    #    """
    #    """
    #    return self._current
    ##
    #@property
    #def combination(self):
    #    """
    #    """
    #    return self._combination    
#
#
# HYDRODYNAMICS Section
#

