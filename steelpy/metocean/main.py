# 
# Copyright (c) 2019-2023 steelpy
#
from __future__ import annotations
# 
# Python stdlib imports
#from typing import Union, Dict

# package imports
#from steelpy.metocean.irregular.main import WaveIrregular
from steelpy.metocean.regular.main import RegularWaves
#from steelpy.metocean.irregular.spectrum import Sprectrum
from steelpy.process.units.main import Units
#from steelpy.metocean.interface.wave import RegularWaves, IregularWaves
#from steelpy.metocean.interface.wind import Winds
#from steelpy.metocean.interface.current import Currents
#from steelpy.f2uModel.load.combination import MetoceanCombination
from steelpy.metocean.regular.process.bsotm import BSOTM
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
    def regular_wave(self, values:None|list=None,
                     df=None):
        """
        """
        if values:
            print('-->')
            1/0
        else:
            try:
                df.columns            
                self._regular_wave.df(df)
            except AttributeError:
                pass
        return self._regular_wave
    #  
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
    #
    def get_load(self, mesh, kinematic,
                 condition:int = 1,
                 rho:float = 1025):
        """
        condition :
            1 - Linear wave (dafault)
            2 - Non-linear wave
        """
        wforce = BSOTM(kinematic, condition, rho)
        wforce.wave_force(mesh=mesh)
        print('-->')
    #
    #
    def pile_response(self, D:float|Units,
                      L:float|Units,
                      kinematic,
                      condition: int = 1,
                      rho:float = 1025):
        """
        D : Pile diametre
        L : Pile length
        kinematics : kinematic dataframe
        condition :
            1 - Linear wave (dafault)
            2 - Non-linear wave
        """
        Dp = D.convert('metre').value
        Lp = L.convert('metre').value
        #
        #if kinematic._type == 'regular':
        #bs, ovtm = bsotm_reg(kinematic, D_pile, condition)
        #else:
        #    #bs, ovtm = bsvtm(kinematic, D_pile, condition)
        #    raise NotImplemented
        #
        bsotm = BSOTM(kinematic, condition, rho)
        bs, otm = bsotm.solve(D=Dp, L=Lp)
        #
        return bs, otm 
        #return surface    
#
#
# HYDRODYNAMICS Section
#

