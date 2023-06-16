# 
# Copyright (c) 2009-2023 fem2ufo
#
# Python stdlib imports
from __future__ import annotations

#
# package imports
from .marine_growth import MarineGrowth
from .morison import CdCmCoefficients
from .hydro_diametre import HydroDiametre
from .flooding import Flooding


#
#@dataclass
class Hydrodynamic:
    """
    """
    __slots__ = ['flooding', 
                 'conductor_shielding',
                 'element_refining',
                 '_buoyancy_area',
                 '_cdcm', 
                 '_air_drag',
                 '_marine_growth', #'mg_default',
                 '_hydro_diametre',
                 '_non_hydro']
    
    def __init__(self) -> None:
        """
        """
        #
        #global mg_default
        #mg_default = None
        #
        self.flooding = Flooding()
        self._cdcm = CdCmCoefficients()
        self._hydro_diametre = HydroDiametre()
        self._marine_growth = MarineGrowth()
        # TODO: 
        self._non_hydro = {}
        self._buoyancy_area = {}
        self._air_drag = {}
    #  
    #
    #@property
    def CdCm(self):
        """
        """
        return self._cdcm
    #
    @property
    def Morison(self):
        """
        """
        return self._cdcm
    #
    @property
    def flooded(self):
        """
        """
        return self.flooding
    #
    #@property
    def marine_growth(self):
        """
        """
        return self._marine_growth
    #
    @property
    def diameter(self):
        """
        """
        return self._hydro_diametre
    #
    @property
    def buoyancy_area(self):
        """
        """
        return self._buoyancy_area
    #
    @property
    def non_hydro(self):
        """
        """
        return self._non_hydro
    #
    #
    @property
    def air_drag(self):
        """
        """
        return self._air_drag    
#
#
