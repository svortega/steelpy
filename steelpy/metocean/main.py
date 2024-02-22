# 
# Copyright (c) 2019 steelpy
# 
# Python stdlib imports
from __future__ import annotations
#from collections.abc import Mapping
from datetime import datetime as dt
#from dataclasses import dataclass
#from typing import NamedTuple
#import os
import re


# package imports
from steelpy.utils.units.main import Units
from steelpy.utils.sqlite.main import ClassMainSQL
#
from steelpy.metocean.wave.regular.process.bsotm import BSOTM
from steelpy.metocean.hydrodynamic.main import HydroProperty
from steelpy.metocean.process.criteria import HydroCriteria


#
#
class Metocean(ClassMainSQL):
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
    
    **Parameter**:  
      :number:  integer internal number 
      :name:  string node external name
      
      : seastate : metocean combination
    """
    #
    __slots__ = ['_name', '_criteria', '_rho_w',
                 '_properties', '_build', 'db_file']
    #
    def __init__(self, name:str|None = None,
                 sql_file:str|None = None):
        """
        """  
        super().__init__(name=name, sql_file=sql_file)
        #
        self._rho_w: float = 1032.0  # kg/m^3  
        #
        self._properties = HydroProperty(rho_w=self._rho_w,
                                         db_file=self.db_file)
        #
        self._criteria = HydroCriteria(rho_w= self._rho_w,
                                       properties=self._properties, 
                                       db_file=self.db_file)
    #
    #
    #
    # ------------------------------------------
    #
    @property
    def rho_w(self):
        """ """
        units = Units()
        return self._rho_w * units.kg / units.m**3

    @rho_w.setter
    def rho_w(self, value:Units):
        """ """
        self._rho_w = value.convert('kilogram/metre^3').value
    #
    # ------------------------------------------
    #@property
    def property(self):
        """
        """
        return self._properties
    #
    #
    def criteria(self):
        """Metocean criteria"""
        return self._criteria
    #
    # ------------------------------------------
    #
    def get_load(self, mesh, kinematic,
                 condition:int|None = None,
                 rho:float = 1025):
        """
        condition :
            1 - Linear wave (dafault)
            2 - Non-linear wave
        """
        #if not condition:
            
        wforce = BSOTM(kinematic, condition, rho)
        wforce.wave_force(mesh=mesh)
        #print('-->')
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
        bs, otm = bsotm.solveBSOTM(D=Dp, L=Lp)
        #
        return bs, otm 
        #return surface
#

