#
# Copyright (c) 2009-2023 steelpy
#
from __future__ import annotations
#
# Python stdlib imports
#from typing import NamedTuple

# package imports
#from steelpy.metocean.regular.current.main_current import MeanCurrent
from steelpy.metocean.regular.stokes.Stokes import StokesModule
from steelpy.metocean.regular.fourier.Fourier import FourierModule
from steelpy.metocean.regular.cnoidal.Cnoidal import CnoidalModule

class RegularWaves:
    __slots__ = ['_wave']

    def __init__(self)-> None:
        """
        """
        self._wave = FourierModule()
    #
    @property
    def Stokes(self):
        """ """
        return StokesModule
    #
    @property
    def Cnoidal(self):
        """ """
        return CnoidalModule
    #
    @property
    def Fourier(self):
        """ """
        return FourierModule
    #
    #def current(self):
    #    """ """
    #    return self._current
    #
    def __setitem__(self, case_name: int,
                    case_data: list[float]|dict[str, float]) -> None:
        """
        case_name : Wave name
         H : Wave height [unit length]
         T : Wave Period [second]
         d : Water depth LTH + Tide and Surge [unit length]
         Lw : Wave Length [unit length]
         Phase : Wave phase [degree]
         Order : ??
        """
        self._wave[case_name] = case_data
    #
    #
    def __getitem__(self, case_name: int|str) -> tuple:
        """
        case_name : Wave name
        """
        return self._wave[case_name]
    #
    #
    def df(self, df):
        """ """
        1/0
    #
    #