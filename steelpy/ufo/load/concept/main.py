# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
#
#
# package imports
from steelpy.ufo.load.process.main import MasterLoad
from steelpy.ufo.load.concept.wave_load import MetoceanLoadIM
#from steelpy.ufo.load.utils.actions import SelfWeight
#from steelpy.ufo.load.concept.timehistory import TimeHistory
from steelpy.ufo.load.concept.combination import LoadCombConcept
from steelpy.ufo.load.concept.load_case import BasicLoadConcept
#
#
class ConceptLoad(MasterLoad):
    __slots__ = ["_basic", "_combination", '_hydro', '_mass',
                 '_point', '_elements', '_boundaries', '_component']

    def __init__(self, points, elements, boundaries,
                 component: int | str) -> None:
        """
        """
        super().__init__(component)
        #
        self._point = points
        self._elements = elements
        self._boundaries = boundaries
        # self._component = component
        #
        self._basic = BasicLoadConcept(points=self._point,
                                       elements=self._elements,
                                       component=component)
        #
        self._hydro = MetoceanLoadIM(elements=self._elements)
        # self.th = TimeHistory()
        self._combination = LoadCombConcept(basic_load=self._basic,
                                            component=component)
        self._mass = LoadCombConcept(basic_load=self._basic,
                                     component=component)

    #
    #
    def basic(self):
        """
        """
        return self._basic

    #
    #
    def combination(self):
        """
        """
        return self._combination

    #
    #
    def metocean(self, condition: dict | list | None = None,
                 df=None):
        """
        design_load : max_BS
        criterion : select design load based on local (member) or global system
        """
        #
        # if isinstance(condition, dict):
        if condition:
            # cases = []
            for key, item in condition.items():
                self._hydro[key] = item

            # 1 / 0
            # if isinstance(values[0], (list, tuple)):
            #    for value in values:
            #        self._hydro[value[0]] = value[1:]
            # else:
            #    self._hydro[values[0]] = [*values[1:], 'local']
        # elif isinstance(condition, (list, tuple)):
        #    1 / 0

        #
        # dataframe input
        try:
            columns = list(df.columns)
            header = {}
            1 / 0
        except AttributeError:
            pass
            #
        return self._hydro

    #
    #
    # def time_history(self):
    #    """
    #    """
    #    return self.th
    #
    def mass(self):
        """
        """
        return self._mass
    #
    #
    # def __str__(self) -> str:
    #    """ """
    #    output = "\n"
    #    output += self._basic.__str__()
    #    output += self._combination.__str__()
    #    return output
    #
    # --------------------
    # Plotting
    # --------------------
    #
    # @property
    # def plot(self, figsize:tuple = (10, 10)):
    #    """ """
    #    return PlotLoad(cls=self, figsize=figsize)
#