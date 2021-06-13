#
# Copyright (c) 2009-2019 fem2ufo
#


# Python stdlib imports
import math
from typing import NamedTuple, Tuple, List


# package imports
#import numpy as np
#import matplotlib.pyplot as plt
from steelpy.process.units.units import Units



class PlotModel:

    def __init__(self):
        """
        """
        self._units = Units()

    #
    @property
    def mesh(self):
        """"""
        return self._mesh

    @mesh.setter
    def mesh(self, value):
        """"""
        self._mesh = value
    #
    #
    def plot_model(self, scale:float=1,
                   offset:List[float]=[0,0,0]):
        """
        """
        nodes = self._mesh.nodes
        max_coord, min_coord = nodes._get_maxmin
        centre_x = (max_coord[0] - min_coord[0]) / 2 + min_coord[0] + - offset[0]
        centre_y = (max_coord[1] - min_coord[1]) / 2 + min_coord[1] + - offset[1]
        centre_z = (max_coord[2] - min_coord[2]) / 2 + min_coord[2] + - offset[2]
        max_plot_range = max(max_coord)
        max_plot_range



