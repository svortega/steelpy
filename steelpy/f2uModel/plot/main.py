#
# Copyright (c) 2009-2023 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
#

# package imports
from steelpy.f2uModel.plot.frame import PlotFrame
from steelpy.f2uModel.plot.load import PlotLoad
#from steelpy.f2uModel.plot.frame import (init_frame_plot,
#                                          add_beams, add_supports,
#                                          add_materials,
#                                          add_beams2, add_nodes)
#                                          #add_global_axes)
#import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.widgets import RadioButtons, CheckButtons
#
#
#
#
#
class PlotConcept:
    __slots__ = ['_mesh']
    
    def __init__(self, cls):
        """
        """
        self._mesh = cls
    #
    def concept(self, verbosity:bool=False, show=True):
        """ """
        plot_num = 2
        points = self._concept.points
        lims = points._get_maxmin
        ax = init_3D_plot([lims[0][0], lims[1][0]],
                          [lims[0][1], lims[1][1]],
                          [lims[0][2], lims[1][2]])        
        add_items_concept(self._concept, ax, plot_num, 
                          verbosity=False)
        #
        basic = self._concept.load.basic[1]
        elements =  self._concept.beam
        add_concept_load(ax, elements, points, basic)        
        #
        if show:
            #plt.draw()
            plt.show()
        else:
            return ax
#
#
#
class PlotMesh:
    __slots__ = ['_mesh', '_frame', '_load']
    
    def __init__(self, mesh, figsize):
        """
        """
        self._mesh = mesh
        self._frame = PlotFrame(mesh=self._mesh, figsize=figsize)
        self._load = PlotLoad(mesh=self._mesh, figsize=figsize)
    #
    def frame(self, show=True):
        """ plot mesh element, nodes & boundary"""
        #print('--')
        #
        ax = self._frame.frame()
        #
        if show:
            plt.show()
        else:
            return ax        
    #
    def material(self):
        """plot material"""
        self._frame.materials()
    #
    #
    def section(self):
        """plot material"""
        self._frame.sections()    
    #
    def load(self):
        """ """
        return self._load
    #
    #
#
#
#
