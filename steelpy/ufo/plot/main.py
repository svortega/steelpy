#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
#

# package imports
from steelpy.ufo.plot.frame import PlotFrame
from steelpy.ufo.plot.load import PlotLoad
#from steelpy.ufo.plot.frame import (init_frame_plot,
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
# Frame
#
class PlotMain:
    __slots__ = ['_cls', '_frame', '_load']
    
    def __init__(self, cls, figsize:tuple)->None:
        """
        """
        self._cls = cls
        self._frame = PlotFrame(figsize=figsize)
    #
    #
    def material(self):
        """plot material"""
        nodes = self._cls._nodes
        lims = nodes.get_maxmin()
        materials = self._cls._materials
        elements = self._cls._elements
        self._frame.materials(materials, elements, lims)
    #
    #
    def section(self):
        """plot material"""
        nodes = self._cls._nodes
        lims = nodes.get_maxmin()
        sections = self._cls._sections
        elements = self._cls._elements
        self._frame.sections(sections, elements, lims)    
    #    
#
#
class PlotConcept(PlotMain):
    __slots__ = ['_cls', '_frame', '_load']
    
    def __init__(self, cls, figsize:tuple)->None:
        """
        """
        super().__init__(cls, figsize)
    #
    #
    #
    def frame(self, show:bool =True):
        """ plot mesh element, nodes & boundary"""
        #print('--')
        elements = self._cls._elements
        nodes = self._cls._nodes
        #
        mbc = self._cls._boundaries.supports()
        supports = mbc._nodes        
        #
        ax = self._frame.frame(nodes, elements, supports)
        #
        if show:
            plt.show()
        else:
            return ax    
#
#
class PlotMesh(PlotMain):
    __slots__ = ['_cls', '_frame', '_load']
    
    def __init__(self, cls, figsize:tuple)->None:
        """
        """
        super().__init__(cls, figsize)
    #
    def frame(self, show:bool=True):
        """ plot mesh element, nodes & boundary"""
        #print('--')
        elements = self._cls._elements
        nodes = self._cls._nodes
        #lims = nodes.get_maxmin()
        #
        mbc = self._cls._boundaries
        supports = mbc.support()        
        #
        ax = self._frame.frame(nodes, elements, supports)
        #
        if show:
            plt.show()
        else:
            return ax
    #
    #
    #
#
#
#
#
#