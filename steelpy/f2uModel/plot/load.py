#
# Copyright (c) 2009-2023 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
#
#
#
# package imports
#
from steelpy.f2uModel.plot.frame import PlotBasic
from steelpy.f2uModel.plot.plot3D import (get_vnorm,
                                          plot_circle,
                                          get_line_coord,
                                          plot_lload)
#
#import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#import matplotlib.cm as cm
#from matplotlib.widgets import RadioButtons
#
#
#
#
# --------------------
# Concept
# --------------------
#
#
def add_concept_load(ax, elements, nodes, basic):
    """ """
    # nodal load
    #nodeload(ax, nodes, 
    #         nload=basic.point_node)
    # beam point load
    try:
        beampointload(ax, elements, nodes,
                      beamload=basic.beams)
    except AttributeError:
        pass
    #
    #
    # beam line load
    beamlineload(ax, elements, nodes, 
                 beamload=basic.beams)
#
#
#
# --------------------
# Mesh
# --------------------
#
#
class PlotLoad(PlotBasic):
    __slots__ = ['_mesh', 'figsize']
    
    def __init__(self, mesh, figsize:tuple = (10, 10)):
        """
        """
        super().__init__(mesh, figsize)
    #
    def basic(self):
        """ """
        basic = self._mesh._load._basic
        bname = [key for key in basic.keys()]
        #
        #
        fig, ax = self.init_frame()
        #
        #
        self.get_load(ax, load=basic)
        #
        #axcolor = 'lightgoldenrodyellow'
        #rax = plt.axes([0.05, 0.7, 0.15, 0.15], facecolor=axcolor)
        #lines_by_label = {'Nodes': s2, 'Elements': s3, 'Supports': s4}
        #
        plt.show()
        #else:
        return ax
    #
    def get_load(self, ax, load):
        """ """
        nodes = self._mesh.nodes()
        beams = self._mesh._elements.beams()
        #
        for key, item in load.items():
            #key, item
            # nodal load
            nodeload(ax, nodes, nload=item._node)
            # beam point load
            #beampointload(ax, beams, nodes, 
            #              pointload=basic._beam._load._point)
            # beam line load
            beamlineload(ax, beams, nodes, 
                         lineload=item._beam._load._line)            
#
#
def add_mesh_load(ax, beams, nodes, basic):
    """ """
    # nodal load
    nodeload(ax, nodes, 
             nload=basic._node)
    # beam point load
    beampointload(ax, beams, nodes, 
                  pointload=basic._beam._load._point)
    # beam line load
    beamlineload(ax, beams, nodes, 
                 lineload=basic._beam._load._line)
#
def nodeload(ax, nodes, nload):
    """ """
    scale = 0.000001
    
    for key, nloads in nload._node._load.items():
        node = nodes[key]
        x, y, z = node[:3]
        for nload in  nloads:
            # nodal force
            force = [ item * scale for item in nload[:6]]
            if any(force[:3]):
                Fx, Fy, Fz = force[:3]
                ax.quiver(x, y, z, Fx, Fy, Fz, color='r')
            # nodal moment
            if any(force[3:6]):
                Mx, My, Mz = [item*20 for item in force[3:6]]
                plot_circle(ax, [Mx, My, Mz], [x, y, z])            
        
    #for lname, nloading in nload.items():
    #    node = nodes[lname]
    #    x, y, z = node[:3]
    #    1 / 0
    #    for key, nloading in nloads._load.items():
    #        for nload in nloads:
    #            # nodal force
    #            force = [ item * scale for item in nload[:6]]
    #            if any(force[:3]):
    #                Fx, Fy, Fz = force[:3]
    #                ax.quiver(x, y, z, Fx, Fy, Fz, color='r')
    #            # nodal moment
    #            if any(force[3:6]):
    #                Mx, My, Mz = [item*20 for item in force[3:6]]
    #                plot_circle(ax, [Mx, My, Mz], [x, y, z])
    #
    #
#
def beamlineload(ax, beams, nodes, lineload,
                 scale:float=1.0):
    """ """
    scale = 0.00001
    for bname, line in lineload.items():
        try:
            beam = beams[bname]
        except IndexError:
            continue
        n1, n2 = beam.nodes
        normalized = get_vnorm(n1, n2)
        #
        for lload in line:
            force = [item * scale for item in lload[:6]]
            L0 = lload.L0
            L1 = lload.L1
            # qy
            if any([force[1], force[4]]):
                q0, q1 = force[1], force[4]
                x, y, z = get_line_coord('y', n1, n2, q0, q1, L0, L1, normalized)
                ax = plot_lload(ax, x,y,z)
            # qz
            if any([force[2], force[5]]):
                q0, q1 = force[2], force[5]
                x, y, z = get_line_coord('z', n1, n2, q0, q1, L0, L1, normalized)
                ax = plot_lload(ax, x,y,z)    
#
def beampointload(ax, beams, nodes, pointload,
                  scale:float=1.0):
    """ """
    for bname, point in pointload.items():
        try:
            beam = beams[bname]
        except IndexError:
            continue        
        n1, n2 = beam.nodes
        normalized = get_vnorm(n1, n2)
        #
        for pload in point:
            force = [ item * scale for item in pload[:6]]
            x = n1.x + normalized[ 0 ] * pload.L0
            y = n1.y + normalized[ 1 ] * pload.L0
            z = n1.z + normalized[ 2 ] * pload.L0
            # point force
            if any(force[:3]):
                Fx, Fy, Fz = force[:3]
                ax.quiver(x, y, z, Fx, Fy, Fz, color='r')
            # point moment
            if any(force[3:6]):
                Mx, My, Mz = [item*100 for item in force[3:6]]
                ax.quiver(x, y, z, Mx, My, Mz, color='b')    
#