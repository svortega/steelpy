#
# Copyright (c) 2009-2023 fem2ufo
#
from __future__ import annotations
# Python stdlib imports
import math
from enum import Enum
from typing import NamedTuple, Tuple, List


# package imports
#import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.cm as cm
import numpy as np

from steelpy.process.units.main import Units
from steelpy.process.geometry.euclid import Point3

#
class PlotItems(Enum):
    bc = 'bc'
    bc_id = 'bc_id'
    beam_index = 'beam_index'
    deformed = 'deformed'
    forces = 'forces'
    global_axes = 'global_axes'
    moments = 'moments'
    node_uids = 'node_uids'
    nodes = 'nodes'
    undeformed = 'undeformed'

    @classmethod
    def to_list(cls):
        return [e.value for e in cls]
#
#
class C:
    BOX_BC = 'black'
    TXT_BC = 'white'

    BOX_NODE_UID = 'orange'
    TXT_NODE_UID = 'black'

    BOX_BEAM_IDX = 'navy'
    TXT_BEAM_IDX = 'white'

    UNDEFORMED = 'grey'
    DEFORMED = 'red'

    FORCE = 'steelblue'
    MOMENT = 'purple'

    GLOBAL_SYS = 'blue'
    LOCAL_SYS = 'green'

    MASS = 'maroon'
#
#
class PlotModel:
    __slots__ = ['_mesh', '_load', '_concept']
    
    def __init__(self, mesh):
        """
        """
        self._mesh = mesh

    #
    #
    #
    #def __call__(self, verbosity:bool=False, show=True):
    #    """ """
    #    self.concept()
    #    print('here')
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
    def mesh(self, verbosity:bool=False , show=True):
        """ """
        plot_num = 2
        nodes = self._mesh._nodes
        lims = nodes.get_maxmin()
        ax = init_3D_plot([lims[0][0], lims[1][0]],
                          [lims[0][1], lims[1][1]],
                          [lims[0][2], lims[1][2]])
        add_items_per_beam(self._mesh, ax, plot_num, verbosity)
        add_boundary_conditions(self._mesh, ax, plot_num, verbosity)
        if verbosity:
            add_global_axes(ax, plot_num)
        #
        #rfiles.set('plots', file_list)
        if show:
            plt.show()
        else:
            return ax
    #
    def basic_load(self, name, verbosity:bool=False):
        """ """
        ax = self.mesh(verbosity, show=False)
        load = self._mesh._load
        basic = load._basic[name]
        nodes = self._mesh._nodes
        #
        elements =  self._mesh.elements()
        beams = elements._elements
        add_mesh_load(ax, beams, nodes, basic)
        plt.show()
    #
    #
    def load_combination(self, combination):
        """ """
        pass    
#
#
def add_items_concept(m, ax, plot_num, 
                      verbosity:bool=False):
    """ """
    points = m.points
    #named_nodes = []
    #
    #ax.scatter(points._x, points._y, points._z,
    #            marker='o', color='red', s=50, alpha=0.5) # , zdir="y"
    #
    for key, point in points.items():
        ax.scatter(point.x, point.y, point.z, edgecolor='indigo',
                    marker='o', color='purple', s=50, alpha=0.5,
                    linewidth=1)
        #
        ax.text3D(point.x, point.y, point.z, str(key))
    #
    for key, beam in m.beam.items():
        n1,n2 = beam.connectivity
        #n1 = nodes[conn[0]]
        #n2 = nodes[conn[1]]
        #named_nodes = [n1, n2]
        #named_nodes.extend(conn[:])
        x = [n1.x, n2.x]
        y = [n1.y, n2.y]
        z = [n1.z, n2.z]
        #xyz = abm.get_all_points(beam_idx)
        #x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]

        # ----- Undeformed mesh -----
        ax.plot(x, y, z, **args_plot(color='teal', l_width=2))
#
#
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

        return np.min(zs)
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
    for lname, nloads in nload.items():
        node = nodes[lname]
        x, y, z = node[:3]
        for key, item in nloads._load.items():
            for nload in item:
                # nodal force
                force = [ item * scale for item in nload[:6]]
                if any(force[:3]):
                    Fx, Fy, Fz = force[:3]
                    ax.quiver(x, y, z, Fx, Fy, Fz, color='r')
                # nodal moment
                if any(force[3:6]):
                    Mx, My, Mz = [item*100 for item in force[3:6]]
                    plot_circle(ax, [Mx, My, Mz], [x, y, z])    
#
def beamlineload(ax, beams, nodes, lineload,
                 scale:float=1.0):
    """ """
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
                plot_lload(ax, x,y,z)
            # qz
            if any([force[2], force[5]]):
                q0, q1 = force[2], force[5]
                x, y, z = get_line_coord('z', n1, n2, q0, q1, L0, L1, normalized)
                plot_lload(ax, x,y,z)    
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
#
def get_vnorm(n1, n2):
    """get vector normalized"""
    A = Point3(*n1[:3])
    B = Point3(*n2[:3])
    vector = B - A
    return vector.normalized()    
#
#
def plot_circle(ax, moments, coord):
    """ """
    x, y, z = coord
    Mx, My, Mz = moments
    # Theta varies only between pi/2 and 3pi/2. to have a half-circle
    theta = np.linspace(np.pi/2., 4*np.pi/2.)
    # Mx
    if Mx:
        r = Mx*0.10
        xc = np.zeros_like(theta) + x # x=0
        yc = r*np.cos(theta) + y # y - y0 = r*cos(theta)
        zc = r*np.sin(theta) + z # z - z0 = r*sin(theta)
        ax.plot(xc, yc, zc)
        # arrow heads
        get_arrow(ax, xc, yc, zc, color="b")
    # My
    if My:
        r = My*0.10    
        yc = np.zeros_like(theta) + y
        zc = r*np.cos(theta) + z
        xc = r*np.sin(theta) + x
        ax.plot(xc, yc, zc, color="b")
        # arrow heads
        get_arrow(ax, xc, yc, zc)    
    # Mz
    if Mz:
        r = Mz*0.10    
        zc = np.zeros_like(theta) + z
        yc = r*np.cos(theta) + y 
        xc = r*np.sin(theta) + x
        ax.plot(xc, yc, zc, color="b")
        # arrow heads
        get_arrow(ax, xc, yc, zc)
    #
#
def get_arrow(ax, xc, yc, zc):
    """ """
    a = Arrow3D([xc[-2], xc[-1]], [yc[-2], yc[-1]], 
                [zc[-2], zc[-1]], mutation_scale=10, 
                lw=2, arrowstyle="->", color="b")
    ax.add_artist(a)     
#
#
#
def get_line_coord(axis:str, n1, n2, q0, q1, L0, L1, normalized):
    """ """
    x = [ 0, 0, 0, 0, 0]
    y = [ 0, 0, 0, 0, 0]
    z = [ 0, 0, 0, 0, 0]
    # corner 1
    x[0] = n1.x + normalized[ 0 ] * L0
    y[0] = n1.y + normalized[ 1 ] * L0
    z[0] = n1.z + normalized[ 2 ] * L0
    # corner 2
    x[1] = n2.x - normalized[ 0 ] * L1
    y[1] = n2.y - normalized[ 1 ] * L1
    z[1] = n2.z - normalized[ 2 ] * L1
    if axis.lower() == 'z':
        # corner 3
        x[2] = x[1]
        y[2] = y[1]
        z[2] = z[1] + q1
        # corner 4
        x[3] = x[0]
        y[3] = y[0]
        z[3] = z[0] + q0
    else:
        # corner 3
        x[2] = x[1] + q1
        y[2] = y[1] 
        z[2] = z[1] 
        # corner 4
        x[3] = x[0] + q0
        y[3] = y[0]   
        z[3] = z[0]       
    # corner 1
    x[4] = x[0]
    y[4] = y[0]
    z[4] = z[0]
    return x, y, z
#
def plot_lload(ax, x, y, z):
    """ """
    verts = [list(zip( x, y, z))]
    tri = Poly3DCollection(verts)
    tri.set_color('palegreen')
    tri.set_alpha(0.50)
    tri.set_edgecolor('lightgrey')
    ax.add_collection3d(tri)    
#
#
def add_items_per_beam(m, ax, plot_num, 
                       verbosity:bool=False):
    """ """
    nodes = m._nodes
    beams = m._elements.beams()
    #
    named_nodes = []
    for key, beam in beams.items():
        conn = beam.connectivity
        n1 = nodes[conn[0]]
        n2 = nodes[conn[1]]
        #named_nodes = [n1, n2]
        named_nodes.extend(conn[:])
        x = [n1.x, n2.x]
        y = [n1.y, n2.y]
        z = [n1.z, n2.z]
        #xyz = abm.get_all_points(beam_idx)
        #x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]

        # ----- Undeformed mesh -----
        #if PlotItems.undeformed.value in to_show:
        ax.plot(x, y, z, **args_plot(color='k', 
                                     marker = '.', l_width=1.5))

        # ----- Deformed mesh -----
        #if PlotItems.deformed.value in to_show:
        #    d = m.results.get('tensors').get('comp:U')
        #    scale = ps.get('scale_deformation', 1)
        #    xd = x + scale*abm.gbv(d['ux'], beam_idx)
        #    yd = y + scale*abm.gbv(d['uy'], beam_idx)
        #    zd = z + scale*abm.gbv(d['uz'], beam_idx)
        #    ax.plot(xd, yd, zd, **args_plot(m, C.DEFORMED, marker=marker))
        #
        # ----- Forces -----
        #if PlotItems.forces.value in to_show:
        #    d = m.results.get('tensors').get('comp:F')
        #    scale = ps.get('scale_forces', 1)
        #    Fx = scale*abm.gbv(d['Fx'], beam_idx)
        #    Fy = scale*abm.gbv(d['Fy'], beam_idx)
        #    Fz = scale*abm.gbv(d['Fz'], beam_idx)
        #    if ps.get('deform_loads', True):
        #        ax.quiver(xd, yd, zd, Fx, Fy, Fz, color=C.FORCE)
        #    else:
        #        ax.quiver(x, y, z, Fx, Fy, Fz, color=C.FORCE)
        #
        # ----- Moments -----
        #if PlotItems.moments.value in to_show:
        #    d = m.results.get('tensors').get('comp:F')
        #    scale = ps.get('scale_moments', 1)
        #    Fx = scale*abm.gbv(d['Mx'], beam_idx)
        #    Fy = scale*abm.gbv(d['My'], beam_idx)
        #    Fz = scale*abm.gbv(d['Mz'], beam_idx)
        #    if ps.get('deform_loads', True):
        #        ax.quiver(xd, yd, zd, Fx, Fy, Fz, color=C.MOMENT)
        #    else:
        #        ax.quiver(x, y, z, Fx, Fy, Fz, color=C.MOMENT)
        #
        # ----- Beam index -----
        if verbosity:
            coord = [0,0,0]
            mid_length = beam.length * 0.50
            normalized = get_vnorm(n1, n2)
            #
            coord[0] = n1.x + normalized[0] * mid_length
            coord[1] = n1.y + normalized[1] * mid_length
            coord[2] = n1.z + normalized[2] * mid_length
            #
            ax.text(*coord, str(key), **args_text(f_size=8, c_txt='red', c_box='b'))
    #
    # ----- Named nodes -----
    if verbosity:
        named_nodes = set(named_nodes)
        for node_id in named_nodes:
            node = nodes[node_id]
            ax.text(*node[:3], node.number, **args_text(f_size=8, c_txt='g', c_box='b'))
#
#
def args_plot(color, l_width:int=2, marker=None, m_size:int=5):
    args = {'linewidth': l_width,
            #'markersize': m_size,
            'color': color,
            #'edgecolor': color,
    }
    if marker is not None:
        args['marker'] = marker
        args['markersize'] = m_size
    return args
#
def args_scatter(color, marker, l_width:int=4):
    args = {'linewidth': l_width, 
            'color': color,
            'marker':marker,}
    return args

#
def args_text(f_size, c_txt, c_box):
    args = {'fontsize': f_size,
            'color': c_txt,
            #'bbox': dict(facecolor=c_box, alpha=0.5),
            'horizontalalignment': 'center',
            'verticalalignment': 'bottom',}
    return args
#
#
#
def init_3D_plot(x_lims, y_lims, z_lims):
    """
    Inititalize the 3D plot

    Args:
        :x_lims: (tuple) min and max x-value
        :y_lims: (tuple) min and max y-value
        :z_lims: (tuple) min and max z-value
    """

    plt.figure(figsize=(10, 10))
    ax = plt.axes(projection='3d')

    # Avoid setting same min and max value by adding diff
    diff = (-1e-6, 1e-6)
    ax.set_xlim(*(x+d for x, d in zip(x_lims, diff)))
    ax.set_ylim(*(y+d for y, d in zip(y_lims, diff)))
    ax.set_zlim(*(z+d for z, d in zip(z_lims, diff)))

    set_equal_aspect_3D(ax)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    return ax
#
def set_equal_aspect_3D(ax):
    """
    Set aspect ratio of plot correctly

    Args:
        :ax: (obj) axis object
    """

    # See https://stackoverflow.com/a/19248731
    # ax.set_aspect('equal') --> raises a NotImplementedError
    # See https://github.com/matplotlib/matplotlib/issues/1077/

    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:, 1] - extents[:, 0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
#
#
#
class GlobalSystem:
    Origin = np.array([0, 0, 0])
    X = np.array([1, 0, 0])
    Y = np.array([0, 1, 0])
    Z = np.array([0, 0, 1])
#
def _coordinate_system(plot, origin, axes, axes_names, color, scale=1):
    axes = [scale*np.array(axis) for axis in axes]
    x_axis, y_axis, z_axis = axes

    for axis, axis_name in zip(axes, axes_names):
        x, y, z = origin
        u, v, w = axis
        plot.scatter(x, y, z)
        plot.quiver(x, y, z, u, v, w, length=1)
        plot.text(x+u, y+v, z+w, axis_name)

    # Plot xy-plane
    p1 = np.array(origin)
    p2 = np.array(origin + y_axis)
    p3 = np.array(origin + z_axis)
    p4 = np.array(origin + y_axis + z_axis)
    points = np.array([p1, p2, p3, p4])
    xx = points[:, 0].reshape(2, 2)
    yy = points[:, 1].reshape(2, 2)
    z = points[:, 2].reshape(2, 2)
    plot.plot_surface(xx, yy, z, alpha=0.4, color=color)
#
def add_global_axes(ax, plot_num):
    #to_show = m.get('post_proc').get('plot')[plot_num]
    #to_show = 'global_axes'
    #if 'global_axes' in to_show:
    #if to_show:
    orig = GlobalSystem.Origin
    X = GlobalSystem.X
    Y = GlobalSystem.Y
    Z = GlobalSystem.Z
    ax = _coordinate_system(ax, orig, (X, Y, Z), ('X', 'Y', 'Z'), color=C.GLOBAL_SYS)
#
#
def add_boundary_conditions(m, ax, plot_num, verbosity:bool=False):
    """ """
    mbc = m.boundaries()
    supports = mbc.supports()
    nodes = m.nodes()
    #
    for key, bc in supports.items():
        node = nodes[key]
        xyz = node[:3]
        bcn = bc[:6]
        if all(bcn):
            if sum(bcn) == 0:
                continue
            marker ='s'
        else:
            marker ='^'
        #
        ax.scatter(*xyz, **args_scatter(color='m', marker=marker))
        if verbosity:
            bc_id = [int(item) for item in bcn]
            ax.text(*xyz, f'{bc_id}', **args_text(f_size=8, c_txt='m', c_box='b'))
#
#



