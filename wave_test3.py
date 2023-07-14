#!/usr/bin/env python
# coding: utf-8
#
import matplotlib.pyplot as plt

from steelpy import f2uModel
from steelpy import Metocean
from steelpy import Trave2D
from steelpy import Units


units = Units()
#
# ----------------------------------------------------
#
f2umodel = f2uModel(component="wave_test3")
#
#
# ----------------------------------------------------
# Material input
# ----------------------------------------------------
#
material = f2umodel.materials()
matlinear = material.elastic()
matlinear[1] = [345.0 * units.MPa]
print(material)
#
#
# ----------------------------------------------------
# Section Input
# ----------------------------------------------------
#
section = f2umodel.sections()
section[2] = ['Tubular', 500 * units.mm, 25 * units.mm]
print(section)
#
# ----------------------------------------------------
# Mesh Model
# ----------------------------------------------------
#
mesh = f2umodel.mesh()
#
# Node input
nodes = mesh.nodes()
nodes[10] = [0, -100, 0 ]
nodes[20] = [0, -75 , 0 ]
nodes[30] = [0, -50 , 0 ]
nodes[40] = [0, -25 , 0 ]
nodes[50] = [0,  0, 0 ]
nodes[60] = [0,  25, 0 ]
nodes[70] = [0, 50, 0 ]
print(nodes)
#
#
# boundary Input
boundary = mesh.boundaries()
#
supports = boundary.supports()
supports[10] = 'fixed'
supports[20] = 'pinned'
supports[30] = 'pinned'
supports[40] = 'pinned'
#supports[50] = 'pinned'
supports[60] = 'pinned'
supports[70] = 'fixed'
print(boundary)
#
#
# Element input
#
#
elements = mesh.elements()
#
beams = elements.beams()
# beam[number] = [material, section, node1, node2, roll_angle]
beams[12] = [10, 20, 1, 2, 0]
beams[23] = [20, 30, 1, 2, 0]
beams[34] = [30, 40, 1, 2, 0]
beams[45] = [40, 50, 1, 2, 0]
beams[56] = [50, 60, 1, 2, 0]
beams[67] = [60, 70, 1, 2, 0]
#
print(elements)
#
#
# ----------------------------------------------------
# Metocean 
# ----------------------------------------------------
#
meto = Metocean()
#
# ----------------------------------------------------
# hydrodynamic parameters input
# ----------------------------------------------------
#
#
#
mg = meto.marine_growth()
mg['MG_1'] = "profile"
mg['MG_1'].level[1] = [ 10 * units.m, 60 * units.mm]
mg['MG_1'].level[2] = [ 00 * units.m, 60 * units.mm]
mg['MG_1'].level[3] = [-10 * units.m, 30 * units.mm]
mg['MG_1'].level[4] = [-50 * units.m, 10 * units.mm]
mg['MG_1'].level[5] = [-100 * units.m, 0 * units.mm]
mg['MG_1'].set_default()
#
mg['MG_2'] = "constant"
mg['MG_2'].thickness = 100 * units.mm
mg['MG_2'].elements = 12, 23
#
#
#cdcm = hydro.CdCm
#
#
# ----------------------------------------------------
# Regular wave [Stokes/Fourier/Cnoidal]
# ----------------------------------------------------
#
wave = meto.regular_wave()
#
# -------------------------------------------------------------------
#
# Stokes
#
#wave = waveReg.Stokes()
#wave = waveReg.Fourier()
#wave.mean_current.Euler = 0.31 * units.m/units.sec
#wave.infinite_water_depth
# [H, T, d]
#wave['100yrs'] = {'Hw':15.0 * units.m, 'Tw':12.0 * units.sec, 'd':100*units.m}
# 
# [Hw, Tw, d, wave_type(Fourier/Stokes)]
wave['100yrs'] = [15 * units.m, 12.0 * units.sec, 100*units.m, 'Stokes']
Ls = wave['100yrs'].L
print(f'Wave length = {Ls: 1.4e} m')
#surface = wave['100yrs'].surface(surface_points=18)
#surface.plot(phase=True)
#print(surface)
#
#kinematic = wave['100yrs'].kinematics()
#print(kinematic)
#kinematic.plot()
#kinematic.plot_vectorfield()
#pressure = kinematic.pressure
#phi = kinematic.phi
#
#BS, OTM = meto.pile_response(D=1*units.m, L=100*units.m, 
#                             kinematic=kinematic)
# get maximum
#BSgrp = BS.groupby(['length'])['BS'].sum()
#BSgrp.plot(kind="line",
#           xlabel='Wave length [m]', ylabel='BS [N]')
#plt.show()
#
#meto.get_load(mesh=mesh, kinematic=kinematic,
#              condition=2)
#
metload = meto.load()
#
metload['sea_1'] = 'storm_0deg'
# wave =[wave_name, Direction(deg), Kinematics, title]
metload['sea_1'].wave = ['100yrs', 0.0, 0.95]
# current [current_name,  Direction(deg), Blockage, Stretching, title]
metload['sea_1'].current = ['current_1', 0.0, 0.85, True]
# wind [wind_name, Direction(deg), title]
metload['sea_1'].wind = ['wind_1', 0.0]
#
#
#
#
# ----------------------------------------------------
# Basic Load
#
# loading
load = mesh.load()
#
# ----------------------------------------------------
# Basic Load
basic = load.basic()
#
#
basic[10] = 'wave load'
basic[10].wave(wave_load = metload['sea_1'],
               design_load = 'max_BS')
#
# ----------------------------------------------------
# Meshing input
# ----------------------------------------------------
#
f2umodel.build()
#
#
# ----------------------------------------------------
# Structural Analysis
# ----------------------------------------------------
#
frame = Trave2D()
frame.mesh = mesh
results = frame.run_static()
results.print()
#
print('-->')