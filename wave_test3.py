#!/usr/bin/env python
# coding: utf-8
#
#import matplotlib.pyplot as plt

from steelpy import UFOmodel
from steelpy import Metocean
from steelpy import Trave2D, Trave3D
from steelpy import Units
#
#
units = Units()
#
# ----------------------------------------------------
#
f2umodel = UFOmodel(name="wave_test3")
#
mesh = f2umodel.mesh()
#
# ----------------------------------------------------
# Material input
# ----------------------------------------------------
#
material = mesh.material()
matlinear = material.linear()
matlinear[1] = [345.0 * units.MPa]
print(material)
#
#
# ----------------------------------------------------
# Section Input
# ----------------------------------------------------
#
section = mesh.section()
section[2] = ['Tubular', 500 * units.mm, 25 * units.mm]
print(section)
#
# ----------------------------------------------------
# Mesh Model
# ----------------------------------------------------
#
#
# Node input
nodes = mesh.node()
nodes[10] = [0 * units.m, -100* units.m, 0 * units.m]
nodes[20] = [0 * units.m, -75 * units.m, 0 * units.m]
nodes[30] = [0 * units.m, -50 * units.m, 0 * units.m]
nodes[40] = [0 * units.m, -25 * units.m, 0 * units.m]
nodes[50] = [0 * units.m,  0 * units.m, 0 * units.m]
nodes[60] = [0 * units.m,  25* units.m, 0 * units.m]
nodes[70] = [0 * units.m, 50* units.m, 0 * units.m]
print(nodes)
#
#
# boundary Input
boundary = mesh.boundary()
#
supports = boundary.nodes()
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
elements = mesh.element()
#
beams = elements.beam()
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
meto = Metocean(name="metocean_1")
#
# ----------------------------------------------------
# hydrodynamic parameters input
# ----------------------------------------------------
#
proph = meto.property()
#
cdcm = proph.CdCm()
# ['coefficients', Cd, Cm]
#cdcm['cdcm_1'] = ['cdcm_default', 0.70, 2.0]
#cdcm['cdcm_1'].elements = [12, 23]
cdcm['cdcm_1'] = [[ 10 * units.m, 0.70, 2.0],
                  [ 00 * units.m, 0.70, 2.0],
                  [-10 * units.m, 0.70, 2.0],
                  [-50 * units.m, 0.70, 2.0],
                  [-100 * units.m, 0.70, 2.0]]
#
#
#
mg = proph.MG()
#       [title, density, thickness]
#mg['MG_1'] = ['MG_1' , meto.rho_w * 1.35]
# profile [elevation, thickness]
mg['MG_1'] = [[ 10 * units.m, 60 * units.mm],
              [ 00 * units.m, 60 * units.mm],
              [-10 * units.m, 30 * units.mm],
              [-50 * units.m, 10 * units.mm],
              [-100 * units.m, 0 * units.mm]]
mg['MG_1'].rho = meto.rho_w * 1.35
#
# Wave kinematics factor
wkrf = proph.WKF()
#  [text, type (constant/profile), factor]
#wkrf['wkf1'] = [ 'wkf1', 0.95]
wkrf['wkf1'] = [[ 10 * units.m, 0.95],
                [ 00 * units.m, 0.95],
                [-10 * units.m, 0.95],
                [-50 * units.m, 0.95],
                [-100 * units.m, 0.95]]
#
#
#
cbf = proph.CBF()
#  [text, type (constant/profile), factor]
#cbf['cbf1'] = [ 'cbf1', 0.85]
cbf['cbf1'] = [[ 10 * units.m, 0.85],
               [ 00 * units.m, 0.85],
               [-10 * units.m, 0.85],
               [-50 * units.m, 0.85],
               [-100 * units.m, 0.85]]
#
#
# ----------------------------------------------------
# Current
#
#
MetCriteria = meto.criteria()
MetCriteria[1] = 'test_1'
#
current = MetCriteria[1].current()
#
# [title, profile (linear/exponential/user), velocity_top, velocity_bottom]
#current['curr_1'] = 'current_1' # [1.5 * units.m/units.sec, 0.10 * units.m/units.sec]
# profile [elevation, velocity]
cvel = 1.54 * units.m/units.sec
current['curr_1'] = [[  5 * units.m, 1.0 * cvel],
                     [-10 * units.m, 1.0 * cvel],
                     [-30 * units.m, 0.70 * cvel],
                     [-60 * units.m, 0.50 * cvel],
                     [-80 * units.m, 0.20 * cvel]]
#
#
# ----------------------------------------------------
# Wind
#wind = meto.wind()
#wind['wind_1'] =  1.293 * units.kg / units.m**3
#
#
# ----------------------------------------------------
# Regular wave [Stokes/Fourier/Cnoidal]
# ----------------------------------------------------
#
wave = MetCriteria[1].wave()
#
# -------------------------------------------------------------------
#
# Stokes
regwave = wave.regular()
#
#
#waveCnoidal = regwave.Cnoidal(H=0.30, T=6.4, d=1)
#print(waveCnoidal.Lw)
#surface = waveCnoidal.surface()
#surface.plot()
#kin = waveCnoidal.kinematics()
#kin.plot()
#
#waveStoke = regwave.Stokes(H=15, T=12.0, d=100)
#print(waveStoke.Lw)
#surface = waveStoke.surface()
#surface.plot()
#kin = waveStoke.kinematics()
#kin.plot()
#wave = regwave.Fourier()
#wave.mean_current.Euler = 0.31 * units.m/units.sec
#wave.infinite_water_depth
# [H, T, d]
#wave['100yrs'] = {'Hw':15.0 * units.m, 'Tw':12.0 * units.sec, 'd':100*units.m}
# 
# [Hw, Tw, d, wave_type(Fourier/Stokes)]
regwave['100yrs_1'] = [15*units.m, 12.0*units.sec, 100*units.m, 'Stokes']
#Ls = regwave['100yrs'].L
#print(f'Wave length = {Ls: 1.4e} m')
#surface = regwave['100yrs'].surface(surface_points=18)
#surface.plot(phase=True)
#print(surface)
#
regwave['100yrs_2'] = [15*units.m, 12.0*units.sec, 100*units.m]
#regwave['100yrs_3'] = [0.30 * units.m, 6.4 * units.sec, 1*units.m, 'Cnoidal']
#regwave['100yrs_4'] = [0.30 * units.m, 6.4 * units.sec, 1*units.m]
#
#kinematic = regwave['100yrs'].kinematics()
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
metcond = MetCriteria[1].condition()
# [title]
metcond[10] = 'storm_0deg'    #['storm_0deg', 'MG_1', False, 0.85]
# wave =[wave_id, Direction(deg), Kinematics, crest_elevation]
metcond[10].wave = ['100yrs_1', 10.0, 'wkf1']
# current [current_id,  Direction(deg), Blockage, Stretching]
metcond[10].current = ['curr_1', 20.0, 'cbf1', True]
# wind [wind_id, Direction(deg)]
#metload[1].wind = [1, 30.0]
#
# [ marine_growth, CdCm, Flooding, conductor_shielding, element_refinament, airCdCm]
metcond[10].properties = ['MG_1', 'cdcm_1']
#
# [Design load, Buoyancy(False/True), 
#  Criterion (local/global), Scale factor, title]
metcond[10].parameters = ['max_BS', None, 'local']
#
#
# ----------------------------------------------------
#
#wave.solve(surface_points=36,
#           depth_points=100)
#
MetCriteria.solve(surface_points=36,
                  depth_points=100)
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
loadCase = load
hydro = loadCase.metocean(metcond)
#
#hydro['MET_10'] = metcond[10]
#
#
# ----------------------------------------------------
# Meshing input
# ----------------------------------------------------
#
mesh.build()
#
print("Load")
loadm = mesh.load()
basicLoad = loadm.basic()
print(basicLoad)
#
# ----------------------------------------------------
# Structural Analysis
# ----------------------------------------------------
#
frame = Trave2D(mesh=mesh)
frame.static()
results = frame.results()
print(results)
#
print('-->')