#!/usr/bin/env python
# coding: utf-8
#
#import numpy as np

from steelpy import Metocean
from steelpy import Units


meto = Metocean()
units = Units()
#
# Regular wave
#
waveReg = meto.regular_wave()
#
waveReg['100yrs'] = {'Hw':15.0 * units.m, 'Tw':12.0 * units.sec, 'd':100*units.m}
#waveReg['100yrs'] = [15.0 * units.m, 12.0 * units.sec, 100*units.m]
surface = waveReg['100yrs'].surface(surface_points=18)
surface.printout()
#surface.plot()
#
kinematic = waveReg['100yrs'].kinematics(depth_points=10)
print(kinematic['u'].max())
kinematic.printout()
kinematic.plot()
kinematic.plot_vectorfield()
#
#
# -------------------------------------------------------------------
#
# Stokes
#
stokes5 = waveReg.Stokes()
#stokes5.mean_current.Euler = 0.31 * units.m/units.sec
#stokes5.infinite_water_depth
# [H, T, d]
stokes5['100yrs'] = {'Hw':15.0 * units.m, 'Tw':12.0 * units.sec, 'd':100*units.m}
Ls = stokes5['100yrs'].L
surface = stokes5['100yrs'].surface(surface_points=18)
surface.plot()
surface.printout()
#
kinematic = stokes5['100yrs'].kinematics(depth_points=10)
kinematic.printout()
kinematic.plot()
#kinematic.plot_vectorfield()
#pressure = kinematic.pressure
#
# Fourier
#
fourier = waveReg.Fourier()
#fourier.mean_current.Euler = 0.31 * units.m/units.sec
# [H, T, d]
fourier['100yrs'] = {'Hw':15.0 * units.m, 'Tw':12.0 * units.sec, 'd':100*units.m}
Lf = fourier['100yrs'].L
surface = fourier['100yrs'].surface(surface_points=18)
#surface.plot()
surface.printout()
#
kinematic = fourier['100yrs'].kinematics(depth_points=10)
#kinematic.plot()
#kinematic.printout()
#
#
#
print('-->')