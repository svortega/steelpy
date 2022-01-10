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
waveReg = meto.regular_waves
#
waveReg['100yrs'] = {'Hw':15.0 * units.m, 'Tw':12.0 * units.sec, 'd':100*units.m}
#waveReg['100yrs'] = [15.0 * units.m, 12.0 * units.sec, 100*units.m]
surface = waveReg['100yrs'].surface(surface_points=18)
surface.printout()
surface.plot()
#
kinematic = waveReg['100yrs'].kinematic(depth_points=10)
kinematic.printout()
kinematic.plot()
#
#
# -------------------------------------------------------------------
#
# Stokes
#
stokes5 = waveReg.Stokes()
#stokes5.infinite_water_depth
# [H, T, d]
stokes5['100yrs'] = [15.0 * units.m, 12.0 * units.sec, 100*units.m]
#L = stokes5['100yrs'].L
surface = stokes5['100yrs'].surface(surface_points=18)
#surface.plot()
#surface.printout()
#
kinematic = stokes5['100yrs'].kinematic(depth_points=10)
#kinematic.printout()
#kinematic.plot()
kinematic.plot_vectorfield()
#pressure = kinematic.pressure
#
# Fourier
#
fourier = waveReg.Fourier()
# [H, T, d]
fourier['100yrs'] = [15.0 * units.m, 12.0 * units.sec, 100*units.m]
#L = stokes5['100yrs'].L
surface = fourier['100yrs'].surface()
surface.plot()
surface.printout()
#
kinematic = fourier['100yrs'].kinematic(surface_points=18)
kinematic.plot()
kinematic.printout()
#
#
#
print('-->')