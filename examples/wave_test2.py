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
# ===========================
# Irregular wave
#
wave = meto.irregular_wave()
sptm = wave.spectrum()
# set spectrum 
jonswap = sptm.Jonswap(Hs=12.0 * units.m, Tp=14.0 * units.sec)
jonswap.plot(harmonic=1)
#tsim = 1*units.hour
#print(f'time = {tsim.convert("second").value}')
surface = sptm.simulation(t=120*units.sec,  #1800
                          dt=0.50*units.sec,
                          d=100*units.m)
                          #wave_order=2)
surface.plot()
#
kinematic = wave.kinematics(surface)
kinematic.plot()
#
Dp = 1.0 * units.m
#BS, OTM = wave.pile_response(D, kinematic)
BS, OTM = meto.pile_response(D=Dp, L=100*units.m,
                             kinematic=kinematic)
#
print('-->')
#
#
#
# -------------------------------------------------------------------
#
# Cnoidal
#
cnoidal = waveReg.Cnoidal()
cnoidal.mean_current.Euler = 0.31 * units.m/units.sec
# [H, T, d]
cnoidal['test'] = [0.30 * units.m, 6.4 * units.sec, 1*units.m]
surface = cnoidal['test'].surface()
#surface.plot()
#
# Fourier
#
fourier = waveReg.Fourier()
fourier.mean_current.Euler = 0.31 * units.m/units.sec
#fourier.infinite_water_depth
#fourier.mean_current.Stokes = 0.31 * units.m/units.sec
# [H, T, d]
fourier['test'] = [0.30 * units.m, 6.4 * units.sec, 1*units.m]
surface = fourier['test'].surface()
surface.printout()
#surface.plot()
kinematic = fourier['test'].kinematics()
kinematic.printout()
#kinematic.plot()
#
# -------------------------------------------------------------------
#
# Stokes
#
stokes5 = waveReg.Stokes()
# [H, T, d]
stokes5['100yrs'] = [15.0 * units.m, 12.0 * units.sec, 100*units.m]
#L = stokes5['100yrs'].L
surface = stokes5['100yrs'].surface()
#surface.plot()
#surface.printout()
kinematic = stokes5['100yrs'].kinematics(surface_points=18)
#kinematic.plot()
#kinematic.printout()
#pressure = kinematic.pressure
#
# Fourier
#
fourier = waveReg.Fourier()
# [H, T, d]
fourier['100yrs'] = [15.0 * units.m, 12.0 * units.sec, 100*units.m]
surface = fourier['100yrs'].surface()
kinematic = fourier['100yrs'].kinematics(surface_points=18)
#kinematic.plot()
#

#
print('-->')