from steelpy import Units
from steelpy import Vibration
import numpy as np

#
vib = Vibration()
units = Units()
#
#
mdof = vib.mdof()
#
# Ciria Example
#
m1 = 1522 # units.kg
m2 = 635  # units.kg
M = [[m1, 0], [0, m2]]
#
K = [[531 * 10**3, -166 * 10**3],
     [-166 * 10**3, 66 * 10**3]]
#
mdof['ciria'] = [M, K]
mdof['ciria'].free_response()
#
#
t = np.linspace(0,20,1001)
f = np.sin(2.5*t)
#
th = [[t, f], None]
mdof['ciria'].forced_response(th=th, dur=1000)
#
print('--> end')