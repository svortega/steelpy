from steelpy import  Beam
from steelpy import Units

# set units
units = Units()
#
# beam class
beam = Beam(beam_length=5*units.m)
#
# set section
beam.section = 'Rectangle'
beam.section.d = 350 * units.mm
beam.section.w = 350 * units.mm
#
# set supports : pinned/fixed/free/guided
beam.support[1] = "pinned"
beam.support[2] = "pinned"
#
# set loading
beam.load[1].point = {'L1':2.5*units.m, 'fy':100*units.kN}
beam.load[2].moment = {'L1':3.0*units.m, 'my':200 * units.kN*units.m}
beam.load[3].line = {'qy1':50*units.kN/units.m, 'qy2':100*units.kN/units.m,
                     'L1':1*units.m, 'L2':1*units.m}
#
# Get support reactions : R/M/theta/w
R1 = beam.support[1].in_plane.R.value
R2 = beam.support[2].in_plane.R.value
#
# plot beam results
beam.bending_moment.plot("in_plane")
beam.shear.plot("in_plane")
beam.deflection.plot("in_plane")
beam.slope.plot("in_plane")
#
