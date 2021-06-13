from steelpy import Beam
from steelpy import Units

# set units
units = Units()
#
# beam class
beam = Beam()
beam.length = 6*units.m
#
# set section
beam.section = ['Tubular', 500 * units.mm, 25 * units.mm]
#beam.section = ['Rectangle', 350 * units.mm, 200 * units.mm]
#beam.section.d = 350 * units.mm
#beam.section.w = 350 * units.mm
#
print(beam.section)
#
# set supports : pinned/fixed/free/guided/spring [k=F/x]
beam.support[1] = "pinned"
#beam.support[1] = ["spring", 200 * units.kN/units.m]
#
beam.support[2] = "pinned"
#beam.support[2] = ["spring", 200 * units.kN/units.m]
#
# set loading
#     point = [L1, Fy, Fz]
#beam.load[1].point = [3*units.m, -1000*units.N, 000*units.N]
#beam.load[1].point = [2*units.m, 18*units.kN, 18*units.kN]
#beam.load[1].point = {'L1':2.5*units.m, 'fy':100*units.kN}
beam.load[1].line = [-500*units.N/units.m, 0*units.N/units.m,
                     -1000*units.N/units.m, 0*units.N/units.m]
#     moment = [L1, My, Mz]
#beam.load[2].moment = [3.0*units.m, 200*units.kN*units.m, 300*units.kN*units.m]
#beam.load[2].moment = {'L1':3.0*units.m, 'my':200 * units.kN*units.m}
#     line = [qy1,qz1, qy2,qz2, L1,L2]
#beam.load[3].line = [50*units.kN/units.m, 50*units.kN/units.m,
#                     100*units.kN/units.m, 100*units.kN/units.m,
#                     1*units.m, 1*units.m]
#beam.load[3].line = {'qy1':50*units.kN/units.m, 'qy2':50*units.kN/units.m,
#                     #'qz1':100*units.kN/units.m, 'qz2':100*units.kN/units.m,
#                     'L1':1*units.m, 'L2':1*units.m}
#
print(beam.load)
#
#
#beam.load_combination[1] = [1, 0.50]
#
#
print(beam.support)
#
#for key, support in beam.support.items():
    #print(key, support.in_plane.R.value)
    #print(support)
    #print(support.in_plane)
#
# Get support reactions : R/M/theta/w
#R1 = beam.support[1].in_plane.R.value
#R2 = beam.support[2].in_plane.R.value
#
# plot beam results
data = beam.shear()
#data = beam.bending_moment()
#beam.bending_moment.plot("in_plane")
#beam.shear.plot("in_plane")
beam.deflection.plot("in_plane")
#beam.slope.plot("in_plane")
#
#beam.bending_moment.plot("out_plane")
#beam.shear.plot("out_plane")
#beam.deflection.plot("out_plane")
#beam.slope.plot("out_plane")
#
# Design
#stress = beam.response.stress()
#design = beam.design
#design.API()
#design.print_results()
#print(design)
#
print("End")
