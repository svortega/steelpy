from steelpy import Beam #, SimpleBeam
from steelpy import Units

# set units
units = Units()
#
#
#
#
#test = SimpleBeam(L=5, E=Emat)
#test.supports(end1="pinned", end2="pinned")
#test.load['test'] = ['point', -1000, 2.5]
#R = test.reactions
#resy = test.response(x=2.5, I=Iy, load_list=['test'])
#
#
# beam class
beam = Beam(name='beam_test')
beam.L = 5*units.m
#
# set section
#beam.section = section['beam']
#beam.section = ['Tubular', 500 * units.mm, 25 * units.mm]
#beam.section = ['Rectangle', 350 * units.mm, 200 * units.mm]
beam.section = ['ub', 303.4*units.mm, 6*units.mm, 165*units.mm, 10.20*units.mm]
#beam.section.d = 350 * units.mm
#beam.section.w = 350 * units.mm
#
#print(beam.section)
#
#
#beam.material = material['beam']
beam.material = ['elastic', 345.0 * units.MPa]
print(beam.material)
#
#
# set supports : pinned/fixed/free/guided/spring [k=F/x]
#
beam.support = ["pinned", "pinned"]
#
#beam.support[1] = "pinned"
#beam.support[1] = ["spring", 200 * units.kN/units.m]
#
#beam.support[2] = "pinned"
#beam.support[2] = ["spring", 200 * units.kN/units.m]
#
# set loading
#     point = [L1, Fx, Fy, Fz]
beam.P = [2.5*units.m, 60*units.N, 5000*units.N, 500*units.N,
          1*units.kN*units.m, 0*units.kN*units.m,  0*units.kN*units.m, 
          'point_1']

#beam.P = [3.0*units.m, 0*units.N, 0*units.N, 0*units.N,
#          0*units.kN*units.m, 200*units.kN*units.m, 300*units.kN*units.m,
#          'point_2']
#
#beam.P = {'L1':2.5*units.m, 'fy':100*units.N, 'name': 'point_3'}
#beam.load.point = [2.5*units.m, 1000*units.N, 2000*units.N]
#
#beam.load.point = [2*units.m, 18*units.kN, 18*units.kN]
beam.q = [10*units.N/units.m,  -500*units.N/units.m, 0*units.N/units.m,
          10*units.N/units.m, -1000*units.N/units.m, 0*units.N/units.m,
          1*units.m, 1*units.m, 'udl_1']
#
#     line = [qy1,qz1, qy2,qz2, L1,L2]
#beam.load.line = [50*units.kN/units.m, 50*units.kN/units.m,
#                     100*units.kN/units.m, 100*units.kN/units.m,
#                     1*units.m, 1*units.m]
#beam.q = {'qy1':50*units.kN/units.m, 'qy2':50*units.kN/units.m,
#          #'qz1':100*units.kN/units.m, 'qz2':100*units.kN/units.m,
#          'L1':1*units.m, 'L2':1*units.m, 'name': 'udl_2'}
#
#beam.load[1].line = {'qy1':9*units.kN/units.m, 'qz1':-9*units.kN/units.m}
#
print(beam.load)
#
#
beam.load_combination = ['udl_1', 0.80]
beam.load_combination = ['point_1', 0.85]
#
beam.load_combination = {'point_1': 0.50, 'point_2': 0.75, 'point_3': 1.0}
#
#
# beam results
#
reactions = beam.reactions()
print(reactions)
#
#
forces = beam.response()
print(forces)
#
stress = beam.stress()
print(stress)
#
#data = beam.shear()
#data = beam.bending_moment()
#beam.bending_moment.plot("in_plane")
#beam.shear.plot("in_plane")
#beam.deflection.plot("in_plane")
#beam.slope.plot("in_plane")
#
#beam.bending_moment.plot("out_plane")
#beam.shear.plot("out_plane")
#beam.deflection.plot("out_plane")
#beam.slope.plot("out_plane")
#
# Design
#design = beam.design
#design.API()
#design.print_results()
#print(design)
#
#
#
#
print("End")
