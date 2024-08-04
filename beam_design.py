from steelpy import Beam
from steelpy import Units
#
from steelpy import UFOmodel
from steelpy import Trave3D, Trave2D
#
import matplotlib.pyplot as plt

# set units
units = Units()
#
Fx_factor = 1.0
#
#
# beam class
beam = Beam(name='beam_test', steps=4)
#beam.L = 15*units.ft
beam.L = 28*units.ft
#beam.L = 5*units.m
#
# set section
#beam.section = ['Rectangle', 200 * units.mm, 100 * units.mm]
#beam.section = ['trapeziod', 200 * units.mm, 100 * units.mm, 150 * units.mm]
#beam.section = ['Circular', 500 * units.mm]
#beam.section = ['Tubular', 500 * units.mm, 25 * units.mm]
#beam.section = ['ub', 253.5*units.mm, 8.60*units.mm,
#                254.0*units.mm, 14.22*units.mm,           # top flange
#                1. * 254.0*units.mm, 1. * 14.22*units.mm, # bottom flange
#                0*units.mm]                              # radious
#beam.section = ['box', 500 * units.mm, 6 * units.mm, 250 * units.mm, 12 * units.mm]
#beam.section = ['tee', 200 * units.mm, 6 * units.mm, 150 * units.mm, 8 * units.mm]
#beam.section = ['channel', 250 * units.mm, 6 * units.mm, 150 * units.mm, 12 * units.mm]
#beam.section = ['Angle', 150 * units.mm, 6 * units.mm, 150 * units.mm, 6 * units.mm]
#
beam.section = ['ub', 351*units.mm, 8.64*units.mm,
                204.0*units.mm, 15.10*units.mm]
#
#print(beam.section)
#
#
beam.material = ['elastic', 275.0 * units.MPa]
print(beam.material)
#
#
# set supports : pinned/fixed/free/guided/spring [k=F/x]
#
#beam.support = ["fixed", "free"]
beam.support = ["fixed", "fixed"]
#
#beam.support[1] = "pinned"
#beam.support[1] = ["spring", 200 * units.kN/units.m]
#
#beam.support[2] = "pinned"
#beam.support[2] = ["spring", 200 * units.kN/units.m]
#
#
#beam.selfweight = [0, -1* units.gravity, 0] #* units.gravity
#
# set loading
#     point = [L1, Fx, Fy, Fz]
#beam.P = [2.0*units.m, 0*units.N, 100*units.kN, 0*units.N,
#          7.50*units.kN*units.m, 0*units.kN*units.m,  0*units.kN*units.m, 
#          'example_1']
#
#beam.P = [3.0*units.m, 0*units.N, 0*units.N, 0*units.N,
#          0*units.kN*units.m, 200*units.kN*units.m, 300*units.kN*units.m,
#          'point_2']
#
beam.P = {#'L1':4.5*units.m,
#          'Fy':-75 * units.kN,
#          #'Mz': 33.75  * units.kN * units.m,
#          #'Fx': -0.1 * units.N,
          'L1':14*units.ft,
#          #'Fy':15*units.kips,
          'Mx':-90*units.kip*units.inch,
          'name': 'selfweight_y'}
#
#beam.P = {'L1':28*units.ft,
#          'Fx': -445 * Fx_factor * units.kN,
#          'Fy': -4.45 * units.kN,
#          'name': 'AISC_C-C2.3'}
#
#beam.P = {'L1':90*units.inch, 'Mx':90*units.kip*units.inch, 'name': 'point_3'}
#beam.load.point = [2.5*units.m, 1000*units.N, 2000*units.N]
#
#beam.load.point = [2*units.m, 18*units.kN, 18*units.kN]
#beam.q = [10*units.N/units.m,  -500*units.N/units.m, 0*units.N/units.m,
#          10*units.N/units.m, -1000*units.N/units.m, 0*units.N/units.m,
#          1*units.m, 1*units.m, 'udl_1']
#
#     line = [qy1,qz1, qy2,qz2, L1,L2]
#beam.load.line = [50*units.kN/units.m, 50*units.kN/units.m,
#                     100*units.kN/units.m, 100*units.kN/units.m,
#                     1*units.m, 1*units.m]
#beam.q = {'qy1':50*units.kN/units.m, 'qy2':50*units.kN/units.m,
#          #'qz1':100*units.kN/units.m, 'qz2':100*units.kN/units.m,
#          'L1':1*units.m, 'L2':1*units.m, 'name': 'udl_2'}
#
#
#beam.q = {'qy':-15*units.kN/units.m,
#          #'qx': -0.001*units.N/units.m,
#          'name': 'selfweight_y'}
#
#beam.q = {'qz1':15*units.kN/units.m, 'qz2':15*units.kN/units.m,
#          'name': 'selfweight_z'}
#
#beam.q = {'qy':0.20 * units.kip / units.ft,
#          'name': 'AISC_C-C2.2'}
#
#beam.load[1].line = {'qy1':9*units.kN/units.m, 'qz1':-9*units.kN/units.m}
#
#print(beam.load)
#
# TODO: fix combination
#beam.load_combination = ['AISC_C-C2.2', 1.0]
#beam.load_combination = ['AISC_C-C2.3', 1.0]
#
#beam.load_combination = {'point_1': 0.50, 'point_2': 0.75, 'point_3': 1.0}
#
#
# beam results
#Mx =  -90*units.kip*units.inch
#print(Mx.convert('newton*metre'))
#
#reactions = beam.reactions()
#print(reactions)
#
#
#forces = beam.response()
#print(forces)
#
#
#force_grp = (forces.groupby(['load_name', 'component_name',
#                             'load_title', 'load_type',
#                             'load_level', 'load_system',
#                             'element_name' , 'node_end'])
#             [['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']].sum())
#
#force_grp = force_grp.reset_index(names=['load_name', 'component_name',
#                                         'load_title', 'load_type',
#                                         'load_level', 'load_system',
#                                         'element_name' ,'node_end'])
#
#force_grp.plot.line(x='node_end', y='Fx')
#force_grp.plot.line(x='node_end', y='Fy')
#force_grp.plot.line(x='node_end', y='Fz')
#force_grp.plot.line(x='node_end', y='Mx')
#force_grp.plot.line(x='node_end', y='My')
#force_grp.plot.line(x='node_end', y='Mz')
#force_grp.plot.line(x='node_end', y='B')
#force_grp.plot.line(x='node_end', y='Tw')
#plt.show()
#
#stress = beam.stress()
#print(stress)
#
#
#stress.plot.line(x='node_end', y='tau_y')
#stress.plot.line(x='node_end', y='sigma_x')
#stress.plot.line(x='node_end', y='sigma_z')
#plt.show()
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
#1 / 0
#
f2umodel = UFOmodel("BeamDesign_0")
mesh = f2umodel.mesh()

mesh[10] = 'beam design example 0'

mesh[10].material([[10, 'linear', 275.0 * units.MPa]])

mesh[10].section([[155, 'ub', 253.5*units.mm, 8.60*units.mm,
               254.0*units.mm, 14.22*units.mm,  
               1. * 254.0*units.mm, 1. * 14.22*units.mm,
               0*units.mm],
              [150, 'Tubular', 500 * units.mm, 25 * units.mm],
              [15, 'ub', 351*units.mm, 8.64*units.mm,
               204.0*units.mm, 15.10*units.mm]])

mesh[10].node([(1, 0*units.ft,   0*units.ft),
           ##(2, 15*units.ft,  0*units.ft),
           #(2, 1*units.m,  0*units.m),
           #(3, 2*units.m,  0*units.m),
           #(4, 3*units.m,  0*units.m),
           #(5, 4*units.m,  0*units.m),
           #(6, 5*units.m,  0*units.m)])
           (6, 28*units.ft,  0*units.ft)])

mesh[10].boundary([[1, 'support', 'fixed'],
               #[6, 'support', 'free']])
                [6, 'support', (0,1,1,1,0,0)]])
#
# [element_id, node1, node2, material, section, roll_angle]
mesh[10].element([#(1,  'beam',  1, 2, 10, 15, 0),
              #(2,  'beam',  2, 3, 10, 15, 0),
              #(3,  'beam',  3, 4, 10, 15, 0),
              #(4,  'beam',  4, 5, 10, 15, 0),
              #(5,  'beam',  5, 6, 10, 15, 0)])
              (5,  'beam',  1, 6, 10, 15, 0)])
#
#
# ----------------------------------------------------
# Load input
# ----------------------------------------------------
#
load = mesh[10].load()
#
# ----------------------------------------------------
# Basic Load
basic = load.basic()
#
#
#basic[1] = 'dispExample'
#nodeload = basic[1].node()
#nodeload[2].displacement = {'y': -0.103 * units.m,
#                            'rz': 0.03 * units.rad,}
#
# load_title, 'beam', beam_id, 'point',  L0,x,y,z,mx,my,mz, comment(optional)
#basic.beam([
#             ['snow load', 1, 'point',
#             90 * units.inch,  # L0
#             0 * units.kip,  # x
#             15 * units.kip, # y
#             0 * units.kip,  # z
#             90 * units.kip*units.inch,  # mx
#             'point_3'],
#            ['snow load', 1, 'line',
#             0 * units.N / units.m,   # qx
#             15 * units.kN / units.m, # qy
#             0 * units.N / units.m,   # qz
#             0 * units.N / units.m,   # qt
#             'udly_3']])
#
#
#basic.beam({'load': 'snow load',
#            'beam': 5, #(1, 2, 3, 4, 5),
#            'type': 'line',
#            'qy': -7.5 * units.kN / units.m})
#
#basic.beam({'load': 'snow load',
#            'beam': 5, 
#            'type': 'line',
#            'qy': -7.5 * units.kN / units.m})
#
#basic.beam({'load': 'snow load',
#            'beam': 5, 
#            'type': 'point',
#            'Fy': -75 * units.kN,
#            'L0': 4.50 * units.m})
#
#basic.beam({'load': 'snow load',
#            'beam': 4, 
#            'type': 'line',
#            'qy': -15 * units.kN / units.m})
#
#basic.beam({'load': 'snow load',
#            'beam': 5, 
#            'type': 'line',
#            'qy': -15 * units.kN / units.m})

basic.beam({'load': 'AISC_C-C2.2',
            'beam': 5,
            'type': 'line',
            'qy': -0.20 * units.kip / units.ft})
#
basic.beam({'load': 'snow load',
            'beam': 5, 
            'type': 'point',
            'Fy': -75 * units.kN,
            'L0': 0.50 * units.m,})
#
#basic.beam(
           #{'load':'snow load',
           # 'beam': 1,
           # 'type': 'point',
           # 'L0': 90 * units.inch,
           # 'y': 15 * units.kip,
           # 'mx': 90 * units.kip*units.inch,
           # 'title': 'point_3',},
           #{'load':'snow load',
           # 'beam': 1,
           # 'type': 'line',
           # 'qy': 15 * units.kN / units.m,
           # 'title': 'udly_3'})
#
#
# [load_title, 'node', node_id, 'point',  Fx,Fy,Fz,mx,my,mz, comment(optional)],
# [load_title, 'node', node_id, 'displacement',  x,y,z,rx,ry,rz, comment(optional)],
# [load_title, 'node', node_id, 'mass' ,  x,y,z, comment(optional)]
#
#basic.node([['dispExample', 2, 'displacement',
#             0 * units.m,
#             -0.103 * units.m,
#             0 * units.m, # x,y,z,
#             0 * units.rad,
#             0 * units.rad,
#             0.03 * units.rad,  # rx,ry,rz
#             'nodey_2']])
#            ['wind load', 3, 'load', -200 * Punit, 'nodex_4']])
#
#basic.node({'load': 'dispExample',
#            'node': 6,
##            'type': 'force',
##            'Fy': -75 * units.kN ,})
#            'type': 'displacement', 
#            'y': -0.0511813 * units.m,})
#            #'rz': 0.0136483 * units.rad,})
#
basic.node({'load': 'AISC_C-C2.3',
            'node': 6,
            'type': 'force',
            'Fx': -667 * Fx_factor * units.kN,})
#            'Fy': -4.45 * units.kN})
#            'y': -0.0511813 * units.m,})
#            #'rz': 0.0136483 * units.rad,})
#
# ----------------------------------------------------
# Load Combinations
#
#
#comb = load.combination([['factored comb 1', 'basic', 11, 1.20],
#                         ['factored comb 1', 'basic', 22, 1.25],
#                         ['factored comb 1', 'basic', 33, 1.30]])
#
#
comb = load.combination()
comb[100] = 'Factored Comb 1'
comb[100].basic['AISC_C-C2.2'] = 1.00
comb[100].basic['AISC_C-C2.3'] = 1.00
#comb[100].basic['snow load'] = 1.25
#comb[100].basic[33] = 1.30
#
#comb[200] = 'factored comb 2'
#comb[200].basic[11] = 1.35
#comb[200].basic[22] = 1.40
#comb[200].basic[33] = 1.45
#
comb[300] = 'factored comb 3'
comb[300].combination[100] = 1.50
comb[300].basic['snow load'] = 1.60
#
#
# ----------------------------------------------------
# Meshing input
# ----------------------------------------------------
#
mesh.build()
#
# ----------------------------------------------------
# Plot mesh
# ----------------------------------------------------
#
#plot = mesh.plot()
#plot.frame()
#
# ----------------------------------------------------
# Structural Analysis
# ----------------------------------------------------
#
frame = Trave2D(mesh=mesh[10])
frame.static(second_order=True)
results = frame.results(beam_steps=4)
#
# ----------------------------------------------------
# Results
# ----------------------------------------------------
#
beamres = results.beam()
#print(beamres)
bdisp = beamres.displacement()
bdisp2 = bdisp.df
bforce = beamres.force()
forces2 = bforce.df
bstress = beamres.stress()
bstress2 = bstress.df
#
print("End")
