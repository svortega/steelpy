from steelpy import Trave3D
from steelpy import Units
from steelpy import f2uModel

f2umodel = f2uModel(component="test1")
                   #mesh_type="inmemory")
units = Units()
#
# f2umodel.materials = [['elastic', 345.0 * units.MPa],
#                       ['elastic', 345.0 * units.MPa, 490.0 * units.MPa,
#                        205000.0 * units.MPa, 77200.0 * units.MPa]]
#
material = f2umodel.materials
# material[n] = [elastic, Fy, Fu, E, G, Poisson, density, alpha]
#material[1] = 'elastic'
material[1] = ['elastic', 345.0 * units.MPa]
#
material[3] = ['elastic', 345.0 * units.MPa, 490.0 * units.MPa,
               205000.0 * units.MPa, 77200.0 * units.MPa]
#
print(material)
#
section = f2umodel.sections
section[1] = ['Rectangle', 0.40, 0.35]
section[2] = ['Tubular', 500 * units.mm, 25 * units.mm]
#
section[5] = ['Circular', 500 * units.mm]
#section[10] = ['trapeziod', 400 * units.mm, 350 * units.mm, 200 * units.mm] #
section[20] = ['ub', 700 * units.mm, 12 * units.mm, 264 * units.mm, 15.5 * units.mm]
section[30] = ['box', 500 * units.mm, 6 * units.mm, 250 * units.mm, 12 * units.mm]
section[40] = ['channel', 250 * units.mm, 6 * units.mm, 150 * units.mm, 12 * units.mm]
section[50] = ['tee', 150 * units.mm, 6 * units.mm, 150 * units.mm, 6 * units.mm]
#
print(section)
#
mesh = f2umodel.mesh
#
# nodes corrdinates [x, y, z]
#mesh.nodes = [[1, 0, 0, 0],
#              [2, 3, 4, 0],
#              [3, 8, 4, 0]]
nodes = mesh.nodes
nodes[1] = [0 , 0 , 0 ]
nodes[2] = [0 , 6 , 0 ]
#nodes[2] = [0 , 3 , 0 ]
nodes[3] = [4 , 3 , 0 ]
nodes[4] = [4 , 0 , 0 ]
# boundary nodes
boundary = mesh.boundaries
boundary.node[1] = 'fixed'
boundary.node[4] = 'fixed'
# elements [beam, node1, node2, material, section, roll_angle]
elements = mesh.elements
elements[1] = ['beam', 1, 2, 1, 1, 0]
elements[2] = ['beam', 2, 3, 1, 1, 0]
elements[3] = ['beam', 3, 4, 1, 1, 0]
#
#
properties = f2umodel.properties
hydro = properties.hydrodynamic
mg = hydro.marine_growth
mg['MG_1'] = "profile"
mg['MG_1'].level[1] = [   0 * units.m, 60 * units.mm]
mg['MG_1'].level[2] = [ -10 * units.m, 60 * units.mm]
mg['MG_1'].level[3] = [ -50 * units.m, 30 * units.mm]
mg['MG_1'].level[4] = [-100 * units.m, 10 * units.mm]
mg['MG_1'].level[5] = [-120 * units.m,  0 * units.mm]
mg['MG_1'].set_default()
#
mg['MG_2'] = "constant"
mg['MG_2'].thickness = 100 * units.mm
mg['MG_2'].elements = 1, 3
#
#
cdcm = hydro.CdCm
#
#
#
#mesh.renumbering()
#
#
# groups
#groups = mesh.groups
#groups[1] = 'element'
#
# loading
load = f2umodel.load
#
# load.basic.system = 'local'  # This will affect beam load only (global default)
# load.basic = [[load_title, 'node', node_number, 'point',  x,y,z,mx,my,mz, comment(optional)],
#               [load_title, 'node', node_number, 'mass' ,  x,y,z, comment(optional)]
#               [load_title, 'beam', beam_number, 'line' ,  qx0,qy0,qz0, qx1,qy1,qz1, L0,L1, comment(optional)],
#               [load_title, 'beam', beam_number, 'point',  L0,x,y,z,mx,my,mz, comment(optional)]]
#
#load.basic = [['wind load x', 'node', 2, 'point', 0, -4_000_000_000, 0],
#              ['wind load x', 'node', 3, 'point', 0, -2_000_000_000, 0],
#              ['snow load'  , 'beam', 2, 'line'  , -1_000_000, 0, 0]]
#
#
basic = load.basic
basic[11] = 'wind load x'
#
# basic[1].node = [[node_number, 'point', x,y,z,mx,my,mz, comment(optional)],
#                  [node_number, 'mass' , x,y,z, comment(optional)]]
#
# basic[1].beam = [[beam_number, 'point', L0,x,y,z,mx,my,mz, comment(optional)],
#                  [beam_number, 'line' , qx0,qy0,qz0, qx1,qy1,qz1, L0,L1, comment(optional)]]
#
##basic[1].udl_beam.coordinate_system = 'local'
basic[11].point_node[2] = [0, -4_000_000_000, 'wind_1']
basic[11].point_node[3] = [0, -2_000_000_000, 'wind_2']
basic[11].line_beam[1] = [-1_000_000, 0, 0, 'wind_1']
##basic[1].udl_beam[2] = [-10_000_000, 0, 0]
##
basic[22] = 'snow load'
##basic[2].udl_beam.coordinate_system = 'local'
basic[22].line_beam[2] = [0, 0, -480_000, 'snow load on roof']
##basic[2].udl_beam[2] = [-10_000_000, 0, 0]
##basic[2].udl_beam[2] = [-4_800_000, 0, 0]
##basic[2].point_beam[2] = [2.5, 0, 0, -2_400_000]
##basic[2].point_beam[2] = [2.0, 0, 0, -1_920_000, 0, 0, 0]
##basic[2].point_beam[2] = [2.5, 0, 0, 0,  -100_000_000, 0, 0]
##
##
##basic = load.basic
##basic[1] = 'wind load y'
##basic[1].point_node[2] = [ 0, 0, -400_000, 0, 0, 0]
##basic[1].point_node[3] = [ 0, 0, -200_000, 0, 0, 0]
##
basic[33] = 'crane load'
## tilt
##basic[3].point_beam.coordinate_system = 'local'
basic[33].point_beam[2] = [2.5, 0, 0, 0,  0, 0, -100_000_000, 'lifting module A']
##print(basic[3].point_beam.coordinate_system)
##
##basic[3].point_beam[2] = [2.0, -100_000_000]
## flat
##basic[3].point_beam[2] = [2.0, -100_000_000, 0,  0, 0, 0, 0]
##
##basic[2] = 'snow load'
##basic[2].udl_beam.coordinate_system = 'local'
##print(basic[2].udl_beam.coordinate_system)
basic[22].line_beam[2] = [0, -100_000_000, 0, 0, -200_000_000, 0, 1.0, 1.0]
##basic[2].udl_beam[2] = [0, -480_000_000*0.5, 0, 0, -480_000_000*0.5, 0]
##basic[2].udl_beam[2] = [0, -480_000_000*0.5, 0, 0, -480_000_000*0.5, 0]
#
# create new basic load
basic[44] = 'dead load'
basic[44].selfweight.y = -1
#
#print(basic)
#
##
comb = load.combination
comb[100] = 'factored comb 1'
comb[100].basic_load[11] = 1.20
comb[100].basic_load[22] = 1.25
comb[100].basic_load[33] = 1.30
#
comb[200] = 'factored comb 2'
comb[200].basic_load[11] = 1.35
comb[200].basic_load[22] = 1.40
comb[200].basic_load[33] = 1.45
#
comb[300] = 'factored comb 3'
comb[300].load_combination[100] = 1.50
comb[300].load_combination[200] = 1.55
comb[300].basic_load[44] = 1.60
#
#print(comb)
#
print(load)
#
f2umodel.build()
#
#
frame = Trave3D()
frame.f2u_model = f2umodel
frame.run_static()
frame.print_results()
print('-->')