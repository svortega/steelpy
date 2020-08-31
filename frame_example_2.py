from steelpy import Trave
from steelpy import Units
from steelpy import f2uModel

f2umodel = f2uModel(component='test1')
units = Units()

#
material = f2umodel.materials
material[1] = 'elastic'
material[1].Fy = 345.0 * units.MPa
material[1].E = 205000.0 * units.MPa
material[1].G = 77200.0 * units.MPa
#
material[3] = 'elastic'
material[3].Fy = 345.0 * units.MPa
material[3].E = 205000.0 * units.MPa
material[3].G = 77200.0 * units.MPa
#
section = f2umodel.sections
section[1] = 'Rectangle'
section[1].d = 350 * units.mm
section[1].w = 350 * units.mm
#
section[10] = 'Rectangle'
section[10].d = 350 * units.mm
section[10].w = 350 * units.mm
#
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
# elements [type, node1, node2, material, section]
elements = mesh.elements
elements[1] = ['beam', 1, 2, 1, 1]
elements[2] = ['beam', 2, 3, 1, 1]
elements[3] = ['beam', 3, 4, 1, 1]
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
basic = load.basic
basic[1] = 'wind load x'
##basic[1].udl_beam.coordinate_system = 'local'
basic[1].point_node[2] = [0, -4_000_000_000]
basic[1].point_node[3] = [0, -2_000_000_000]
basic[1].udl_beam[1] = [-1_000_000, 0, 0]
##basic[1].udl_beam[2] = [-10_000_000, 0, 0]
##
basic[2] = 'snow load'
##basic[2].udl_beam.coordinate_system = 'local'
basic[2].udl_beam[2] = [0, 0, -480_000]
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
basic[3] = 'crane load'
## tilt
##basic[3].point_beam.coordinate_system = 'local'
basic[3].point_beam[2] = [2.5, 0, 0, 0,  0, 0, -100_000_000]
##print(basic[3].point_beam.coordinate_system)
##
##basic[3].point_beam[2] = [2.0, -100_000_000]
## flat
##basic[3].point_beam[2] = [2.0, -100_000_000, 0,  0, 0, 0, 0]
##
##basic[2] = 'snow load'
##basic[2].udl_beam.coordinate_system = 'local'
##print(basic[2].udl_beam.coordinate_system)
basic[2].udl_beam[2] = [0, -100_000_000, 0, 0, -200_000_000, 0, 1.0, 1.0]
##basic[2].udl_beam[2] = [0, -480_000_000*0.5, 0, 0, -480_000_000*0.5, 0]
##basic[2].udl_beam[2] = [0, -480_000_000*0.5, 0, 0, -480_000_000*0.5, 0]
#
#
# create new basic load
basic[4] = 'dead load'
basic[4].selfweight.y = -1
#
##
comb = load.combination
comb[1] = 'factored comb'
comb[1].basic_load[1] = 1.20
comb[1].basic_load[2] = 1.35
comb[1].basic_load[3] = 1.50
#
#
f2umodel.get_mesh()
#
#
frame = Trave()
frame.f2u_model = f2umodel
frame.solve_static()
frame.print_results()
print('-->')