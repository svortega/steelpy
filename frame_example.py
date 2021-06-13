from steelpy import Trave3D
from steelpy import Units
from steelpy import f2uModel

f2umodel = f2uModel(component='test2', mesh_type="inmemory") #
units = Units()
#
material = f2umodel.materials
material[1] = 'elastic'
material[1].Fy = 345.0 * units.MPa
material[1].E = 205000.0 * units.MPa
material[1].G = 77200.0 * units.MPa
print(material)
#
section = f2umodel.sections
section[1] = 'tubular'
section[1].D = 500 * units.mm
section[1].t = 25 * units.mm
print(section)
#
#
mesh = f2umodel.mesh
#
# nodes corrdinates [x, y, z]
nodes = mesh.nodes
nodes[1] = [0, 0, 0]
nodes[2] = [3, 0, 0]
nodes[3] = [6, 0, 0]
# boundary nodes
boundary = mesh.boundaries
boundary.node[1] = [1,1,1,1,0,0]
boundary.node[3] = [1,1,1,0,0,0]
# elements [type, node1, node2, material, section]
elements = mesh.elements
elements[1] = ['beam', 1, 2, 1, 1]
elements[2] = ['beam', 2, 3, 1, 1]
#
#elements[1] = ['beam', 1, 3, 1, 1]
#
# loading
#
load = f2umodel.load

basic = load.basic
basic[111] = 'wind load'
#
#basic[111].point_node[2] = [0, -1000, 0, 'test']
#basic[111].point_beam[1] = [3.0, 0, -1000, 0, 'test']
basic[111].line_beam[1] = [0, -1_000, 0, 'wind_1']
basic[111].line_beam[2] = [0, -1_000, 0, 'wind_1']
#
print(load)
#
f2umodel.build()
#
frame = Trave3D()
frame.f2u_model = f2umodel
frame.run_static()
frame.print_results()
print('-->')