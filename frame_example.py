from steelpy import Trave3D
from steelpy import Units
from steelpy import f2uModel

f2umodel = f2uModel(component='test2', mesh_type="inmemory") #
units = Units()
#
material = f2umodel.materials()
material[1] = ['elastic', 345.0 * units.MPa]
print(material)
#
section = f2umodel.sections()
section[1] = ['tubular', 500 * units.mm, 25 * units.mm]
print(section)
#
#
mesh = f2umodel.mesh()
#
# nodes corrdinates [x, y, z]
nodes = mesh.nodes()
nodes[1] = [0, 0, 0]
nodes[2] = [3, 0, 0]
nodes[3] = [6, 0, 0]
# boundary nodes
boundary = mesh.boundaries()
boundary.node[1] = [1,1,1,1,1,1]
#boundary.node[3] = [1,1,1,0,0,0]
print(boundary)
# elements [type, node1, node2, material, section]
elements = mesh.elements()
elements[1] = ['beam', 1, 2, 1, 1]
elements[2] = ['beam', 2, 3, 1, 1]
#
#elements[1] = ['beam', 1, 3, 1, 1]
#
# loading
#
load = mesh.load()

#
#load.basic = [['wind load x', 'node', 2, 'point', 0, -4_000_000_000, 0],
#                ['wind load x', 'node', 3, 'point', 0, -2_000_000_000, 0]]
#         ['snow load'  , 'beam', 2, 'line'  , -1_000_000, 0, 0]]
#
basic = load.basic()
#
basic[111] = 'wind load'
#
#basic[111].node = [[2, 'load', 0, -1000, 0, 'test1'],
#                   [3, 'load', 0, -1000, 0, 'test2']]
#
basic[111].node[3].load = [0, -1000, 0, 'ytest']
basic[111].node[3].load = [0, 0, -1000, 'ztest']
#
node2 = basic[111].node[2]
node2.load = [-1000, 0, 0, 'xtest']
node2.mass = [0, -9.81, 0, 'ymass_test']
node2.mass = [0, 0, -9.81, 'zmass_test']
#
node3 = basic[111].node[3]
node3.load = [-1000, 0, 0, 'xtest3']
#
#
basic[111].beam[1].point = [3.0, 0, -1000, 0, 'test']
basic[111].beam[1].line = [0, -1_000, 0, 'wind_1']
basic[111].beam[2].line = [0, 0, -1_000, 'wind_2']
#
beam_load = basic[111].beam[2]
beam_load.line = [0, -2_000, 0, 'wind_3']
#
#
#basic[111].beam = [[1, 'point', 3.0, 0, -1000, 0, 'test'],
#                   [1, 'line', 0, -1_000, 0, 'wind_1'],
#                   [2, 'line', 0, 0, -1_000, 'wind_2']]
#
print(load)
#
#for load_name, lcase in basic.items():
#    for name, loads in lcase.node.items():
#        for key, items in loads.point.items():
#            for item in items:
#                print(item)
#    #
#    for name, loads in lcase.beam.items():
#        for key, items in loads.point.items():
#            for item in items:
#                print( item )
#        #
#        for key, items in loads.line.items():
#            for item in items:
#                print( item )
#
f2umodel.build()
#
frame = Trave3D()
frame.mesh = mesh
results = frame.run_static()
results.print()
print('-->')