# Import steelpy modules
from steelpy import Units
from steelpy import f2uModel
from steelpy import Trave2D
from steelpy import Trave3D
#
#
units = Units()
#
f2umodel = f2uModel(component="test1")
#
# ----------------------------------------------------
# Material input
# ----------------------------------------------------
# [elastic, Fy, Fu, E, G, Poisson, density, alpha]
f2umodel.materials([[10, 'linear', 345.0 * units.MPa, 490.0 * units.MPa, 200 * units.GPa],
                    [15, 'linear', 245.0 * units.MPa, 490.0 * units.MPa, 200 * units.GPa]])
#
print(f2umodel.materials())
#
#
# ----------------------------------------------------
# Section Input
# ----------------------------------------------------
#
f2umodel.sections([[20, 'ub', 240*units.mm, 6.2*units.mm, 120*units.mm, 9.8*units.mm],
                   [25, 'Tubular', 300 * units.mm, 10 * units.mm]])
#
print(f2umodel.sections())
#
#
# ----------------------------------------------------
# ----------------------------------------------------
# Mesh 
# ----------------------------------------------------
# ----------------------------------------------------
#
mesh = f2umodel.mesh()
#
#
# ----------------------------------------------------
# Node input
# ----------------------------------------------------
#
storeyBase = 0.0 * units.m
storeyHeight1 = 6.0 * units.m
storeyHeight2 = 3.0 * units.m
bayWidth = 4.0 * units.m
#
# nodes corrdinates [node_id, x, y, z=0]
#
mesh.nodes([(1, storeyBase,   storeyBase),
            (2, storeyBase, storeyHeight1),
            (3, bayWidth, storeyHeight2),
            (4, bayWidth,   storeyBase)])
#
print(mesh.nodes())
#
#
# ----------------------------------------------------
# boundary Input
# ----------------------------------------------------
#
# [id, type, fixity]
mesh.boundaries([[1, 'node', 'fixed'],
                 [4, 'node', 'fixed']])
#
print(mesh.boundaries())
#
# ----------------------------------------------------
# Element input
# ----------------------------------------------------
#
# Example:
# Elements[number] = [beam, material, section, node1, node2, roll_angle]
# Elements[number] = [plate, material, section, node1, node2, node3, node4]
#
#
mesh.elements([(1,  'beam',  1, 2, 10, 20, 0),
               (2,  'beam',  2, 3, 15, 25, 0),
               (3,  'beam',  3, 4, 10, 20, 0)])
#
#
#
# ----------------------------------------------------
# mesh data
# ----------------------------------------------------
#
print(mesh.elements().beams())
#
#
# ----------------------------------------------------
# Load input
# ----------------------------------------------------
#
#
# ----------------------------------------------------
# Basic Load (global system default)
#
# loading
load = mesh.load()
#
# load.basic.system = 'local'  # This will affect beam load only (global default)
#
#
# load numbering is automatic (consecutive)
# load.basic([[load_title, 'node', node_number, 'point',  x,y,z,mx,my,mz, comment(optional)],
#             [load_title, 'node', node_number, 'mass' ,  x,y,z, comment(optional)]
#             [load_title, 'beam', beam_number, 'line' ,  qx0,qy0,qz0, qx1,qy1,qz1, L0,L1, comment(optional)],
#             [load_title, 'beam', beam_number, 'point',  L0,x,y,z,mx,my,mz, comment(optional)]])
#
#basic = load.basic([['wind load x', 'node', 2, 'point', 0, -4_000_000_000, 0],
#                    ['wind load x', 'node', 3, 'point', 0, -2_000_000_000, 0],
#                    ['snow load'  , 'beam', 2, 'line' , -1_000_000, 0, 0]])
#
#print(basic)
#
basic = load.basic()
#
#
# basic[1].node([[node_number, 'point', x,y,z,mx,my,mz, comment(optional)],
#                [node_number, 'mass' , x,y,z, comment(optional)]])
#
# basic[1].beam([[beam_number, 'point', L0,x,y,z,mx,my,mz, comment(optional)],
#                [beam_number, 'line' , qx0,qy0,qz0, qx1,qy1,qz1, L0,L1, comment(optional)]])
#
#basic[11].node([[2, 'load', 0, -4_000_000_000, 'wind_1'],
#                [3, 'load', 0, -2_000_000_000, 'wind_2']])
#
#
#basic[10] = 'Beam load'
#basic[10].beam([[8, 'line', 0, -1000, 'udl_1'],
#                [9, 'line', 0, -1000, 'udl_2'],
#                [10, 'line', 0, -1000, 'udl_2']])
#
#
nullLoad = 0 * units.N
pointLoad = -1_000 * units.N
basic[11] = 'example'
basic[11].node([[2, 'load',
                 400_000 * units.N, nullLoad, nullLoad,
                 100_000 * units.N * units.m, 100_000 * units.N * units.m,
                 'nodex_1'],
                [3, 'load',
                 -1 * 200_000 * units.N, nullLoad, nullLoad, 'nodex_2']])

#
nullUDL= 0 * units.N / units.m
basic[11].beam([[2, 'line',
                nullUDL, -50_000 * units.N/units.m, nullUDL,
                nullUDL, -50_000 * units.N/units.m, 20_000 * units.N/units.m,
                'udly_1'],
                [1, 'point', 3* units.m, 10 * units.kN, 'point_1']])
#
#
print(basic)
#
#for key, items in basic.items():
#    key, items
#    for key2, items2 in items.node().items():
#        key2, items2
#        for items3 in items2.load:
#            print(items3)
#
#
# ----------------------------------------------------
# Meshing
# ----------------------------------------------------
#
#
mesh.build()
#
# ----------------------------------------------------
# Plotting
# ----------------------------------------------------
#
# Structure
#
plot = mesh.plot()
#plot.frame()
#plot.material()
#plot.section()
#
# Loading
#
plotload = load.plot()
plotload.basic()
#
#
#
# ----------------------------------------------------
# Structural Analysis
# ----------------------------------------------------
#
frame = Trave2D()
frame.mesh = mesh
frame.static()
results = frame.solve()
print(results)
print('-->')