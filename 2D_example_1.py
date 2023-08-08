# Import steelpy modules
from steelpy import Units
from steelpy import f2uModel
from steelpy import Trave2D
#
#
units = Units()
#
f2umodel = f2uModel(component="test1")
#
# ----------------------------------------------------
# Material input
# ----------------------------------------------------
#
f2umodel.materials([10, 'linear', 345.0 * units.MPa])
print(f2umodel.materials())
#
#
# ----------------------------------------------------
# Section Input
# ----------------------------------------------------
#
f2umodel.sections([20, 'ub', 240*units.mm, 6.2*units.mm, 120*units.mm, 9.8*units.mm])
#f2umodel.sections([20, 'Tubular', 300 * units.mm, 10 * units.mm])
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
basePoint = 0.0 * units.m
storeyHeight = 4.0 * units.m
bayWidth = 6.0 * units.m
#
# nodes corrdinates [node_id, x, y, z=0]
#
mesh.nodes([(1, basePoint,    basePoint),
            (2, bayWidth,     basePoint),
            (3, 2.0*bayWidth, basePoint),
            (4, basePoint,    storeyHeight),
            (5, bayWidth,     storeyHeight),
            (6, 2.0*bayWidth, storeyHeight),
            (7, basePoint,    2.0*storeyHeight),
            (8, bayWidth,     2.0*storeyHeight),
            (9, 2.0*bayWidth, 2.0*storeyHeight)])
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
                 [2, 'node', 'fixed'],
                 [3, 'node', 'fixed']])
#
print(mesh.boundaries())
#
# ----------------------------------------------------
# Element input
# ----------------------------------------------------
#
# Example:
# Elements[number] = [id, beam, material, section, node1, node2, roll_angle]
# Elements[number] = [id, plate, material, section, node1, node2, node3, node4]
#
#
mesh.elements([(1,  'beam',  1, 4, 10, 20, 0),
               (2,  'beam',  2, 5, 10, 20, 0),
               (3,  'beam',  3, 6, 10, 20, 0),
               (4,  'beam',  4, 5, 10, 20, 0),
               (5,  'beam',  5, 6, 10, 20, 0),
               (6,  'beam',  4, 7, 10, 20, 0),
               (7,  'beam',  5, 8, 10, 20, 0),
               (8,  'beam',  6, 9, 10, 20, 0),
               (9, 'beam',  7, 8, 10, 20, 0),
               (10, 'beam',  8, 9, 10, 20, 0)])
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
# Basic Load
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
nullLoad = 0 * units.N
pointLoad = -1_000 * units.N
basic[11] = 'Buckling Example'
basic[11].node([[4, 'load', nullLoad, pointLoad, 'buckling_1'],
                [5, 'load', nullLoad, pointLoad, 'buckling_2'],
                [6, 'load', nullLoad, pointLoad, 'buckling_3'],
                [7, 'load', nullLoad, pointLoad, 'buckling_4'],
                [8, 'load', nullLoad, pointLoad, 'buckling_5'],
                [9, 'load', nullLoad, pointLoad, 'buckling_6']])
#
print(basic)
#
#
#
# ----------------------------------------------------
# Meshing input
# ----------------------------------------------------
#
#
f2umodel.build()
#
# ----------------------------------------------------
# Plot mesh
# ----------------------------------------------------
#
#plot = f2umodel.plot()
#plot.mesh()
#plot.basic_load(name=22)
#
#
#from steelpy.beam.frame2D.process.assembleK import assembleKa, get_Kf
#nodes = mesh.nodes()
#beams = mesh.elements().beams()
#boundaries = mesh.boundaries()
#Ka, KlStorageFrame, TlgStorageFrame = assembleKa(nodes, frameElements=beams)
#kf = get_Kf(Ka, nodes, boundaries)
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