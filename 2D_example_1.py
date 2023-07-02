# Import steelpy modules
from steelpy import Units
from steelpy import f2uModel
#from steelpy import Trave3D
#from steelpy.beam.frame2D.process.assembleK import assembleKa, get_Kf
#from steelpy.beam.frame2D.process.assemble_f2u import solve_deflections
#from steelpy.beam.frame2D.process.postprocess_f2u import beam_end_force
# 
from steelpy.trave3D.postprocessor.operations import beam_end_force,  beam_force
from steelpy.f2uModel.mesh.process.matrix.Kassemble import assemble_Kmatrix_np
from steelpy.trave3D.processor.static_solver import solve_deflections
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
f2umodel.materials([10, 'elastic', 345.0 * units.MPa])
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
storeyHeight = 4.0
bayWidth = 6.0
#
# nodes corrdinates [node_id, x, y, z=0]
#
mesh.nodes([(1, 0.0,          0.0),
            (2, bayWidth,     0.0),
            (3, 2.0*bayWidth, 0.0),
            (4, 0.0,          storeyHeight),
            (5, bayWidth,     storeyHeight),
            (6, 2.0*bayWidth, storeyHeight),
            (7, 0.0,          2.0*storeyHeight),
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
mesh.boundaries([[1, 'support', 'fixed'],
                 [2, 'support', 'fixed'],
                 [3, 'support', 'fixed']])
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
mesh.elements([(1,  'beam',  1, 4, 10, 20, 0),
               (2,  'beam',  2, 5, 10, 20, 0),
               (3,  'beam',  3, 6, 10, 20, 0),
               (4,  'beam',  4, 5, 10, 20, 0),
               (5,  'beam',  5, 6, 10, 20, 0),
               (6,  'beam',  4, 7, 10, 20, 0),
               (7,  'beam',  5, 8, 10, 20, 0),
               (9,  'beam',  6, 9, 10, 20, 0),
               (10, 'beam',  7, 8, 10, 20, 0),
               (11, 'beam',  8, 9, 10, 20, 0)])
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
#
#
basic[11] = 'Buckling Example'
basic[11].node([[4, 'load', 0, -1_000, 'buckling_1'],
                [5, 'load', 0, -1_000, 'buckling_2'],
                [6, 'load', 0, -1_000, 'buckling_3'],
                [7, 'load', 0, -1_000, 'buckling_4'],
                [8, 'load', 0, -1_000, 'buckling_5'],
                [9, 'load', 0, -1_000, 'buckling_6']])
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
# ----------------------------------------------------
# Structural Analysis
# ----------------------------------------------------
#
nodes = mesh.nodes()
boundaries = mesh.boundaries()
elements = mesh.elements()
beams = elements.beams()
load = mesh.load()
basic = load.basic()
df_nload = basic.node_df()
#
#
neq = nodes.neq(supports=boundaries._nodes)
jbc = nodes.jbc(supports=boundaries._nodes)
#Kf = assemble_Kmatrix_np(elements=elements, jbc=jbc, neq=neq, m2D=True)
jbc, K = mesh.K(solver=False, log=False, m2D=True)
#
df_ndisp, df_nload = solve_deflections(df_nload, method=None, m2D=True)
#
df_nforce = beam_end_force(elements=elements, df_ndisp=df_ndisp, m2D=True)
#
# get beam end node forces
# -----------------------------------
# get beam force along lenght 
df_membf = beam_force(elements, 
                      basic_load=basic,
                      df_ndisp=df_ndisp, 
                      df_nforce=df_nforce,
                      m2D=True)
#
#frame = Trave3D()
#frame.mesh = mesh
#results = frame.run_static(method = 'banded')
#results.print()
print('-->')