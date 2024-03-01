# Import steelpy modules
from steelpy import Units
from steelpy import UFOmodel
from steelpy import Trave2D
from steelpy import Trave3D
#
#
units = Units()
#
f2umodel = UFOmodel("2D_example_0")
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
# Material input
# ----------------------------------------------------
# [elastic, Fy, Fu, E, G, Poisson, density, alpha]
mesh.material([[10, 'linear', 345.0 * units.MPa, 490.0 * units.MPa, 200 * units.GPa],
               [15, 'linear', 245.0 * units.MPa, 490.0 * units.MPa, 200 * units.GPa]])
#
material = mesh.material()
print(material)
#

#
# ----------------------------------------------------
# Section Input
# ----------------------------------------------------
#
mesh.section([[15, 'rectangle', 8 * units.inch, 4 * units.inch],
              [20, 'ub', 240*units.mm, 6.2*units.mm, 120*units.mm, 9.8*units.mm],
              [25, 'Tubular', 300 * units.mm, 10 * units.mm]])
#
print(mesh.section())
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
# nodes coordinates [node_id, x, y, z=0]
#
mesh.node([(1, storeyBase,   storeyBase),
           (2, storeyBase, storeyHeight1),
           (3, bayWidth, storeyHeight2),
           (4, bayWidth,   storeyBase)])
#
print(mesh.node())
#
#
# ----------------------------------------------------
# boundary Input
# ----------------------------------------------------
#
# [node_id, type, fixity]
mesh.boundary([[1, 'support', 'pinned'],
               [4, 'support', 'fixed']])
#
print(mesh.boundary())
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
mesh.element([(1,  'beam',  1, 2, 10, 15, 0),
              (2,  'beam',  2, 3, 15, 20, 0),
              (3,  'beam',  3, 4, 10, 25, 0)])
#
print(mesh.element().beam())
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
# load.basic([[load_title, 'node', node_id, 'point',  x,y,z,mx,my,mz, comment(optional)],
#             [load_title, 'node', node_id, 'mass' ,  x,y,z, comment(optional)]
#             [load_title, 'beam', beam_id, 'line' ,  qx0,qy0,qz0, qx1,qy1,qz1, L0,L1, comment(optional)],
#             [load_title, 'beam', beam_id, 'point',  L0,x,y,z,mx,my,mz, comment(optional)]])
#
#
Pnull = 0 * units.N
Lnull= 0 * units.N / units.m
Punit = units.kN
Munit = units.kN * units.m
Lunit = units.N/units.m
#
#load.basic([['wind load', 'node', 2, 'load', 400 * Punit, Pnull, Pnull, 100 * Munit, 100 * Munit, 'nodex_1'],
#            ['wind load', 'node', 3, 'load', -200 * Punit, 'nodex_2'],
#            ['snow load', 'beam', 2, 'line', Lnull, -50_000 * Lunit, 20_000 * Lunit, 'udly_1'],
#            ['snow load', 'beam', 3, 'line', Lnull, -50_000 * Lunit, 20_000 * Lunit, 'udly_2'],
#            ['snow load', 'beam', 1, 'point', 3* units.m, 10 * Punit, 'point_1']])
#
#print(basic)
#
basic = load.basic()
#
#basic.node([['wind load', 2, 'load', 400 * Punit, Pnull, Pnull, 100 * Munit, 100 * Munit, 'nodex_3'],
#            ['wind load', 3, 'load', -200 * Punit, 'nodex_4']])
#
#
basic.beam([['snow load', 2, 'line', 1.5* units.m, Lnull, -50_000 * Lunit, 20_000 * Lunit, 'udly_3']])
#            ['snow load', 3, 'line', Lnull, -50_000 * Lunit, 20_000 * Lunit, 'udly_4'],
#            ['snow load', 1, 'point', 3* units.m, 10 * Punit, 'point_2']])
#
#
#basic.nodal([['gust load', 2, 'load', 400 * Punit, Pnull, Pnull, 100 * Munit, 100 * Munit, 'nodex_5'],
#             ['gust load', 3, 'load', -200 * Punit, 'nodex_6']])
#
#
#basic.beam([['car load', 2, 'line', Lnull, -50_000 * Lunit, 20_000 * Lunit, 'udly_5'],
#            ['car load', 3, 'line', Lnull, -50_000 * Lunit, 20_000 * Lunit, 'udly_6'],
#            ['car load', 1, 'point', 3* units.m, 10 * Punit, 'point_3']])
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
print("Nodes")
nodes = mesh.node()
nodedf = nodes.df
print(nodedf)
#
print("boundaries")
bds = mesh.boundary()
bdsdf = bds.df
print(bdsdf)
#
print("Elements")
elements = mesh.element()
elementsdf = elements.df
print(elements)
#
print("Load")
loadm = mesh.load()
basicLoad = loadm.basic()
bl_nodedf = basicLoad.node().df
print(bl_nodedf)
#
beamload = basicLoad.beam()
bl_bpointdf = beamload.point.df
print(bl_bpointdf)
bl_blinedf = beamload.line.df
print(bl_blinedf)
#
#
#bload = loadm.beams().df
#
# mesh.to_excel()
#
# ----------------------------------------------------
# Plotting
# ----------------------------------------------------
#
# Structure
#
#plot = mesh.plot()
#plot.frame()
#plot.material()
#plot.section()
#
# Loading
#
#plotload = load.plot()
#plotload.basic()
#
#
#
# ----------------------------------------------------
# Structural Analysis Implicit
# ----------------------------------------------------
#
frame = Trave3D(mesh=mesh)
frame.static()
results = frame.results()
#
noderes = results.nodes()
#print(noderes)
ndisp = noderes.displacement()
# get pandas df
ndisp_df = ndisp.df
nreacc = noderes.reaction()
# get individual node results
ndisp = noderes[1].displacement()
#
beamres = results.beam()
#print(beamres)
bdisp = beamres.displacement()
bforce = beamres.force()
bstress = beamres.stress()
# get individual beam results
bforce = beamres[1].displacement()
bdisp = beamres[2].force()
bstress = beamres[3].stress()
# get pandas df
bstress_df = bstress.df
#print(bstress)
#
print(results)
print('-->')