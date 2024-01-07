from steelpy import f2uModel
from steelpy import Units
from steelpy import Trave3D
#
#
units = Units()
#
# -----------------------------------
# Start conceptual modelling
# -----------------------------------
f2u_model = f2uModel(component="EA_PDQ_Caisson_C1")
concept = f2u_model.concept()
#
# -----------------------------------
# define material
material = concept.materials()
material["MAT345"] = ['elastic', 345.0 * units.MPa]
#
#
# -----------------------------------
# Define sections
section = concept.sections()
# Caisson
section["T1350x40"] = ['Tubular', 1350 * units.mm, 40 * units.mm]
section["T1350x25"] = ['Tubular', 1350 * units.mm, 25 * units.mm]
section["T1350x15"] = ['Tubular', 1350 * units.mm, 15 * units.mm]
# stub member
section["T559x16"] = ['Tubular', 559 * units.mm, 15.9 * units.mm]
section["T610x16"] = ['Tubular', 610 * units.mm, 15.9 * units.mm]
# jacket
section["T800x30"] = ['Tubular', 800 * units.mm, 30 * units.mm]
section["T900x35"] = ['Tubular', 900 * units.mm, 35 * units.mm]
section["T100x25"] = ['Tubular', 100 * units.mm, 25 * units.mm]
#
# -----------------------------------
# define Elevations
elevation = concept.points()
elevation[1] = [0*units.m, 14.2*units.m, 0*units.m]
elevation[2] = [0*units.m, 8.0*units.m, 0*units.m]
elevation[3] = [0*units.m, -14.0*units.m, 0*units.m]
elevation[4] = [0*units.m, -39.0*units.m, 0*units.m]
elevation[5] = [0*units.m, -64.0*units.m, 0*units.m]
elevation[6] = [0*units.m, -74.0*units.m, 0*units.m]
#
point = concept.points()
#point[22] = [0*units.m, 8.0*units.m, -1.5*units.m]
#point[33] = [0*units.m, -14.0*units.m, -1.5*units.m]
#point[44] = [2.7*units.m, -39.0*units.m, 0*units.m]
#point[55] = [2.7*units.m, -64.0*units.m, 0*units.m]
#
point[22] = [-1.5*units.m, 8.0*units.m, 0*units.m]
point[33] = [-1.5*units.m, -14.0*units.m, 0*units.m]
point[44] = [0*units.m, -39.0*units.m, 2.7*units.m]
point[55] = [0*units.m, -64.0*units.m, 2.7*units.m]
#
# -----------------------------------
elements = concept.elements()
#
# Start beam modelling
beam = elements.beams()
# set material & section default
material.default = "MAT345"
#
# Define Caisson from bottom to top
#
section.default = "T1350x40"
beam["bm12"] = elevation[6], elevation[5]
beam["bm12"].step[1].length = 7.0 * units.m
beam["bm12"].step[1].section = section["T1350x40"]
#
section.default = "T1350x40"
beam["bm9"] = elevation[5], elevation[4]
beam["bm9"].step[1].length = 1.5 * units.m
beam["bm9"].step[1].section = section["T1350x40"]
beam["bm9"].step[2].length = 23.5 * units.m
beam["bm9"].step[2].section = section["T1350x40"]
#
beam["bm6"] = elevation[4], elevation[3]
beam["bm6"].step[1].length = 1.5 * units.m
beam["bm6"].step[1].section = section["T1350x40"]
beam["bm6"].step[2].length = 23.5 * units.m
beam["bm6"].step[2].section = section["T1350x40"]
#
beam["bm3"] = elevation[3], elevation[2]
beam["bm3"].step[1].length = 1.5 * units.m
beam["bm3"].step[1].section = section["T1350x40"]
beam["bm3"].step[2].length = 9.0 * units.m
beam["bm3"].step[2].section = section["T1350x40"]
beam["bm3"].step[3].length = 20.5 * units.m
beam["bm3"].step[3].section = section["T1350x40"]
#
section.default = "T1350x40"
beam["bm27"] = elevation[2], elevation[1]
beam["bm27"].step[1].length = 1.5 * units.m
beam["bm27"].step[1].section = section["T1350x40"]
#
# Stub members
section.default = "T1350x40"
beam["bm14"] = elevation[2], point[22]
section.default = "T1350x40"
beam["bm15"] = elevation[3], point[33]
beam["bm16"] = elevation[4], point[44]
beam["bm17"] = elevation[5], point[55]
#
# -----------------------------------
# Define boundary conditions
boundary = concept.boundaries()
supports = boundary.supports()
#
supports['pinned'] = 'pinned'
#
supports['pinned'].points = elevation[1]
supports['pinned'].points = point[22]
supports['pinned'].points = point[33]
supports['pinned'].points = point[44]
supports['pinned'].points = point[55] 
#
#
# -----------------------------------
# Start concept loading
#
load = concept.load()
# define basic load
basic = load.basic()
#
# create new basic load
basic[1] = 'wave_loading'
##
##basic[1].point = elevation[2]
##basic[1].point.load["point_22"] = {'fz': 5000 * units.kN}
##
#basic[1].beam = beam["bm12"]
#basic[1].beam.line_load["wave_1"] = {'qz': 10 * units.kN / units.m}
##
#basic[1].beam = beam["bm9"]
#basic[1].beam.line_load["wave_2"] = {'qz': 10 * units.kN / units.m}
##
#basic[1].beam = beam["bm6"]
#basic[1].beam.line_load["wave_3"] = {'qz': 10 * units.kN / units.m}
##
basic[1].beam = beam["bm3"]
basic[1].beam.line = {'qz': 10 * units.kN / units.m, 'name': "wind_4"}
##
#basic[1].beam = beam["bm27"]
#basic[1].beam.line_load["wave_5"] = {'qz': 10 * units.kN / units.m}
##
##basic[1].beam.line_load["snow_2"] = {'qy1': 0 * units.kN / units.m, # in plane triangular load
##                                     'qy2':-2 * units.kN / units.m} # from node to node
## trapezoidal out plane load
#basic[1].beam.global_system # reset load coord system global
#basic[1].beam.line_load["snow_3"] = {'qz1': 2 * units.kN / units.m, # start load value
#                                     'qz2': 4 * units.kN / units.m, # end load value
#                                     'd1': 0.5 * units.m, # load start 0.5m from node 1
#                                     'd2': 1.0 * units.m} # load end 1m from node 2
##
#basic[2] = 'point_loading'
##
#basic[2].point = elevation[1]
#basic[2].point.load["point_1"] = {'fz': 1000 * units.kN} # nodal load in plane
##
#basic[2].point = elevation[2]
#basic[2].point.load["point_2"] = {'fz': 1000 * units.kN}
##
#basic[2].point = elevation[3]
#basic[2].point.load["point_3"] = {'fz': 1000 * units.kN}
##
#basic[2].point = elevation[4]
#basic[2].point.load["point_4"] = {'fz': 1000 * units.kN}
##
#basic[2].point = elevation[5]
#basic[2].point.load["point_5"] = {'fz': 1000 * units.kN}
##
#basic[2].point = elevation[6]
#basic[2].point.load["point_6"] = {'fz': 1000 * units.kN}
##
##
#
##point[100] = [0*units.m, 9.5*units.m, 0*units.m]
##basic[2].point = point[100]
##basic[2].point.load["point_100"] = {'fz': 1000 * units.kN}
#
# ----------------------------------------------------
# Plotting
# ----------------------------------------------------
#
# Structure
#
plot = concept.plot()
plot.frame()
#plot.material()
#plot.section()
#
# Loading
#
#plotload = load.plot()
#plotload.basic()
#
#
# ----------------------------------------------------
# Meshing
# ----------------------------------------------------
#
mesh = concept.mesh()
#
print("Materials")
print(material)
#
print("Sections")
print(section)
#
nodes = mesh.nodes
print(nodes)
#
bds = mesh.boundaries
print("boundaries")
print(bds)
#
print("")
elements = mesh.elements
print(elements)
#
print("Load basic")
print(basic)
#
#
# ----------------------------------------------------
# Structural Analysis
# ----------------------------------------------------
#
#
frame = Trave3D(mesh=mesh)
frame.static()
results = frame.results()
print(results)
print('-->')
