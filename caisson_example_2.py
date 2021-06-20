from steelpy import f2uModel
from steelpy import Units
from steelpy import Trave3D
#
#
#
# -----------------------------------
# Start conceptual modelling
# -----------------------------------
f2u_model = f2uModel(component="EA_PDQ_Caisson_C1", mesh_type="inmemory")
concept = f2u_model.concept
units = Units()
#
# -----------------------------------
# define material
material = f2u_model.materials
material["MAT345"] = 'elastic'
material["MAT345"].Fy = 345.0 * units.MPa
#
# -----------------------------------
# Define sections
section = f2u_model.sections
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
elevation = concept.points
elevation[1] = [0*units.m, 14.2*units.m, 0*units.m]
elevation[2] = [0*units.m, 8.0*units.m, 0*units.m]
elevation[3] = [0*units.m, -14.0*units.m, 0*units.m]
elevation[4] = [0*units.m, -39.0*units.m, 0*units.m]
elevation[5] = [0*units.m, -64.0*units.m, 0*units.m]
elevation[6] = [0*units.m, -74.0*units.m, 0*units.m]
#
# -----------------------------------
# Start beam modelling
beam = concept.beam
# set material & section default
material.default = "MAT345"
#
# Define Caisson from bottom to top
#
section.default = "T1350x40"
beam["bm12"] = elevation[6], elevation[5]
#
beam["bm9"] = elevation[5], elevation[4]
#
beam["bm6"] = elevation[4], elevation[3]
#
beam["bm3"] = elevation[3], elevation[2]
#
section.default = "T1350x25"
beam["bm27"] = elevation[2], elevation[1]
#
# -----------------------------------
# Define boundary conditions
boundary = concept.boundary
#
boundary["sp2"] = elevation[2]
boundary["sp2"].support = "fixed"
#
boundary["sp3"] = elevation[3]
boundary["sp3"].support = 'pinned'
#
boundary["sp4"] = elevation[4]
boundary["sp4"].support = 'pinned'
#
boundary["sp5"] = elevation[5] 
boundary["sp5"].support = 'fixed'
#
#
# -----------------------------------
# Start concept loading
#
load = concept.load
# define basic load
basic = load.basic
#
# create new basic load
basic[1] = 'wave_loading'
basic[1].beam = beam["bm12"]
basic[1].beam.line_load["wave_1"] = {'qz': 9 * units.kN / units.m}
#
basic[1].beam = beam["bm9"]
basic[1].beam.line_load["wave_2"] = {'qz': 9 * units.kN / units.m}
#
basic[1].beam = beam["bm6"]
basic[1].beam.line_load["wave_3"] = {'qz': 9 * units.kN / units.m}
#
basic[1].beam = beam["bm3"]
basic[1].beam.line_load["wave_4"] = {'qz': 9 * units.kN / units.m}
#
basic[1].beam = beam["bm27"]
basic[1].beam.line_load["wave_5"] = {'qz': 9 * units.kN / units.m}
#
#basic[1].beam.line_load["snow_2"] = {'qy1': 0 * units.kN / units.m, # in plane triangular load
#                                     'qy2':-2 * units.kN / units.m} # from node to node
## trapezoidal out plane load
#basic[1].beam.global_system # reset load coord system global
#basic[1].beam.line_load["snow_3"] = {'qz1': 2 * units.kN / units.m, # start load value
#                                     'qz2': 4 * units.kN / units.m, # end load value
#                                     'd1': 0.5 * units.m, # load start 0.5m from node 1
#                                     'd2': 1.0 * units.m} # load end 1m from node 2
#
#
#
#
f2u_model.build()
#
print("Materials")
print(material)
#
print("Sections")
print(section)
#
nodes = f2u_model.mesh.nodes
print(nodes)
#
bds = f2u_model.mesh.boundaries
print("boundaries")
print(bds)
#
print("")
elements = f2u_model.mesh.elements
print(elements)
#
print("Load")
print(f2u_model.load.basic)
#
#
#
frame = Trave3D()
frame.f2u_model = f2u_model
frame.run_static()
frame.print_results()
print('-->')
