from steelpy import f2uModel
from steelpy import Units
from steelpy import Trave3D
#
#
#
# -----------------------------------
# Start conceptual modelling
# -----------------------------------
f2u_model = f2uModel(component="caisson_1", mesh_type="inmemory") # set name component
#                    , mesh_type="inmemory"
concept = f2u_model.concept # call conceptual model
# call units module
units = Units()
#
# -----------------------------------
# define material
material = f2u_model.materials
material["MAT45"] = 'elastic'
material["MAT45"].Fy = 345.0 * units.MPa
#print(material["MAT45"].Fy)
material["MAT45"].E = 205000.0 * units.MPa # optional
#print(material["MAT45"].E)
material["MAT45"].G = 77200.0 * units.MPa  # optional
#print(material["MAT45"].G)
#
# -----------------------------------
# Define sections
section = f2u_model.sections
#
section["TUB500"] = 'Tubular'
section["TUB500"].d = 500 * units.mm
section["TUB500"].t = 25 * units.mm
#
section["TUB400"] = ['Tubular', 400 * units.mm, 15 * units.mm]
#section["TUB400"].d = 400 * units.mm
#section["TUB400"].t = 15 * units.mm
#
# -----------------------------------
# define points
point = concept.points
point[3] = [4*units.m, 3*units.m]
point[4] = {"x":4*units.m, "y":0*units.m}
#print(point[4])
#
# -----------------------------------
# Start beam modelling
beam = concept.beam
# set material & section default
material.default = "MAT45"
section.default = "TUB500"
#
# define beam via coordinades with list = [x, y, z=0]start, [x, y, z=0]end
beam["bm1"] = [0*units.m, 0*units.m], [0*units.m, 6*units.m]
#
# define beam via coordinades with dict and concept Point = {x, y, z=0}start, point_end
beam["bm2"] = {"x":0*units.m, "y":6*units.m}, point[3]
# segmented beam
#beam["bm2"].step[1] = [1*units.m, 6*units.m]
beam["bm2"].step[1].length = 1.0 * units.m
beam["bm2"].step[1].section = section["TUB400"]
#
beam["bm2"].step[2] = [3*units.m, 6*units.m]
#beam["bm2"].step[2].length = 3.0 * units.m
#
# define beam via concept Points = point_start, point_end
beam["bm3"] = point[3], point[4]
#
# -----------------------------------
# Define boundary conditions
boundary = concept.boundary
# define boundary via coordinades with list = [x, y, z=0]
boundary["sp1"] = [0*units.m, 0*units.m]
boundary["sp1"].support = "fixed"
# define beam via concept Points = point
boundary["sp2"] = point[4]
boundary["sp2"].support = 'fixed'
#
#
# -----------------------------------
# Start concept loading
#
load = concept.load
# define basic load
basic = load.basic
# create new basic load
basic[1] = 'wind load'
basic[1].point = [0*units.m, 6*units.m]
basic[1].point.load["wind_1"] = {'fy': -2 * units.MN} # nodal load in plane
basic[1].point.load["wind_2"] = {'mz': -3 * units.MN*units.m} # bending moment in plane
#
# create new basic load
basic[1].point = point[3]
basic[1].point.load["wind_3"] = {'fz': -4 * units.MN} # nodal load out plane
#
# create new basic load
basic[2] = 'snow basic'
basic[2].beam = beam["bm2"]
basic[2].beam.local_system # set load coord system local
basic[2].beam.line_load["snow_1"] = {'qy': -1 * units.kN / units.m} # in plane udl from node to node

basic[2].beam.line_load["snow_2"] = {'qy1': 0 * units.kN / units.m, # in plane triangular load
                                     'qy2':-2 * units.kN / units.m} # from node to node
# trapezoidal out plane load
basic[2].beam.global_system # reset load coord system global
basic[2].beam.line_load["snow_3"] = {'qz1': 2 * units.kN / units.m, # start load value
                                     'qz2': 4 * units.kN / units.m, # end load value
                                     'd1': 0.5 * units.m, # load start 0.5m from node 1
                                     'd2': 1.0 * units.m} # load end 1m from node 2
#
# create new basic load
basic[3] = 'crane load'
basic[3].beam = beam["bm2"]
basic[3].beam.point_load["crane_1"] = {'fx':-100 * units.kN,  # beam point axial load
                                       'd1': 2.5 * units.m}   # 2.5m from node 1
#
# create new basic load
basic[4] = 'dead load'
basic[4].selfweight.y = -1
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
#print("Nodes")
print(nodes)
#for key, node in nodes.items():
#    print(node)
#
print("")
elements = f2u_model.mesh.elements
print("Elements")
print(elements)
#for key, element in elements.items():
#    print(element)
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
