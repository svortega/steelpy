from steelpy import UFOmodel
from steelpy import Units
from steelpy import Trave3D
#
#
units = Units()
#
# -----------------------------------
# Start conceptual modelling
# -----------------------------------
ufo_model = UFOmodel(name="FrameExample_1")
#
concept = ufo_model.concept()
concept[10] = 'Frame_Concept_1'
#
# -----------------------------------
# define material
material = concept[10].material({'name':'MAT250',
                                 'type':'linear',
                                 'Fy': 250.0 * units.MPa})
material["MAT345"] = ['elastic', 345.0 * units.MPa]
# set material & section default
material.default = "MAT345"
#
print(material.df)
#
# -----------------------------------
# Define sections
# Caisson
#section = concept[10].section()
concept[10].section({'type':'tubular',
                     'name':["T1350x40","T1350x25","T1350x15"],
                     'd':1350 * units.mm,
                     'tw':[40 * units.mm, 25 * units.mm, 15 * units.mm]})
# Caisson
#section["T1350x40"] = ['Tubular', 1350 * units.mm, 40 * units.mm]
#section["T1350x25"] = ['Tubular', 1350 * units.mm, 25 * units.mm]
#section["T1350x15"] = ['Tubular', 1350 * units.mm, 15 * units.mm]
# stub member
concept[10].section({'type':'tubular',
                     'name':["T559x16", "T610x16"],
                     'd':[559 * units.mm, 610 * units.mm],
                     'tw':[15.9 * units.mm, 15.9 * units.mm]})
#section["T559x16"] = ['Tubular', 559 * units.mm, 15.9 * units.mm]
#section["T610x16"] = ['Tubular', 610 * units.mm, 15.9 * units.mm]
#
concept[10].section({'type':'Circular bar',
                     'name': [20, 30],
                     'd': [100*units.mm, 200*units.mm]})
#
# jacket
section = concept[10].section()
section["T800x30"] = ['Tubular', 800 * units.mm, 30 * units.mm]
section["T900x35"] = ['Tubular', 900 * units.mm, 35 * units.mm]
section["T1000x25"] = ['Tubular', 1000 * units.mm, 25 * units.mm]
#
#section = [["T800x30", 'Tubular', 800 * units.mm, 30 * units.mm],
#          ["T900x35", 'Tubular', 900 * units.mm, 35 * units.mm],
#          ["T100x25", 'Tubular', 100 * units.mm, 25 * units.mm]]
#
# set section default
section.default = "T1350x40"
#
print(section.df)
#
# -----------------------------------
# define Elevations
elevation = concept[10].point()
#elevation[1] = [0*units.m, 14.2*units.m, 0*units.m]
elevation[2] = [0*units.m, 8.0*units.m, 0*units.m]
elevation[3] = [0*units.m, -14.0*units.m, 0*units.m]
#elevation[4] = [0*units.m, -39.0*units.m, 0*units.m]
#elevation[5] = [0*units.m, -64.0*units.m, 0*units.m]
#elevation[6] = [0*units.m, -74.0*units.m, 0*units.m]
#
point = concept[10].point()
point[22] = [0*units.m, 8.0*units.m, -1.5*units.m]
point[33] = [0*units.m, -14.0*units.m, -1.5*units.m]
#point[44] = [2.7*units.m, -39.0*units.m, 0*units.m]
#point[55] = [2.7*units.m, -64.0*units.m, 0*units.m]
#
#point[22] = [-1.5*units.m, 8.0*units.m, 0*units.m]
#point[33] = [-1.5*units.m, -14.0*units.m, 0*units.m]
#point[44] = [0*units.m, -39.0*units.m, 2.7*units.m]
#point[55] = [0*units.m, -64.0*units.m, 2.7*units.m]
#
#
# -----------------------------------
# Define boundary conditions
boundary = concept[10].boundary()
#supports = boundary.support()
spoint = boundary.point()
#
#supports['fixed'] = 'fixed'
#supports['pinned'] = 'pinned'
#
spoint["sp1"] = point[22]
spoint["sp1"].restrain = "fixed"
#
spoint["sp4"] = point[33]
spoint["sp4"].restrain = 'pinned'
#
# Test
#spoint["sp5"] = [0*units.m, 4.0*units.m, 0*units.m]
#spoint["sp5"].restrain = 'pinned'
#
#
# -----------------------------------
# Start beam modelling
elements = concept[10].element()
beam = elements.beam()
#
# Define Caisson from bottom to top
#
#beam["bm12"] = elevation[6], elevation[5]
#beam["bm12"].step[1].length = 7.0 * units.m
#beam["bm12"].step[1].section = section["T1350x40"]
##
#section.default = "T1350x40"
#beam["bm9"] = elevation[5], elevation[4]
#beam["bm9"].step[1].length = 1.5 * units.m
#beam["bm9"].step[1].section = section["T1350x40"]
#beam["bm9"].step[2].length = 23.5 * units.m
#beam["bm9"].step[2].section = section["T1350x40"]
##
#beam["bm6"] = elevation[4], elevation[3]
#beam["bm6"].step[1].length = 1.5 * units.m
#beam["bm6"].step[1].section = section["T1350x40"]
#beam["bm6"].step[2].length = 23.5 * units.m
#beam["bm6"].step[2].section = section["T1350x40"]
##
beam["bm3"] = elevation[3], elevation[2]
beam["bm3"].step[1].length = 1.5 * units.m
beam["bm3"].step[1].section = section["T1000x25"]
beam["bm3"].step[2].length = 9.0 * units.m
beam["bm3"].step[2].section = section["T900x35"]
beam["bm3"].step[3].length = 20.5 * units.m
beam["bm3"].step[3].section = section["T800x30"]
##
#section.default = "T1350x40"
#beam["bm27"] = elevation[2], elevation[1]
#beam["bm27"].step[1].length = 1.5 * units.m
#beam["bm27"].step[1].section = section["T1350x40"]
#
# Stub members
#section.default = "T559x16"
beam["bm14"] = elevation[2], point[22], material["MAT250"] , "T559x16" 
#
#section.default = "T610x16"
beam["bm15"] = elevation[3], point[33], section["T610x16"]
#beam["bm16"] = elevation[4], point[44]
#beam["bm17"] = elevation[5], point[55]
#
#
# -----------------------------------
# Start concept loading
#
load = concept[10].load()
# define basic load
basic = load.basic()
#
# create new basic load
basic[1] = 'wave_loading'
#
# Noda load
basic[1].point = elevation[2]
basic[1].point.load = {'fz': 5000 * units.kN, 'title': "point_F22"}
basic[1].point.load = {'fz': 5000 * units.kN, 'title': "point_F22"}
basic[1].point.load = {'mz': 100 * units.kN*units.m, 'title': "point_M22"}
#
# beam load 
basic[1].beam = beam["bm3"]
basic[1].beam.global_system()
# beam linear trapezoidal load
basic[1].beam.line = {'qy1': 2 * units.kN / units.m, # start load value
                      'qy2': 4 * units.kN / units.m, # end load value
                      'd1': 0.5 * units.m, # load start 0.5m from node 1
                      'd2': 1.0 * units.m, # load end 1m from node 2
                      'title': "snow_3"}
# beam point load
basic[1].beam.point = {'Fx': 10 * units.kN,
                       'L': 0.0 * units.m,
                       'title': 'snow_4'}
#
# new basic load
basic[3] = 'crane load'
basic[3].beam = beam["bm3"]
basic[3].beam.point = {'fx': -5000 * units.kN,  # beam point axial load
                       'd1': (1.5+3.75) * units.m,  # 2.5m from node 1
                       'title': "crane_1"}
#
print(basic)
#
#
#
#
#f2u_model.plot()
#
#
mesh = concept[10].mesh()
#
print("Materials")
print(material)
#
print("Sections")
print(section)
#
nodes = mesh.node()
print(nodes)
#
bds = mesh.boundary()
print("boundaries")
print(bds)
#
print("")
elements = mesh.element()
print(elements)
#
print("Load")
mload = mesh.load()
bload = mload.basic()
print(bload)
#
#
#mesh.plot()
#
#
#for load_name, lcase in bload.items():
#    for bname, loads in lcase.line_beam.items():
#        for item in loads:
#            item
#
#
#
#from steelpy.f2uModel.plot.main import PlotModel
#plotm = PlotModel()
#plotm.mesh = f2u_model.mesh
#plotm.load = f2u_model.load
#plotm.plot_mesh()
#f2u_model.plot.plot_mesh(verbosity=True)
#f2u_model.plot.basic_load(name=1)
#
frame = Trave3D(mesh=mesh)
frame.static()
results = frame.results()
print(results)
print('--> xx')
