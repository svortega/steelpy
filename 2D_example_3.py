# Import steelpy modules
from steelpy import Units
from steelpy import f2uModel
from steelpy import Trave2D
from steelpy import Trave3D
#
#
units = Units()
#
f2umodel = f2uModel()
#
#
# ----------------------------------------------------
# ----------------------------------------------------
# Mesh 
# ----------------------------------------------------
# ----------------------------------------------------
#
mesh = f2umodel.mesh(sql_file="2D_example_2")
#
#
#
# ----------------------------------------------------
# mesh data
# ----------------------------------------------------
#
#
mesh.build()
#
print(mesh.materials())
#
print(mesh.sections())
#
nodes = mesh.nodes()
print(nodes)
#
bds = mesh.boundaries()
print("boundaries")
print(bds)
#
print("")
elements = mesh.elements()
print(elements)
#
loadm = mesh.load().basic()
print("Load")
print(loadm)
#
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
# Structural Analysis
# ----------------------------------------------------
#
frame = Trave2D(mesh=mesh)
# FIXME: calcs 
frame.static()
results = frame.results()
print(results)
print('-->')