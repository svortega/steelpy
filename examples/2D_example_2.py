# Import steelpy modules
from steelpy import Units
from steelpy import f2uModel
from steelpy import Trave3D
from steelpy import Spreadsheet
#
#
units = Units()
#
f2umodel = f2uModel()
#
# ----------------------------------------------------
# Data
#
ss = Spreadsheet()
wb = ss.read_book("trave3D.xlsx")
sheets = wb.sheet_names
print(sheets)
#
# ----------------------------------------------------
# Mesh 
# ----------------------------------------------------
# ----------------------------------------------------
#
mesh = f2umodel.mesh(name="2D_example_2")
#
# ----------------------------------------------------
# Material input
# ----------------------------------------------------
#
data = wb.sheets["Material"]
matdata = data.to_df()
mesh.materials(df=matdata)
print(mesh.materials())
#
#
# ----------------------------------------------------
# Section Input
# ----------------------------------------------------
#
data = wb.sheets["Sections"]
setdata = data.to_df()
mesh.sections(df=setdata)
print(mesh.sections())
#
#
#
#
# ----------------------------------------------------
# Node input
# ----------------------------------------------------
#
data = wb.sheets["Nodes"]
nodedata = data.to_df()
mesh.nodes(df=nodedata)
print(mesh.nodes())
#
#
# ----------------------------------------------------
# boundary Input
# ----------------------------------------------------
#
support =  nodedata.drop(["x", "y", "z"], axis=1)
support["type"] = "supports"
#
mesh.boundaries(df=support)
print(mesh.boundaries())
#
# ----------------------------------------------------
# Element input
# ----------------------------------------------------
#
data = wb.sheets["Elements"]
membdata = data.to_df()
mesh.elements(df=membdata)
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
basic = load.basic()
#
# Nodal load
data = wb.sheets["Node_Load"]
nodeldf = data.to_df()
basic.node(df=nodeldf)
nodaldf = basic.node().df
#
# beam load
data = wb.sheets["Beam_Load"]
beamldf = data.to_df()
basic.beam(df=beamldf)
#beamsdf = basic.beam().df
#
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
mesh.build()
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
# ----------------------------------------------------
# Plot mesh
# ----------------------------------------------------
#
plot = mesh.plot()
plot.frame()
#plot.material()
#
# Loading
#
#plotload = load.plot()
#plotload.basic()
#
# ----------------------------------------------------
# Structural Analysis
# ----------------------------------------------------
#
frame = Trave3D(mesh=mesh)
frame.static()
results = frame.results()
print(results)
print('-->')