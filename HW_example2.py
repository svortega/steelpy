print("start")
#import datetime
from steelpy import Spreadsheet
from steelpy import f2uModel
from steelpy import Metocean
from steelpy import Trave3D
from steelpy import Units
#
units = Units()
#
#
#
###########################################################
#
# Pass Equipment ID from user interface (python trigger)
equipmentId = "143"
#
###########################################################
# Start conceptual modelling
###########################################################
f2u_model = f2uModel(component=equipmentId)
concept = f2u_model.concept()
#
#
# -----------------------------------
# Read data from spreadsheet
ss = Spreadsheet()
wb = ss.read_book("Clair Caisson C3 Test Data HW v1 26-June-2023.xlsx")
sheets = wb.sheet_names
print(sheets)
#
# -----------------------------------
# Asset Data
#
ws = wb.sheets["Asset"]
asset = ws.to_df()
#print(asset)
asset['WaterDepth_m'] *= units.m
asset['DeckElevation_m'] *= units.m
AssetID = dict(zip(asset['AssetID'], asset['WaterDepth_m']))
#
# -----------------------------------
# Define boundary conditions
#
ws = wb.sheets["Caisson Supports"]
# get data as dataframe
data = ws.to_df()
#print(data.tabulate())
print(data)
# update data name
bc = data[["NodeNo", "x", "y", "z", "Support Fixity"]].copy()
bc.rename(columns={"NodeNo": "name", "Support Fixity": "support"},
          inplace=True)
#bc['units'] = "metre"
bc["x"] *= units.m
bc["y"] *= units.m
bc["z"] *= units.m
#print(bc.tabulate())
#
# -----------------------------------
# Define boundary conditions
concept.boundaries(df=bc)
#boundary = concept.boundary()
#boundary.df(bc)
#
#
# -----------------------------------
# Define Material
#
ws = wb.sheets["Caisson Sections"]
data = ws.to_df()
print(data)
#print(data.tabulate())
#
mat = data[["Yield"]].copy()
mat["type"] = "elastic"
mat["name"] = mat["Yield"].apply(lambda x: f"mat_{str(x)}")
mat["Yield"] *= units.MPa
mat.rename(columns={"Yield": "Fy"}, inplace=True)
#
# -----------------------------------
concept.materials(df=mat)
#
#
# -----------------------------------
# Define sections
#
sect = data[["OD", "WT"]].copy()
sect["type"] = 'tubular'
#print(sect.tabulate())
#def naming(row):
sect["name"] = sect.apply(lambda row: f"TUB_{str(row.OD)}x{str(row.WT)}", axis=1)
print(sect)
sect["OD"] *= units.mm
sect["WT"] *= units.mm
sect.rename(columns={"OD": "d", "WT": "tw"}, inplace=True)
#
concept.sections(df=sect)
#f2u_model.sections(df=sect)
#section = f2u_model.sections()
#section.df = sect
sect = concept.sections().df
#
# -----------------------------------
# beam concept modelling
#
memb = data[["NodeNo", "x_Node_Start", "y_Node_Start", "z_Node_Start",
             "x_Node_End", "y_Node_End", "z_Node_End"]].copy()
# rename
memb.columns = ["name", "coordx_1", "coordy_1", "coordz_1",
                "coordx_2", "coordy_2", "coordz_2"]
#memb['units'] = "metre"
#memb.join(sect["name"])
memb.insert(6, column="material", value=mat["name"])
memb.insert(6, column="section", value=sect["name"])
#
print(memb)
#
concept.beams(df=memb)
#beams = concept.beams()
#beams.df = memb
#
#
#
# -----------------------------------
# Start concept loading
#
load = concept.load()
# define basic load
basic = load.basic()
# create new basic load
#basic[1] = 'Gravity load example'
#basic[1].gravity = [0, -1* units.gravity, 0] #* units.gravity
#
# create new basic load
basic[2] = 'Point load example'
basic[2].point = [0* units.m, 5.0* units.m] #* units.m
basic[2].point.load = {'fx': -1 * units.MN, 'name': "deck_2"} # nodal load in plane
#
#
basic[3] = 'Beam load example'
basic[3].beam = 'SECT143-1'

basic[3].beam.point = {'fz': 2 * units.MN, 'L1': 3 * units.m,
                       'name': "deck_3"}

basic[3].beam.line = {'qx': -3 * units.kN/units.m, 'name': "wind_3"}
#
#
#basic[4] = 'wave load'
#1 / 0
# box [point1, point2, width, height]
#basic[4].box = [0* units.m, -42 * units.m], [0 * units.m, 2.0 * units.m]
#basic[4].box.load = {'qx': -1 * units.kN/units.m, 'name': "SW"}
#
# -----------------------------------
#
# define th load
#TH = load.time_history()
#TH['WiD'] = 'Wave in Deck'
#TH['WiD'].basic_load[2] = 1.0
#TH['WiD'].points = [[0*units.sec, 4*units.sec, 6*units.sec, 8*units.sec, 10*units.sec],
#                    [0, 0, 1, 0, 0]]
#
#
# -----------------------------------
#
#mass = load.mass()
#mass['sw'] = 'selfweight'
#mass['sw'].basic_load[1] = 1.0
#
#----------------------------------------------------------
# Metocean criteria data
#
#
Csheets = wb.sheets["Current Profile"]
current = Csheets.to_df()
current['HeightFromSeaBed_m'] *= units.m
current['CurrentSpeed_ms'] *= units.m / units.sec

#
MGsheet = wb.sheets["Marine Growth"]
mg = MGsheet.to_df()
mg['TopElevation_m'] *= units.m
mg['BottomElevation_m'] *= units.m
mg['Thickness_mm'] *= units.m
#
#
Msheet = wb.sheets["Metocean Criteria"]
met = Msheet.to_df()
print(met)
met['StormSurge_m'] *= units.m
met['StormTide_m'] *= units.m
met['WaveHeightHmax_m'] *= units.m
met['StormTide_m'] *= units.m
met['CrestElevation_m'] *= units.m
#
met['WindSpeed_ms'] *= units.m / units.sec
#
#meto = Metocean()
# Regular wave
#meto.regular_wave()
#
# -----------------------------------
#
mesh = f2u_model.build()
#
nodes = mesh.nodes()
print(nodes)
#
boundary = mesh.boundaries()
print(boundary)
#
elements = mesh.elements()
print(elements)
#
load = mesh.load()
print(load)
#
frame = Trave3D()
frame.mesh = mesh
results = frame.run_static()
#frame.nfreq()
#results = frame.run_dynamic(end_time=10, delta_t=0.10)
results.print()
print('--')

