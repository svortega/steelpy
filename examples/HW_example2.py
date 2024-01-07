print("start")
#import datetime
from steelpy import Spreadsheet
from steelpy import UFOmodel
from steelpy import Metocean
from steelpy import Trave2D
from steelpy import Units
#
units = Units()
#
#
#
###########################################################
#
# Pass Equipment ID from user interface (python trigger)
equipId = "143"
#
###########################################################
# Start conceptual modelling
###########################################################
ufo_model = UFOmodel(name='HW_example2')
concept = ufo_model.concept()
concept[equipId] = 'HW_example2_143'
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
bc.rename(columns={"NodeNo": "name", "Support Fixity": "boundary"},
          inplace=True)
bc['type'] = "support"
bc["x"] *= units.m
bc["y"] *= units.m
bc["z"] *= units.m
#print(bc.tabulate())
#
# -----------------------------------
# Define boundary conditions
concept[equipId].boundaries(df=bc)
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
concept[equipId].materials(df=mat)
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
concept[equipId].sections(df=sect)
sect = concept[equipId].sections().df
#
# -----------------------------------
# beam concept modelling
#
memb = data[["NodeNo", "x_Node_Start", "y_Node_Start", "z_Node_Start",
             "x_Node_End", "y_Node_End", "z_Node_End"]].copy()
# rename
memb.columns = ["name", "coordx_1", "coordy_1", "coordz_1",
                "coordx_2", "coordy_2", "coordz_2"]
#
memb.insert(6, column="material", value=mat["name"])
memb.insert(6, column="section", value=sect["name"])
#
print(memb)
elements = concept[equipId].elements()
elements.beams(df=memb)

#
#
#
# ----------------------------------------------------
# Metocean 
# ----------------------------------------------------
#
# ----------------------------------------------------
# hydrodynamic parameters input
# 
meto = Metocean(name="HW_example2_met")
#
MetCriteria = meto.criteria()
MetCriteria[5] = 'test_1'
#
#
# ----------------------------------------------------
#
proph = meto.properties()
#
# Marine Growth
#
MGsheet = wb.sheets["Marine Growth"]
mg = MGsheet.to_df()
mg['TopElevation_m'] *= units.m
mg['BottomElevation_m'] *= units.m
mg['Thickness_mm'] *= units.mm
#
mg.rename(columns={"AssetID": "name",
                   "TopElevation_m": "elevation",
                   "Thickness_mm": "thickness"}, inplace=True)
#
MG = proph.MG(df=mg)
#
# ----------------------------------------------------
# Current
#
Csheets = wb.sheets["Current Profile"]
current = Csheets.to_df()
current['HeightFromSeaBed_m'] -= current['HeightFromSeaBed_m'].max()
current['HeightFromSeaBed_m'] *= units.m
current['CurrentSpeed_ms'] *= units.m / units.sec
#
current.rename(columns={"CriteriaID": "name",
                        "HeightFromSeaBed_m" : "elevation",
                        "CurrentSpeed_ms": "velocity"}, inplace=True)
current = current[['name', 'elevation', 'velocity']]
#
MetCriteria[5].current(df=current)
#
#
# ----------------------------------------------------
# Wave setup
#
Msheet = wb.sheets["Metocean Criteria"]
metsetup = Msheet.to_df()
print(metsetup)
metsetup['StormSurge_m'] *= units.m
metsetup['StormTide_m'] *= units.m
metsetup['WaveHeightHmax_m'] *= units.m
metsetup['StormTide_m'] *= units.m
metsetup['CrestElevation_m'] *= units.m
metsetup['WavePeriod_Sec'] *= units.sec
#
metsetup['WindSpeed_ms'] *= units.m / units.sec
#
#
wave_type = "regular"
wave_theory = 'Stokes'
grpmet = metsetup.groupby(['AssetID'])
for key, item in AssetID.items():
    waveitem = grpmet.get_group(name=key)[['CriteriaID',
                                           'WaveHeightHmax_m',
                                           'WavePeriod_Sec']]
    waveitem.rename(columns={"CriteriaID": "name",
                             "WaveHeightHmax_m": "Hw",
                             "WavePeriod_Sec": "Tw"}, inplace=True)
    waveitem['d'] = item
    waveitem['type'] = wave_type
    waveitem['theory'] = wave_theory
    MetCriteria[5].wave(df=waveitem)
#
# Metocean combination
#
metcomb = Msheet.to_df()
metcomb.drop(columns=['WaveHeightHmax_m', 'WavePeriod_Sec',
                      'ReturnPeriod_Yrs', 'SurgeTideReference',
                      'EWLCrestElevation_m', 'WindSpeed_ms'],
             inplace=True, axis=1)
#
grpcomb = metcomb.groupby(['AssetID'])
grpcurr = current.groupby(['name'])
grpmg = mg.groupby(['name'])
for row in asset.itertuples():
    waveitem = grpcomb.get_group(name=row.AssetID)
    criteriaID = waveitem['CriteriaID'].values.tolist()[0]
    #
    waveitem['wave_name'] = criteriaID
    waveitem['wave_direction'] = 0
    waveitem['CrestElevation_m'] *= units.m
    #
    try:
        WKF = waveitem['WaveKinematicsFactor'].values.tolist()
        # Wave kinematics factor
        wkrf = proph.WKF()
        wkrf[criteriaID] = ['constant', row.WaterDepth_m, WKF[0]]
        waveitem['WaveKinematicsFactor'] = criteriaID
    except KeyError:
        waveitem['WaveKinematicsFactor'] = None    
    #
    #
    # Current
    try:
        curritem = grpcurr.get_group(name=criteriaID)
        waveitem['current_name'] = criteriaID
        waveitem['current_direction'] = 0
        #waveitem['current_blockage'] = None
        waveitem['current_stretching'] = True
        #
        CBF = waveitem['CurrentBlockage'].values.tolist()
        cbf = proph.CBF()
        cbf[criteriaID] = ['constant', row.WaterDepth_m, CBF[0]]
        waveitem['CurrentBlockage'] = criteriaID
        #
    except KeyError:
        waveitem['current_name'] = None
    #
    # Wind
    # TODO: wind not implemented
    #waveitem['wind_name'] = criteriaID
    #waveitem['wind_direction'] = 0
    #
    # properties
    #
    # Marine Growth
    try:
        mgitem = grpmg.get_group(name=row.AssetID)
        #MG[row.AssetID]
        waveitem['MG'] = row.AssetID
    except KeyError:
        waveitem['MG'] = None
    #
    #
    cdcm = proph.CdCm()
    cdcm[criteriaID] = ['constant',row.WaterDepth_m, 0.70, 2.0]
    waveitem['CdCm'] = criteriaID
    #
    #
    # parameters
    #
    waveitem['design_load'] = 'max_BS'
    waveitem['buoyancy'] = False
    waveitem['criterion'] = 'local'
    #
    #
    # Setup input
    #
    waveitem.rename(columns={'CriteriaID': 'name',
                             'MetoceanEvent': 'title',
                             'CrestElevation_m': 'crest_elevation',
                             'WaveKinematicsFactor': 'wave_kinematics',
                             'CurrentBlockage': 'current_blockage',
                             'ConductorShielding': 'conductor_shielding'}, inplace=True)
    #
    # Input load
    #
    MetCriteria[5].condition(df=waveitem)
#
metcond = MetCriteria[5].condition()
#
MetCriteria.solve(surface_points=36,
                  depth_points=100)
#
#
#
#
# -----------------------------------
# Start concept loading
# -----------------------------------
#
load = concept.load()
#
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
#
# -----------------------------------
#
hydro = load.metocean(metcond)
#

#
#
#
#
# ----------------------------------------------------
# Meshing
# ----------------------------------------------------
#
#
mesh = concept.mesh()
#
#
# ----------------------------------------------------
# Get mesh in dataframe pandas
#
print("Nodes")
nodes = mesh.nodes()
nodedf = nodes.df
print(nodedf)
#
print("boundaries")
supports = mesh.boundaries()
supportsdf = supports.df
print(supportsdf)
#
print("Elements")
elements = mesh.elements()
elementsdf = elements.df
print(elements)
#
print("Load")
loadmesh = mesh.load()
basiclm = loadmesh.basic()
#
nodalblm = basiclm.node()
print("Nodal Load")
bl_nodedf = nodalblm.df
print(bl_nodedf)
#
beamblm = basiclm.beam()
print("Beam Point Load")
bl_bpointdf = beamblm.point.df
print(bl_bpointdf)
#
print("Beam Line Load")
bl_blinedf = beamblm.line.df
print(bl_blinedf)
#
# ----------------------------------------------------
# Export mesh to excel 
#mesh.to_excel()
#
#
# ----------------------------------------------------
# Plotting
# ----------------------------------------------------
#
#plot = mesh.plot()
#plot.frame()
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
frame = Trave2D(mesh=mesh)
frame.static()
results = frame.results()
#
#
noderes = results.nodes()
#print(noderes)
ndisp = noderes.displacement()
# get pandas df
ndisp_df = ndisp.df
nreacc = noderes.reaction()
# get indidual node results
ndisp = noderes[1].displacement()
#
beamres = results.beam()
#print(beamres)
bdisp = beamres.displacement()
bforce = beamres.force()
bstress = beamres.stress()
# get indidual beam results
bforce = beamres[1].displacement()
bdisp = beamres[2].force()
bstress = beamres[3].stress()
# get pandas df
bstress_df = bstress.df
#print(bstress)
#
#
print(results)
print('--')

