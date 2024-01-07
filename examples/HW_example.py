print("start")


#import datetime

from steelpy import f2uModel
from steelpy import Trave3D
from steelpy import Units
#
#
# Pass Asset ID from user interface (python trigger)
_assetId = 17
# print('')
# print (_assetId)

# Pass Equipment ID from user interface (python trigger)
_equipmentId = "143"
# print('')
# print (_equipmentId)

# Pass metocean Criteria ID from user interface (python trigger)
_metoceancriteriaId = 112
# print('')
# print (_metoceancriteriaId)
#
###########################################################
# Start conceptual modelling
###########################################################
f2u_model = f2uModel(component=_equipmentId) # set name component , mesh_type="inmemory"
concept = f2u_model.concept # call conceptual model
# call units module
units = Units()
#
#
###########################################################
# Get data from database
###########################################################
#
#
#----------------------------------------------------------
# Equipment support data
df_support_data = [{"NodeNo": "SP143-1", "SupportNo": 209, "EquipmentID": 143, "x": 0, "y":  26.44  , "z":  0, "SupportType":   "DWS_Flange"  , "SupportFixity": "Fixed"},
                    {"NodeNo":  "SP143-2", "SupportNo": 210, "EquipmentID": 143, "x": 0, "y":  10	  , "z":  0, "SupportType":   "Guide"       , "SupportFixity": "Pinned"},
                    {"NodeNo":  "SP143-3", "SupportNo": 211, "EquipmentID": 143, "x": 0, "y":  -15.5  , "z":  0, "SupportType":   "Guide"       , "SupportFixity": "Pinned"},
                    {"NodeNo":  "SP143-4", "SupportNo": 212, "EquipmentID": 143, "x": 0, "y":  -43.1  , "z":  0, "SupportType":   "Guide"       , "SupportFixity": "Pinned"},
                    {"NodeNo":  "SP143-5", "SupportNo": 213, "EquipmentID": 143, "x": 0, "y":  -44.188, "z":  0, "SupportType":   "Termination" , "SupportFixity": "Free"}]

# Data read column designations - supports (nodes)
# [0]   df['NodeNo'] - generated for python
# [1]   df['SupportNo'] - unique database number
# [2]   df['EquipmentID']
# [3]   df['x'] - in-plane offset
# [4]   df['y'] - vertical (support) elevation
# [5]   df['z'] - out-of-plane offset
# [6]   df['SupportType']
# [7]   df['SupportFixity']
# [8]   df['SupportSNCurve']

# print (df_support_data)

# -----------------------------------
# Define boundary conditions
boundary = concept.boundary
for index, row in enumerate(df_support_data):
    boundary[row['NodeNo']] = [row['x']*units.m, row['y']*units.m, row['z']*units.m]
    boundary[row['NodeNo']].support = row['SupportFixity']

#print("Node boundary list")
#print(_boundary_list)
#print("")

#
#----------------------------------------------------------
# Equipment beam span data
# NOT USED
#df_beam_span_data = []

# Data read column designations - elements
# [0]   df['SpanNo'] - generated for python
# [1]   df['EquipmentID']
# [2]   df['x_Node_Start'] - node start element in-plane offset
# [2]   df['y_Node_Start'] - end 1 element elevation
# [2]   df['z_Node_Start'] - end 1 element out-of-plane offset
# [3]   df['x_Node_End'] - end 2 element in-plane offset
# [4]   df['y_Node_End'] - end 2 element elevation
# [5]   df['z_Node_End'] - end 2 element out-of-plane offset

# print (df_beam_span_data)

#_span_list = []
#
#for index, row in df_beam_span_data.iterrows():
#    _span_list.append([row['SpanNo'],
#                       [row['x_Node_Start'], row['y_Node_Start'], row['z_Node_Start']],
#                       [row['x_Node_End'], row['y_Node_End'], row['z_Node_End']]])
#
#print("Span list")
#print(_span_list)
#print("")

#
#----------------------------------------------------------
# Equipment section data - used to create steps
df_section_data = [{"SectionNo": "SECT143-1", "EquipmentID":  143, "x_Node_Start":  0, "y_Node_Start":  28.1 , "z_Node_Start":  0, "x_Node_End":  0, "y_Node_End":    11.5, "z_Node_End":  0, "OD":  762, "WT":  20, "Yield":  345},
                    {"SectionNo": "SECT143-2", "EquipmentID":  143, "x_Node_Start":  0, "y_Node_Start":  11.5 , "z_Node_Start":  0, "x_Node_End":  0, "y_Node_End":     8.5, "z_Node_End":  0, "OD":  782, "WT":  30, "Yield":  345},
                    {"SectionNo": "SECT143-3", "EquipmentID":  143, "x_Node_Start":  0, "y_Node_Start":  8.5  , "z_Node_Start":  0, "x_Node_End":  0, "y_Node_End":     -14, "z_Node_End":  0, "OD":  762, "WT":  20, "Yield":  345},
                    {"SectionNo": "SECT143-4", "EquipmentID":  143, "x_Node_Start":  0, "y_Node_Start":  -14  , "z_Node_Start":  0, "x_Node_End":  0, "y_Node_End":     -17, "z_Node_End":  0, "OD":  782, "WT":  30, "Yield":  345},
                    {"SectionNo": "SECT143-5", "EquipmentID":  143, "x_Node_Start":  0, "y_Node_Start":  -17  , "z_Node_Start":  0, "x_Node_End":  0, "y_Node_End":   -41.6, "z_Node_End":  0, "OD":  762, "WT":  20, "Yield":  345},
                    {"SectionNo": "SECT143-6", "EquipmentID":  143, "x_Node_Start":  0, "y_Node_Start":  -41.6, "z_Node_Start":  0, "x_Node_End":  0, "y_Node_End":   -44.2, "z_Node_End":  0, "OD":  782, "WT":  30, "Yield":  345}]

# Data read column designations - sections
# [0]   df['SectionNo']
# [1]   df['EquipmentID']
# [2]   df['x_Node_Start']
# [3]   df['y_Node_Start']
# [4]   df['z_Node_Start']
# [5]   df['x_Node_End']
# [6]   df['y_Node_End']
# [7]   df['z_Node_End']
# [8]   df['OD']
# [9]   df['WT']
# [10]   df['Yield']

#
# -----------------------------------
# define material
material = f2u_model.materials
# -----------------------------------
# Define sections
section = f2u_model.sections
# -----------------------------------
# beam modelling
beam = concept.beam
# -----------------------------------
# Start concept loading
load = concept.load
basic = load.basic # define basic load
# create new basic load
basic[1] = 'dummy load'
#
#_section_list = []
for index, row in enumerate(df_section_data):
    # set material
    mat_name = "MAT" + str(row['Yield'])
    try:
        material[mat_name] = ['elastic', row['Yield'] * units.MPa]
    except IOError:
        pass
    #
    # set section
    sec_name = "TUB" + str(row['OD']) + "x" + str(row['WT'])
    try:
        #section[sec_name] = ['Tubular', row['OD'] * units.mm, row['WT'] * units.mm]
        section[sec_name] = ['Tubular', row['OD'] * units.mm, row['WT'] * units.mm]
        #print('set :',sec_name )
    except IOError:
        #print('exist', sec_name)
        pass
    #
    # create beam concepts
    beam[row['SectionNo']] = [row['x_Node_Start'] * units.m, row['y_Node_Start'] * units.m, row['z_Node_Start'] * units.m], \
                             [row['x_Node_End'] * units.m, row['y_Node_End'] * units.m, row['z_Node_End'] * units.m]
    beam[row['SectionNo']].material = mat_name
    beam[row['SectionNo']].section = sec_name
    #
    # example basic load UDL X-axis global system
    basic[1].beam = beam[row['SectionNo']]
    basic[1].beam.line_load["UDL"] = {'qy': 1 * units.kN / units.m}
#
#print("Section list")
#print(_section_list)
#print("")
#
###########################################################
# Start Metocean modelling
###########################################################
#
#----------------------------------------------------------
# Metocean criteria data
df_metocean_criteria_data = []
# print(df_metocean_criteria_data)

# Data read column designations - metocean criteria
# [0]   df['CriteriaID']
# [1]   df['AssetID']
# [2]   df['MetoceanEvent']
# [3]   df['ReturnPeriod_Yrs']
# [4]   df['StormSurge_m']
# [5]   df['StormTide_m']
# [6]   df['WaveHeightHmax_m']
# [7]   df['WavePeriod_Sec']
# [8]   df['CrestElevation_m']
# [9]   df['SurgeTideReference']
# [10]   df['EWLCrestElevation_m']

_metocean_criteria_list = [{"CriteriaID": 112, "AssetID": 17, "MetoceanEvent": "Extreme Storm", "ReturnPeriod_Yrs": 100,
                            "StormSurge_m": 0.88,  "StormTide_m": 1.11, "WaveHeightHmax_m": 31.7, "WavePeriod_Sec": 19.3,
                            "CrestElevation_m": 19.85, "SurgeTideReference": "MSL", "EWLCrestElevation_m": "NULL"}]

for index, row in enumerate(df_metocean_criteria_data):
    _metocean_criteria_list.append([row['StormSurge_m'], row['StormTide_m'], row['WaveHeightHmax_m'], row['WaveHeightHmax_m'], row['WavePeriod_Sec']])

print("Metocean criteria")
#print(_metocean_criteria_list)
print("")
#
#
#
###########################################################
# Start Properties modelling
###########################################################
#
properties = f2u_model.properties
hydro = properties.hydrodynamic
mg = hydro.marine_growth
#
#----------------------------------------------------------
# Equipment marine growth profile data
df_marine_growth_profile_data = [{"MarineGrowthID": 89, "AssetID":  17, "TopElevation_m":    3, "BottomElevation_m":   -14, "Thickness_mm":  40},
                                {"MarineGrowthID": 90, "AssetID":  17, "TopElevation_m":  -14, "BottomElevation_m":   -39, "Thickness_mm":  50},
                                {"MarineGrowthID": 91, "AssetID":  17, "TopElevation_m":  -39, "BottomElevation_m":   -79, "Thickness_mm":  40},
                                {"MarineGrowthID": 92, "AssetID":  17, "TopElevation_m":  -79, "BottomElevation_m":  -140, "Thickness_mm":  20}]
#print(df_marine_growth_profile_data)

# Data read column designations - metocean criteria
# [0]   df['MarineGrowthID']
# [1]   df['AssetID']
# [2]   df['TopElevation_m']
# [3]   df['BottomElevation_m']
# [4]   df['Thickness_mm']
mg_list = []
_marine_growth_profile_list = []
mg['MG_1'] = "profile"
for index, row in enumerate(df_marine_growth_profile_data):
    #_marine_growth_profile_list.append([row['MarineGrowthID'],
    #                                    row['TopElevation_m'], 
    #                                    row['BottomElevation_m'],
    #                                    row['Thickness_mm']])
    #mg_list.append([row['MarineGrowthID'], row['TopElevation_m'], row['Thickness_mm']])
    mg_list.append([row['MarineGrowthID'], row['BottomElevation_m'], row['Thickness_mm']])
    #
    mg['MG_1'].level[row['MarineGrowthID']] = [row['TopElevation_m'] * units.m, 
                                               row['Thickness_mm'] * units.mm]
#
#mg_list = [item for item in _marine_growth_profile_list]
#
print("Marine growth")
#print(_marine_growth_profile_list)
print("")
#
#
#from vpython import vector, cylinder, box, ring
#member = {}
#for key, item in beam.items():
#    print(key)
#    rad = item.section.d.value * 0.5
#    length = item.length.value
#    points = item.connectivity
#    unitv = item.unit_vector
#    pre = vector(points[0].x.value, points[0].y.value, points[0].z.value)
#    #post = vector(points[1].x.value, points[1].y.value, points[1].z.value)
#    post = vector(*unitv)
#    member[key] = cylinder(pos=pre, axis=post, radius=rad, length=length)
#
#supp = {}
#for key, item in boundary.items():
#    point = vector(*item.point)
#    suptype = item.support
#    supp[key] = box(pos=point, length=1.0, height=1.0, width=1.0)
#    #supp[key] = ring(pos=point, axis=vector(0,1,0),
#    #                radius=0.5, thickness=0.1)
#
#
#
#
mesh = f2u_model.build()
#
#
print("Materials")
for key, mat in material.items():
    print(key, mat)
#
print("Sections")
for key, sec in section.items():
    print(key, sec)
#
nodes = mesh.nodes
print("Nodes")
for key, node in nodes.items():
    print(node)
#
print("")
elements = mesh.elements
print("Elements")
for key, element in elements.items():
    print(element)
#
frame = Trave3D()
frame.static(mesh=mesh)
results = frame.results()
print(results)
print('-->')
