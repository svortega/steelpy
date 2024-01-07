#
import steelpy as steelpy
#
#
# ========================================================
#                    Assigment 3 
#            Based on ISO design method
# ========================================================
# 
Clamp1 = steelpy.design.clamp
units = Clamp1.units
#
# --------------------
# 1.- Selecting Desing 
#     Code/Method
#
# Select Design Method :
# OTH88283 - Joint Industry Repairs Research Project,
#            Grouted & Mechanical Strengthening and 
#            Repair of Tubular Steel Offshore Structures
#            Department of Energy 1988 [OTH 88283]
# ISO - ISO/FDIS 19902 2007 - Appendix A
#
Clamp1.design_method = "ISO"
#Clamp1.design_method= "OTH88283"
#
# --------------------
# 2.- General Data
#
# Clamp Name/ID
Clamp1.name = "Assigment_3_A"
#
# Select Clamp Type
# G - Grouted
# L - Lined
# M - Mechanical
Clamp1.clamp_type = "MECHANICAL"
#
# Select Clamp's Flage Plate Arrangement:
# Continuous (Flage)
# Discontinuous (Flange)
Clamp1.flange_type = "Discontinuous"
#
# Select Design Condition:
# Operating/Extreme
Clamp1.design_condition = "extreme"
#
# ISO Specific Data
#
# Select Application of the Clamp
# end-to-end/member-addition/tubular-joint
Clamp1.clamp_application = "member-addition"
#
#
# --------------------
# 3.- Clamp Geometrical Data
#
# Existing Tubular
#
# External diameter (mm)
D = 530.0 * units.mm
# Wall Thickness (mm)
tp = 12.0 * units.mm
#
# --------------------
# Clamp
# Length (mm)
Lc = 1500.0 * units.mm
# pull-up gap
Hgap = 20.0 * units.mm
# number of bolts per row
Nbp = 9.0
#Nbp = 0
# Actual Stiffener Spacing
Lst = 200 * units.mm
# Internal diameter (mm)
#Dci = 278
# shell thickness
#tc = 22
#
# --------------------
# 4.- Material Propertiies
#
# Existing Tubular
# Yield Strength Material (N/mm2)
Fytub = 345.0 * units.MPa
# Young Modulus (N/mm2)
Etub = 205000.0 * units.MPa
#
# Clamp
# Yield Strength Tubular (N/mm2)
Fy = 345.0 * units.MPa
# Young Modulus (N/mm2)
E = 205000.0 * units.MPa
#
# --------------------
# 5.- Applied Loads OTH88283
# Axial [-ve compression](N)
F1 = 250000.0 * units.N
# Shear In Plane (N)
F2 = 55000.0 * units.N
# Shear Out of Plane
F3 = 0 * units.N
# Torsion (N.mm)
M1 = 83000000. * units.N * units.mm
# Bending In Plane (N.mm)
M2 = 500000000.0 * units.N * units.mm
# Bending Out of Plane (N.mm)
M3 = 0 * units.N * units.mm
#
# --------------------
# 5.- Applied Loads ISO
# Axial [-ve compression](N)
F1 = 0 * units.N
# Shear In Plane (N)
F2 = 55000.0 * units.N
# Shear Out of Plane
F3 = 250000.0 * units.N
# Torsion (N.mm)
M1 = 0 * units.N * units.mm
# Bending In Plane (N.mm)
M2 = 500000000.0  * units.N * units.mm
# Bending Out of Plane (N.mm)
M3 =  83000000. * units.N * units.mm
#
#
# --------------------
# 6.- Bolt Data
#
shell = Clamp1.shell
bolt = shell.bolt
bolt.name = "M24"
bolt.fp = 0.80
#bolt.Lsb = 256.0 * units.mm
#
#
# --------------------
# 8.- Stiffener Data
# edge distance from centreline of bolts 
# to edge of flange plate (mm)
Cf = 80.0 * units.mm
# Bolt spacing (mm)
Bsp = 160.0 * units.mm
# e is the distance between the
# end bolt and the end of the
# clamp (mm)
e = 110.0 * units.mm
# Bolt eccentricity
eBolt = 27.0 * units.mm
# Height of stiffener
h = 150.0 * units.mm
# Total width of flange plate
bf = 172 * units.mm
# stiffener thickness (mm)
tst = 20.0 * units.mm
#
#
# --------------------
# 9.- Flange Plate Data
# flange thickness
tf = 35.0 * units.mm
tf = 0 * units.mm
#
#
# ----------------------
# 11.- Hydrostatic Pressure
#     based on Wave Data
#     (Optional)
#
# Wave Height (mm)
Hw = 11000 * units.mm
# Period (seconds)
T = 11.40
# Water Depth (mm)
d = 61000 * units.mm
# Depth of the member relative to still water level
# Measure positive upwards (mm)
z = -44000 * units.mm
# Wave Length according to Genie (mm)
#WaveLength = 238200
#
#
#
# ====> CALCULATIONS
#
#
# Geometric Data
# Existing Tubular
# Clamp1.tubular.geometry(D, tp)
# Basic Clamp Data
#Clamp1.clamp_geometric_data(Lc, Hgap, Nbp, Lst)
#Clamp1.SaddlePlate(Dci, tc)
#
# Material Properies
# Existing Tubular
chord = Clamp1.substrate
chord.section.D = D
chord.section.t = tp
chord.material.Fy = Fytub
chord.material.E = Etub
#
# Clamp
#Clamp1.clamp_material_data(Fy, E)
Clamp1.shell.material.Fy = Fy
Clamp1.shell.material.E = E
#
# Bolt Data
#Clamp1.stud_bolt.data(BoltID, fp, Lsb=0 * units.mm)
#
# Neoprene Data
# N/A
#
# Stiffener Data
#Clamp1.StiffenerData(Cf, Bsp, e, eBolt)
#Clamp1.FlangePlate(tf, bf)
#Clamp1.StiffenerPlate(h, bf)
#
# Include Hydrostatic Pressure
wave_data = Clamp1.hydrostatic_pressure_data
wave_data.Hw = Hw
wave_data.T = T
wave_data.d = d
wave_data.z = z
#Clamp1.hydrostatic_pressure_data(Hw, T, d, z)
#
# Applied Loads
Clamp1.applied_force.member_actions(F1=F1, F2=F2, F3=F3, M1=M1, M2=M2, M3=M3)
#
# Print Results
Clamp1.print_results()
#
#