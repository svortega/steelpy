#
from steelpy import Clamp
from steelpy import Units
#
#
# ========================================================
#                    Assigment 1 
#            Based on ISO design method
# ========================================================
# 
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
Clamp1 = Clamp()
units = Units()

Clamp1.name = "914Brace"
Clamp1.design_method = "ISO"
#Clamp1.design_method = "OTH88283"
# G - Grouted
# L - Lined
# M - Mechanical
Clamp1.clamp_type = "LINED"
# Operating/Extreme
Clamp1.design_condition = "OPERATING"
# end-to-end/member-addition/tubular-joint
Clamp1.clamp_application = "member-addition"
#
# --------------------
# 2.- General Data
#
#
# --------------------
# 3.- Clamp Geometrical Data
#
# Existing Tubular
#
chord = Clamp1.substrate
chord.section.D = 914.0 * units.mm
chord.section.t = 35.0 * units.mm
chord.theta = 90 * units.degrees
chord.material.Fy = 345.0 * units.MPa
chord.material.E = 205000.0 * units.MPa
#
#chord_load = chord.actions
#chord_load.Fx = 10 * units.kN
#chord_load.Fy =  0 * units.kN
#chord_load.Fz =  0 * units.kN
#chord_load.Mx = 64 * units.kN * units.m
#chord_load.Mx = 0.0 * units.kN * units.m
#chord_load.Mz = 0.0 * units.kN * units.m
#
# Clamp
# geometry
shell = Clamp1.shell
shell.Lc = 800.0 * units.mm
# Internal diameter (mm)
shell.Dci = 914.0 * units.mm
# Clamp Shell Enclosed Angle
shell.theta_C = 180 * units.degrees
# number of bolts rows
shell.Nbp = 4.0
#
# material
shell.material.Fy = 345.0 * units.MPa
shell.material.E = 205000.0 * units.MPa
#
# --------------------
# add-on member
riser_support = Clamp1.add_on_member
riser_support.Length = 1 * units.mm
riser_support.section.D = 250 * units.mm
riser_support.section.t = 15 * units.mm
# material
riser_support.material.Fy = 345.0 * units.MPa
riser_support.material.E = 205000.0 * units.MPa
# Angle between the added member and the substrate
riser_support.theta = 90 * units.degrees
#
riser_load = riser_support.actions
riser_load.Fx = 53. * units.kN
riser_load.Fy = 10. * units.kN
riser_load.Fz = 82. * units.kN
riser_load.Mx = 6.0 * units.kN * units.m
riser_load.My = 64. * units.kN * units.m
riser_load.Mz = 0.0 * units.kN * units.m
#
# --------------------
# 5.- Applied Loads
# Axial [-ve compression](N)
#F1 = 82.0 * units.kN
# Shear In Plane (N)
#F2 = 10.0 * units.kN
# Shear Out of Plane
#F3 = 53.0 * units.kN
# Torsion (N.mm)
#M1 = 0.0 * units.kN * units.m
# Bending In Plane (N.mm)
#M2 = 6.0 * units.kN * units.m
# Bending Out of Plane (N.mm)
#M3 = 64.0 * units.kN * units.m
#
# --------------------
# 6.- Bolt Data
#
# Bolt Data
bolt = Clamp1.bolt
bolt.name = "M20"
bolt.fp = 0.15
bolt.Lsb = 913.0 * units.mm
#Clamp1.stud_bolt.data(BoltID, fp, Lsb)
#
# --------------------
# 7.- Neoprene Data
#
neoprene = Clamp1.neoprene
neoprene.name = "IRHD65"
neoprene.Ln = 800.0 * units.mm
neoprene.tn = 20 * units.mm
neoprene.Mu = 0.20
neoprene.material.E = 5.850 * units.MPa
neoprene.material.G = 1.3700 * units.MPa
neoprene.Kn = 0.54
neoprene.E_alpha = 1210.0 * units.MPa
#Clamp1.neoprene.data(NeopreneID, Ln, tn, Mu, Eo, Gn, Kn, Ealpha)
#
# --------------------
# 8.- Stiffener Data
# pull-up gap
Hgap = 20.0 * units.mm
# --------------------
# 9.- Flange Plate Data
# flange thickness
#tf = 35.0 * units.mm
#
# Stiffener Data
stiffener = shell.stiffener
stiffener.Cf = 60.0 * units.mm
stiffener.Bsp = 180.0 * units.mm
stiffener.e = 130.0 * units.mm
stiffener.eBolt = 22.0 * units.mm
stiffener.h = 275.0 * units.mm
stiffener.tsp = 20.0 * units.mm
#stiffener.Bbp = 0 * units.mm
#
#
#
# ====> CALCULATIONS
#
# Print Results
Clamp1.print_results()
#
#