#
import steelpy as steelpy
#
#
# ========================================================
#                    Assigment 1 
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
#
# --------------------
# 2.- General Data
#
# Clamp Name/ID
Clamp1.name = "273_Brace"
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
Clamp1.design_condition = "OPERATING"
#
# ISO Specific Data
#
# Select Application of the Clamp
# end-to-end/member-addition/tubular-joint
Clamp1.clamp_application = "member-addition"
#
# Angle between the added member and the substrate
# member (degrees)
#Theta = 90.0
#
# --------------------
# 3.- Clamp Geometrical Data
#
# Existing Tubular
#
#
chord = Clamp1.substrate
chord.section.D = 273.0 * units.mm
chord.section.t = 16.0 * units.mm
chord.material.Fy = 345.0 * units.MPa
chord.material.E = 205000.0 * units.MPa
#
# --------------------
# Clamp
#
shell = Clamp1.shell
shell.Lc = 700.0 * units.mm
# Internal diameter (mm)
shell.Dci = 278 * units.mm
shell.tc = 22 * units.mm
# Clamp Shell Enclosed Angle
#shell.theta_C = 180 * units.degrees
#shell.Lst = 160 * units.mm
#shell.Hgap = 10.0 * units.mm
# number of bolts rows
shell.Nbp = 4.0
#
# material
shell.material.Fy = 345.0 * units.MPa
shell.material.E = 205000.0 * units.MPa
#
#
Lst = 160 * units.mm
#Hgap = 10.0 * units.mm
#
# --------------------
# 6.- Bolt Data
#
bolt = shell.bolt
bolt.name = "M24"
bolt.fp = 0.80
bolt.Lsb = 256.0 * units.mm
#
#
# --------------------
# 8.- Stiffener Data
# edge distance from centreline of bolts 
# to edge of flange plate (mm)
#Cf = 80.0 * units.mm
# Bolt spacing (mm)
#Bsp = 160.0 * units.mm
# e is the distance between the
# end bolt and the end of the
# clamp (mm)
#e = 110.0 * units.mm
# Bolt eccentricity
#eBolt = 27.0 * units.mm
# Height of stiffener
#h = 150.0 * units.mm
# Total width of flange plate
#bf = 172 * units.mm
# stiffener thickness (mm)
#tst = 20.0 * units.mm
#
#
# --------------------
# 9.- Flange Plate Data
# flange thickness
#tf = 35.0 * units.mm
#tf = 0 * units.mm
#
# Stiffener Data
stiffener = shell.stiffener
stiffener.Cf = 80.0 * units.mm
stiffener.Bsp = 160.0 * units.mm
stiffener.e = 110.0 * units.mm
stiffener.eBolt = 27.0 * units.mm
stiffener.h = 150.0 * units.mm
stiffener.tsp = 20.0 * units.mm
stiffener.Bbp = 172 * units.mm
#
# Print Results
Clamp1.print_results()
#
#