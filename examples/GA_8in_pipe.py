#
#import math
from steelpy import CodeCheck
from steelpy import Units
from steelpy import Materials
from steelpy import Sections
#
units = Units()
code = CodeCheck()
#
# ---------------------------
# Pipe section
#
Dp = 219.10 * units.mm
tp = 14.3 * units.mm
#
section = Sections()
section["TUB"] = ['Tubular', Dp, tp]
#
#
# -----------------------------------
# define material
SMYS = 448 * units.MPa
SMTS = 517 * units.MPa
E = 207000.0 * units.MPa
alpha = (1.170E-05 + 1/273.15 ) / units.K
#
material = Materials()
material["X65"] = ['elastic', SMYS, SMTS]
#
# --------------------------
#
# Op Pressure
Pe = 67.50 * units.bar
Pi = 93 * units.bar
#
#
PipeTemp = 0*(50 + 273.15) * units.K
#
#
# =====================
#
pipegl = code.pipe
pipegl.design_code("PD8010 Part 2")
pipegl.design_data(design_condition = "operating", 
                  design_method = 'stress')
#
# ----------------------
# Pipe details
#
pipegl.pipe_data(pipe_type = "Seabed", 
                 pipe_name="Exp-Flowline",
                 pipe_restraint=True,
                 pipe_history = False)
#
pipegl.section = section["TUB"]
pipegl.pipe_section(tcorr = 3*units.mm, 
                    tol = 0.08, fo = 0.025)
#
#
#
# ----------------------
# 4.- Material Data (N/mm2)
#
pipegl.material = material["X65"]
#
#
#pipegl.material_derate(material_type="CMn",
#                      derate_code="DNV",
#                      Tmax=60.0)
#
# ----------------------
# 6.- Design Pressure (N/mm2)
#
pipegl.design_pressure(Pi = Pi, Po = Pe)
#
# ----------------------
# 7.- Design Temperature in degrees Celsius (C)
#
pipegl.design_temperature(T1 = 0 * units.K, 
                          T2 = PipeTemp)
# Oil and petroleum products
pipegl.set_substance_category("B")
#
# ----------------------
# Perform Calcs & Print Resuls
#
pipegl.get_results()
#Pipe2.print_results()
#
#
#
#
print('end')