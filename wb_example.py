from steelpy import Spreadsheet

ss = Spreadsheet()
#
# basic functionality
#
from math import sin, pi
ss.tools.update(sin=sin, pi=pi, len=len)
ss['A1'] = 5
ss['A2'] = '=A1*6'
ss['A3'] = '=A2*7'
print(ss['A1'], ss['A2'], ss['A3'])
#
ss['B1'] = '=sin(pi/4)'
print(ss['B1'])
#
# read existing spreadsheet
#
wb = ss.read_book("Clair Caisson C3 Test Data HW v0 25-Aug-2020.xlsx")
sheets = wb.sheet_names
print(sheets)
#
ws = wb.sheets["Caisson Supports"]
#
# get data as dataframe
#
data = ws.dataframe()
print(data)
item = data['y']
print(item.max())
print(data['y'])
index = item.idxmax()
print('--')