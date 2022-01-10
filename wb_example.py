from steelpy import Spreadsheet

sp = Spreadsheet()
wb = sp.workbook("Clair Caisson C3 Test Data HW v0 25-Aug-2020.xlsx")
sheets = wb.sheet_names
#print(sheets)
#
ws = wb.sheets["Caisson Supports"]
data = ws.head()
print(data)
print(data['x'])
print('--')