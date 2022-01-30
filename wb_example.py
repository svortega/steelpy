from steelpy import Spreadsheet

sp = Spreadsheet()
wb = sp.read_book("Clair Caisson C3 Test Data HW v0 25-Aug-2020.xlsx")
sheets = wb.sheet_names
print(sheets)
#
ws = wb.sheets["Caisson Supports"]
data = ws.dataframe()
print(data)
item = data['y']
print(item.max())
print(data['y'])
index = item.idxmax()
print('--')