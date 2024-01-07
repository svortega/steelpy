from steelpy import Spreadsheet
import pandas as pd


# Pandas

mydict = [{'a': 1, 'b': 2, 'c': 3, 'd': 4},
          {'a': 100, 'b': 200, 'c': 300, 'd': 400},
          {'a': 1000, 'b': 2000, 'c': 3000, 'd': 4000 }]

#df = pd.DataFrame(mydict)

df = pd.read_excel("Clair Caisson C3 Test Data HW v0 25-Aug-2020.xlsx",
                   sheet_name="Caisson Supports")

print(df)

test1 = df.iloc[0]
print(test1)

test2 = df.iloc[[0]]
#print(test2)
#
for row in df.itertuples():
    print( f"{row}")
#
print('----')
for key, row in df.iterrows():
    print(row)
#
#
#gby = df.groupby(['a', 'd'])
gby = df.groupby(["Support Type", "SupportNo", "Support Fixity"])
gbyg = gby.groups
print(gbyg)

#
#
#ss = Spreadsheet()
##
## basic functionality
##
#from math import sin, pi
#ss.tools.update(sin=sin, pi=pi, len=len)
#ss['A1'] = 5
#ss['A2'] = '=A1*6'
#ss['A3'] = '=A2*7'
#print(ss['A1'], ss['A2'], ss['A3'])
##
#ss['B1'] = '=sin(pi/4)'
#print(ss['B1'])
##
## read existing spreadsheet
##
#wb = ss.read_book("Clair Caisson C3 Test Data HW v0 25-Aug-2020.xlsx")
#sheets = wb.sheet_names
#print(sheets)
##
#ws = wb.sheets["Caisson Supports"]
##
## get data as dataframe
##
#data = ws.dataframe()
#print(data)
#item = data['y']
#print(item.max())
#print(data['y'])
#index = item.idxmax()
##for index, row in data.iterrows():
##   print( f"{index}: {row['NodeNo']}")
##
#for index, row in data.itertuples():
#   print( f"{index}: {row}")
##
##print(data.columns)
#idx = len(data.columns)
##
#print(data.iloc[0])
#print(data.iloc[[0]])
##
##print(data.iloc[:3])
##print(data.iloc[:,:idx].astype('str'))
##
#print(data.tableize())
##
print('--')