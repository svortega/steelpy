#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
from typing import NamedTuple, Union
#

# package imports
#from steelpy.f2uModel.material.mechanical import GetMaterial

class GetMaterial(NamedTuple):
    """ Linear Material Data"""
    name:Union[int,str]
    type:str
    E:float
    Fy:float
    Fu:float
    G:float
    nu:float
    rho:float
    alpha:float
#
#
def get_materials(conn, component_name):
    """
    """
    materials = {}
    cur = conn.cursor()
    #cur.execute("SELECT * FROM tb_materials;")
    cur.execute("SELECT tb_Materials.name, tb_Materials.type,\
                tb_MatElastoPlastic.*\
                FROM tb_Materials, tb_MatElastoPlastic\
                WHERE  tb_materials.number = tb_MatElastoPlastic.number;")
    rows = cur.fetchall()
    for row in rows:
        #print(row)
        #if row[3] == "curve":
        #    pass
        #else:
        materials[row[0]] = GetMaterial(name=row[0], type=row[1],
                                        E=row[4], Fy=row[5], Fu=row[6], G=row[7],
                                        nu=row[8], rho=row[9], alpha=row[10])
        #    print("boundary : ",row[9])
    #conn.close()
    #print("--->")
    return materials    