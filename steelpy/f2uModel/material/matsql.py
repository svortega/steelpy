# 
# Copyright (c) 2009-2020 fem2ufo
# 


# Python stdlib imports
from collections.abc import Mapping
from typing import Dict, List, Union, NamedTuple
#

# package imports
import steelpy.f2uModel.material.operations as operations
from steelpy.f2uModel.material.mechanical import MaterialElasticSQL
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table
#
#
class MaterialSQL(Mapping):
    __slots__ = ['db_file', 'db_system', '_labels', #'_default',
                 '_number', '_type', '_elastic', '_curve']

    def __init__(self, db_file: str,
                 db_system:str="sqlite") -> None:
        """
        db_system: sqlite
        """
        self.db_file = db_file
        self.db_system = db_system
        #
        self._labels:List[Union[str,int]] = []
        self._number: List[int] = []
        self._type:List[Union[str,int]] = []
        #self._default:Union[str,None] = None
        # create node table
        self._create_table()
        self._elastic = MaterialElasticSQL(self.db_file)
    #
    def __setitem__(self, material_name: Union[str, int],
                    properties: List[Union[str, float]]) -> None:
        """
        """
        try:
            self._labels.index(material_name)
            raise IOError('   error material {:} already exist'.format(material_name))
        except ValueError:
            material_type = properties[0]
            self._labels.append(material_name)
            self._type.append(material_type)
            #
            conn = create_connection(self.db_file)
            with conn:
                mat_number = self._push_material(conn, material_name, material_type)
                #conn.commit()
                self._number.append(mat_number)
            # set material
            if 'curve' == material_type :
                #self._material[material_name]
                raise Exception('--> No ready')
            elif 'elastic' == material_type :
                self._elastic[mat_number] = [material_name, *properties[1:]]
            else:
                raise IOError(' material type {:} not recognised'
                              .format(material_type))
    #
    def __getitem__(self, material_name: Union[str, int]):
        """
        """
        try:
            index = self._labels.index(material_name)
            mat_number = self._number[index]
            material_type = self._type[index]
        except IndexError:
            raise KeyError('Invalid material name : {:}'.format(material_name))
        #
        #
        if 'curve' == material_type :
            return self._curve[mat_number]
        elif 'elastic' == material_type :
            return self._elastic[mat_number]
        else:
            raise IOError(' material type {:} not recognised'
                          .format(material_type))        
    #
    def _push_material(self, conn, material_name: Union[str, int],
                       material_type:str):
        """
        """
        project = (material_name, None, material_type)
        sql = 'INSERT INTO  tb_Materials(name, title, type) VALUES(?,?,?)'
        cur = conn.cursor ()
        cur.execute (sql, project)
        return cur.lastrowid

    #
    def _create_table(self) -> None:
        """ """
        _table_materials = "CREATE TABLE IF NOT EXISTS tb_Materials (\
                            number INTEGER PRIMARY KEY,\
                            name INTEGER NOT NULL,\
                            title TEXT,\
                            type TEXT NOT NULL);"
        #
        conn = create_connection(self.db_file)
        create_table(conn, _table_materials)
    #
    def __len__(self):
        return len(self._labels)

    def __iter__(self):
        """
        """
        return iter(self._labels)

    def __contains__(self, value):
        return value in self._labels
    #
    #
    #@property
    #def _default(self):
    #    """ """
    #    print('here')
    #
    #
#
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
#
#
#
def get_materialSQL(conn, material_number:int):
    """
    """
    cur = conn.cursor()
    cur.execute("SELECT tb_Materials.name, tb_Materials.type,\
                tb_MatElastoPlastic.*\
                FROM tb_Materials, tb_MatElastoPlastic\
                WHERE  tb_materials.number = tb_MatElastoPlastic.number\
                AND tb_materials.number = {:};".format(material_number))
    row = cur.fetchone()
    materials = GetMaterial(name=row[0], type=row[1],
                            E=row[6], Fy=row[4], Fu=row[5], G=row[7],
                            nu=row[8], rho=row[9], alpha=row[10])
    #print("--->")
    return materials