# 
# Copyright (c) 2009 steelpy
# 
#
# Python stdlib imports
from __future__ import annotations
#from array import array
from collections.abc import Mapping
from dataclasses import dataclass
from typing import NamedTuple
import re
#

# package imports
from steelpy.utils.units.main import Units
from steelpy.utils.sqlite.main import ClassBasicSQL
#from ..process.operations import
from ..process.print_report import print_isomat
from steelpy.utils.sqlite.utils import create_connection, create_table
#
from ..process.operations import get_isomat_prop, get_isomat_prop_df
from ..process.mechanical import MaterialItem 
#
from steelpy.utils.dataframe.main import DBframework
#
#
# ----------------------------------------
#
class MatBasicSQL(ClassBasicSQL):
    __slots__ = ['db_file']
    
    def __init__(self, db_file: str) -> None:
        """
        db_system: sqlite
        """
        super().__init__(db_file)
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._component, )
        table = 'SELECT name FROM Material\
                 WHERE component_id = ? ;'
        #
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]
    #
    #
    #
    #
    #@property
    #def _default(self):
    #    """ """
    #    print('here')
    #
    #def get_number(self, start:int=1)-> Iterable[int]:
    #    """
    #    """
    #    try:
    #        n = max(self._number) + 1
    #    except ValueError:
    #        n = start
    #    #
    #    while True:
    #        yield n
    #        n += 1
    #     
#
#
# ----------------------------------------
#
class MaterialSQL(MatBasicSQL):
    
    __slots__ = ['db_file', 'db_system', '_component', 
                 '_default', '_elastic', '_curve']

    def __init__(self, db_file: str,
                 component: str|int = None,
                 db_system:str="sqlite") -> None:
        """
        db_system: sqlite
        """
        super().__init__(db_file)
        #
        self._component = component
        #self.db_system = db_system
        self._default: str|None = None
        self._elastic = MaterialElasticSQL(component=self._component,
                                           db_file=self.db_file)
    #    
    #
    def __setitem__(self, material_name: str|int,
                    properties: list[str|float]) -> None:
        """
        """
        try:
            self._labels.index(material_name)
            raise IOError('   error material {:} already exist'.format(material_name))
        except ValueError:
            material_type = properties[0]
            #
            conn = create_connection(self.db_file)
            with conn:
                push_material(conn, material_name, material_type, self._component)
            #
            # set material
            if re.match(r"\b(curve)\b", material_type, re.IGNORECASE):
                #self._material[material_name]
                raise NotImplementedError('material type not implemented')
            elif re.match(r"\b(elastic|linear|isotropic)\b", material_type, re.IGNORECASE):
                self._elastic[material_name] = properties[1:] # [mat_number, *properties[1:]]
            else:
                raise IOError(f' material type {material_type} not valid')
        #
        self._default = material_name
    #
    def __getitem__(self, material_name: str|int):
        """
        """
        try:
            index = self._labels.index(material_name)
            self._default = material_name
        except (IndexError, ValueError):
            raise KeyError(f'Invalid material: {material_name}')
        #
        conn = create_connection(self.db_file)
        with conn:        
            mat_number, mat_type = get_mat_specs(conn, material_name,
                                                 self._component)
        #
        if re.match(r"\b(curve)\b", mat_type, re.IGNORECASE):
            return self._curve[material_name]
        
        elif re.match(r"\b(elastic|linear|isotropic)\b", mat_type, re.IGNORECASE):
            return self._elastic[material_name]
        
        else:
            raise IOError(' material type {:} not recognised'
                          .format(mat_type))        
    #
    #
    # ------------------
    # SQL ops
    # ------------------    
    #  
    #
    def _new_table(self, conn) -> None:
        """ """
        table = "CREATE TABLE IF NOT EXISTS Material (\
                    number INTEGER PRIMARY KEY,\
                    name NOT NULL,\
                    component_id INTEGER NOT NULL REFERENCES Component(number), \
                    type TEXT NOT NULL,\
                    title TEXT);"
        #
        create_table(conn, table)
    #
    # ------------------
    # Opretations
    # ------------------     
    #
    @property
    def default(self):
        """ """
        return self._default

    @default.setter
    def default(self, material_name):
        """ """
        if not material_name in self._labels:
            raise IOError ( f'material {material_name} missing' )
        self._default = material_name
    #
    #
    #@elastic.setter
    def elastic(self, values:None|list=None,
                df=None):
        """ """
        if values:
            if isinstance(values, list):
                for item in values:
                    material_name = item[0]
                    #material_type = "elastic"
                    #self._labels.append(material_name)
                    #self._type.append(material_type)
                    #mat_number = next(self.get_number())
                    #self._number.append(mat_number)
                    #self.__setitem__(item[0], item[1:])
                    properties = get_isomat_prop(item[1:])
                    self._elastic[material_name] = properties # [material_name, *properties]
            
            else:
                raise IOError('material input not valid')
        # dataframe input
        try:
            df.columns
            df = df.drop_duplicates(['name'], keep='first')
            #
            group = df.groupby("type")
            elastic = group.get_group("elastic")
            elastic = get_isomat_prop_df(elastic)
            #elastic = elastic.drop_duplicates(['name'])
            #
            #self._labels.extend(elastic.name)
            #self._type.extend(elastic.type)
            #material_id = [next(self.get_number()) for _ in elastic.name]
            #self._number.extend(material_id)
            #
            self._elastic.df = elastic
        except AttributeError:
            pass
        #print('--')
        #return self._elastic
        return MaterialType(mat_type="elastic",
                            cls_type=self._elastic, cls=self)
#
#
#
#
#
@dataclass
class MaterialType:
    __slots__ = ['_cls_type', '_item_type',
                 'db_file', '_component']
    
    def __init__(self, mat_type: str, cls_type, cls):
        """
        """
        self._cls_type = cls_type
        self._item_type = mat_type
        #self._labels = cls._labels
        #self._type = cls._type
        self._component = cls._component
        self.db_file = cls.db_file
        
    def __setitem__(self, material_name:str|int,
                    properties:list[str|float]) -> None:
        """
        """
        #self._labels.append(material_name)
        #self._type.append(self._item_type)
        # set material master
        conn = create_connection(self.db_file)
        with conn:
            push_material(conn, material_name,
                          self._item_type,
                          component=self._component)
        #
        #self._number.append(mat_number)
        # set material type
        prop = get_isomat_prop(properties)
        self._cls_type[material_name] = prop # [mat_number, *prop]
    
    def __getitem__(self, material_name:str|int):
        """
        """       
        return self._cls_type[material_name]
    #
    @property
    def df(self):
        """ """
        return self._cls_type.df
    #
    #def get_number(self, start:int=1)-> Iterable[int]:
    #    """
    #    """
    #    try:
    #        n = max(self._number) + 1
    #    except ValueError:
    #        n = start
    #    #
    #    while True:
    #        yield n
    #        n += 1
    #   
#
#  
#
# 
#
def push_material(conn, material_name: str|int, 
                  material_type:str, component: int):
    """
    """
    project = (material_name, component, material_type, None, )
    sql = 'INSERT INTO  Material(name, component_id, type, title)\
           VALUES(?,?,?,?)'
    cur = conn.cursor ()
    cur.execute (sql, project)
    #return cur.lastrowid
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
def get_materials(conn, component: int):
    """
    """
    query = (component, )
    table = "SELECT Material.name, Material.type,\
                        MatElastoPlastic.*\
                FROM Material, MatElastoPlastic, Component \
                WHERE  materials.number = MatElastoPlastic.number \
                AND Component.number = ?;"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    rows = cur.fetchall()
    #
    materials = {}
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
def get_materialSQL(conn, material_name:int|str,
                    component: int):
    """
    """
    query = (material_name, component, )
    table = "SELECT Material.name, Material.type, \
                    MaterialElastic.*, Component.name \
                FROM Component, Material, MaterialElastic\
                WHERE  Material.number = MaterialElastic.material_id\
                AND Material.name = ? \
                AND Component.number = ? ;"
    #
    with conn: 
        cur = conn.cursor()
        cur.execute(table, query)
        row = list(cur.fetchone())
    #
    #component = row.pop(0)
    #
    return MaterialItem(name=row[0], number=row[3], 
                        Fy=row[4], Fu=row[5],
                        E=row[6], G=row[7],
                        poisson=row[8], density=row[9],
                        alpha=row[10])
    #print("--->")
    #return materials
#
#
#
@dataclass
class GetMaterialSQL:
    """ Linear Material Data"""

    __slots__ = ['number', 'index', 'cls', 'type', 'db_file', 'units']

    def __init__(self, cls, material_id: int) -> None:
        """
        """
        1 / 0
        self.index: int = cls._labels.index(material_id)
        self.cls = cls
        self.number: int = material_id
        self.db_file: str = cls.db_file
        # get material name
        self.type: str = "elastic"
        self.units = Units()
    #
    #
    def __str__(self) -> str:
        return print_isomat(self)    
    #
    # -------------------------------------------
    #
    @property
    def Fy(self):
        return self.get_item(item="Fy") * self.units.Pa

    @Fy.setter
    def Fy(self, value) -> None:
        self.update_item(item='Fy', value=value.convert("pascal").value)

    #
    @property
    def E(self):
        return self.get_item(item="E") * self.units.Pa

    @E.setter
    def E(self, value):
        self.update_item(item='E', value=value.convert("pascal").value)

    #
    @property
    def G(self):
        return self.get_item(item="G") * self.units.Pa

    @G.setter
    def G(self, value):
        self.update_item(item='G', value=value.convert("pascal").value)

    #
    @property
    def Fu(self):
        Fu = self.get_item(item="Fu")
        try:
            1 / Fu
            return Fu * self.units.Pa
        except ZeroDivisionError:
            return self.Fy / 0.75

    @Fu.setter
    def Fu(self, value):
        """
        Fu :
        """
        self.update_item(item='Fu', value=value.convert("pascal").value)

    #
    @property
    def density(self):
        return self.get_item(item="density") * self.units.kg / self.units.m ** 3

    @density.setter
    def density(self, value):
        self.update_item(item='density', value=value.value)

    #
    #
    @property
    def alpha(self):
        """"""
        return self.get_item(item="alpha") * self.units.K

    #
    @property
    def poisson(self):
        """"""
        return self.get_item(item="poisson")

    #
    @property
    def name(self):
        """ """
        return self.cls._title[self.index]

    @name.setter
    def name(self, name: Union[str, int]):
        """ """
        self.cls._title[self.index] = name
    #
    # -------------------------------------------
    #
    def update_item(self, item: str, value: float):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            self._update_item(conn, self.number, item, value)
            conn.commit()

    #
    def _update_item(self, conn, number: int,
                     item: str, value: float):
        """ """
        query = (value, number)
        table = f'UPDATE MatElastoPlastic SET {item} = ? \
                  WHERE material_id = ? ;'
        #
        cur = conn.cursor()
        cur.execute(table, query)

    #
    #
    def get_item(self, item: str):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            value = self._get_item(conn, self.number, item)
        return value

    #
    def _get_item(self, conn, number: int, item: str):
        """ """
        query = (number,)
        table = f'SELECT {item} FROM MatElastoPlastic \
                  WHERE material_id = ?'
        #
        cur = conn.cursor()
        cur.execute(table, query)
        record = cur.fetchone()
        return record[0]

    #
    #
    # def set_default(self):
    #    """ """
    #    #name = self.cls._cls._labels[]
    #    index = self.cls._number.index(self.number)
    #    self.cls._default = self.cls._labels[index]
#
#
#
#
#
class MaterialElasticSQL(MatBasicSQL):
    __slots__ = ['db_file', '_default', '_component']

    def __init__(self, component: int, db_file: str):
        """
        """
        super().__init__(db_file)
        self._component = component
    #
    #
    def __setitem__(self, name: int|str,
                    properties: list[float]) -> None:
        """
        """
        try:
            self._labels.index(name)
            conn = create_connection(self.db_file)
            with conn:
                self._push_material(conn,
                                    material_name=name, 
                                    properties=properties)
                # conn.commit()            
        except ValueError:
            raise Exception(f' *** warning material {name} missing')
    #
    def __getitem__(self, material_name: int|str) -> tuple:
        """
        """
        try:
            index = self._labels.index(material_name)
            conn = create_connection(self.db_file)
            return get_materialSQL(conn, material_name,
                                   component=self._component)
        except ValueError:
            raise IndexError(f' *** material {material_name} not valid')

    #
    #
    def _new_table(self, conn) -> None:
        """ """
        table = "CREATE TABLE IF NOT EXISTS MaterialElastic(\
                    number INTEGER PRIMARY KEY,\
                    material_id INTEGER NOT NULL REFERENCES Material(number),\
                    Fy DECIMAL NOT NULL,\
                    Fu DECIMAL NOT NULL,\
                    E DECIMAL NOT NULL,\
                    G DECIMAL NOT NULL,\
                    poisson DECIMAL NOT NULL,\
                    density DECIMAL NOT NULL, \
                    alpha DECIMAL NOT NULL);"
        #
        conn = create_connection(self.db_file)
        create_table(conn, table)

    #
    def _push_material(self, conn, material_name: int|str,
                       properties):
        """
        """
        mat_number, mat_type = get_mat_specs(conn, material_name,
                                             self._component)
        query = (mat_number, *properties)
        table = 'INSERT INTO  MaterialElastic(material_id,\
                                            Fy, Fu, E, G, poisson, density , alpha)\
                                VALUES(?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(table, query)
        #return cur.lastrowid
    #
    #
    @property
    def df(self):
        """ raw data for dataframe"""
        query = (self._component, )
        table = "SELECT Material.name, Material.type, Material.title, \
                        MaterialElastic.*\
                        FROM Material, MaterialElastic, Component\
                        WHERE  Material.number = MaterialElastic.material_id \
                        AND Component.number = ? ;"
        #
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            rows = cur.fetchall()            
        
        #
        db = DBframework()
        header = ['name', 'type', 'title',
                  'number', 'material_id',
                  'Fy', 'Fu', 'E', 'G', 'poisson', 'density', 'alpha']
        matdf = db.DataFrame(data=rows, columns=header)        
        #
        header = ['name', 'type', 
                  'Fy', 'Fu', 'E', 'G', 'poisson', 'density', 'alpha', 'title']        
        return matdf[header]
    
    @df.setter
    def df(self, df):
        """ """
        #
        conn = create_connection(self.db_file)
        #
        # push material header
        dfmat = df[['name', 'type']].copy()
        dfmat['title'] = None
        dfmat['component_id'] = self._component
        with conn:
            dfmat.to_sql('Material', conn,
                         index_label=['name', 'component_id',
                                      'type', 'title'], 
                         if_exists='append', index=False)
        #
        #
        query = (self._component, )
        table = "SELECT * from Material \
                 WHERE Material.component_id = ?"
        #
        cur = conn.cursor()
        cur.execute(table, query)
        rows = cur.fetchall()
        mat_name = {item[1]: item[0] for item in rows}
        #
        # push elastic material
        dfel = df[['Fy', 'Fu', 'E', 'G',
                   'poisson', 'density' , 'alpha']].copy()
        dfel['material_id'] = [mat_name[item] for item in df.name]
        #dfel['material_id'] = dfel['material_id'].astype(int)
        with conn:
            dfel.to_sql('MatElastoPlastic', conn,
                        index_label=['material_id',
                                     'Fy', 'Fu', 'E', 'G',
                                     'poisson', 'density' , 'alpha'], 
                        if_exists='append', index=False)
        #
        # TODO: update matIM
        #self._number.extend(dfel['material_id'].tolist())
        #self._labels.extend(dfmat['name'].tolist())
        #print('--')
#
#
#
def get_mat_specs(conn, material_name:int|str, component: int):
    """
    """
    query = (material_name, component, )
    table = "SELECT number, type FROM Material \
             WHERE  name = ? AND component_id =?;"
    
    cur = conn.cursor()
    cur.execute(table, query)
    row = cur.fetchone()
    return list(row)
#
#