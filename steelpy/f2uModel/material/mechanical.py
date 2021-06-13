# 
# Copyright (c) 2009-2019 fem2ufo
# 
"""
defines mechanical material type: linear or non-linear material
"""

# Python stdlib imports
from array import array
from collections.abc import Mapping
from dataclasses import dataclass
from typing import NamedTuple, ClassVar, Dict, List, Tuple, Union, Iterable

# package imports
from steelpy.process.units.main import Units
import steelpy.f2uModel.material.print_report as print_material
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table
 
#
@dataclass
class GetMaterialSQL:
    """ Linear Material Data"""

    __slots__ = ['number', 'index', 'cls', 'type', 'db_file']

    def __init__(self, cls, material_number:int) -> None:
        """
        """
        self.index: int = cls._labels.index(material_number)
        self.cls: ClassVar = cls
        self.number: int = material_number
        #self.type = cls.type
        self.db_file = cls.db_file
        # get material name
        self.type = "elastic"
    #
    #
    @property
    def Fy(self)-> Units:
        return self.get_item(item="Fy")  * self.cls.f2u_units.Pa
    
    @Fy.setter
    def Fy(self, value:Units)-> None:
        self.update_item(item='Fy', value=value.convert("pascal").value)
    #
    @property
    def E(self):
        return self.get_item(item="E")  * self.cls.f2u_units.Pa

    @E.setter
    def E(self, value):
        self.update_item(item='E', value=value.convert("pascal").value)
    #
    @property
    def G(self):
        return self.get_item(item="G")  * self.cls.f2u_units.Pa

    @G.setter
    def G(self, value):
        self.update_item(item='G', value=value.convert("pascal").value)
    #
    @property
    def Fu(self):
        Fu = self.get_item(item="Fu")
        try:
            1 / Fu
            return Fu * self.cls.f2u_units.Pa
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
        return self.get_item(item="density") * self.cls.f2u_units.kg / self.cls.f2u_units.m**3

    @density.setter
    def density(self, value):
        self.update_item(item='density', value=value.value)
    #
    #
    @property
    def alpha(self):
        """"""
        return self.get_item(item="alpha") * f2u_units.K
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
    def name(self, name:Union[str, int]):
        """ """
        self.cls._title[self.index] = name    
    #
    def update_item(self, item:str, value:float):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            self._update_item(conn, self.number, item, value)
            conn.commit()
    #
    def _update_item(self, conn, name, item, value):
        """ """
        project = (value, name)
        sql = 'UPDATE tb_MatElastoPlastic SET {:} = ? WHERE material_number = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    #
    def get_item(self, item:str):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            value = self._get_item(conn, self.number, item)
        return value
    #
    def _get_item(self, conn, name, item):
        """ """
        project = (name,)
        sql = 'SELECT {:} FROM tb_MatElastoPlastic WHERE material_number = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
        record = cur.fetchone()
        return record[0]
    #
    #
    #def set_default(self):
    #    """ """
    #    #name = self.cls._cls._labels[]
    #    index = self.cls._number.index(self.number)
    #    self.cls._default = self.cls._labels[index]
    #
    def __str__(self) -> str:
        return print_material.print_plastic_material(self)    
#
class MaterialElasticSQL(Mapping):
    
    __slots__ = ['db_file', '_labels', '_default', 
                 'f2u_units', '_title', '_number']
    
    def __init__(self, db_file:str):
        """
        """
        #global f2u_units
        self.f2u_units = Units()
        #self._cls = cls
        self.db_file = db_file
        self._title: List[Union[str, int]] = []
        self._labels: array = array('I', [])
        self._number: array = array('I', [])
        # create node table
        self._create_table()
    #
    #@property
    #def type(self) -> str:
    #    """ Material type classification"""
    #    return  "elastic"
    #
    def __setitem__(self, material_number: int,
                    properties: List[float]) -> None:
        """
        """
        try:
            self._labels.index(material_number)
            raise Exception('    *** warning material {:} already exist'
                            .format(material_number))
        except ValueError:
            self._labels.append(material_number)
            self._number.append(material_number)
            self._title.append(properties.pop(0))
            #
            conn = create_connection(self.db_file)
            with conn:
                self._push_material(conn, material_number, properties)
                #conn.commit()
    #
    def __getitem__(self, material_number: int) -> Tuple:
        """
        """
        try:
            self._labels.index(material_number)
            #index = self._cls._number.index(material_number)
            #material_name = self._cls._labels[index]
            return GetMaterialSQL(self, material_number)
        except ValueError:
            raise IndexError('   *** material {:} does not exist'.format(material_number))
    #
    def _create_table(self) -> None:
        """ """
        _table_material = "CREATE TABLE IF NOT EXISTS tb_MatElastoPlastic(\
                            number INTEGER PRIMARY KEY,\
                            material_number INTEGER NOT NULL,\
                            Fy DECIMAL NOT NULL,\
                            Fu DECIMAL NOT NULL,\
                            E DECIMAL NOT NULL,\
                            G DECIMAL NOT NULL,\
                            poisson DECIMAL NOT NULL,\
                            density DECIMAL NOT NULL, \
                            alpha DECIMAL NOT NULL);"
        #
        conn = create_connection(self.db_file)
        create_table(conn, _table_material)
    #
    def _push_material(self, conn, material_number, properties):
        """
        """
        project = (material_number, *properties)
        sql = 'INSERT INTO  tb_MatElastoPlastic(material_number,\
                                                Fy, Fu, E, G, poisson, density , alpha)\
                                                VALUES(?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    def __len__(self) -> float:
        return len(self._labels)

    def __iter__(self):
        """
        """
        return iter(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._labels
    #
#
#
#
@dataclass
class GetMaterial:
    """ Linear Material Data"""

    __slots__ = ['number','index', 'cls', 'type']

    def __init__(self, cls, material_number:int) -> None:
        """
        """
        self.index: int = cls._labels.index(material_number)
        self.cls: ClassVar = cls
        self.number = material_number
        # get material name
        #self.name = material_name
        self.type = "elastic"

    #
    @property
    def Fy(self):
        return self.cls._Fy[self.index] * f2u_units.Pa

    
    @Fy.setter
    def Fy(self, item):
        self.cls._Fy[ self.index ] = item.convert("pascal").value

    #
    @property
    def E(self):
        return self.cls._E[ self.index ] * f2u_units.Pa

    @E.setter
    def E(self, item):
        self.cls._E[ self.index ] = item.convert("pascal").value

    #
    @property
    def G(self):
        return self.cls._G[ self.index ] * f2u_units.Pa

    @G.setter
    def G(self, item):
        self.cls._G[ self.index ] = item.convert("pascal").value
        #

    @property
    def Fu(self):
        try:
            1 / self.cls._Fu[ self.index ]
            return self.cls._Fu[ self.index ] * f2u_units.Pa
        except ZeroDivisionError:
            return self.cls._Fy[ self.index ] / 0.75 * f2u_units.Pa

    @Fu.setter
    def Fu(self, item):
        """
        Fu :
        """
        self.cls._Fu[ self.index ] = item.convert("pascal").value

    #
    @property
    def density(self):
        return self.cls._density[ self.index ] * f2u_units.kg / f2u_units.m ** 3

    @density.setter
    def density(self, item):
        self.cls._density[ self.index ] = item.value

    #
    #
    @property
    def alpha(self):
        """"""
        return self.cls._alpha[ self.index ] * f2u_units.K

    #
    @property
    def poisson(self):
        """"""
        return self.cls._poisson[ self.index ]

    #
    @property
    def name(self):
        """ """
        return self.cls._title[self.index]
    
    @name.setter
    def name(self, name:Union[str, int]):
        """ """
        self.cls._title[self.index] = name
    #
    def get_name(self, number, _grade=False):
        """
        """
        if not _grade:
            try:
                _grade = round ( self.cls._grade[ self.index ] )
            except TypeError:
                _grade = 'Zero'
        # in case material already named (i.e. sacs case)
        _name = '_' + self.cls._title[ self.index ]
        if '_MT' in _name:
            _name = ''
        #
        self.cls._title[ self.index ] = 'MT' + str ( number ).zfill ( 3 ) + '_' + str ( _grade ) + _name
    #
    def equal(self, other, grade=None):
        """
        """
        if not grade:
            grade = other.grade
        # TODO : check if material type is needed
        if self.cls._type[ self.index ] == other.type \
                and self.cls._E[ self.index ] == other.E \
                and self.cls._grade[ self.index ] == grade \
                and self.cls._poisson[ self.index ] == other.poisson \
                and self.cls._density[ self.index ] == other.density:
            return True
        else:
            return False
    #
    def print_properties(self):
        """ """
        return print_material.print_plastic_material(self)
        # for line in text:
        #    print(line.rstrip())
    #
    #def set_default(self):
    #    """ """
    #    self.cls._cls._default = self.cls._title[self.index]
    #
    def __str__(self) -> str:
        return print_material.print_plastic_material(self)
#
#
class MaterialElastic(Mapping):
    """
    Represents an istotrop linear elastic material used for FEM simulations
    
    Material
        |_ name
        |_ number
        |_ type
        |_ Emodulus
        |_ poisson
        |_ Fy
        |_ density
        |_ alpha
        |_ damping
        |_ stiffness [k1, k2, k3,..., kn]
        |_ spring    [sp[0],...]
        |_ Gmodulus (shear modulus)
    
    **Parameters**:  
      :number:  integer internal number 
      :name:  string node external name
      :coordinate:  list node coordinates 
    """
    __slots__ = ['_E', '_poisson', '_density', '_default',
                 '_Fy', '_Fu', '_damping', '_alpha', 'f2u_units',
                  '_grade', '_G', '_labels', '_number', '_title']
    
    def __init__(self)-> None:
        """
        """
        global f2u_units
        f2u_units = Units()
        #
        #self._cls = cls
        #self.type = material_type
        #
        self._title: List[Union[str, int]] = []
        self._labels: array = array('I', [])
        self._number: array = array('I', [])
        self._grade: List[str] = []
        self._Fy: array = array('f', [])
        self._Fu: array = array('f', [])
        self._poisson: array = array('f', [])
        self._density: array = array('f', [])
        self._E: array = array('f', [])
        self._G: array = array('f', [])
        self._alpha: array = array('f', [])
    #
    def __setitem__(self, material_number: int,
                    properties: List[Union[str, float]]) -> None:
        """
        """
        try:
            self._labels.index(material_number)
            raise Exception('    *** warning material {:} already exist'
                            .format(material_number))
        except ValueError:
            #material_number = next(self.get_number())
            self._labels.append(material_number)
            self._number.append(material_number)
            self._grade.append(-1)
            #
            self._title.append(properties.pop(0))
            # set properties with default when value missing
            self._Fy.append(properties[0]) # Pa
            self._Fu.append(properties[1]) # Pa
            self._E.append(properties[2])  # Pa
            self._G.append(properties[3])  # Pa
            self._poisson.append(properties[4])
            self._density.append(properties[5]) # kg/m^3
            self._alpha.append(properties[6])   # K
            #
            #index = self._labels.index(material_number)
            #self._Fy[index] = properties[0]
    #
    def __getitem__(self, material_number: int) -> Tuple:
        """
        """
        try:
            index = self._labels.index(material_number)
            #material_number = self._number[index]
            #index = self._cls._number.index(material_number)
            #material_name = self._cls._labels[index]            
            #_index = self._labels.index(material_number)
            return GetMaterial(self, material_number)
        except ValueError:
            raise IndexError('   *** material {:} does not exist'.format(material_number))        
    #
    def get_number(self, start:int=1)-> Iterable[int]:
        """
        """
        try:
            n = max(self._number) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
    #
    def __len__(self) -> float:
        return len(self._labels)

    def __iter__(self):
        """
        """
        return iter(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._labels
    #    
#
#
class Spring(NamedTuple):
    """
    """
    #number: int
    force: List
    displacement: List
#
@dataclass
class Curve:
    """
    """
    __slots__ = ['type', '_spring', 'cls']
    
    def __init__(self, cls):
        """
        """
        self.cls = cls
    #
    @property
    def spring(self):
        """
        """
        _force = [_spring[1] for _spring in self._spring]
        _disp = [_spring[0] for _spring in self._spring]
        
        return Spring(force=_force, 
                      displacement=_disp)
    
    @spring.setter
    def spring(self, data):
        """
        """
        #self._force.append(data[0])
        #self._displacement.append(data[1])
        self._spring = data
#