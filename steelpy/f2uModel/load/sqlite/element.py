#
# Copyright (c) 2009-2021 fem2ufo
# 

# Python stdlib imports
from array import array
from collections.abc import Mapping
#from collections import defaultdict
#from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, Union, Iterable, Dict  


# package imports
from steelpy.f2uModel.load.operations.operations import get_beam_load
from steelpy.f2uModel.load.operations.element import (LineBeam, PointBeam,
                                                      get_beam_point_load,
                                                      BeamDistMaster)
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table

#
# ---------------------------------
#
class BeamDistributedSQL(BeamDistMaster):
    
    __slots__ = ['_cls', '_labels', '_system_flag']
    
    def __init__(self, cls) -> None:
        """
        """
        super().__init__()
        self._cls = cls
        self._load_number: array = array("I", [])
        # create node table
        self._create_table()
    #
    def __setitem__(self, element_number: Union[int, str], 
                    udl: Union[List[float], Dict[str,float]]) -> None:
        """
        """
        # get load data
        index = self._cls._index
        #load_tile = self._cls._title[index]
        load_name = self._cls._labels[index]
        #load_number = self._cls._number[index]
        # set element load
        self._labels.append(element_number)
        self._load_number.append(load_name)
        if isinstance(udl[-1], str):
            load_title = udl[-1]
            udl.pop()
        else:
            load_title = 'NULL'
        # update load
        udl = get_beam_load(udl)
        # push to SQL
        bd_file = self._cls.bd_file
        conn = create_connection(bd_file)
        with conn:
            self._push_beam_load(conn, element_number, load_title,
                                 self._system_flag, udl)
        #print("-->")
    #
    def __getitem__(self, element_number:int) -> List[Tuple]:
        """
        """
        #1/0
        #_index_list: List = [x for x, _item in enumerate(self._labels)
        #                     if _item == element_number]
        index = self._cls._index
        load_number = self._cls._labels[index]
        bd_file = self._cls.bd_file
        #try:
        #self._load_number.index(load_number)
        # get beam load
        conn = create_connection(bd_file)
        with conn:
            udl = self._get_beam_load(conn, element_number, load_number)
        return udl
        #except ValueError:
        #    return []
    #
    #
    def __len__(self) -> float:
        """ """
        index = self._cls._index
        load_name = self._cls._labels[index]
        return self._load_number.count(load_name)
        #return len(self._labels)

    def __contains__(self, value) -> bool:
        """ """
        index = self._cls._index
        load_name = self._cls._labels[index]
        indexes = [x for x, item in enumerate(self._load_number)
                   if item == load_name]
        items = [self._labels[x] for x in indexes]        
        return value in items

    def __iter__(self) -> Iterable:
        """
        """
        index = self._cls._index
        load_name = self._cls._labels[index]
        load_list = list(set(self._load_number))
        indexes = [x for x, item in enumerate(load_list)
                   if item == load_name]
        items = [self._labels[x] for x in indexes]        
        #items = list(set(self._labels))
        return iter(items)
    #
    #@property
    #def name(self) -> str:
    #    """
    #    """
    #    return self._title[self._index]
    #
    #@name.setter
    #def name(self, load_name:str) -> None:
    #    """
    #    """
    #    index = self._cls._index
    #    load_name = self._cls._labels[index]
    #    bd_file = self._cls.bd_file
    #    conn = create_connection(bd_file)
    #    load_number = get_basic_load_number(conn, load_name)
    #    cur = conn.cursor()
    #    cur.execute("UPDATE tb_LoadBeamLine\
    #                 SET title = {:}\
    #                 WHERE load_number = {:}"
    #                .format(load_name, load_number))
    #    #try:
    #    #    self._title[self._index] = load_name
    #    #except AttributeError:
    #    #    #self.load_name = load_name
    #   #     raise IndexError("load name not found")
    #
    #
    @property
    def coordinate_system(self):
        if self._system_flag != 0:
            return "local"
        return "global"
    
    @coordinate_system.setter
    def coordinate_system(self, system:Union[str,int]):
        """
        Coordinate system for load : global or local (member)
        """
        self._system_flag = 0
        if system in ['local', 'member', 1]:
            self._system_flag = 1
        
    #    
    #
    def _create_table(self) -> None:
        """ """
        _table_element_line = "CREATE TABLE IF NOT EXISTS tb_LoadBeamLine(\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                load_number INTEGER NOT NULL REFERENCES tb_Load(number),\
                                element_name INTEGER NOT NULL REFERENCES tb_Elements(name),\
                                title TEXT,\
                                system INTEGER NOT NULL,\
                                Lnode1 DECIMAL,\
                                qx1 DECIMAL,\
                                qy1 DECIMAL,\
                                qz1 DECIMAL,\
                                qx1i DECIMAL,\
                                qy1i DECIMAL,\
                                qz1i DECIMAL,\
                                Lnode2 DECIMAL,\
                                qx2 DECIMAL,\
                                qy2 DECIMAL,\
                                qz2 DECIMAL,\
                                qx2i DECIMAL,\
                                qy2i DECIMAL,\
                                qz2i DECIMAL);"
        #
        bd_file = self._cls.bd_file
        conn = create_connection(bd_file)
        create_table(conn, _table_element_line)
    #
    def _push_beam_load(self, conn, element_number:int,
                        load_title:str, load_system:int,
                        udl:List[float]):
        """ """
        #print("-->")
        index = self._cls._index
        load_name = self._cls._labels[index]        
        #load_name = self._load_number[self._index]
        load_number = get_basic_load_number(conn, load_name)
        project = (load_number, element_number,
                   load_title, load_system,
                   udl[6], *udl[:3],
                   'NULL', 'NULL', 'NULL',
                   udl[7], *udl[3:6],
                   'NULL', 'NULL', 'NULL')
        #
        sql = 'INSERT INTO tb_LoadBeamLine(load_number, element_name,\
                                            title, system,\
                                            Lnode1, qx1, qy1, qz1, qx1i, qy1i, qz1i,\
                                            Lnode2, qx2, qy2, qz2, qx2i, qy2i, qz2i)\
                                            VALUES(?,?,?,?,\
                                                   ?,?,?,?,?,?,?,\
                                                   ?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)      
    #
    #
    def _get_beam_load(self, conn, beam_number:int, load_name:int):
        """ """
        cur = conn.cursor()
        # beam line load
        load_number = get_basic_load_number(conn, load_name)
        #
        cur.execute("SELECT tb_Load.title, tb_LoadBeamLine.*\
                    FROM tb_LoadBeamLine, tb_Load\
                    WHERE tb_LoadBeamLine.load_number = {:}\
                    AND tb_LoadBeamLine.load_number = tb_Load.number\
                    AND tb_LoadBeamLine.element_name = {:};"
                    .format(load_number, beam_number))
        rows = cur.fetchall()
        beam_line = [] # defaultdict(list)
        for row in rows:
            data = [*row[7:10], *row[14:17], row[6], row[13],
                    beam_number, row[4], row[5], 0]
            beam_line.append(LineBeam._make(data))
        return beam_line
    #
    #
    #
    #def get_nodal_load(self, elements, materials, sections) -> List:
    #    """
    #    """
    #    items = line2node(self, elements, materials, sections)
    #    return items
#
#
def get_basic_load_number(conn, load_name:int):
    """ """
    cur = conn.cursor()
    cur.execute("SELECT * FROM tb_Load\
                 WHERE name = {:} \
                 AND type = 'basic'".format(load_name))
    loads = cur.fetchone()
    return loads[0]
#
#
class BeamPointSQL(Mapping):
    __slots__ = ['_cls', '_labels', '_system_flag', '_load_number']

    def __init__(self, cls) -> None:
        """
        """
        self._cls = cls
        self._labels: List[Union[str, int]] = []
        self._load_number: array = array ("I", [ ])
        # 0-global/ 1-local
        self._system_flag: int = 0
        # create node table
        self._create_table()
    #
    def __setitem__(self, element_number:int,
                    point_load: Union[List[float], Dict[str,float]]) -> None:
        """
        """
        # get load data
        index = self._cls._index
        #load_title = self._cls._title[index]
        load_name = self._cls._labels[index]
        #load_number = self._cls._number[index]
        # set element load
        self._labels.append(element_number)
        if isinstance(point_load[-1], str):
            load_title = point_load[-1]
            point_load.pop()
        else:
            load_title = 'NULL'
        #self._title.append(load_title)
        self._load_number.append(load_name)
        # update load
        point_load = get_beam_point_load(point_load)
        # push to SQL
        bd_file = self._cls.bd_file
        conn = create_connection(bd_file)
        with conn:
            self._push_beam_load(conn, element_number,
                                 load_title, self._system_flag, point_load)
        # print("-->")
    #
    def __getitem__(self, element_number:int)-> List[Tuple]:
        """
        """
        index = self._cls._index
        load_number = self._cls._labels[index]
        bd_file = self._cls.bd_file
        # get beam load
        conn = create_connection(bd_file)
        with conn:
            udl = self._get_beam_load(conn, element_number, load_number)
        return udl
    #
    #
    def __len__(self) -> float:
        """ """
        index = self._cls._index
        load_name = self._cls._labels[index]
        return self._load_number.count(load_name)

    def __contains__(self, value) -> bool:
        """ """
        index = self._cls._index
        load_name = self._cls._labels[index]
        indexes = [x for x, item in enumerate(self._load_number)
                   if item == load_name]
        items = [self._labels[x] for x in indexes]        
        return value in items

    def __iter__(self) -> Iterable:
        """
        """
        index = self._cls._index
        load_name = self._cls._labels[index]
        load_list = list(set(self._load_number))
        indexes = [x for x, item in enumerate(load_list)
                   if item == load_name]
        items = [self._labels[x] for x in indexes]
        return iter(items)
    #    
    #
    def _create_table(self) -> None:
        """ """
        _table_element_point = "CREATE TABLE IF NOT EXISTS tb_LoadBeamPoint(\
                                    number INTEGER PRIMARY KEY NOT NULL,\
                                    load_number INTEGER NOT NULL REFERENCES tb_Load(number),\
                                    element_name INTEGER NOT NULL REFERENCES tb_Elements(name),\
                                    title TEXT,\
                                    system INTEGER NOT NULL,\
                                    Lnode1 DECIMAL,\
                                    fx DECIMAL,\
                                    fy DECIMAL,\
                                    fz DECIMAL,\
                                    mx DECIMAL,\
                                    my DECIMAL,\
                                    mz DECIMAL,\
                                    fxi DECIMAL,\
                                    fyi DECIMAL,\
                                    fzi DECIMAL,\
                                    mxi DECIMAL,\
                                    myi DECIMAL,\
                                    mzi DECIMAL);"
        #
        bd_file = self._cls.bd_file
        conn = create_connection(bd_file)
        create_table(conn, _table_element_point)
    #
    def _push_beam_load(self, conn, element_number:int,
                        load_title: str, load_system: int,
                        point_load:List[float]):
        """ """
        #print("-->")
        index = self._cls._index
        load_name = self._cls._labels[index]
        #load_name = self._load_number[self._index]
        load_number = get_basic_load_number(conn, load_name)
        project = (load_number,element_number,
                   load_title, load_system,
                   point_load[6], *point_load[:6],
                   'NULL', 'NULL', 'NULL',
                   'NULL', 'NULL', 'NULL')
        #
        sql = 'INSERT INTO tb_LoadBeamPoint(load_number, element_name,\
                                            title, system, \
                                            Lnode1, fx, fy, fz, mx, my, mz,\
                                            fxi, fyi, fzi, mxi, myi, mzi)\
                                            VALUES(?,?,?,?,?,\
                                                   ?,?,?,?,?,?,?,\
                                                   ?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    #
    def _get_beam_load(self, conn, beam_number:int, load_name:int):
        """ """
        cur = conn.cursor()
        # beam line load
        load_number = get_basic_load_number(conn, load_name)
        cur.execute("SELECT tb_Load.title, tb_LoadBeamPoint.* \
                    FROM tb_LoadBeamPoint, tb_Load\
                    WHERE tb_LoadBeamPoint.load_number = {:}\
                    AND tb_LoadBeamPoint.load_number = tb_Load.number\
                    AND tb_LoadBeamPoint.element_name = {:};"
                    .format(load_number, beam_number))
        rows = cur.fetchall()
        beam_line = []
        for row in rows:
            data = [*row[7:13], row[6], beam_number,
                    row[4], row[5], 0]
            beam_line.append(PointBeam._make(data))
        return beam_line
    #    
    #
    #def get_nodal_load(self, elements, materials, sections) -> List:
    #    """
    #    """
    #    items = point2node(self, elements, materials, sections)
    #    return  items
#   # 
#

