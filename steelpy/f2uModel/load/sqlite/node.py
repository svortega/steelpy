#
# Copyright (c) 2009-2021 fem2ufo
# 

# Python stdlib imports
from array import array
#from collections.abc import Mapping
#from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, Iterator, Dict, Iterable, ClassVar, Union
#import re

# package imports
from steelpy.f2uModel.load.operations.nodes import NodeLoadBasic, PointNode, get_nodal_load
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table

# ---------------------------------
#
#
class NodeLoadSQL(NodeLoadBasic):
    __slots__ = ['_cls', '_labels', '_load_number']
                 #'_system', '_title',
                 #'_system_flag', '_load_name', '_index',
                 #'_complex']
    
    def __init__(self, cls) -> None:
        """
        """
        super().__init__()
        self._cls = cls
        self._load_number: array = array("I", [])
        # create node table
        self._create_table()
    #
    def __setitem__(self, node_number:int,
                    point_load:List) -> None:
        """
        """
        # get load data
        index = self._cls._index
        #load_tile = self._cls._title[index]
        load_name = self._cls._labels[index]
        #load_number = self._cls._number[index]
        # set element load        
        self._labels.append(node_number)
        self._load_number.append(load_name)
        if isinstance(point_load[-1], str):
            load_title = point_load[-1]
            point_load.pop()
        else:
            load_title = 'NULL'
        #self._system.append(self._system_flag)
        #self._index = len(self._labels)-1
        #self._complex.append(0)
        #
        point_load = get_nodal_load(point_load)
        # push to SQL
        bd_file = self._cls.bd_file
        conn = create_connection(bd_file)
        with conn:
            self._push_node_load(conn, node_number, load_title, point_load)
            conn.commit()
    #
    def __getitem__(self, node_number:int) -> Union[List[Tuple],None]:
        """
        """
        try:
            index = self._cls._index
            load_name = self._cls._labels[index]
            #load_title = self._title[index]
            bd_file = self._cls.bd_file
            #load_number = self._load_number.index(load_name)
            conn = create_connection(bd_file)
            #load_number = get_basic_load_number(conn, load_name)
            # get beam load
            with conn:
                node_load = self._get_node_load(conn, node_number, load_name)
            return node_load
        except ValueError:
            return None
    #
    def _push_node_load(self, conn, node_number: int,
                        load_title:str, node_load: List[float]):
        """
        """
        index = self._cls._index
        load_name = self._cls._labels[index]
        #load_name = self._load_number[self._index]
        load_number = get_basic_load_number(conn, load_name)
        system = "global"
        project = (load_number, node_number, load_title,
                   system, *node_load,
                   'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL')
        #
        # print('-->')
        #
        sql = 'INSERT INTO tb_LoadNode(load_number, node_name, \
                                        title, system, \
                                        fx, fy, fz, mx, my, mz,\
                                        fxi, fyi, fzi, mxi, myi, mzi)\
                                    VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    def _get_node_load(self, conn, node_number:int, load_name:int):
        """ """
        #print('--')
        cur = conn.cursor()
        # beam line load
        load_number = get_basic_load_number(conn, load_name)
        cur.execute("SELECT * FROM tb_LoadNode\
                    WHERE load_number = {:}\
                    AND node_name = {:};"
                    .format(load_number, node_number))
        rows = cur.fetchall()
        node_load = []
        for row in rows:
            #data = [*row[7:10], *row[14:17], row[6], row[13],
            #        node_number, row[2], *row[4:6]]
            data = [*row[5:11], *row[2:4], 0, 0]
            node_load.append(PointNode._make(data))
        return node_load        
    #
    def _create_table(self) -> None:
        """ """
        _table_node_load = "CREATE TABLE IF NOT EXISTS tb_LoadNode(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            load_number INTEGER NOT NULL REFERENCES tb_Load(number),\
                            node_name INTEGER NOT NULL REFERENCES tb_Nodes(name),\
                            title TEXT,\
                            system TEXT,\
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
        create_table(conn, _table_node_load)
    #
    #
    #@property
    def get_group_nodes(self):
        """
        """
        items = []
        set_items = set(self._labels)
        index = self._cls._index
        load_tile = self._cls._title[index]
        load_name = self._cls._labels[index]
        try:
            self._load_number.index(load_name)
            bd_file = self._cls.bd_file
            conn = create_connection(bd_file)
            load_number = get_basic_load_number(conn, load_name)
            for node_number in set_items:
                with conn:
                    sum_load = self.get_sum_column(conn, node_number, load_number)
                    items.append(PointNode(*sum_load, node_number,
                                            load_tile, 0, 0))
            return items
        except ValueError:
            return items
    #
    #
    def get_sum_column(self, conn, node_name, load_number):
        """ """
        project = (node_name, load_number,)
        sql = 'SELECT SUM(fx), SUM(fy), SUM(fz), SUM(mx), SUM(my), SUM(mz)\
               FROM tb_LoadNode WHERE node_name = ? AND load_number = ?'
        cur = conn.cursor()
        cur.execute( sql, project )
        record = cur.fetchone()
        return record
    #
    #
    def update(self, other) -> None:
        """
        """
        pl = other._point
        try:
            title = pl.name
            for index, node_name in enumerate(pl._labels):
                point_load = [pl._fx[index], pl._fy[index], pl._fz[index],
                              pl._mx[index], pl._my[index], pl._mz[index],
                              title]
                self.__setitem__(node_name, point_load)
        except AttributeError:
            pass        
        #print('----')
    #
    def __len__(self) -> float:
        """ """
        index = self._cls._index
        load_name = self._cls._labels[index]
        return self._load_number.count(load_name)
    #
    def __contains__(self, value) -> bool:
        """ """
        index = self._cls._index
        load_name = self._cls._labels[index]
        indexes = [x for x, item in enumerate(self._load_number)
                   if item == load_name]
        items = [self._labels[x] for x in indexes]        
        return value in items
    #
    def __iter__(self)-> Iterable:
        """
        """
        index = self._cls._index
        load_name = self._cls._labels[index]
        indexes = [x for x, item in enumerate(self._load_number)
                 if item == load_name]
        items = [self._labels[x] for x in indexes]
        return iter(items)
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