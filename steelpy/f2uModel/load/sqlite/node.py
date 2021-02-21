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
from steelpy.f2uModel.load.operations.actions import PointNode
from steelpy.f2uModel.load.operations.operations import NodeLoadMaster, get_nodal_load
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table

# ---------------------------------
#
#
class NodeLoadSQL:
    __slots__ = ['_cls', '_labels', '_system', '_title',
                 '_system_flag', '_load_number', '_index',
                 '_complex']
    
    def __init__(self, cls) -> None:
        """
        """
        self._cls = cls
        self._labels: List[Union[str, int]] = [ ]
        self._title: List[str] = []
        # 0-global/ 1-local
        self._system_flag:int = 0
        self._load_number: array = array("I", [])
        self._complex: array = array("I", [])
        self._system: array = array("I", [])
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
        load_number = self._cls._number[index]
        # set element load        
        self._labels.append(node_number)
        #self._title.append(load_tile)
        self._load_number.append(load_number)        
        self._system.append(self._system_flag)
        self._index = len(self._labels)-1
        self._complex.append(0)
        #
        point_load = get_nodal_load(point_load)
        # push to SQL
        bd_file = self._cls.bd_file
        conn = create_connection(bd_file)
        with conn:
            self._push_node_load(conn, node_number, point_load)
            conn.commit()
    #
    def __getitem__(self, node_number:int) -> List[Tuple]:
        """
        """
        1/0
        try:
            index = self._cls._index
            load_name = self._cls._labels[index]
            bd_file = self._cls.bd_file
            load_number = self._load_number.index(load_name)
            # get beam load
            conn = create_connection(bd_file)
            with conn:
                udl = self._get_beam_load(conn, node_number, load_number)
            return udl
        except ValueError:
            return []    
    #
    def _push_node_load(self, conn, node_number: int,
                        node_load: List[float]):
        """
        """
        load_name = self._load_number[self._index]
        load_number = get_basic_load_number(conn, load_name)        
        system = "global"
        project = (load_number, node_number, system, 
                   self._complex[self._index], 
                   *node_load,
                   'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL')
        #
        sql = 'INSERT INTO tb_LoadNode(load_number, node_name,\
                                         system, load_complex, \
                                        fx, fy, fz, mx, my, mz,\
                                        fxi, fyi, fzi, mxi, myi, mzi)\
                                    VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    def _get_beam_load(self, conn, node_number:int, load_name:int):
        """ """
        #print('--')
        cur = conn.cursor()
        # beam line load
        load_number = get_basic_load_number(conn, load_name)
        cur.execute("SELECT * FROM tb_LoadNode\
                    WHERE load_number = {:}\
                    AND element_name = {:};"
                    .format(load_number, node_number))
        rows = cur.fetchall()
        beam_line = []
        for row in rows:
            data = [*row[7:10], *row[14:17], row[6], row[13],
                    node_number, row[2], *row[4:6]]
            beam_line.append(PointNode._make(data))
        return beam_line        
    #
    def _create_table(self) -> None:
        """ """
        _table_node_load = "CREATE TABLE IF NOT EXISTS tb_LoadNode(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            load_number INTEGER NOT NULL REFERENCES tb_Load(number),\
                            node_name INTEGER NOT NULL REFERENCES tb_Nodes(name),\
                            system TEXT,\
                            load_complex INTEGER NOT NULL,\
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
        load_number = self._cls._labels[index]
        try:
            self._load_number.index(load_number)
            bd_file = self._cls.bd_file
            conn = create_connection(bd_file)
            for node_number in set_items:
                with conn:
                    sum_load = self.get_sum_column(conn, node_number, load_number)
                    items.append(PointNode(*sum_load, node_number,
                                            load_tile, 'global', 0))
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
    #
    def update(self, other) -> None:
        """
        """
        pl = other._point
        try:
            _name = pl.name
            _name
            1/0
        except AttributeError:
            pass        
        #print('----')
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