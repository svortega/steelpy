#
# Copyright (c) 2009-2023 fem2ufo
# 

# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
#from collections import defaultdict
#from dataclasses import dataclass
#from typing import NamedTuple, Tuple, List, Union, Iterable, Dict
import re


# package imports
#from ....process.units.buckingham import Number
#from steelpy.f2uModel.load.operations.operations import get_beam_load
#from ..process.operations import (check_point_dic, check_list_units,
#                                  check_beam_dic,  
#                                  get_beam_node_load, get_beam_udl_load)
from steelpy.f2uModel.mesh.process.process_sql import  check_nodes
# steelpy.f2uModel.load
from ..process.beam import (LineBeam, PointBeam,
                            BeamDistMaster,
                            BeamLoadItem, BeamLoad)
from ..process.nodes import NodeLoadBasic, PointNode
# steelpy.
from steelpy.f2uModel.mesh.process.process_sql import create_connection, create_table, get_load_data
from steelpy.f2uModel.mesh.sqlite.beam import BeamItemSQL
#
from steelpy.process.dataframe.main import DBframework
#
# ---------------------------------
#
#
def check_element(conn, element_name):
    """ """
    cur = conn.cursor()
    cur.execute ("SELECT * FROM tb_Elements\
                WHERE tb_Elements.name = {:};".format(element_name))
    row = cur.fetchone()
    return row
#
#
class BeamLoadItemSQL(BeamLoadItem):
    __slots__ = ['_labels', '_load', '_bd_file', ]

    def __init__(self, load_name: int|float, bd_file: str) -> None:
        """
        """
        super().__init__()
        self._bd_file =  bd_file
        self._load = BeamLoadSQL(load_name=load_name,
                                     bd_file=self._bd_file)
        #
        self._node_eq = BeamToNodeSQL(load_name=load_name, 
                                        bd_file=bd_file)
    #
    def __setitem__(self, beam_name: int|str,
                    beam_load: list) -> None:
        """
        """
        conn = create_connection(self._bd_file)
        with conn:  
            beam = check_element(conn, beam_name)        
        try:
            memb_type = beam[3]
        except TypeError:
            raise IOError(f"beam {beam_name} not found")        
        #
        #
        self._beam_id = beam_name
        self._labels.append(beam_name)
        #
        if re.match(r"\b(point|node)\b", str(beam_load[0]), re.IGNORECASE):
            self._load._point[beam_name] = beam_load[1:]
            
        elif re.match(r"\b(line|udl)\b", str(beam_load[0]), re.IGNORECASE):
            self._load._line[beam_name] = beam_load[1:]
            
        else:
            raise IOError(f'Beam lod type {beam_load[0]} not implemented')

    #
    def __getitem__(self, beam_name: int | str):
        """
        """
        conn = create_connection(self._bd_file)
        with conn:  
            #beam = check_element(conn, beam_name)
            beam =  BeamItemSQL(beam_name, self._bd_file)
        try:
            memb_type = beam.type # beam[3]
            if memb_type != 'beam':
                raise ValueError(f"element {beam_name} type {memb_type} not valid")
        except TypeError:
            raise IOError(f"beam {beam_name} not found")
        #
        if not beam_name in self._labels:
            self._labels.append(beam_name)
        
        return self._load(beam=beam)
    #
    #
    def fer(self):
        """ Return Fix End Reactions (FER) global system"""
        #beams = self._f2u_beams
        for key in self._labels:
            #beam = beams[key]
            beam =  BeamItemSQL(key, self._bd_file)
            end_nodes = beam.connectivity
            res = self._load(beam=beam).fer()
            for gnload in res:
                self._node_eq[key] = [[end_nodes[0], *gnload[4], gnload[1], gnload[2]],
                                      [end_nodes[1], *gnload[5], gnload[1], gnload[2]]]
        #print('--> get_end_forces')
        #1 / 0    
#
#
class BeamLoadSQL(BeamLoad):
    __slots__ = ['_system_flag', #'_beam_id',
                 '_line', '_point', '_beam']

    def __init__(self, load_name: int|float, bd_file: str): #, beams
        """
        """
        super().__init__()
        self._line = BeamDistributedSQL(load_name=load_name, bd_file=bd_file)
        self._point = BeamPointSQL(load_name=load_name, bd_file=bd_file)
        #self._bd_file = bd_file
    #
    #
    def __call__(self, beam):
        """ """
        #self._beam_id = beam_name
        self._beam = beam # =  BeamItemSQL(beam_name, self._bd_file)
        return self    
    #
    #
#
#    
#
#
#
#
# ---------------------------------
#
#
class BeamDistributedSQL(BeamDistMaster):
    __slots__ = ['_labels', '_title', '_index', '_complex',
                 '_system', '_bd_file'] # '_system_flag', 
    
    def __init__(self, load_name: int|float, bd_file: str) -> None:
        """
        """
        super().__init__()
        self._name = load_name
        self._bd_file =  bd_file
        # create node table
        conn = create_connection(self._bd_file)
        with conn:        
            self._create_table(conn)
    #
    def __setitem__(self, beam_name: int|str, line_load: list) -> None:
        """
        """
        conn = create_connection(self._bd_file)
        with conn:  
            beam = check_element(conn, beam_name)        
        try:
            beam_number = beam[0] 
        except TypeError:
            raise IOError(f"beam {beam_name} not found")        
        # get load data
        # set element load
        self._labels.append(beam_name)
        title = line_load.pop()
        self._title.append(title)
        system = line_load.pop() #line_load[8]
        #
        # push to SQL
        bd_file = self._bd_file
        conn = create_connection(bd_file)
        with conn:
            self._push_beam_load(conn, beam_number, title,
                                 system, line_load)
        #print("-->")
    #
    def __getitem__(self, beam_name:int|str) -> list:
        """
        """
        conn = create_connection(self._bd_file)      
        with conn:
            udl = self._get_beam_load(conn, beam_name=beam_name, load_name=self._name)
        return udl
    #
    #
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
    #    
    #
    def _create_table(self, conn) -> None:
        """ """
        _table_element_line = "CREATE TABLE IF NOT EXISTS tb_LoadBeamLine(\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                load_number INTEGER NOT NULL REFERENCES tb_Load(number),\
                                element_number INTEGER NOT NULL REFERENCES tb_Elements(number),\
                                title TEXT,\
                                system INTEGER NOT NULL,\
                                L_end1 DECIMAL,\
                                qx1 DECIMAL,\
                                qy1 DECIMAL,\
                                qz1 DECIMAL,\
                                qx1i DECIMAL,\
                                qy1i DECIMAL,\
                                qz1i DECIMAL,\
                                L_end2 DECIMAL,\
                                qx2 DECIMAL,\
                                qy2 DECIMAL,\
                                qz2 DECIMAL,\
                                qx2i DECIMAL,\
                                qy2i DECIMAL,\
                                qz2i DECIMAL);"
        #
        #bd_file = self._bd_file
        #conn = create_connection(bd_file)
        create_table(conn, _table_element_line)
    #
    def _push_beam_load(self, conn, beam_number:int,
                        load_title:str, load_system:int,
                        udl:List[float]):
        """ """
        #beam = check_element(conn, beam_name)
        #beam_number = beam[0]        
        #print("-->")
        load_data = get_load_data(conn, self._name, load_type='basic')
        load_number = load_data[0]
        #
        project = (load_number, beam_number,
                   load_title, load_system,
                   udl[6], *udl[:3],
                   'NULL', 'NULL', 'NULL',
                   udl[7], *udl[3:6],
                   'NULL', 'NULL', 'NULL')
        #
        sql = 'INSERT INTO tb_LoadBeamLine(load_number, element_number,\
                                            title, system,\
                                            L_end1, qx1, qy1, qz1, qx1i, qy1i, qz1i,\
                                            L_end2, qx2, qy2, qz2, qx2i, qy2i, qz2i)\
                                            VALUES(?,?,?,?,\
                                                   ?,?,?,?,?,?,?,\
                                                   ?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)      
    #
    #
    def _get_beam_load(self, conn, beam_name:int, load_name:int):
        """ """
        # get beam data
        beam = check_element(conn, beam_name)
        beam_number = beam[0]
        # beam line load
        load_data = get_load_data(conn, self._name, load_type='basic')
        load_number = load_data[0]
        #
        cur = conn.cursor()
        cur.execute("SELECT tb_Load.name, tb_LoadBeamLine.*\
                    FROM tb_LoadBeamLine, tb_Load\
                    WHERE tb_LoadBeamLine.load_number = {:}\
                    AND tb_LoadBeamLine.load_number = tb_Load.number\
                    AND tb_LoadBeamLine.element_number = {:};"
                    .format(load_number, beam_number))
        rows = cur.fetchall()
        beam_line = [] # defaultdict(list)
        for row in rows:
            data = [*row[7:10], *row[14:17], row[6], row[13],
                    beam_name, row[4], row[0], # name, title, load_name,
                    row[5], 0, "Line Load"]    # system, load_complex, load_type
            beam_line.append(LineBeam._make(data))
        return beam_line
    #
    #
#
#
class BeamPointSQL(NodeLoadBasic):
    __slots__ = ['_labels', '_title', '_complex', 
                 '_system_flag', '_system', 
                 '_bd_file', '_name']

    def __init__(self, load_name: int|float, bd_file: str) -> None:
        """
        """
        super().__init__()
        self._name = load_name
        self._bd_file =  bd_file
        # create node table
        conn = create_connection(self._bd_file)
        with conn:        
            self._create_table(conn)
    #
    def __setitem__(self, beam_name: int|str, point_load: list) -> None:
        """
        """
        # get load data
        # set element load
        self._labels.append(beam_name)
        title = point_load.pop()
        self._title.append(title)
        system = point_load.pop() #line_load[8]        
        # 
        # push to SQL
        bd_file = self._bd_file
        conn = create_connection(bd_file)
        with conn:
            self._push_beam_load(conn, beam_name, title,
                                 system, point_load)
        # print("-->")
    #
    def __getitem__(self, beam_name:int|str)-> list:
        """
        """
        bd_file = self._bd_file
        # get beam load
        conn = create_connection(bd_file)
        with conn:
            pl = self._get_beam_load(conn, beam_name=beam_name, load_name=self._name)
        return pl
    #
    #
    def _create_table(self, conn) -> None:
        """ """
        _table_element_point = "CREATE TABLE IF NOT EXISTS tb_LoadBeamPoint(\
                                    number INTEGER PRIMARY KEY NOT NULL,\
                                    load_number INTEGER NOT NULL REFERENCES tb_Load(number),\
                                    element_number INTEGER NOT NULL REFERENCES tb_Elements(number),\
                                    title TEXT,\
                                    system INTEGER NOT NULL,\
                                    L_end1 DECIMAL,\
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
        #bd_file = self._bd_file
        #conn = create_connection(bd_file)
        create_table(conn, _table_element_point)
    #
    def _push_beam_load(self, conn, beam_name:int|str,
                        load_title: str, load_system: int,
                        point_load:list[float]):
        """ """
        #print("-->")
        beam = check_element(conn, beam_name)
        beam_number = beam[0]
        #
        load_data = get_load_data(conn, self._name, load_type='basic')
        load_number = load_data[0]
        #
        #
        project = (load_number, beam_number,
                   load_title, load_system,
                   point_load[6], *point_load[:6],
                   'NULL', 'NULL', 'NULL',
                   'NULL', 'NULL', 'NULL')
        #
        sql = 'INSERT INTO tb_LoadBeamPoint(load_number, element_number,\
                                            title, system, \
                                            L_end1, fx, fy, fz, mx, my, mz,\
                                            fxi, fyi, fzi, mxi, myi, mzi)\
                                            VALUES(?,?,?,?,?,\
                                                   ?,?,?,?,?,?,?,\
                                                   ?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    #
    def _get_beam_load(self, conn, beam_name:int, load_name:int):
        """ """
        # get beam data
        beam = check_element(conn, beam_name)
        beam_number = beam[0]
        # beam line load
        load_data = get_load_data(conn, self._name, load_type='basic')
        load_number = load_data[0]
        #
        # beam line load
        cur = conn.cursor()
        cur.execute("SELECT tb_Load.name, tb_LoadBeamPoint.* \
                    FROM tb_LoadBeamPoint, tb_Load\
                    WHERE tb_LoadBeamPoint.load_number = {:}\
                    AND tb_LoadBeamPoint.load_number = tb_Load.number\
                    AND tb_LoadBeamPoint.element_number = {:};"
                    .format(load_number, beam_number))
        rows = cur.fetchall()
        beam_line = []
        for row in rows:
            data = [*row[7:13], row[6],
                    beam_name, row[4], row[0],  # name, title, load_name,
                    row[5], 0, "Point Load"]    # system, load_complex, load_type
            beam_line.append(PointBeam._make(data))
        return beam_line
#    
#
class BeamToNodeSQL(NodeLoadBasic):
    __slots__ = ['_labels', '_title', '_complex', 
                '_system_flag', '_system', 
                '_bd_file', '_name']
    def __init__(self, load_name: int|float, bd_file: str) -> None:
        """
        """
        super().__init__()
        self._name = load_name
        self._bd_file =  bd_file
        # create node table
        conn = create_connection(self._bd_file)
        with conn:        
            self._create_table(conn)
    #
        #
    def __setitem__(self, beam_name: int|str,
                    node_load: list[float]|dict[str,float]) -> None:
        """
        """
        # get load data
        # set element load
        self._labels.append(beam_name)
        #
        #
        # push to SQL
        bd_file = self._bd_file
        conn = create_connection(bd_file)
        for node in node_load:
            title = node.pop()
            system = node.pop()
            with conn:
                self._push_beam_load(conn, beam_name, title,
                                     system, node)
        # print("-->")
    #
    def __getitem__(self, beam_name:int|str)-> list:
        """
        """
        bd_file = self._bd_file
        # get beam load
        conn = create_connection(bd_file)
        with conn:
            pl = self._get_beam_load(conn, beam_name=beam_name, load_name=self._name)
        return pl
    #
    def _create_table(self, conn) -> None:
        """ """
        _table_element_point = "CREATE TABLE IF NOT EXISTS tb_LoadBeamToNode(\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                load_number INTEGER NOT NULL REFERENCES tb_Load(number),\
                                title TEXT,\
                                system TEXT,\
                                element_number INTEGER NOT NULL REFERENCES tb_Elements(number),\
                                node_number INTEGER NOT NULL REFERENCES tb_Nodes(number),\
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
        create_table(conn, _table_element_point)
    #
    def _push_beam_load(self, conn, beam_name:int|str,
                        load_title: str, load_system: int,
                        node_load:list[float]):
        """ """
        #print("-->")
        beam = check_element(conn, beam_name)
        beam_number = beam[0]
        #
        node_name = node_load.pop(0)
        node = check_nodes(conn, node_name)
        node_number = node[0]
        #
        load_data = get_load_data(conn, self._name, load_type='basic')
        load_number = load_data[0]
        #
        try:
            1 / load_system
            load_system = 'global'
        except ZeroDivisionError:
            raise RuntimeError('node load in global system')
        #
        project = (load_number, load_title, load_system,
                   beam_number, node_number,
                   *node_load,
                   'NULL', 'NULL', 'NULL',
                   'NULL', 'NULL', 'NULL')
        #
        sql = 'INSERT INTO tb_LoadBeamToNode(load_number, title, system, \
                                            element_number,\
                                            node_number, fx, fy, fz, mx, my, mz,\
                                            fxi, fyi, fzi, mxi, myi, mzi)\
                                            VALUES(?,?,?,?,?,\
                                                   ?,?,?,?,?,?,?,\
                                                   ?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    def _get_beam_load(self, conn, beam_name:int, load_name:int):
        """ """
        # get beam data
        beam = check_element(conn, beam_name)
        beam_number = beam[0]
        # beam line load
        load_data = get_load_data(conn, self._name, load_type='basic')
        load_number = load_data[0]
        #
        # beam line load
        cur = conn.cursor()
        cur.execute("SELECT tb_Load.title, tb_LoadBeamToNode.* \
                    FROM tb_LoadBeamPoint, tb_LoadBeamLine, tb_LoadBeamToNode, tb_Load\
                    WHERE tb_LoadBeamToNode.load_number = {:}\
                    AND tb_LoadBeamToNode.load_number = tb_Load.number\
                    AND tb_LoadBeamToNode.element_number = {:};"
                    .format(load_number, beam_number))
        rows = cur.fetchall()
        node_load = []
        1/0
        for row in rows:
            #data = [*row[7:10], *row[14:17], row[6], row[13],
            #        node_number, row[2], *row[4:6]]
            data = [*row[5:11],
                    #*row[2:4],
                    load_number, self._name, 
                    0, 0, self._type]
            node_load.append(PointNode._make(data))
        return node_load
    #
    @property
    def df(self):
        """nodes in dataframe format"""
        db = DBframework()
        conn = create_connection(self._bd_file)
        #
        with conn:
            cur = conn.cursor()
            cur.execute("SELECT tb_Load.name, tb_Load.type, \
                        tb_Nodes.name, tb_Elements.name, \
                        tb_LoadBeamToNode.title, tb_LoadBeamToNode.system, \
                        tb_LoadBeamToNode.fx, tb_LoadBeamToNode.fy, tb_LoadBeamToNode.fz, \
                        tb_LoadBeamToNode.mx, tb_LoadBeamToNode.my, tb_LoadBeamToNode.mz \
                        FROM tb_Load, tb_Nodes, tb_Elements, tb_LoadBeamToNode\
                        WHERE tb_LoadBeamToNode.load_number = tb_Load.number\
                        AND tb_LoadBeamToNode.node_number = tb_Nodes.number \
                        AND tb_LoadBeamToNode.element_number = tb_Elements.number \
                        AND tb_Load.name = {:};".format(self._name))
        #
        rows = cur.fetchall()
        #
        cols = ['load_name', 'load_type', 'node_name', 'element_name', 
                'load_comment', 'load_system',
                'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        df = db.DataFrame(data=rows, columns=cols)
        #df = db.read_sql_query("SELECT * FROM tb_LoadNode", conn)
        df = df[['load_name', 'load_type', 'load_comment', 'load_system',
                 'element_name', 'node_name',
                 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']]
        return df    
#
#

