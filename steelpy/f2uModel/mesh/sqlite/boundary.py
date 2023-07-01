# Copyright (c) 2009-2023 fem2ufo

# Python stdlib imports
from __future__ import annotations
from array import array
#from collections.abc import Mapping
import re
from typing import NamedTuple


# package imports
from steelpy.f2uModel.mesh.process.boundary import BoundaryItem, BoundaryNode
from steelpy.f2uModel.mesh.sqlite.process_sql import create_connection, create_table, check_nodes

#
#
#
#
class BoundaryNodeSQL(BoundaryNode):
    
    __slots__ = ['_db_file', '_labels']
    
    def __init__(self, db_file:str):
        """
        """
        super().__init__()
        self._db_file = db_file
        # create node table
        self._create_table()
    #
    def __setitem__(self, node_name: int,
                    fixity:list|tuple|dict|str) -> None:
        """
        """
        conn = create_connection(self._db_file)
        with conn:  
            node = check_nodes(conn, node_name)       
        try:
            node_number = node[0]
        except TypeError:
            raise IOError(f"Node {node_name} not found")
        
        try:
            # TODO : update data option if needed?
            self._labels.index(node_name)
            raise Warning('    *** warning node {:} boundary already exist'.format(node_name))
        except ValueError:
            self._labels.append(node_name)
            fixity = self._get_fixity(fixity)
            conn = create_connection(self._db_file)
            with conn:
                self._push_boundary(conn, node_number, fixity)
        #
    #
    def __getitem__(self, node_name: int) -> tuple | bool:
        """
        """
        conn = create_connection(self._db_file)
        with conn:  
            node = check_nodes(conn, node_name)       
        try:
            node_number = node[0]
        except TypeError:
            raise IOError(f"Node {node_name} not found")
        
        try:
            self._labels.index(node_name)
            conn = create_connection (self._db_file)
            data = get_boundary(conn, node_number)
            return BoundaryItem(*data[3:], number=data[0],
                                name=data[1], node=node_name)
        except ValueError:
            return False
    #
    def _create_table(self) -> None:
        """ """
        _table_boundary = "CREATE TABLE IF NOT EXISTS tb_Boundaries(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            title TEXT,\
                            node_number INTEGER NOT NULL REFERENCES tb_Nodes(name),\
                            x DECIMAL, y DECIMAL, z DECIMAL,\
                            rx DECIMAL, ry DECIMAL, rz DECIMAL);"
        #
        conn = create_connection(self._db_file)
        create_table(conn, _table_boundary)
    #
    #
    def _push_boundary(self, conn, node_number: int, fixity: list):
        """
        """
        try:
            title = fixity[6]
        except IndexError:
            title = None

        project = (title, node_number, *fixity[:6])
        sql = 'INSERT INTO tb_Boundaries(title, node_number,\
                                         x, y, z, rx, ry, rz)\
                                         VALUES(?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    def transposed(self):
        """??? """
        1 / 0
#
#
def get_boundary(conn, node_number, item:str="*"):
    """
    """
    #
    project = (node_number,)
    sql = 'SELECT {:} FROM tb_Boundaries WHERE node_number = ?'.format(item)
    cur = conn.cursor()
    cur.execute(sql, project)
    record = cur.fetchone()
    return record
#
#
#
class BoundarySQL:
    
    def __init__(self, db_file:str,
                 db_system:str="sqlite") -> None:
        """
        """
        self._nodes = BoundaryNodeSQL(db_file)
    #
    @property
    def node(self):
        """"""
        return self._nodes
    
    @node.setter
    def node(self, values):
        """"""
        for value in values:
            self._nodes[value[0]] = value[1:]
#
