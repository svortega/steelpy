# Copyright (c) 2009-2020 fem2ufo


# Python stdlib imports
from array import array
from collections.abc import Mapping
import re
from typing import NamedTuple, Tuple, List, Iterator, Iterable, Union, Dict


# package imports
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table

#
class BoundaryItem(NamedTuple):
    """
    """
    x: float
    y: float
    z: float
    rx: float
    ry: float
    rz: float
    number: int
    name: Union[str,None]
    node:int
    
    def __str__(self) -> str:
        if (name := self.name) == 'NULL':
            name = ""
        return "{:12d} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f} {: 8.0f} {:>12s}\n"\
            .format(self.node, self.x, self.y, self.z, self.rx, self.ry, self.rz, name)


#
class BoundaryNodeSQL(Mapping):
    
    __slots__ = ['db_file', '_labels']
    
    def __init__(self, db_file:str):
        """
        """
        self.db_file = db_file
        self._labels : array = array('I', [])
        # create node table
        self._create_table()
    #
    def __setitem__(self, node_number: int,
                    fixity:Union[List, Tuple, Dict, str]) -> None:
        """
        """
        try:
            # TODO : update data option if needed?
            self._labels.index(node_number)
            raise Warning('    *** warning node {:} boundary already exist'.format(node_number))
        except ValueError:
            self._labels.append(node_number)
            fixity = self._get_fixity(fixity)
            conn = create_connection(self.db_file)
            with conn:
                self._push_boundary(conn, node_number, fixity)
                conn.commit()
    #
    def __getitem__(self, node_number: int) -> Union[Tuple,bool]:
        """
        """
        try:
            index = self._labels.index(node_number)
            conn = create_connection (self.db_file)
            data = get_boundary(conn, node_number)
            return BoundaryItem(*data[3:], number=data[0],
                                name=data[1], node=node_number)
        except ValueError:
            return False
    #
    def _create_table(self) -> None:
        """ """
        _table_boundary = "CREATE TABLE IF NOT EXISTS tb_Boundaries(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            title TEXT,\
                            node_name INTEGER NOT NULL REFERENCES tb_Nodes(name),\
                            x DECIMAL, y DECIMAL, z DECIMAL,\
                            rx DECIMAL, ry DECIMAL, rz DECIMAL);"
        #
        conn = create_connection(self.db_file)
        create_table(conn, _table_boundary)
    #
    def __len__(self) -> float:
        return len(self._labels)

    def __iter__(self) -> Iterator:
        """
        """
        return iter(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    #
    def _push_boundary(self, conn, node_name, fixity):
        """
        """
        try:
            title = fixity[6]
        except IndexError:
            title = None

        project = (title, node_name, *fixity[:6])
        sql = 'INSERT INTO tb_Boundaries(title, node_name,\
                                         x, y, z, rx, ry, rz)\
                                         VALUES(?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    def _get_fixity(self, fixity):
        """ """
        if isinstance(fixity, str):
            if re.match(r"\b(fix(ed)?)\b", fixity, re.IGNORECASE):
                return [1,1,1,1,1,1]
            elif re.match(r"\b(pinn(ed)?|roll)\b", fixity, re.IGNORECASE):
                return [1,1,1,1,0,0]
            elif re.match(r"\b(free)\b", fixity, re.IGNORECASE):
                return None
            else:
                raise IOError("boundary type {:} not implemented".format(fixity))
        elif isinstance(fixity, (list, tuple)):
            return fixity
        elif isinstance(fixity, dict):
            return [fixity['x'], fixity['y'], fixity['z'], 
                    fixity['rx'], fixity['ry'], fixity['rz']]
        else:
            raise Exception('   *** Boundary input format not recognized')        
#
#
#
def get_boundary(conn, node_number, item:str="*"):
    """
    """
    #
    project = (node_number,)
    sql = 'SELECT {:} FROM tb_Boundaries WHERE node_name = ?'.format(item)
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
