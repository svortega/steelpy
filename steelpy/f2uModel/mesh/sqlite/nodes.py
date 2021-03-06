# Copyright (c) 2009-2021 fem2ufo

# Python stdlib imports
from array import array
from collections.abc import Mapping
from math import isclose, dist
from typing import NamedTuple, Tuple, List, Iterator, Iterable, Union, Dict

# package imports
#from steelpy.process.units.main import Units
from steelpy.f2uModel.mesh.operations.nodes import (check_point_list, check_point_dic, 
                                                    get_coordinate_system)
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table



#
class NodeSQL(Mapping):
    """
    This is a fem2ufo model node class


    Parameters
    ----------
    boundaries: object
        f2u boundary object


    Attributes
    ----------
    _labels : array
        node internal number
    _x : array
        coordinate x
    _y : array
        coordinate y
    _z : array
        coordinate y
    _sets : List[tuple]
        set with node/element
    """
    __slots__ = ['_system', 'db_file', '_labels']

    def __init__(self, db_file: str,
                 db_system:str="sqlite") -> None:
        """
        """
        #self._units = Units()
        self.db_file = db_file
        self._labels : array = array('I', [])
        self._system = 'cartesian'
        # create node table
        conn = create_connection(self.db_file)
        with conn:
            new_node_table(conn)
    #
    # ---------------------------------
    #
    def __setitem__(self, node_name: int,
                    coordinates: Union[List[float], Dict[str, float]]) -> None:
        """
        """
        try:
            self._labels.index(node_name)
            raise Exception('    *** warning point {:} already exist'
                            .format(node_name))
        except ValueError:
            coordinates = self._get_coordinates(coordinates)
            self._labels.append(node_name)
            #
            conn = create_connection(self.db_file)
            with conn:
                self._push_node(conn, node_name, coordinates)
                conn.commit()
    #
    def __getitem__(self, node_name: int) -> Tuple:
        """
        node_name : node number
        """
        try:
            index = self._labels.index(node_name)
            conn = create_connection(self.db_file)
            data = get_node(conn, node_name)
            system = get_coordinate_system(data[2])
            # FIXME : include rest system
            #number = data[0] - 1
            return system(x=data[3], y=data[4], z=data[5],
                          name=node_name, number=data[0], 
                          index=data[0]-1)
        except ValueError:
            raise IndexError('   *** node {:} does not exist'.format(node_name))    
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
    @property
    def system(self) -> str:
        """
        """
        return self._system

    @system.setter
    def system(self, value:str) -> None:
        """
        """
        self._system = value
    #
    #
    def _get_coordinates(self, coordinates):
        """ """
        if isinstance(coordinates, (list, tuple)):
            coordinates = check_point_list(coordinates, steps=3)
        elif isinstance(coordinates, dict):
            coordinates = check_point_dic(coordinates)
        else:
            raise Exception('   *** Node input format not recognized')
        return coordinates
    #
    def _push_node(self, conn, node_name, coordinates):
        """
        Create a new project into the projects table
        """
        #number = len(self._labels) - 1
        if self.system == 'cylindrical': # z, r, theta,
            project = (node_name, self.system,
                       None, None, *coordinates, None)
        
        elif self.system == 'spherical': # r, theta, phi
            project = (node_name, self.system,
                       None, None, None, *coordinates)
        
        else:
            project = (node_name, self.system,
                       *coordinates, None, None, None)
        #
        sql = 'INSERT INTO tb_Nodes(name, type,\
                                    x, y, z, r, theta, phi)\
                                    VALUES(?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    #
    def get_new_point(self, coordinates):
        """ """
        #create a new point
        while True:
            #node_name = "pnt_{:}".format(str(next(self.get_number())))
            node_name = next(self.get_number())
            try:
                self._labels.index(node_name)
            except ValueError:
                break
        self.__setitem__(node_name, coordinates)
        return node_name
    #
    def get_point_name(self, coordinates,
                       tol:float=0.01, rel_tol:float=1e-6) -> int:
        """
        tol: absolte tolerance in metres (0.010 m default)
        """
        # check if point already exist
        system = get_coordinate_system(self.system)
        if isinstance(coordinates, system):
            return coordinates.name
        # get index of x coord location in existing database
        coord = self._get_coordinates(coordinates)
        #
        items = self._iscloseSQL(key='*', item='x', value=coord[0],
                                 abs_tol=tol, rel_tol=rel_tol)
        # check if y and z coord match
        if items:
            for item in items:
                if isclose(coord[1], item[4], abs_tol=tol, rel_tol=rel_tol):
                    if isclose(coord[2], item[5], abs_tol=tol, rel_tol=rel_tol):
                        return item[1]
        raise IOError('   error coordinate not found')
    #
    def _iscloseSQL(self, key:str, item:str, value:float,
                    rel_tol:float=1e-9, abs_tol:float=0.0)-> tuple:
        """ """
        sql = 'SELECT {:} FROM tb_Nodes WHERE ABS({:} - {:}) <= {:} * MAX(ABS({:}), ABS({:}), {:})'\
              .format(key, item, value, rel_tol, item, value, abs_tol)
        conn = create_connection(self.db_file)
        cur = conn.cursor()
        cur.execute(sql)
        records = cur.fetchall()
        return records
    #
    def get_number(self, start:int=1)-> Iterable[int]:
        """
        """
        try:
            n = max(self._labels) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
    #
    #
    def renumbering(self, new_numbers:List[int]):
        """ """
        indexes = [self._labels.index(node_name) 
                   for node_name in new_numbers]
        conn = create_connection(self.db_file)
        with conn:
            nodes = get_nodes(conn)
            #nindex = [item[1] for item in nodes]
            nodes = [nodes[indx][1:] for indx in indexes]
            update_table(conn, nodes)
            conn.commit()
    #
    #def update_number(self, node_name:int, value:Union[float,int]):
    #    """ """
    #    conn = create_connection(self.db_file)
    #    with conn:
    #        nodes = get_nodes(conn)
    #        nindex = [item[1] for item in nodes]
    #        row = nodes.pop(nindex.index(node_name))
    #        nodes.insert(value-1, row)
    #        update_table(conn, nodes)
    #        conn.commit()
    #
    def _update_item(self, conn, name, item, value):
        """ """
        project = (value, name)
        sql = 'UPDATE tb_Nodes SET {:} = ? WHERE name = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    def _orderby(self):
        """
        """
        #print('--')
        conn = create_connection(self.db_file)
        sql = 'SELECT * FROM tb_Nodes ORDER BY number ASC'
        cur = conn.cursor()
        cur.execute(sql)
        conn.commit()
        conn.close()
#
#
#
def get_node(conn, node_name:int, item:str='*'):
    """ """
    with conn:
        value = get_item_table(conn, node_name, item)
    return value
#
def get_item_table(conn, node_name, item):
    """ """
    project = (node_name,)
    sql = 'SELECT {:} FROM tb_Nodes WHERE name = ?'.format(item)
    cur = conn.cursor()
    cur.execute(sql, project)
    record = cur.fetchone()
    return record
#
def get_nodes(conn):
    """ """
    sql = 'SELECT * FROM tb_Nodes ORDER BY number ASC'
    cur = conn.cursor()
    cur.execute(sql)
    record = cur.fetchall()
    return record
#
def update_table(conn, nodes):
    """ """
    # drop table
    sql = 'DROP TABLE tb_Nodes'
    cur = conn.cursor()
    cur.execute(sql)
    #
    new_node_table(conn)
    push_nodes(conn, nodes)
#
def new_node_table(conn) -> None:
    """ """
    # conn = create_connection(self.db_file)
    _table_nodes = "CREATE TABLE IF NOT EXISTS tb_Nodes (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name INTEGER NOT NULL,\
                    type TEXT NOT NULL,\
                    x DECIMAL,\
                    y DECIMAL,\
                    z DECIMAL,\
                    r DECIMAL,\
                    theta DECIMAL,\
                    phi DECIMAL);"
    #
    create_table(conn, _table_nodes)
#
def push_nodes(conn, nodes):
    """
    Create a new project into the projects table
    :param conn:
    :param project:

    :return: project id
    """
    project = nodes
    #project = [item[1:] for item in nodes]
    # number = len(self._labels) - 1
    #if csystem == 'cylindrical':  # z, r, theta,
    #    project = (node_name, csystem,
    #               None, None, *coordinates, None)
    #
    #elif csystem == 'spherical':  # r, theta, phi
    #    project = (node_name, csystem,
    #               None, None, None, *coordinates)
    #
    #else:
    #    project = (node_name, csystem,
    #               *coordinates, None, None, None)
    #
    sql = 'INSERT INTO tb_Nodes(name, type,\
                                x, y, z, r, theta, phi)\
                                VALUES(?,?,?,?,?,?,?,?)'
    cur = conn.cursor()
    cur.executemany(sql, project)
#
#
#
#