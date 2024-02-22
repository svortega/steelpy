# Copyright (c) 2009 steelpy

# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
from math import isclose, dist
from typing import NamedTuple
import re

# package imports
from steelpy.ufo.mesh.elements.nodes import NodePoint, NodeBasic
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework


#
class NodeSQL(NodeBasic):
    """
    steelpy model node class


    Parameters
    ----------
    boundaries: object
        UFOmodel Node object

    Attributes
    ----------
    _labels : array
        node internal number
    x : array
        coordinate x
    y : array
        coordinate y
    z : array
        coordinate y
    sets : list[tuple]
        set with node/element
    """
    __slots__ = ['_system', 'db_file', '_component']

    def __init__(self,
                 db_file: str,
                 plane: NamedTuple,
                 component: str|int, 
                 db_system:str="sqlite",
                 system:str = 'cartesian') -> None:
        """
        """
        super().__init__(system)
        #
        self.db_file = db_file
        self._plane = plane
        self._component =  component
        # create node table
        conn = create_connection(self.db_file)
        with conn:
            self._create_table(conn)
        #print('--> update labels?')
        #1 / 0
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._component, )
        table = 'SELECT Node.name FROM Node \
                 WHERE component_id = ? \
                 ORDER BY number ASC ;'
        #
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]
    #
    # ---------------------------------
    #
    def __setitem__(self, node_number: int,
                    coordinates: list[float]|dict[str, float]) -> None:
        """
        """
        if not isinstance(node_number, int):
            raise IOError('node id must me an integer')
        #
        try:
            self._labels.index(node_number)
            raise Exception(f' warning node {node_number} already exist')
        except ValueError:
            coordinates = self._get_coordinates(coordinates)
            #self._labels.append(node_number)
            #
            conn = create_connection(self.db_file)
            with conn:
                self._push_node(conn, node_number, coordinates)
                #conn.commit()
    #
    def __getitem__(self, node_number: int) -> tuple:
        """
        node_number : node number
        """
        try:
            self._labels.index(node_number)
            conn = create_connection(self.db_file)
            with conn:
                node = get_node(conn, node_number, self._component)
            return node
        except ValueError:
            raise IndexError(f' node id : {node_number} not valid')
    #
    # ---------------------------------
    # SQL ops
    #
    def _create_table(self, conn) -> None:
        """ """
        # Node main
        table = "CREATE TABLE IF NOT EXISTS Node (\
                number INTEGER NOT NULL,\
                name INTEGER NOT NULL,\
                component_id INTEGER NOT NULL REFERENCES Component(number), \
                title TEXT,\
                mesh_idx INTEGER NOT NULL, \
                PRIMARY KEY (number));"
        #
        create_table(conn, table)
        #
        # Node coordinate
        table = "CREATE TABLE IF NOT EXISTS NodeCoordinate (\
                number INTEGER NOT NULL,\
                node_id INTEGER NOT NULL REFERENCES Node(number), \
                system TEXT NOT NULL,\
                x DECIMAL NOT NULL,\
                y DECIMAL NOT NULL,\
                z DECIMAL NOT NULL,\
                r DECIMAL, \
                theta DECIMAL, \
                phi DECIMAL, \
                PRIMARY KEY (number));"
        #
        create_table(conn, table)
    #
    #
    def _push_node(self, conn, node_number: int,
                   coordinates: list):
        """
        Create a new project into the projects table
        """
        # get row number
        cur = conn.cursor()
        sql = 'SELECT max(number) from Node;'
        if (idx := max(cur.execute(sql))[0]) == None:
            idx = 0
        #
        # -------------------------------------------
        #
        query = (node_number, self._component, None, idx)
        table = 'INSERT INTO Node(name, component_id, \
                                  title, mesh_idx) \
                                  VALUES(?,?,?,?);' 
        #
        #cur = conn.cursor()
        node_id = cur.execute(table, query).lastrowid
        #node_id = node_id.lastrowid
        #
        # -------------------------------------------
        #
        if self.system == 'cylindrical': # z, r, theta,
            project = (node_id, self.system,
                       None, None, *coordinates, None, )
        
        elif self.system == 'spherical': # r, theta, phi
            project = (node_id, self.system,
                       None, None, None, *coordinates, )
        
        else:
            project = (node_id, self.system,
                       *coordinates, None, None, None,)
        #
        sql = 'INSERT INTO NodeCoordinate(node_id, system,\
                                          x, y, z, r, theta, phi)\
                            VALUES(?,?,?,?,?,?,?,?)'
        # push
        #cur = conn.cursor()
        cur.execute(sql, project)
        #print('-->')
    #
    def _isclose(self, key:str, item:str, value:float,
                    rel_tol:float=1e-9, abs_tol:float=0.0)-> tuple:
        """ """
        query = (self._component, )
        table = 'SELECT NodeCoordinate.{:} \
                 FROM Node, NodeCoordinate \
                 WHERE Node.component_id = ? \
                 AND ABS({:} - {:}) <= {:} * MAX(ABS({:}), ABS({:}), {:})'\
            .format(key, item, value, rel_tol, item, value, abs_tol)
        #
        conn = create_connection(self.db_file)
        cur = conn.cursor()
        cur.execute(table, query)
        records = cur.fetchall()
        return records
    #
    def _update_item(self, conn, name, item, value):
        """ """
        1 / 0
        project = (value, name, self._component, )
        sql = f'UPDATE Node SET {item} = ? \
                WHERE name = ? AND component_id = ? '
        
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    def _orderby(self):
        """
        """
        1 / 0
        query = (self._component, )
        conn = create_connection(self.db_file)
        table = 'SELECT * FROM Node \
                WHERE component_id = ? \
                ORDER BY number ASC'
        cur = conn.cursor()
        cur.execute(table, query)
        conn.commit()
        conn.close()
    #
    # ---------------------------------
    # ops
    #
    def get_point_name(self, coordinates,
                       tol:float=0.01, rel_tol:float=1e-6) -> int:
        """
        tol: absolte tolerance in metres (0.010 m default)
        """
        # check if point already exist
        #try:
        #    #system = self.system # get_coordinate_system(self.system)
        #    #if isinstance(coordinates, system):
        #    return coordinates.name
        #except AttributeError:
        # get index of x coord location in existing database
        coord = self._get_coordinates(coordinates)
        #
        items = self._isclose(key='*', item='x', value=coord[0],
                              abs_tol=tol, rel_tol=rel_tol)
        # check if y and z coord match
        if items:
            for item in items:
                if isclose(coord[1], item[4], abs_tol=tol, rel_tol=rel_tol):
                    if isclose(coord[2], item[5], abs_tol=tol, rel_tol=rel_tol):
                        return item[1]
        raise IOError('   error coordinate not found')
    #
    #
    #
    #
    #
    def renumbering(self, new_numbers:list[int]):
        """ """
        indexes = [self._labels.index(node_name) 
                   for node_name in new_numbers]
        indexes = [(val, j + 1, self._component, )
                   for j, val in enumerate(indexes)]
        #
        conn = create_connection(self.db_file)
        with conn:
            update_colum(conn, colname='mesh_idx',
                         newcol=indexes)
            #
            #nodes = get_nodes(conn)
            #nindex = [item[1] for item in nodes]
            #nodes = [nodes[indx][1:] for indx in indexes]
            #update_table(conn, nodes)
            #conn.commit()
        #print('-->?')
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
    def get_maxmin(self):
        """
        """
        def maxmin(head: str, col: str):
            #
            query = (self._component, )
            table = f'SELECT {head.upper()}(NodeCoordinate.{col}) \
                      FROM Node, NodeCoordinate \
            WHERE Node.component_id = ? \
            AND Node.number = NodeCoordinate.node_id'
            
            cur = conn.cursor()
            cur.execute(table, query)
            record = cur.fetchone()
            return record[0]
        #
        conn = create_connection(self.db_file)
        with conn:
            max_x = maxmin(head='max', col='x')
            min_x = maxmin(head='min', col='x')
            #
            max_y = maxmin(head='max', col='y')
            min_y = maxmin(head='min', col='y')
            #
            max_z = maxmin(head='max', col='z')
            min_z = maxmin(head='min', col='z')
        return [max_x, max_y, max_z], [min_x, min_y, min_z]
    #
    #
    #def get_number(self, start:int=1):
    #    """
    #    """
    #    try:
    #        n = max(self._labels) + 1
    #    except ValueError:
    #        n = start
    #    #
    #    while True:
    #        yield n
    #        n += 1    
    #
    # ---------------------------------
    # dataframe
    #
    @property
    def df(self):
        """get node dataframe"""
        db = DBframework()
        conn = create_connection(self.db_file)
        table = f'SELECT Node.name, Node.component_id, NodeCoordinate.* \
                  FROM Node, NodeCoordinate \
                  WHERE Node.component_id = {self._component} \
                  AND Node.number = NodeCoordinate.node_id ;'
        nodes = db.read_sql_query(table,conn)
        nodes.drop(columns=['node_id'], inplace=True)
        return nodes
    
    @df.setter
    def df(self, df):
        """ """
        1 / 0
        columns = list(df.columns)
        header = {}
        for key in columns:
            if re.match(r"\b(id|name|node(s)?)\b", key, re.IGNORECASE):
                header[key] = 'name'
            
            elif re.match(r"\b(x(\_|\-)?(coord(inate)?(s)?)?)\b", key, re.IGNORECASE):
                header[key] = 'x'
                try:
                    df[key] = df[key].apply(lambda x: x.value)
                except AttributeError:
                    pass
            
            elif re.match(r"\b(y(\_|\-)?(coord(inate)?(s)?)?)\b", key, re.IGNORECASE):
                header[key] = 'y'
                try:
                    df[key] = df[key].apply(lambda x: x.value)
                except AttributeError:
                    pass                
            
            elif re.match(r"\b(z(\_|\-)?(coord(inate)?(s)?)?)\b", key, re.IGNORECASE):
                header[key] = 'z'
                try:
                    df[key] = df[key].apply(lambda x: x.value)
                except AttributeError:
                    pass                
            
            #elif re.match(r"\b(title)\b", key, re.IGNORECASE):
            #    #header.append(key)
            #    header[key] = 'title'
        #
        #df.rename(columns=header, inplace=True)
        #
        nodes = df[header.keys()].copy()
        #nodes[["x", "y", "z"]] = df[coord].values.tolist()
        nodes['component_id'] = self._component
        nodes['type'] = self._system
        nodes[['r', 'theta', 'phi']] = None
        #
        nodes.rename(columns=header, inplace=True)
        #
        #for row in nodes.itertuples():
        #    coordinates = self._get_coordinates([])
        #    self._labels.append(node_name)
        #
        conn = create_connection(self.db_file)
        #
        # get row number
        query = (self._component, )
        table = 'SELECT max(number) from Node WHERE component_id = ?;'
        cur = conn.cursor()
        if (idx := max(cur.execute(table, query))[0]) == None:
            idx = 0
        #
        nodes['idx'] = [item + idx for item in list(nodes.index)]
        #
        with conn:
            nodes.to_sql('Node', conn,
                         index_label=['name', 'component_id', 'type',
                                      'x', 'y', 'z', 'r',
                                      'theta', 'phi', 'idx'], 
                         if_exists='append', index=False)
        #
        #
        #self._labels.extend(nodes['name'].tolist())
        #
#
#
#
def get_node(conn, node_name:int, component: int, item:str='*'):
    """ """
    #with conn:
    data = get_item_table(conn, node_name, component, item)
    #system = get_coordinate_system(data[3])
    #return system(x=data[4], y=data[5], z=data[6],
    #              name=node_name, number=data[0], 
    #              index=data[10])
    #
    boundary = get_boundary(conn, node_id=data[2])
    #
    node = NodePoint(*data, boundary=boundary)
    return node.system()
#
def get_item_table(conn, node_name, component, item):
    """ """
    project = (node_name, component)
    sql = f'SELECT NodeCoordinate.{item}, \
                   Node.title, Node.mesh_idx \
            FROM Node, NodeCoordinate \
            WHERE Node.number =  NodeCoordinate.node_id\
            AND Node.name = ? \
            AND Node.component_id = ?'
    cur = conn.cursor()
    cur.execute(sql, project)
    record = cur.fetchone()
    return [*project, *record[1:]]
#
def get_nodes(conn, component):
    """ """
    query = (component, )
    table = 'SELECT * FROM Node \
            WHERE component_id = ?\
            ORDER BY number ASC'
    cur = conn.cursor()
    cur.execute(table, query)
    record = cur.fetchall()
    return record
#
#
#
def get_boundary(conn, node_id, item:str="*"):
    """
    """
    #
    project = (node_id,)
    sql = 'SELECT {:} FROM NodeBoundary WHERE node_id = ?'.format(item)
    cur = conn.cursor()
    cur.execute(sql, project)
    record = cur.fetchone()
    try:
        boundary = record[2:8]
    except TypeError:
        boundary = [1, 1, 1, 1, 1, 1]
    return boundary
#
#
#def update_table(conn, nodes):
#    """ """
#    # drop table
#    sql = 'DROP TABLE Node'
#    cur = conn.cursor()
#    cur.execute(sql)
#    #
#    new_node_table(conn)
#    push_nodes(conn, nodes)
#
#
#def push_nodes(conn, nodes):
#    """
#    Create a new project into the projects table
#    :param conn:
#    :param project:
#
#    :return: project id
#    """
#    project = nodes
#    #project = [item[1:] for item in nodes]
#    # number = len(self._labels) - 1
#    #if csystem == 'cylindrical':  # z, r, theta,
#    #    project = (node_name, csystem,
#    #               None, None, *coordinates, None)
#    #
#    #elif csystem == 'spherical':  # r, theta, phi
#    #    project = (node_name, csystem,
#    #               None, None, None, *coordinates)
#    #
#    #else:
#    #    project = (node_name, csystem,
#    #               *coordinates, None, None, None)
#    #
#    sql = 'INSERT INTO Node(name, type,\
#                                x, y, z, r, theta, phi)\
#                                VALUES(?,?,?,?,?,?,?,?)'
#    cur = conn.cursor()
#    cur.executemany(sql, project)
#
#
def update_colum(conn, colname: str, newcol: list):
    """update entire column values"""
    #query = (component, )
    table = f'UPDATE Node SET {colname} = ? \
            WHERE rowId = ? \
            AND component_id = ?;'
    cur = conn.cursor()
    cur.executemany(table, newcol)
#    
#