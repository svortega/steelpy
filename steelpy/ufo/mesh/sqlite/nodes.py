# Copyright (c) 2009 steelpy

# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
from itertools import chain, count
#from math import isclose, dist
from typing import NamedTuple
import re

# package imports
from steelpy.utils.math.operations import zeros, to_matrix
from steelpy.utils.sqlite.utils import create_connection, create_table
#
from steelpy.ufo.process.elements.nodes import NodePoint, NodeBasic
from steelpy.ufo.process.elements.boundary import BoundaryItem
#
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
        self.db_file = db_file
        self._plane = plane
        self._component =  component
        # create node table
        conn = create_connection(self.db_file)
        with conn:
            self._create_table(conn)
        #print('--> update labels?')
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
            coord = self.get_coordinates(coordinates)
            #
            try:
                node_id = coordinates.name
            except AttributeError:
                node_id = node_number
            #
            conn = create_connection(self.db_file)
            with conn:
                self._push_node(conn,
                                node_number=node_number,
                                node_name=node_id, 
                                coordinates=coord)
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
                node = self._pull_node(conn, node_number)
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
                title NOT NULL,\
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
                   node_name: int | str,
                   coordinates: list|tuple):
        """
        Create a new project into the projects table
        node_number: int
        node_name: int|str
        coordinate: list|tuple [x,y,z,r,theta,phi]
        """
        # get row number
        cur = conn.cursor()
        sql = 'SELECT max(number) from Node;'
        if (idx := max(cur.execute(sql))[0]) == None:
            idx = 0
        #
        # -------------------------------------------
        #
        query = (node_number, self._component, node_name, idx)
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
    def _pull_node(self, conn, node_number: int,
                   item: str = '*'):
        """ """
        return pull_node(conn, node_name=node_number,
                         component=self._component,
                         item=item)
    #
    def _isclose(self, key:str, item:str, value:float,
                    rel_tol:float=1e-6, abs_tol:float=0.0)-> tuple:
        """ """
        query = (self._component, )
        table = f'SELECT NodeCoordinate.{key}, Node.title \
                  FROM Node, NodeCoordinate \
                  WHERE Node.component_id = ? \
                  AND Node.name = NodeCoordinate.node_id \
                  AND ABS({item} - {value}) <= MAX({rel_tol} * MAX(ABS({item}), ABS({value})), {abs_tol})'
            #.format(key, item, value, rel_tol, item, value, abs_tol)
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
    def renumbering(self, new_numbers:list[int]):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            # [row number, node name]
            rows = pull_node_rows(conn, self._component)
            rows = {item[1]: item[0] for item in rows}

            indexes = [(x, rows[node_name], self._component,)
                       for x, node_name in enumerate(new_numbers)]

            update_colum(conn, colname='mesh_idx',
                         newcol=indexes)
    #
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
    def get_name(self, title:int|str):
        """get node name based on title"""
        conn = create_connection(self.db_file)
        with conn:
            name = pull_node_name(conn, title,
                                  self._component)
        return name
    #
    #
    # ---------------------------------
    #
    #@property
    def jbc(self):
        """ joints with boundary"""
        nnp = len(self._labels)
        jbc = zeros(nnp, 6, code='I')
        #
        #for item in self._labels:
        #    node = self.__getitem__(item)
        #    if node.boundary:
        #        ind = node.index
        #        jbc[ind] = node.boundary[:6]
        #         
        #
        conn = create_connection(self.db_file)
        with conn:        
            fixity = pull_boundary(conn,
                                   component=self._component)
            #
            nidx = pull_node_index(conn, self._component)
            nidx = {item[0]: item[1] for item in nidx}
        #
        # update jbc
        #
        for item in fixity:
            # [node_idx] = [x,y,z,rx,ry,rz]
            jbc[item[0]] = item[4:10]
            #node_name.append(item[1])
        #      
        #
        db = DBframework()
        #
        #jbc = to_matrix(jbc, 6)
        df_jbc = db.DataFrame(data=jbc,
                              columns=['x', 'y', 'z', 'rx', 'ry', 'rz'])
        jbc = df_jbc[self._plane.dof].values.tolist()
        jbc = list(chain.from_iterable(jbc))
        #
        counter = count(start=1)
        jbc = [next(counter) if item == 0 else 0
               for item in jbc]
        # update jbc
        #self._jbc = array('I', jbc)
        #
        #jbc = to_matrix(self._jbc, self._plane.ndof)
        jbc = to_matrix(jbc, self._plane.ndof)
        df_jbc = db.DataFrame(data=jbc, columns=self._plane.dof)
        #node_name = list(self._labels)
        df_jbc['node_name'] = [nidx[idx]
                               for idx, item in enumerate(jbc)]
        df_jbc = df_jbc.set_index('node_name', drop=True)    
        # remove rows with zeros
        #df_jbc = df_jbc[df_jbc.any(axis=1)]
        #df_jbc.replace(0, self._db.nan, inplace=True)
        #df_jbc = df_jbc.notnull()        
        return df_jbc
    #
    def jwbc(self):
        """ joints with boundary condition"""
        project = (self._component, )
        table = 'SELECT Node.mesh_idx, \
                 IIF(Node.number = NodeBoundary.node_id,\
                     (NodeBoundary.x, NodeBoundary.y, NodeBoundary.z, \
                     NodeBoundary.rx, NodeBoundary.ry, NodeBoundary.rz) ,\
                     (0,0,0,0,0,0)) AS JBC \
                FROM Node, NodeBoundary \
                WHERE Node.component_id = ?'
        #table = 'SELECT Node, NodeBoundary \
        #            CASE \
        #                WHEN Node.number = NodeBoundary.node_id \
        #                    THEN Node.mesh_idx,  NodeBoundary.*\
        #                ELSE: Node.mesh_idx \
        #            END \
        #        FROM Node, NodeBoundary \
        #        WHERE Node.component_id = ?'
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, project)
            data = cur.fetchall()
        data
        return data        
    #
    def neq(self):
        """ Number the equations  in jbc from 1 up to the order.
           Start assigning equation numbers for zero dof's
           from 1 up;  only zero given a number. """
        #
        if not self._jbc:
            self.jbc()
        neq = max(self._jbc)
        #jbc = to_matrix(self._jbc, 6)
        return neq
    #
    #
    def newXXX(self):
        """ """
        jbcc = self.jbc().stack().values
        index = list(reversed([i for i, item in enumerate(jbcc)
                               if item == 0]))
        return index    
    #
    def DOF_unreleased(self):
        """ list of the indices for the unreleased DOFs """
        project = (self._component, )
        table = 'SELECT Node.mesh_idx, Node.name, \
                        NodeBoundary.* \
                 FROM Node, NodeBoundary \
                 WHERE Node.number = NodeBoundary.node_id \
                 AND Node.component_id = ? \
                 ORDER BY Node.mesh_idx ASC'
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, project)
            data = cur.fetchall()
        #
        #dof = 6
        #new = []
        #for row in data:
        #    idx = row[0]
        #    ndof = idx * dof
        #    new.extend([ndof + x + 1
        #                for x, item in enumerate(row[4:-1])
        #                if item != 0])
        #
        dof =  self._plane.ndof
        new = [[row[0] * dof + x # + 1
                for x, item in enumerate(row[4:-1]) if item == 0]
                for row in data]
        new = list(chain(*new))
        return new
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
def pull_node(conn, node_name:int, component: int, item:str='*'):
    """ """
    #with conn:
    data = pull_node_item(conn, node_name, component, item)
    #system = get_coordinate_system(data[3])
    #return system(x=data[4], y=data[5], z=data[6],
    #              name=node_name, number=data[0],
    #              index=data[10])
    #
    boundary = pull_node_boundary(conn,
                                  node_name=node_name,
                                  component=component)
    #
    node = NodePoint(*data, boundary=boundary)
    return node.system()
#
def pull_node_item(conn, node_name, component, item):
    """ """
    project = (node_name, component)
    table = f'SELECT NodeCoordinate.{item}, \
                   Node.title, Node.mesh_idx \
            FROM Node, NodeCoordinate \
            WHERE Node.number =  NodeCoordinate.node_id\
            AND Node.name = ? \
            AND Node.component_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    record = cur.fetchone()
    return [*project, *record[1:]]
#
def pull_node_boundary(conn, node_name: int|str,
                       component: int, item:str="*"):
    """
    """
    #
    #project = (node_name, component, )
    #table = 'SELECT NodeBoundary.{:} \
    #         FROM Node, NodeBoundary \
    #         WHERE Node.number = NodeBoundary.node_id \
    #         AND Node.name = ? \
    #         AND Node.component_id = ?'.format(item)
    #cur = conn.cursor()
    #cur.execute(table, project)
    #data = cur.fetchone()
    #
    data = pull_boundary(conn, component,
                         node_name, item)
    try:
        data = data[0]
        boundary = BoundaryItem(*data[4:10],
                                number=data[2],
                                name=data[3],
                                node=node_name)
    except IndexError:
        boundary = None
    #
    return boundary
#
def pull_boundary(conn, component: int,
                  node_name: int|str|None = None,
                  item:str="*"):
    """pull all boundary data"""
    #
    project = [component]
    #
    table = 'SELECT Node.mesh_idx, Node.name, \
             NodeBoundary.{:} \
             FROM Node, NodeBoundary \
             WHERE Node.number = NodeBoundary.node_id \
             AND Node.component_id = ?'.format(item)
    #
    if node_name:
        table += 'AND Node.name = ?'
        project.extend([node_name])
    # 
    cur = conn.cursor()
    cur.execute(table, tuple(project))
    data = cur.fetchall()
    return data
#
def pull_node_number(conn, node_name:int,
                     component: int):
    """ """
    project = (node_name, component)
    table = f'SELECT number FROM Node \
              WHERE name = ? \
              AND component_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    record = cur.fetchone()
    return record[0]
#
def pull_node_rows(conn, component: int):
    """ """
    project = (component,)
    table = f'SELECT number, name FROM Node \
              WHERE component_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    records = cur.fetchall()
    return records
#
def pull_node_index(conn, component: int):
    """ """
    project = (component,)
    table = f'SELECT mesh_idx, name FROM Node \
              WHERE component_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    records = cur.fetchall()
    return records
#
def pull_node_name(conn, node_title: int|str,
                   component: int,):
    """ """
    project = (node_title, component)
    table = f'SELECT name FROM Node \
              WHERE title = ? \
              AND component_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    record = cur.fetchone()
    return record[0]
#
#
# ---------------------------------
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
# ---------------------------------
#
#