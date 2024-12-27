# Copyright (c) 2009 steelpy

# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
from itertools import chain, count
#from math import isclose, dist
#from typing import NamedTuple
import re

# package imports
from steelpy.utils.math.operations import zeros, to_matrix
from steelpy.utils.sqlite.utils import create_connection, create_table
#
from steelpy.ufo.utils.node import NodePoint, NodeBasic, get_node_df
from steelpy.ufo.utils.boundary import BoundaryItem, get_node_boundary

from steelpy.ufo.mesh.sqlite.boundary import push_boundary, push_boundary_node
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
    __slots__ = ['_system', 'db_file', '_mesh_id', '_name']

    def __init__(self,
                 db_file: str,
                 name: str|int,
                 mesh_id: str|int, 
                 db_system:str="sqlite",
                 system:str = 'cartesian') -> None:
        """
        """
        super().__init__(name=name, system=system)
        self.db_file = db_file
        #self._plane = plane
        self._mesh_id =  mesh_id
        # create node table
        conn = create_connection(self.db_file)
        with conn:
            self._new_table(conn)
        #print('--> update labels?')
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._mesh_id, )
        table = 'SELECT Node.name FROM Node \
                 WHERE mesh_id = ? \
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
    def __setitem__(self, name: int|str,
                    coordinates: tuple|list|dict) -> None:
        """
        """
        #
        if not isinstance(name, int):
            raise IOError('node id must me an integer')
        #
        try:
            self._labels.index(name)
            raise Exception(f' warning node {node_number} already exist')
        except ValueError:
            data = self.get_coordinates(coordinates)
            coord = data[:3]
            fixity = data[3]
            if (node_title := data[4]) is None:
                node_title = name
            #node_number = self.get_number()
            #
            conn = create_connection(self.db_file)
            with conn:
                bid = None
                if fixity:
                    mesh_id = self._mesh_id
                    bid = push_node_boundary(conn, name,
                                             mesh_id, fixity)
                
                self._push_node(conn,
                                node_name=name, 
                                coordinates=coord,
                                boundary=bid,
                                node_title=node_title)
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
    def _new_table(self, conn) -> None:
        """ """
        # Node main
        table = "CREATE TABLE IF NOT EXISTS Node (\
                number INTEGER NOT NULL,\
                name INTEGER NOT NULL,\
                mesh_id INTEGER NOT NULL REFERENCES Mesh(number), \
                title NOT NULL,\
                idx INTEGER NOT NULL, \
                boundary_id INTEGER REFERENCES Boundary(number), \
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
    def _push_node(self, conn, 
                   node_name: int | str,
                   coordinates: list|tuple,
                   boundary: int|None = None,
                   node_title: str|None = None):
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
        query = (node_name, self._mesh_id, node_title,
                 idx, boundary)
        table = 'INSERT INTO Node(name, mesh_id, \
                                  title, idx, boundary_id) \
                                  VALUES(?,?,?,?,?);' 
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
                         mesh_id=self._mesh_id,
                         item=item)
    #
    def _isclose(self, item:str, value:float, key:str = '*', 
                 rel_tol:float=1e-6, abs_tol:float=0.0)-> tuple:
        """ """
        query = (self._mesh_id, )
        table = f'SELECT NodeCoordinate.{key}, Node.title \
                  FROM Node, NodeCoordinate \
                  WHERE Node.mesh_id = ? \
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
        project = (value, name, self._mesh_id, )
        sql = f'UPDATE Node SET {item} = ? \
                WHERE name = ? AND mesh_id = ? '

        cur = conn.cursor()
        cur.execute(sql, project)
    #
    def _orderby(self):
        """
        """
        1 / 0
        query = (self._mesh_id, )
        conn = create_connection(self.db_file)
        table = 'SELECT * FROM Node \
                WHERE mesh_id = ? \
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
            rows = pull_node_rows(conn, self._mesh_id)
            rows = {item[1]: item[0] for item in rows}

            indexes = [(x, rows[node_name], self._mesh_id,)
                       for x, node_name in enumerate(new_numbers)]

            update_colum(conn, colname='idx',
                         newcol=indexes)
    #
    #
    def get_maxmin(self):
        """
        """
        def maxmin(head: str, col: str):
            #
            query = (self._mesh_id, )
            table = f'SELECT {head.upper()}(NodeCoordinate.{col}) \
                      FROM Node, NodeCoordinate \
            WHERE Node.mesh_id = ? \
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
            name = pull_node_by_title(conn, title,
                                      self._mesh_id)
        return name
    #
    def get_number(self, start: int = 1):
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            #number = pull_node_rows(conn, self._mesh_id)
            #
            query = (self._mesh_id, )
            table = 'SELECT max(number) from Node WHERE mesh_id = ?;'
            cur = conn.cursor()
            cur.execute(table, query)
            records = cur.fetchone()
        #
        if (idx := max(records)) == None:
            idx = 0
        #
        #if not number:
        #    return 1
        #else:
        #    number = [item[0] for item in number]
        #    return max(number) + 1
        #return idx + 1
        #
        n = idx + 1
        while True:
            yield n
            n += 1        
    #
    # ---------------------------------
    #
    #@property
    def jbc(self, dof: list[str]):
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
                                   mesh_id=self._mesh_id)
            #
            nidx = pull_node_index(conn, self._mesh_id)
            nidx = {item[0]: item[1] for item in nidx}
        #
        # update jbc
        #
        for item in fixity:
            # [node_idx] = [x,y,z,rx,ry,rz]
            jbc[item[0]] = item[5:11]
            #node_name.append(item[1])
        #      
        #
        db = DBframework()
        #
        #jbc = to_matrix(jbc, 6)
        df_jbc = db.DataFrame(data=jbc,
                              columns=['x', 'y', 'z', 'rx', 'ry', 'rz'])
        jbc = df_jbc[dof].values.tolist()
        jbc = list(chain.from_iterable(jbc))
        #
        counter = count(start=1)
        jbc = [next(counter) if item == 0 else 0
               for item in jbc]
        # update jbc
        #self._jbc = array('I', jbc)
        #
        #jbc = to_matrix(self._jbc, self._plane.ndof)
        ndof = len(dof)
        jbc = to_matrix(jbc, ndof)
        df_jbc = db.DataFrame(data=jbc, columns=dof)
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
        project = (self._mesh_id, )
        table = 'SELECT Node.idx, \
                 IIF(Node.number = NodeBoundary.node_id,\
                     (NodeBoundary.x, NodeBoundary.y, NodeBoundary.z, \
                     NodeBoundary.rx, NodeBoundary.ry, NodeBoundary.rz) ,\
                     (0,0,0,0,0,0)) AS JBC \
                FROM Node, NodeBoundary \
                WHERE Node.mesh_id = ?'
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
        project = (self._mesh_id, )
        table = 'SELECT Node.idx, Node.name, \
                        NodeBoundary.* \
                 FROM Node, NodeBoundary \
                 WHERE Node.number = NodeBoundary.node_id \
                 AND Node.mesh_id = ? \
                 ORDER BY Node.idx ASC'
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
        table = f'SELECT Node.name, Node.mesh_id, NodeCoordinate.* \
                  FROM Node, NodeCoordinate \
                  WHERE Node.mesh_id = {self._mesh_id} \
                  AND Node.number = NodeCoordinate.node_id ;'
        nodes = db.read_sql_query(table,conn)
        nodes.drop(columns=['node_id'], inplace=True)
        return nodes
    
    @df.setter
    def df(self, df):
        """ """
        nodes = get_node_df(df, system=self._system)
        nodes['mesh_id'] = self._mesh_id
        # TODO : fix boundary
        if 'boundary' in nodes:
            fixities = {nodes['name'][x]: get_node_boundary(item)
                        if item else None
                        for x, item in enumerate(nodes['boundary'])}
            conn = create_connection(self.db_file)
            with conn:
                blist = []
                for key, item in fixities.items():
                    if not item:
                        blist.append(None)
                        continue
                    blist.append(push_node_boundary(conn, node_id=int(key),
                                                     mesh_id=self._mesh_id,
                                                     fixity=item))
            nodes['boundary_id'] = blist
        else:
            nodes['boundary_id'] = None
        #
        #
        idx = next(self.get_number()) - 1
        nodes['idx'] = [item + idx for item in list(nodes.index)]
        #
        conn = create_connection(self.db_file)
        with conn:
            cols = ['name', 'mesh_id',
                    'title', 'idx', 'boundary_id']
            nodes[cols].to_sql('Node', conn,
                               index_label=cols, 
                               if_exists='append', index=False)
            #
            node_id = pull_node_rows(conn, mesh_id=self._mesh_id)
            node_id = {item[1]: item[0] for item in node_id}
            #
            cols = ['node_id', 'system',
                    'x', 'y', 'z',
                    'r', 'theta', 'phi']
            nodes['node_id'] = [node_id[item] for item in list(nodes.name)]
            nodes[cols].to_sql('NodeCoordinate', conn,
                               index_label=cols, 
                               if_exists='append', index=False)
        #
        #print('-->')
        #
#
#
#
def pull_node(conn, node_name:int|str, mesh_id: int, item:str='*'):
    """ """
    data = pull_node_item(conn, node_name, mesh_id, item)
    boundary = pull_node_boundary(conn,
                                  node_name=node_name,
                                  mesh_id=mesh_id)
    #
    node = NodePoint(*data, boundary=boundary)
    return node.system()
#
def pull_node_item(conn, node_name:int|str, mesh_id: int, item:str='*'):
    """ """
    project = (node_name, mesh_id)
    table = f'SELECT NodeCoordinate.{item}, \
                     Node.title, Node.idx \
            FROM Node, NodeCoordinate \
            WHERE Node.number =  NodeCoordinate.node_id\
            AND Node.name = ? \
            AND Node.mesh_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    record = cur.fetchone()
    return [*project, *record[1:]]
#
def pull_node_boundary(conn, node_name: int|str,
                       mesh_id: int, item:str="*"):
    """
    """
    data = pull_boundary(conn, mesh_id,
                         node_name, item)
    
    try:
        data = data[0]
        boundary = BoundaryItem(*data[5:11],
                                number=data[4],
                                name=data[2],
                                node=node_name)
    except IndexError:
        boundary = None
    #
    return boundary
#
def pull_boundary(conn, mesh_id: int,
                  node_name: int|str|None = None,
                  item:str="*"):
    """pull all boundary data"""
    #
    project = [mesh_id]
    #
    table = f'SELECT Node.idx, Node.name, \
             Boundary.name, \
             BoundaryFixity.{item} \
             FROM Node, BoundaryFixity, Boundary, Mesh \
             WHERE Node.boundary_id = Boundary.number \
             AND BoundaryFixity.boundary_id = Boundary.number \
             AND Node.mesh_id = Boundary.mesh_id \
             AND Mesh.number = ?'
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
def pull_node_number(conn, node_name:int|str,
                     mesh_id: int):
    """ """
    project = (node_name, mesh_id)
    table = f'SELECT number FROM Node \
              WHERE name = ? \
              AND mesh_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    record = cur.fetchone()
    if not record:
        raise IOError(f'node {node_name} not found')
    return record[0]
#
def pull_node_rows(conn, mesh_id: int):
    """ """
    project = (mesh_id,)
    table = f'SELECT number, name FROM Node \
              WHERE mesh_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    records = cur.fetchall()
    return records
#
def pull_node_index(conn, mesh_id: int):
    """ """
    project = (mesh_id,)
    table = f'SELECT idx, name FROM Node \
              WHERE mesh_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    records = cur.fetchall()
    return records
#
def pull_node_by_title(conn, node_title: int|str,
                       mesh_id: int,):
    """ """
    project = (node_title, mesh_id)
    table = f'SELECT name FROM Node \
              WHERE title = ? \
              AND mesh_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    record = cur.fetchone()
    return record[0]
#
#
# ---------------------------------
#
def push_node_boundary(conn, node_id: int, mesh_id: int,
                       fixity: list, btype: str = 'restrain',
                       title: str|None=None):
    """ """
    boundary_id = push_boundary(conn, name=node_id,
                                mesh_id=mesh_id,
                                btype=btype,
                                title=None)
    
    push_boundary_node(conn, boundary_id,
                       fixity, title)
    #1 / 0
    return boundary_id
#
# ---------------------------------
#
def update_colum(conn, colname: str, newcol: list):
    """update entire column values"""
    #query = (mesh_id, )
    table = f'UPDATE Node SET {colname} = ? \
              WHERE rowId = ? \
              AND mesh_id = ?;'
    cur = conn.cursor()
    cur.executemany(table, newcol)
#
# ---------------------------------
#
#