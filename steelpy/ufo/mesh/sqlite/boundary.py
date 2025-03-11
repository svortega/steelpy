# Copyright (c) 2009 steelpy

# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
import re
#from typing import NamedTuple


# package imports
from steelpy.ufo.utils.boundary import (BoundaryItem, BoundaryNode, get_support_df,
                                        get_node_boundary, #get_boundary_beam, find_boundary_type, 
                                        get_boundary_item,
                                        get_boundary_node,
                                        find_support_basic)

from steelpy.ufo.mesh.sqlite.utils import check_nodes, push_dircos
from steelpy.ufo.mesh.sqlite.element import ElementsSQL
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.ufo.mesh.sqlite.utils import pull_boundary # pull_node, 
from steelpy.utils.dataframe.main import DBframework
#
#
#
class BoundaryNodeSQL(BoundaryNode):
    
    __slots__ = ['_db_file', '_mesh_id','_name', '_elements']
    
    def __init__(self, name:str|int, mesh_id: int, db_file:str):
        """
        """
        super().__init__(name=name)
        self._db_file = db_file
        self._mesh_id = mesh_id
        self._elements = ElementsSQL(db_file=self._db_file,
                                     mesh_id=self._mesh_id,
                                     name=self._name)
    #
    def __setitem__(self, name: int,
                    values:list|tuple|dict) -> None:
        """
        name : boundary name int
        values = [node_id, boundary_type, fixity, title]
        """
        #1 / 0
        values = get_boundary_node(values, self._elements._beams)
        node_name = values[0]
        btype = values[1]
        fixity = values[2]
        unit_vector = values[3]
        title = values[4]
        #
        dircos_id = None
        # get nodes
        conn = create_connection(self._db_file)
        with conn:  
            node = check_nodes(conn, node_name,
                               mesh_id=self._mesh_id)
            node_id = node[0]
            #
            if unit_vector:          
                dircos_id = push_dircos(conn, unitvec=unit_vector,
                                        roll_angle=0, 
                                        mesh_id=self._mesh_id)             
        #
        try:
            # TODO : update data option if needed?
            self._labels.index(name)
            raise Warning(f' warning boundary {name} already exist')
        except ValueError:
            # input data
            conn = create_connection(self._db_file)
            with conn:
                boundary_id = push_boundary(conn, name,
                                            self._mesh_id, 
                                            btype=btype,
                                            dircosine_id=dircos_id,
                                            title=title)

                push_boundary_node(conn, boundary_id, fixity, title=title)
                update_node(conn, colname='boundary_id', item=boundary_id,
                            node_id=node_id, mesh_id=self._mesh_id)
        #
        #
        #print('-->')
    #
    def __getitem__(self, node_name: int) -> tuple | bool:
        """
        """
        #1 / 0
        conn = create_connection(self._db_file)
        with conn:  
            node = check_nodes(conn, node_name,
                               mesh_id=self._mesh_id)
            boundary = node[-1]
        #try:
        #    node_id = node[0]
        #except TypeError:
        #    raise IOError(f"Node {node_name} not found")
        if boundary:
        #try:
            self._labels.index(node_name)
            conn = create_connection (self._db_file)
            with conn:
                data = pull_nboundary(conn,
                                      mesh_id=self._mesh_id, 
                                      boundary_id=boundary)
                data = data[0]
                #data = pull_node_boundary(conn, node_id, self._mesh_id)

            return BoundaryItem(*data[4:10],
                                number=data[3],
                                name=data[0],
                                node=node_name,
                                boundary_type=data[1])
        #except ValueError:
        return False
    #
    #
    @property
    def _labels(self):
        """ """
        conn = create_connection(self._db_file)
        with conn:
            items = pull_node_boundaries(conn,
                                         mesh_id=self._mesh_id)
        return [item[1] for item in items]
    #
    #
    # -----------------------------------------
    #
    def _push_boundary(self, conn,
                       boundary_id: int, node_id: int,
                       fixity: list):
        """
        """
        try:
            title = fixity[6]
        except IndexError:
            title = None
        #except TypeError:
        #    title = None
        #    fixity = [0, 0, 0, 0, 0, 0] # free
        #
        push_boundary_node(conn, boundary_id, node_id, fixity[:6], title)
        #
        #project = (node_id, *fixity[:6], title, )
        #sql = 'INSERT INTO NodeBoundary(node_id,\
        #                                 x, y, z, rx, ry, rz,\
        #                                 title)\
        #                                 VALUES(?,?,?,?,?,?,?,?)'
        #cur = conn.cursor()
        #cur.execute(sql, project)
    #
    # -----------------------------------------
    #
    #
    # -----------------------------------------
    #
    @property
    def df(self):
        """Boundary df"""
        query = (self._mesh_id, )
        table = 'SELECT Mesh.name, Node.name, \
                 BoundaryNode.* \
                 FROM Node, BoundaryNode, Boundary, Mesh \
                 WHERE BoundaryNode.node_id = Node.number \
                 AND BoundaryNode.boundary_id = Boundary.number \
                 AND Node.mesh_id = Boundary.mesh_id \
                 AND Boundary.mesh_id = Mesh.number \
                 AND Mesh.number = ? ;'
        #
        conn = create_connection(self._db_file)
        with conn:
            cur = conn.cursor()
            cur.execute (table, query)
            rows = cur.fetchall()
        #
        db = DBframework()
        header = ['mesh_name', 'node_name', 'index', 'number', 'node_id',
                  'x', 'y', 'z', 'rx', 'ry', 'rz', 'title']
        boundf = db.DataFrame(data=rows, columns=header)
        header = ['mesh_name', 'node_name',
                  'x', 'y', 'z', 'rx', 'ry', 'rz',  'title']        
        return boundf[header]
    
    @df.setter
    def df(self, df):
        """ """
        df = get_support_df(df)
        for idx, row in df.iterrows():
            self.__setitem__(name=row['name'],
                             values=row.to_dict())
#
#
#
def pull_node_displacement(conn, mesh_id: int):
    """ """
    rows = []
    # node boundary
    fixity = pull_boundary(conn,
                           mesh_id=mesh_id,
                           boundary_type='constraint')
    #
    rows = [[*item[:2], *item[6:12], item[12]]
            for item in fixity]
    #
    #for item in fixity:
    #    temp = [0 if ndp ==1 else None
    #            for ndp in item[6:12]]
    #    rows.append([*item[:2], *temp, item[12]])
    #
    # node displacement
    #disp = pull_boundary(conn,
    #                     mesh_id=mesh_id,
    #                     boundary_type='displacement')
    #
    #rows = []
    #for item in disp:
    #    temp = [ndp if ndp else 0
    #            for ndp in item[6:12]]
    #    rows.append([*item[:2], *temp, item[12]])   
    #
    cols = ['node_index', 'node_name',  
            #'boundary_name', 'boundary_type', 
            #'number', 'boundary_id', 
            'x', 'y', 'z', 'rx', 'ry', 'rz', 'title']
    # dataframe
    db = DBframework()
    df = db.DataFrame(data=rows, columns=cols)
    #df.drop(columns=['boundary_name', 'boundary_type', 
    #                 'number', 'boundary_id'],
    #        inplace=True)
    df['load_level'] = 'basic'
    return df    
    
#
#
#
def pull_node_boundary(conn, node_id: int,
                       mesh_id: int, 
                       item:str="*"):
    """
    """
    #
    project = (node_id, mesh_id, )
    sql = f'SELECT Boundary.name, BoundaryNode.{item} \
            FROM Node, BoundaryFixity, Boundary, Mesh \
            WHERE BoundaryNode.node_id = ? \
            AND Node.number = BoundaryNode.node_id \
            AND BoundaryNode.boundary_id = Boundary.number \
            AND Node.mesh_id = Boundary.mesh_id \
            AND Mesh.number = ?'
    cur = conn.cursor()
    cur.execute(sql, project)
    record = cur.fetchone()
    return record
#
def push_boundary_node(conn, boundary_id: int, # node_id: int, 
                       fixity: list,
                       title: str|None = None):
    """ """
    project = (boundary_id, *fixity, title, )
    sql = 'INSERT INTO BoundaryFixity(boundary_id,\
                                      x, y, z, rx, ry, rz,\
                                      title)\
                                    VALUES(?,?,?,?,?,?,?,?)'
    cur = conn.cursor()
    cur.execute(sql, project)
#
#
class BoundarySQL:
    __slots__ = ['_db_file', '_mesh_id', '_boundary','_name']
    
    def __init__(self, db_file:str,
                 mesh_id: int, 
                 name:str|int,
                 db_system:str="sqlite") -> None:
        """
        """
        self._db_file = db_file
        self._mesh_id = mesh_id
        #
        # create node table
        conn = create_connection(self._db_file)
        with conn:        
            self._new_table(conn)        
        #
        self._boundary = BoundaryNodeSQL(name=name,
                                         mesh_id=self._mesh_id,
                                         db_file=self._db_file)
        #
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._mesh_id, )
        table = 'SELECT Boundary.name \
                 FROM Boundary, Mesh \
                 WHERE Boundary.mesh_id = Mesh.number \
                 AND Mesh.number = ? ;'
        
        conn = create_connection(self._db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]    
    #
    #
    def __setitem__(self, name: int|str,
                    values: tuple|list|dict) -> None:
        """
        """
        btype = get_boundary_item(values)
        try:
            self._labels.index(name)
            raise IOError(f' error boundary {name} already exist')
        except ValueError:
            match btype:
                case 'node'|'beam':
                    self._boundary[name] = values #[1:]
                case _:
                    raise IOError(f'boundary {btype} not valid')
    #
    def __getitem__(self, name: int|str):
        """
        """
        1/0
        try:
            index = self._labels.index(name)
            boundary_type = self._type[index]
            #b_number = self._number[index]
        except ValueError:
            raise KeyError(f'Invalid boundary name: {name}')
        #
        if re.match(r"\b(node(s)?|coord(inate)?)\b", boundary_type, re.IGNORECASE):
            return self._boundary[name]

        else:
            raise IOError(f' boundary type {boundary_type} not recognised')    
    #
    #
    def __str__(self, units:str="si") -> str:
        """ """
        lenght = ' m'
        space = " "
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += f"{33*space}BOUNDARIES\n"
        output += "\n"
        output += (f"Boundary {11*space} Node {6*space} x {6*space} y {6*space} z {5*space} mx {5*space} my {5*space} mz {5*space}")
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for key, node in self._boundary.items():
            #if sum(node[:6]) == 0:
            #    continue
            output += node.__str__()
        return output
    #    
    # -----------------------------------------
    #
    def _new_table(self, conn) -> None:
        """ """
        #
        table = "CREATE TABLE IF NOT EXISTS Boundary(\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    mesh_id INTEGER NOT NULL REFERENCES Mesh(number), \
                    type TEXT NOT NULL,\
                    dircosine_id INTEGER REFERENCES DirCosines(number),\
                    title TEXT);"
        #
        create_table(conn, table)
        #
        table = "CREATE TABLE IF NOT EXISTS BoundaryFixity(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            boundary_id INTEGER NOT NULL REFERENCES Boundary(number),\
                            x DECIMAL, y DECIMAL, z DECIMAL,\
                            rx DECIMAL, ry DECIMAL, rz DECIMAL, \
                            title TEXT);"
        # node_id INTEGER NOT NULL REFERENCES Node(number),\
        create_table(conn, table)
    #
    def _push_boundary(self, conn, name: int|str, btype: str):
        """
        """
        return push_boundary(conn, name, self._mesh_id, btype)
        
    #
    # -----------------------------------------
    #
    #
    def support(self, values:list|tuple|dict|None = None):
        """
        values: [node_id, 'restrain'|'spring'|'matrix', fixity,
                 'skew'|('beam'|'local'), dircos|beam_id]
        """
        if isinstance(values, (list|tuple)):
            if isinstance(values[0], (list|tuple)):
                for item in values:
                    self._boundary[item[0]] = item[1:]
            elif isinstance(values[0], dict):
                for item in values:
                    bname = item.pop('name')
                    self._boundary[bname] = item

            else:
                self._boundary[values[0]] = values[1:]

        elif isinstance(values, dict):
            bname = values['name']
            if isinstance(bname, (list | tuple)):
                db = DBframework()
                df = db.DataFrame(values)
                self._boundary.df = df
            else:
                sname = values.pop('name')
                self._boundary[sname] = values
        #
        return self._boundary
    #
    #
    @property
    def _nodesX(self):
        """ """
        conn = create_connection(self._db_file)
        with conn:
            nodes =  pull_node_boundaries(conn, mesh_id=self._mesh_id)
            bnid = [row[-1] for row in nodes]            
            items = pull_nboundary(conn, boundary_id=bnid)
                                   #mesh_id=self._mesh_id)
        #
        bnodes = {nodes[x][1]: BoundaryItem(*data[3:9],
                                            number=data[2],
                                            name=data[0],
                                            node=nodes[x][1])
                      for x, data in enumerate(items)}            
        return bnodes
    #
    # -----------------------------------------
    #
    @property
    def df(self):
        """Boundary df"""
        #
        conn = create_connection(self._db_file)
        with conn:        
            nodes =  pull_node_boundaries(conn, mesh_id=self._mesh_id)
        bnid = [row[-1] for row in nodes]
        #
        project = "', '".join(str(x) for x in bnid)
        query = (self._mesh_id, )
        table = f"SELECT Mesh.name, Node.name, Node.idx,\
                  Boundary.type, BoundaryFixity.* \
                  FROM Node, Boundary, BoundaryFixity, Mesh\
                  WHERE BoundaryFixity.boundary_id = Boundary.number \
                  AND Node.boundary_id = Boundary.number \
                  AND Boundary.number IN ('{project}') \
                  AND Node.mesh_id = Boundary.mesh_id \
                  AND Boundary.mesh_id = Mesh.number \
                  AND Mesh.number = ?;"
        #
        conn = create_connection(self._db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            rows = cur.fetchall()
        #
        db = DBframework()
        header = ['mesh_name', 'node_name', 'node_index',
                  'type', 
                  'number', 'boundary_id',
                  'x', 'y', 'z', 'rx', 'ry', 'rz', 'title']
        boundf = db.DataFrame(data=rows, columns=header)
        # TODO: global should not be fixed
        boundf['system'] = 'global'
        header = ['mesh_name', 'type','system', 'node_name', 'node_index',  
                  'x', 'y', 'z', 'rx', 'ry', 'rz',  'title']        
        return boundf[header]

    @df.setter
    def df(self, df):
        """nodes in dataframe format"""
        #beams = self._boundary._elements._beams
        #columns = list(df.columns)
        #header = {key: find_support_basic(key)
        #          for key in columns}
        #df.rename(columns=header, inplace=True)
        #
        df = get_support_df(df) 
        for idx, row in df.iterrows():
            self.__setitem__(name=row['name'],
                             values=row.to_dict())
        #
#
#
def push_boundary(conn, name: int, mesh_id: int,
                  btype: str, dircosine_id: int|None,
                  title: str|None = None):
    """ """
    project = (name, mesh_id, btype, dircosine_id, title, )
    sql = 'INSERT INTO Boundary(name, mesh_id, type, \
                                dircosine_id, title) \
                                VALUES(?,?,?,?,?)'
    cur = conn.cursor()
    cur.execute(sql, project)
    row = cur.lastrowid
    return row
#
#
def pull_nboundary(conn, mesh_id: int,
                   boundary_id: int|list,
                   item:str="*"):
    """
    """
    #nodes =  pull_node_boundaries(conn, mesh_id)
    #bnid = [row[-1] for row in nodes]
    # 
    #if boundary_id in ['*', None]:
    #    project = bnid
    #else:
    #    project = [boundary_id]
    #
    query = (mesh_id, )
    #
    if isinstance(boundary_id, int):
        boundary_id = [boundary_id]
    #
    project = "', '".join(str(x) for x in boundary_id)
    #table = f"SELECT {item} \
    #          FROM Boundary WHERE number in ('{project}')"
    table = f"SELECT Boundary.name, Boundary.type,\
              BoundaryFixity.{item} \
              FROM Boundary, BoundaryFixity\
              WHERE BoundaryFixity.boundary_id = Boundary.number \
              AND mesh_id = ? \
              AND Boundary.number IN ('{project}') ;"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    record = cur.fetchall()
    #
    #bnodes = {nodes[x][1]: BoundaryItem(*data[3:9],
    #                                    number=data[2],
    #                                    name=data[0],
    #                                    node=nodes[x][1])
    #          for x, data in enumerate(record)}
    #
    return record
#
def pull_node_boundaries(conn, mesh_id: int):
    """ """
    project = (mesh_id, )
    table = 'SELECT * FROM Node \
             WHERE boundary_id IS NOT NULL \
             AND mesh_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    nodes = cur.fetchall()
    return nodes
#
def get_boundary_list(items: list):
    """
    [node, fixity]
    """
    node_id = None
    dircos = None
    number = len(items)
    for x in range(0, number, 2):
        print(items[x])
        if re.match(r"\b(node(s)?)\b", items[x], re.IGNORECASE):
            node_id = items[x+1]
            
        elif re.match(r"\b(support|constrain|fixity|suppression)\b", items[x], re.IGNORECASE):
            fixity = get_node_boundary(items[x+1])
            btype = 'restrain'
        
        elif re.match(r"\b(spring(s)?)\b", items[x], re.IGNORECASE):
            btype = 'spring'
            1 / 0
        
        elif re.match(r"\b(beam(s)?)\b", items[x], re.IGNORECASE):
            bname = items[x+1]
            #beam = BeamItemSQL(beam_name=bname,
            #                   mesh_id=self._mesh_id,
            #                   db_file=self.db_file)
            #
            #dircos = beam.dircos
            1 / 0
        
        elif re.match(r"\b((beam_)?end)\b", items[x], re.IGNORECASE):
            end = items[x+1]
            #nodes = beam.nodes
            #node_id = nodes[end].name
            1 / 0
            
        elif re.match(r"\b(translation)\b", items[x], re.IGNORECASE):
            dircos = items[x+1]
    #
    if not node_id:
        raise IOError('node id not found')
    #1 / 0
    return btype, node_id, fixity, dircos
#
#
def update_node(conn, colname: str,
                item: int|str|float,
                node_id: int,
                mesh_id: int):
    """ """
    project = (item, node_id, mesh_id)
    table = f'UPDATE Node SET {colname} = ? \
              WHERE number = ? \
              AND mesh_id = ?;'
    cur = conn.cursor()
    #cur.executemany(table, item, node_id, mesh_id)
    cur.execute(table, project)
    #print('--')
#
#
#
#