# Copyright (c) 2009 steelpy

# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
import re
#from typing import NamedTuple


# package imports
from steelpy.ufo.mesh.elements.boundary import BoundaryItem, BoundaryNode
from steelpy.ufo.mesh.sqlite.utils import check_nodes
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework
#
#
#
#
class BoundaryNodeSQL(BoundaryNode):
    
    __slots__ = ['_db_file', '_component']
    
    def __init__(self, component: int, db_file:str):
        """
        """
        super().__init__()
        self._db_file = db_file
        self._component = component
        # create node table
        conn = create_connection(self._db_file)
        with conn:        
            self._create_table(conn)
    #
    def __setitem__(self, node_name: int,
                    fixity:list|tuple|dict|str) -> None:
        """
        node name : int
        fixity = [x, y, z, rx, ry, rz, title]
        """
        conn = create_connection(self._db_file)
        with conn:  
            node = check_nodes(conn, node_name,
                               component=self._component)       
        try:
            node_id = node[0]
        except TypeError:
            raise IOError(f"Node: {node_name} not valid")
        
        try:
            # TODO : update data option if needed?
            self._labels.index(node_name)
            raise Warning(f' warning node {node_name} boundary already exist')
        except ValueError:
            #self._labels.append(node_name)
            fixity = self._get_fixity(fixity)
            conn = create_connection(self._db_file)
            with conn:
                self._push_boundary(conn, node_id, fixity)
        #
    #
    def __getitem__(self, node_name: int) -> tuple | bool:
        """
        """
        conn = create_connection(self._db_file)
        with conn:  
            node = check_nodes(conn, node_name,
                               component=self._component)       
        try:
            node_id = node[0]
        except TypeError:
            raise IOError(f"Node {node_name} not found")
        
        try:
            self._labels.index(node_name)
            conn = create_connection (self._db_file)
            data = get_boundary(conn, node_id)
            return BoundaryItem(*data[2:8],
                                number=data[0],
                                name=data[1],
                                node=node_name)
        except ValueError:
            return False
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._component, )
        table = 'SELECT Node.name \
                 FROM Node, NodeBoundary \
                 WHERE NodeBoundary.node_id = Node.number \
                 AND Node.component_id = ? ;'
        
        conn = create_connection(self._db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]     
    #
    #
    def _create_table(self, conn) -> None:
        """ """
        table = "CREATE TABLE IF NOT EXISTS NodeBoundary(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            node_id INTEGER NOT NULL REFERENCES Node(name),\
                            x DECIMAL, y DECIMAL, z DECIMAL,\
                            rx DECIMAL, ry DECIMAL, rz DECIMAL, \
                            title TEXT);"
        #
        create_table(conn, table)
    #
    #
    def _push_boundary(self, conn, node_id: int, fixity: list):
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
        #
        project = (node_id, *fixity[:6], title, )
        sql = 'INSERT INTO NodeBoundary(node_id,\
                                         x, y, z, rx, ry, rz,\
                                         title)\
                                         VALUES(?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    #
    #
    @property
    def df(self):
        """Boundary df"""
        query = (self._component, )
        table = "SELECT Node.name, NodeBoundary.* \
                FROM Node, NodeBoundary \
                WHERE NodeBoundary.node_id = Node.number \
                AND Node.component_id = ? ;"
        #
        conn = create_connection(self._db_file)
        with conn:
            cur = conn.cursor()
            cur.execute (table, query)
            rows = cur.fetchall()
        #
        db = DBframework()
        header = ['name', 'number', 'node_id',
                  'ix', 'iy', 'iz', 'rx', 'ry', 'rz', 'title']
        boundf = db.DataFrame(data=rows, columns=header)
        header = ['name', 'ix', 'iy', 'iz', 'rx', 'ry', 'rz',  'title']        
        return boundf[header]
    
    @df.setter
    def df(self, df):
        """ """
        #
        columns = list(df.columns)
        header = {}
        for key in columns:
            if re.match(r"\b(id|name|node(s)?)\b", key, re.IGNORECASE):
                header[key] = 'name'
            
            elif re.match(r"\b(type)\b", key, re.IGNORECASE):
                header[key] = 'type'
            #
            # displacement
            elif re.match(r"\b((i)?(\_|\-|\s*)?x)\b", key, re.IGNORECASE):
                header[key] = 'x'
            
            elif re.match(r"\b((i)?(\_|\-|\s*)?y)\b", key, re.IGNORECASE):
                header[key] = 'y'               
            
            elif re.match(r"\b((i)?(\_|\-|\s*)?z)\b", key, re.IGNORECASE):
                header[key] = 'z'
            #
            # rotation
            elif re.match(r"\b((r)?(\_|\-|\s*)?x)\b", key, re.IGNORECASE):
                header[key] = 'rx'
            
            elif re.match(r"\b((r)?(\_|\-|\s*)?y)\b", key, re.IGNORECASE):
                header[key] = 'ry'               
            
            elif re.match(r"\b((r)?(\_|\-|\s*)?z)\b", key, re.IGNORECASE):
                header[key] = 'rz'
            #
            elif re.match(r"\b(title)\b", key, re.IGNORECASE):
                header[key] = 'title'
        #
        nodes = df[header.keys()].copy()
        nodes.rename(columns=header, inplace=True)
        #support.query("x != '' and y != '' and z != '' and rx != '' and ry != '' and rz != ''",
        #              inplace=True)        
        #
        for row in nodes.itertuples():
            #print(row)
            fixity=[row.x, row.y, row.z,
                    row.rx, row.ry, row.rz]
            if any(fixity):
                #print(fixity)
                self.__setitem__(node_name=row.name,
                                 fixity=fixity)
        #
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
    return record
#
#
#
class BoundarySQL:
    __slots__ = ['_nodes', '_db_file', '_component']
    
    def __init__(self, db_file:str,
                 component: int, 
                 db_system:str="sqlite") -> None:
        """
        """
        self._db_file = db_file
        self._component = component
        self._nodes = BoundaryNodeSQL(component=self._component,
                                      db_file=self._db_file)
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._component, )
        table = 'SELECT Node.name \
                FROM Node, NodeBoundary \
                WHERE NodeBoundary.node_id = Node.number \
                AND Node.component_id = ? ;'
        
        conn = create_connection(self._db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]    
    #
    #
    def __setitem__(self, node_name: int|str,
                    values: list|dict) -> None:
        """
        """
        try:
            self._labels.index(node_name)
            raise IOError(f' error boundary {node_name} already exist')
        except ValueError:        
            boundary_type = values[0]
            #b_number = self._set_item(b_name=node_name, 
            #                          b_type=boundary_type)
            #
            if re.match(r"\b(node(s)?|support(s)?|constrain(s)?)\b", boundary_type, re.IGNORECASE):
                if isinstance(values[1], str):
                    self._nodes[node_name] = values[1]
                else:
                    self._nodes[node_name] = values[1:]
            #elif 'curve' == boundary_type :
            #    raise Exception('--> Mat type No ready')
            else:
                raise IOError(' NodeBoundary type {:} not recognised'.format(boundary_type))


    def __getitem__(self, node_name: int|str):
        """
        """
        try:
            index = self._labels.index(node_name)
            boundary_type = self._type[index]
            #b_number = self._number[index]
        except ValueError:
            raise KeyError('Invalid boundary name : {:}'.format(node_name))
        #
        if re.match(r"\b(node(s)?|coord(inate)?)\b", boundary_type, re.IGNORECASE):
            return self._nodes[node_name]

        else:
            raise IOError(f' boundary type {boundary_type} not recognised')    
    #
    #
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
    #
    def nodes(self, values:list|None=None):
        """"""
        if isinstance(values, list):
            for value in values:
                node_name = value[0]
                boundary_type = "node"
                b_number = self._set_item(b_name=node_name, 
                                          b_type=boundary_type)
                #
                self._nodes[node_name] = value[1:]
        #
        return self._nodes     
    #
    def support(self, values:list|None=None):
        """"""
        return self.nodes(values)
    #
    #
    def get_number(self, start:int=1):
        """
        """
        try:
            n = max(self._number) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
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
        output += (f"Node {14*space} x {6*space} y {6*space} z {5*space} mx {5*space} my {5*space} mz {5*space} title")
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for key, node in self._nodes.items():
            #if sum(node[:6]) == 0:
            #    continue
            output += node.__str__()
        return output
    #
    #
    @property
    def df(self):
        """nodes in dataframe format"""
        #print('nodes df out')
        return self._nodes.df

    @df.setter
    def df(self, df):
        """nodes in dataframe format"""
        try:
            df.columns
            #
            #df[df['typr'].apply(lambda x: 'support'
            #                    #if re.search('^f', x)
            #                    if re.match(r"\b(node(s)?|support(s)?|constrain(s)?)\b", x, re.IGNORECASE)
            #                    else x)]
            #
            df['type'] = df['type'].apply(lambda x: 'support'
                                          if re.match(r"\b(node(s)?|support(s)?|constrain(s)?)\b", x, re.IGNORECASE)
                                          else x)
            #
            grptype = df.groupby('type')
            #
            self._nodes.df = grptype.get_group('support')
        except AttributeError:
            raise IOError('Node df not valid')      
#
