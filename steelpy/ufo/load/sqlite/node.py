#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
#from typing import NamedTuple
import re
#
# package imports
#
from steelpy.ufo.mesh.sqlite.nodes import pull_node, pull_node_number
from steelpy.ufo.load.process.node import (get_nodal_load,
                                           PointNode,
                                           DispNode,
                                           NodeLoadBasic,
                                           get_NodaLoad_df)

from steelpy.ufo.load.sqlite.utils import get_load_data, pull_basic
#
from steelpy.utils.sqlite.main import ClassBasicSQL
from steelpy.utils.sqlite.utils import create_connection, create_table
#
from steelpy.utils.dataframe.main import DBframework
#
#
# ---------------------------------
#
#
#
class NodeLoadGlobalSQL(NodeLoadBasic):
    __slots__ = ['_node', '_db_file', '_component'] 
    
    def __init__(self, component: int, db_file: str) -> None: # 
        """
        """
        super().__init__()
        self._component = component
        self._db_file = db_file
    #
    # -----------------------------------------------
    #
    @property
    def _labels(self):
        """ """
        1 / 0
        query = (self._component, )
        table = "SELECT Node.name \
                FROM Node, LoadNode, Load, LoadBasic \
                WHERE LoadNode.node_id = Node.number \
                AND LoadBasic.number = LoadNode.basic_id \
                AND LoadBasic.load_id = Load.number \
                AND Node.mesh_id = ? ;"
        #
        conn = create_connection(self._db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            rows = cur.fetchall()
        labels = [item[0] for item in rows]
        return labels    
    #
    # -----------------------------------------------
    #
    def __setitem__(self, node_name: int|str,
                    node_load: list|tuple|dict) -> None:
        """
        """
        1 / 0
    
    def __getitem__(self, node_name:int|str) :
        """
        """
        1 / 0
        # get load data
        conn = create_connection(self._db_file)
        with conn:  
            node_id = pull_node_number(conn, node_name,
                                       component=self._component)
        #
        try:
            node_id = node.number
            1 / 0
        except TypeError:
            raise IOError(f"Node {node_name} not found")
    #
    # ----------------------------------------------- 
    #
    @property
    def displacement(self):
        """ nodal displacement"""
        conn = create_connection(self.db_file)
        ndata = pull_NodeLoad_general(conn,
                                      load_name='*',
                                      node_name='*',
                                      component=self._component,
                                      load_type='displacement', 
                                      beam_flag = False)
        return ndata
    #
    # -----------------------------------------------   
    #
    #
    #
    # ----------------------------------------------- 
    #    
    #
    @property
    def df(self):
        """ """
        conn = create_connection(self._db_file)
        df = pull_NodeLoad_df(conn,
                              load_name='*',
                              node_name='*',
                              component=self._component)
                              #beam_flag = True)
        return df
        
    @df.setter
    def df(self, df):
        """ """
        1 / 0 # FIXME: basic load id
        conn = create_connection(self._db_file)
        with conn:
            query = (self._component, )
            table = "SELECT Node.name, Node.number FROM Node \
                     WHERE Node.mesh_id = ?;"
            cur = conn.cursor()
            cur.execute(table, query)
            nodes = cur.fetchall()
            nodes = {item[0]:item[1] for item in nodes}
            #
            #
            query = (self._component, )
            table = "SELECT Element.name, Element.number FROM Element \
                     WHERE Element.mesh_id = ?;"
            cur = conn.cursor()
            cur.execute(table, query)
            elements = cur.fetchall()
            elements = {item[0]:item[1] for item in elements}            
            #
            #
            query = ('basic', self._component, )
            table = "SELECT Load.name, Load.number FROM Load \
                         WHERE Load.level = ? \
                         AND Load.mesh_id = ?;"
            cur = conn.cursor()
            cur.execute(table)
            basic = cur.fetchall()
            basic = {item[0]:item[1] for item in basic}        
        #
        df['load_id'] = df['name'].apply(lambda x: basic[x])
        df['element_id'] = None # df['name'].apply(lambda x: elements[x])
        df['node_id'] = df['node'].apply(lambda x: nodes[x])
        df['system'] = 'global'
        df['type'] = df['type'].apply(lambda x: x.lower())
        #
        #
        header = ['load_id', 'element_id', 'node_id',
                  'title', 'system', 'type',
                  'fx', 'fy', 'fz', 'mx', 'my', 'mz',
                  'x', 'y', 'z', 'rx', 'ry', 'rz']
        #
        nodeconn = df[header].copy()
        nodeconn.replace(to_replace=[''], value=[float(0)], inplace=True)
        nodeconn['element_id'] = df['element_id']
        nodeconn['title'] = df['title']
        #
        with conn:
            nodeconn.to_sql('LoadNode', conn,
                            index_label=header, 
                            if_exists='append', index=False)
        #print('node load in')
        #1 / 0     
#
#
# ---------------------------------
#
class NodeLoadItemSQL(ClassBasicSQL):
    __slots__ = ['_node', '_db_file', '_name', '_component']
    
    def __init__(self, load_name: str|int, db_file: str,
                 component: int) -> None: # 
        """
        """
        #self._db_file = db_file
        self._name = load_name
        self._component = component       
        # create node table
        #conn = create_connection(self._db_file)
        #with conn:        
        #    nodeload_table(conn)
        super().__init__(db_file)
    #
    # -----------------------------------------------
    #
    @property
    def _labels(self):
        """ """
        query = (self._name, self._component, )
        table = "SELECT Node.name \
                FROM Node, LoadNode, Load, LoadBasic \
                WHERE LoadNode.node_id = Node.number \
                AND LoadBasic.number = LoadNode.basic_id \
                AND LoadBasic.load_id = Load.number \
                AND Load.name = ? \
                AND Load.mesh_id = ? ;"
        #
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            rows = cur.fetchall()
        labels = set([item[0] for item in rows])
        return list(labels)       
    #
    # -----------------------------------------------
    #
    def __setitem__(self, node_name: int|str,
                    node_load: list|tuple|dict) -> None:
        """
        """
        point_load = get_nodal_load(node_load)
        #
        conn = create_connection(self.db_file)
        with conn:
            push_node_load(conn,
                           load_name=self._name,
                           node_name=node_name,  
                           node_load=point_load,
                           component=self._component)
        #
    #
    def __getitem__(self, node_name:int|str) :
        """
        """
        # get load data
        conn = create_connection(self.db_file)
        with conn:  
            node = pull_node(conn, node_name,
                             component=self._component)
        #
        try:
            node_id = node.number
            return NodeItemSQL(load_name=self._name,
                               node=node,
                               component=self._component, 
                               db_file=self.db_file)
        except TypeError:
            raise IOError(f"Node {node_name} not valid")
    #
    #
    # -----------------------------------------------
    #
    def __str__(self, units: str = "si") -> str:
        """ """
        conn = create_connection(self.db_file)
        with conn:         
            nodeload = pull_NodeLoad_item(conn,
                                     load_name=self._name,
                                     node_name='*')
        #
        #nodeload = [item if item[]]
        #
        output = ""
        for item in nodeload:
            output += item.__str__()
            # print('---')
        return output
    #
    #
    #@property
    #def load(self):
    #    """
    #    """
    #    return self._node._load    
    #
    # -----------------------------------------------
    #
    def get_group_nodes(self):
        """
        """
        items = []
        set_items = set(self._labels)
        index = self._cls._index
        load_tile = self._cls._title[index]
        load_name = self._cls._labels[index]
        try:
            self._load_id.index(load_name)
            bd_file = self._cls.bd_file
            conn = create_connection(bd_file)
            load_id = pull_basic_load_id(conn, load_name)
            for node_id in set_items:
                with conn:
                    sum_load = self.get_sum_column(conn, node_id, load_id)
                    items.append(PointNode(*sum_load, node_id,
                                            load_tile, 0, 0))
            return items
        except ValueError:
            return items
    #
    def get_sum_column(self, conn, node_name, load_id):
        """ """
        1 / 0
        project = (node_name, load_id,)
        sql = 'SELECT SUM(fx), SUM(fy), SUM(fz), SUM(mx), SUM(my), SUM(mz)\
               FROM LoadNode WHERE node_id = ? AND load_id = ?'
        cur = conn.cursor()
        cur.execute( sql, project )
        record = cur.fetchone()
        return record
    #
    def update(self, other) -> None:
        """
        """
        pl = other._point
        try:
            index = pl._index
            #node_name = pl._labels[index]
            for x, node_name in enumerate(pl._labels):
                point_load = [pl._fx[x], pl._fy[x], pl._fz[x],
                              pl._mx[x], pl._my[x], pl._mz[x],
                              pl._title[x]]
                self.__setitem__(node_name, point_load)
        except AttributeError:
            pass        
        #print('----')
    #
    # -----------------------------------------------
    #
    #
    def _new_table(self, conn):
        """ """
        table = "CREATE TABLE IF NOT EXISTS LoadNode(\
                number INTEGER PRIMARY KEY NOT NULL,\
                basic_id INTEGER NOT NULL REFERENCES Load(number),\
                node_id INTEGER NOT NULL REFERENCES Node(number),\
                title TEXT,\
                system TEXT NOT NULL,\
                type TEXT NOT NULL,\
                fx DECIMAL,\
                fy DECIMAL,\
                fz DECIMAL,\
                mx DECIMAL,\
                my DECIMAL,\
                mz DECIMAL,\
                x DECIMAL,\
                y DECIMAL,\
                z DECIMAL,\
                rx DECIMAL,\
                ry DECIMAL,\
                rz DECIMAL,\
                psi DECIMAL,\
                B DECIMAL,\
                Tw DECIMAL);"
        #
        create_table(conn, table)
    #
    # -----------------------------------------------
    # 
    #
    # -----------------------------------------------   
    #
    @property
    def df(self):
        """ """
        #print('node load out')
        conn = create_connection(self.db_file)
        #ndata = get_SQLdata(conn, load_name=self._name,
        #                    node_name='*',
        #                    beam_flag = False)
        #
        #
        df = pull_NodeLoad_df(conn,
                              load_name=self._name,
                              node_name='*',
                              component=self._component)
        return df
        
    @df.setter
    def df(self, df):
        """ """
        dfdict = get_NodaLoad_df(df)
        #
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            table = "SELECT Node.name, Node.number FROM Node \
                    WHERE Node.mesh_id = ? ;"
            query = (self._component, )
            cur.execute(table, query)
            nodes = cur.fetchall()
            nodes = {item[0]:item[1] for item in nodes}
            #
            table = "SELECT Load.name, Load.number FROM Load \
                     WHERE Load.level = ? \
                     AND Load.mesh_id = ? ;"
            query = ('basic', self._component, )
            cur.execute(table, query)
            basic = cur.fetchall()
            basic = {item[0]:item[1] for item in basic}        
        #
        header = ['basic_id', 'node_id',
                  'title', 'system', 'type',
                  'fx', 'fy', 'fz', 'mx', 'my', 'mz',
                  'x', 'y', 'z', 'rx', 'ry', 'rz',
                  'psi', 'B', 'Tw','step']        
        #
        for key, item in  dfdict.items():
            nodeconn = item.copy()
            nodeconn['basic_id'] = nodeconn['name'].apply(lambda x: basic[x])
            nodeconn['node_id'] = nodeconn['node'].apply(lambda x: nodes[x])
            nodeconn.replace(to_replace=[''], value=[float(0)], inplace=True)
            nodeconn[['psi', 'B', 'Tw']] = None
            #
            if key in ['force']:
                nodeconn[['x', 'y', 'z', 'rx', 'ry', 'rz']] = None
            elif key in ['displacement', 'mass']:
                nodeconn[['fx', 'fy', 'fz', 'mx', 'my', 'mz']] = None
            else:
                raise IOError(f'node load {key} not valid')
            #
            with conn:
                nodeconn[header].to_sql('LoadNode', conn,
                                        index_label=header, 
                                        if_exists='append', index=False)
        #
        #1 / 0     
#
#
@dataclass
class NodeItemSQL:
    __slots__ = ['_title', '_complex', '_component', 
                 '_system_flag', '_system', 
                 '_db_file', '_name', '_node']
    
    def __init__(self, load_name: str|int,
                 node, component: int, 
                 db_file: str) -> None:  
        """
        """
        #super().__init__()
        self._name = load_name
        self._node = node
        self._db_file = db_file
        self._component = component
        #
        # create node table
        #conn = create_connection(self._db_file)
        #with conn:        
        #    nodeload_table(conn)
    #
    # -----------------------------------------------
    #
    @property
    def _labels(self):
        """ """
        1 / 0
        query = (self._component, )
        table = "SELECT Node.name \
                FROM Node, LoadNode, Load, LoadBasic \
                WHERE LoadNode.node_id = Node.number \
                AND LoadBasic.number = LoadNode.basic_id \
                AND LoadBasic.load_id = Load.number \
                AND Node.mesh_id = ? ;"
        #
        conn = create_connection(self._db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            rows = cur.fetchall()
        labels = [item[0] for item in rows]
        return labels
    #
    # -----------------------------------------------
    #
    def __setitem__(self, node_name:int|str,
                    point_load: list) -> None:
        """
        """
        1 / 0
        # get load data
        conn = create_connection(self._db_file)
        #with conn:  
        #    node = check_nodes(conn, node_name)
        #
        #try:
        #    node_id = node[0]           
        #except TypeError:
        #    raise IOError(f"Node {node_name} not found")
        #
        self._labels.append(node_name)
        load_name = point_load.pop(0)
        #load_type = point_load.pop(0)
        #
        #if isinstance(point_load, dict):
        #    point_load = get_nodal_load(point_load)
        #    title = point_load[-1]
        #    point_load.pop()
        #elif isinstance(point_load[-1], str):
        #    title = point_load.pop()
        #    point_load = get_nodal_load(point_load)
        #else:
        #    title = "NULL"
        #    point_load = get_nodal_load(point_load)
        #
        # Push to database
        #
        #if re.match(r"\b(point|load|node)\b", self._type, re.IGNORECASE):
        with conn:
            push_node_load(conn=conn,
                           load_name=load_name,
                           node_name=node_name, 
                           node_load=point_load)
        
        #elif re.match(r"\b(mass)\b", self._type, re.IGNORECASE):
        #    raise NotImplementedError(f'node mass')
        #
        #elif re.match(r"\b(disp(lacement)?)\b", self._type, re.IGNORECASE):
        #    raise NotImplementedError(f'node displacement')
        #
        #else:
        #    raise IOError(f'node load type {self._type} not recognized')    
    #
    def __getitem__(self, node_name:int|str) -> list:
        """
        """
        1 / 0
        # get load data
        conn = create_connection(self._db_file)
        #
        # Get database
        #
        #if re.match(r"\b(point|load|node)\b", self._type, re.IGNORECASE):
        with conn:
            nload = self._get_load(conn, node_name)
        
        #elif re.match(r"\b(mass)\b", self._type, re.IGNORECASE):
        #    raise NotImplementedError(f'node mass')
        #
        #elif re.match(r"\b(disp(lacement)?)\b", self._type, re.IGNORECASE):
        #    raise NotImplementedError(f'node displacement')
        #
        #else:
        #    raise IOError(f'node load type {self._type} not recognized')         
        #
        #
        #_points: list = []
        #for _index in _index_list:
        #    _points.append(PointNode(self._fx[_index], self._fy[_index], self._fz[_index],
        #                             self._mx[_index], self._my[_index], self._mz[_index],
        #                             self._labels[_index], self._title[_index],
        #                             self._system[_index], self._complex[_index], self._type))
        #
        return nload
    #
    # -----------------------------------------------
    #
    @property
    def load(self):
        """ """
        1 / 0
    
    @load.setter
    def load(self, node_load:list|tuple|dict):
        """ """
        # update load with header
        if isinstance(node_load, PointNode):
            node_load = ['force',
                         # Fx, Fy, Fz, Mx, My, Mz,
                         *node_load[:6],
                         node_load.title]
        elif isinstance(node_load, dict):
            node_load.update({'type': 'force',})
        
        elif isinstance(node_load[0], list):
            node_load = [['force', *item] for item in node_load]
        
        else:
            node_load.insert(0, 'force')
        #
        #
        conn = create_connection(self._db_file)
        with conn:
            push_node_load(conn,
                           load_name=self._name,
                           node_name=self._node.name,  
                           node_load=node_load,
                           component=self._component)
    #
    #
    @property
    def displacement(self):
        """ nodal displacement"""
        1 / 0
    
    @displacement.setter
    def displacement(self, node_load: list|tuple|dict):
        """ nodal displacement"""
        # update load with header
        if isinstance(node_load, dict):
            node_load.update({'type': 'displacement',})
        
        elif isinstance(node_load[0], list):
            node_load = [['displacement', *item] for item in node_load]
        
        else:
            node_load.insert(0, 'displacement')
        #
        #
        conn = create_connection(self._db_file)
        with conn:
            push_node_load(conn,
                           load_name=self._name,
                           node_name=self._node.name,  
                           node_load=node_load,
                           component=self._component)
    #
    @property
    def mass(self):
        """nodal mass"""
        1 / 0
    
    @mass.setter
    def mass(self, node_load):
        """ """
        # update load with header
        if isinstance(node_load, dict):
            node_load.update({'type': 'mass',})
        
        elif isinstance(node_load[0], list):
            node_load = [['mass', *item] for item in node_load]
        
        else:
            node_load.insert(0, 'mass')
        #
        #
        conn = create_connection(self._db_file)
        with conn:
            push_node_load(conn,
                           load_name=self._name,
                           node_name=self._node.name,  
                           node_load=node_load,
                           component=self._component)        
    #       
    #
    def _get_load(self, conn, node_name:int|str):
        """ """
        1 / 0
        load_data = get_load_data(conn, self._name, load_type='basic')
        load_id = load_data[0]
        #
        node = check_nodes(conn, node_name)
        #try:
        #    node_id = node[0]
        #except TypeError:
        #    raise IOError(f"Node {node_name} not found")
        # Node load
        cur = conn.cursor()
        cur.execute("SELECT * \
                     FROM LoadNode, LoadBasic, Load \
                     WHERE load_id = {:} \
                     AND node_id = {:}; "
                    .format(load_id, node_id))
        rows = cur.fetchall()
        node_load = []
        for row in rows:
            if row[6] in ['load', 'point', 'node']:
                data = [*row[7:13],                      # [fx, fy, fz, mx, my, mz]
                        node_name, row[4], self._name,   # [name, title, load_name]
                        0, 0, 'load']                    # [system, load_complex, load_type]
                #
                # [fx, fy, fz, mx, my, mz, name, title, load_name, system, load_complex, load_type]
                node_load.append(PointNode._make(data))
            else:
                1 / 0
        return node_load
    #
    #
    # -----------------------------------------
    #
#    
#
#
# ---------------------------------
# Operations
# ---------------------------------
#
#
def push_node_load(conn, load_name: str|int,
                   node_name: int|str,  
                   node_load: list|tuple|dict,
                   component: int,
                   system: str = "global"):
    """
    node_load = [type, Fx, Fy, Fz, Mx, My, Mz, title]
    """
    # clean nodal load input
    #point_load = get_nodal_load(node_load)
    load_type = node_load.pop(0)
    load_title = node_load.pop(-1)
    #
    node_id, load_id = pull_load_specs(conn,
                                      node_name,
                                      load_name,
                                      component=component)
    #
    # TODO: load type
    basic = pull_basic(conn, load_name=load_name)
    try:
        basic_id = basic[0]
    except KeyError:
        raise IOError(f'load {load_name} not found')
    #
    # Load type check
    #
    if re.match(r"\b(point|load|node|force)\b", load_type, re.IGNORECASE):
        project = (basic_id, node_id,                   # basic_id,  node_id
                   load_title, system, 'force',         # title, system, type, 
                   *node_load,                          # fx, fy, fz, mx, my, mz,
                   None, None, None, None, None, None)  # x, y, z, rx, ry, rz
    # TODO: define mass
    elif re.match(r"\b(mass)\b", load_type, re.IGNORECASE):
        project = (basic_id, node_id,                    # basic_id,  node_id
                   load_title, system, 'mass',           # title, system, type, 
                   None, None, None, None, None, None,   # fx, fy, fz, mx, my, mz,
                   *node_load)                           # x, y, z, rx, ry, rz

    elif re.match(r"\b(disp(lacement)?)\b", load_type, re.IGNORECASE):
        project = (basic_id, node_id,                   # basic_id,  node_id
                   load_title, system, 'displacement',  # title, system, type, 
                   None, None, None, None, None, None,  # fx, fy, fz, mx, my, mz,
                   *node_load)                          # x, y, z, rx, ry, rz

    else:
        raise IOError(f'node load type {load_type} not recognized')
    #
    # print('-->')
    sql = 'INSERT INTO LoadNode(basic_id, node_id, \
                                title, system, type, \
                                fx, fy, fz, mx, my, mz,\
                                x, y, z, rx, ry, rz)\
                                VALUES(?,?,?,?,?,?,?,?,?, \
                                       ?,?,?,?,?,?,?,?)'
    cur = conn.cursor()
    cur.execute(sql, project)
    #print('-->')
#
#
# ---------------------------------
# General sql operations
# ---------------------------------
#
#
def pull_load_specs(conn, node_name: str|int,
                   load_name: str|int, component: int):
    """ """
    # Node check
    node_id = pull_node_number(conn, node_name,
                               component=component)
    #
    if not node_id:
        raise IOError(f"Node {node_name} not found")    
    #
    # Load check
    load_data = get_load_data(conn, load_name,
                              load_level='basic',
                              component=component)
    try:
        load_id = load_data[0]
    except TypeError:
        raise IOError(f"Load {load_name} not found")
    #
    #
    return node_id, load_id
#
#
# ---------------------------------
#
def pull_NodeLoad_general(conn, load_name: int|str,
                          node_name: int|str,
                          component: int,
                          load_type: str, 
                          beam_flag: bool = False):
    """ """
    ndata = pull_NodeLoad_SQLdata(conn,
                                  load_name=load_name,
                                  node_name=node_name,
                                  component=component,
                                  load_type=load_type, 
                                  beam_flag = beam_flag)
    #
    dipnode = []
    1 / 0
    if load_type in ['displacement']:
        for row in ndata:
            data = [*row[18:24],                 # [x, y, z, rx, ry, rz]
                    row[3], row[2], row[0],      # [name, title, load_name]
                    row[10], 0, row[11]]         # [system, load_complex, load_type]
            dipnode.append(DispNode._make(data))
    
    elif  load_type in ['load']:
        for row in ndata:
            data = [*row[12:18],                 # [fx, fy, fz, mx, my, mz]
                    row[3], row[2], row[0],      # [name, title, load_name]
                    row[10], 0, row[11]]         # [system, load_complex, load_type]
            dipnode.append(PointNode._make(data))
    
    elif  load_type in ['mass']:
        1 / 0
        
    else:
        raise IOError(f'load type {load_type} not valid')
    #
    return dipnode
#
#
def pull_NodeLoad_SQLdata(conn, load_name: int|str,
                          node_name: int|str,
                          component: int,
                          load_type: str|None = None):
                          #beam_flag: bool = False):
    """ """
    table = ""
    items = []
    #
    # load
    if load_name in ['*', None, '']:
        load_name = None
    else:
        items.extend([load_name]) #(load_name, )
        table += "AND Load.name = ? "
    #
    # node
    if node_name in ['*', None]:
        pass
    else:
        items.extend([node_name])
        table += "AND Node.name = ? "    
    #
    # load type
    if load_type in ['*', None]:
        pass
    else:
        items.extend([load_type])
        table += "AND LoadNode.type = ? "     
    #
    table += "AND Mesh.number = ? ;"
    items.extend([component])
    #
    # Node load
    nodal_load = pull_NodeLoad(conn, items, table)
    #
    return nodal_load
#
def pull_NodeLoad(conn, items: list, utils: str):
    """ """
    # TODO: check it
    table = "SELECT Load.name AS load_name, \
                    Mesh.name AS mesh_name, \
                    Load.title AS load_title, \
                    Node.name AS node_name, \
                    Node.idx as node_index, \
                    LoadNode.* \
            FROM Load, Node, LoadNode, LoadBasic, Mesh \
            WHERE LoadNode.basic_id = LoadBasic.number \
            AND LoadBasic.load_id = Load.number \
            AND LoadNode.node_id = Node.number \
            AND Load.mesh_id = Mesh.number  "
    #
    table += utils
    #
    # Node load
    with conn:
        cur = conn.cursor()
        if items:
            cur.execute(table, tuple(items))
        else:
            cur.execute(table)
        rows = cur.fetchall()
    #
    return rows
#
def pull_BeamNodeLoad(conn, items: list, utils: str):
    """ """
    1 / 0
    table = "SELECT Load.name AS load_name, \
                    Mesh.name AS mesh_name, \
                    Load.title AS load_title, \
                    Node.name AS node_name, \
                    Node.mesh_idx as node_index, \
                    Element.name as element_name, \
                    LoadNode.* \
            FROM Load, Node, Element, LoadNode, LoadBasic, Mesh \
            WHERE LoadNode.basic_id = LoadBasic.number \
            AND LoadBasic.load_id = Load.number \
            AND LoadNode.node_id = Node.number \
            AND LoadNode.element_id is not NULL \
            AND LoadNode.element_id = Element.number "
    #
    table += utils
    #
    # beam load
    with conn:
        cur = conn.cursor()
        if items:
            cur.execute(table, tuple(items))
        else:
            cur.execute(table)
        rows = cur.fetchall()
    #
    return rows
#
# ---------------------------------
#
def pull_NodeLoad_item(conn, load_name: int|str,
                       node_name: int|str):
    """ """
    table = "SELECT Node.name, LoadNode.* \
            FROM Load, Node, LoadNode, LoadBasic \
            WHERE LoadNode.basic_id = LoadBasic.number \
            AND LoadBasic.load_id = Load.number \
            AND LoadNode.node_id = Node.number "
    #
    items = []
    # load
    if load_name in ['*', None, '']:
        load_name = None
    else:
        items.extend([load_name]) #(load_name, )
        table += "AND Load.name = ?"
    #
    # node
    if node_name in ['*', None]:
        pass
    else:
        items.extend([node_name])
        table += "AND Node.name = ? "    
    #
    # beam
    #if beam_flag:
    #    table += "AND LoadNode.element_id is NULL"
    #
    # Node load
    with conn:
        cur = conn.cursor()
        cur.execute(table, tuple(items))
        rows = cur.fetchall()
    #
    #1 / 0
    node_load = []
    for row in rows:
        if row[6] in ['mass']:
            1 / 0
        
        elif row[6] in ['displacement']:
            1 / 0
            data = [*row[13:19],                 # [x, y, z, rx, ry, rz]
                    row[0], row[4], load_name,   # [name, title, load_name]
                    0, 0, 'displacement']         # [system, load_complex, load_type]
            node_load.append(DispNode._make(data))
        else:
            data = [*row[7:13],                      # [fx, fy, fz, mx, my, mz]
                    row[0], row[4], load_name,       # [name, title, load_name]
                    0, 0, 'force']                    # [system, load_complex, load_type]
            #
            # [fx, fy, fz, mx, my, mz, name, title, load_name, system, load_complex, load_type]
            node_load.append(PointNode._make(data))
    #
    return node_load         
#
#
def pull_NodeLoad_df(conn, load_name: int|str,
                     node_name: int|str,
                     component: int):
                     #beam_flag: bool):
    """nodes in dataframe format"""
    ndata = pull_NodeLoad_SQLdata(conn,
                                  load_name=load_name,
                                  node_name=node_name,
                                  component=component)
                                  #beam_flag = beam_flag)
    
    #
    cols = ['load_name', 'mesh_name', 
            'load_title',
            'node_name', 'node_index', 
            #'element_name',
            'number', 
            'load_id',
            #'element_id',
            'node_id',
            'load_comment', 'load_system','load_type',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
            'x', 'y', 'z', 'rx', 'ry', 'rz',
            'psi', 'B', 'Tw', 'step']
    #
    # dataframe
    db = DBframework()
    df = db.DataFrame(data=ndata, columns=cols)
    df['load_level'] = 'basic'
    #df['component'] = component
    #
    df = df[['load_name', 'mesh_name', 
             'load_title', 'load_level',
             'load_id', 'load_system', 'load_comment',
             #'element_name',
             'node_name', 'node_index', 
             'load_type',
             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
             'x', 'y', 'z', 'rx', 'ry', 'rz',
             'psi', 'B', 'Tw', 'step']]
    return df