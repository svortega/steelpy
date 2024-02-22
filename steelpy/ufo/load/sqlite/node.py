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
from steelpy.ufo.mesh.sqlite.nodes import get_node
from steelpy.ufo.load.process.nodes import get_nodal_load, PointNode
from steelpy.ufo.load.sqlite.utils import  get_load_data
from steelpy.utils.sqlite.main import ClassBasicSQL
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework
#
#
# ---------------------------------
#
class NodeMainSQL(Mapping):
    
    # 
    #
    # -----------------------------------------------
    #
    def __contains__(self, value) -> bool:
        return value in self._labels

    def __len__(self) -> int:
        return len(self._labels)

    def __iter__(self):
        """
        """
        items = set(self._labels)
        return iter(items)
    #
#
#
# ---------------------------------
#
class NodeLoadItemSQL(ClassBasicSQL):
    __slots__ = ['_node', '_db_file', '_name', '_component'] # '_type', '_labels', 
    
    def __init__(self, load_name: str|int, db_file: str,
                 component: int) -> None: # 
        """
        """
        #self._db_file = db_file
        self._name = load_name
        self._component = component
        #self._node = NodeItemSQL(self._name, self._db_file)
        #self._node =  NodeLoadSQL(self._db_file)
        #
        #self._labels = []
        #self._type = []        
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
        # FIXME : component
        table = "SELECT Node.name \
                FROM Node, LoadNode, Load \
                WHERE LoadNode.node_id = Node.number \
                AND Load.number = LoadNode.load_id \
                AND Load.name = ? \
                AND Load.component_id = ? ;"
        #
        #if isinstance(self._name, str):
        #    table += f"AND Load.name = '{self._name}' "
        #else:
        #    table += f"AND Load.name = {self._name}"
        #
        conn = create_connection(self.db_file)
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
        conn = create_connection(self.db_file)
        with conn:
            push_node_load(conn,
                           load_name=self._name,
                           node_name=node_name,  
                           node_load=node_load,
                           component=self._component)
        #
    #
    def __getitem__(self, node_name:int|str) :
        """
        """
        # get load data
        conn = create_connection(self.db_file)
        with conn:  
            node = get_node(conn, node_name,
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
            nodeload = get_node_load(conn,
                                     load_name=self._name,
                                     node_name='*', beam_flag=True)
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
            load_id = get_basic_load_id(conn, load_name)
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
    def _create_table(self, conn):
        """ """
        table = "CREATE TABLE IF NOT EXISTS LoadNode(\
                number INTEGER PRIMARY KEY NOT NULL,\
                load_id INTEGER NOT NULL REFERENCES Load(number),\
                element_id INTEGER REFERENCES Element(number),\
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
    #
    # -----------------------------------------------   
    #
    @property
    def df(self):
        """ """
        #print('node load out')
        conn = create_connection(self._db_file)
        #ndata = get_SQLdata(conn, load_name=self._name,
        #                    node_name='*',
        #                    beam_flag = False)
        #
        #
        df = get_nodedf(conn, load_name=self._name,
                        node_name='*',
                        beam_flag = True)
        return df
        
    @df.setter
    def df(self, df):
        """ """
        conn = create_connection(self._db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute("SELECT Node.name, Node.number FROM Node;")
            nodes = cur.fetchall()
            nodes = {item[0]:item[1] for item in nodes}
            #
            cur.execute("SELECT Element.name, Element.number FROM Element;")
            elements = cur.fetchall()
            elements = {item[0]:item[1] for item in elements}            
            #
            cur.execute("SELECT Load.name, Load.number FROM Load \
                         WHERE Load.level = 'basic';")
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
        query = (self._component, )
        table = "SELECT Node.name \
                FROM Node, LoadNode, Load \
                WHERE LoadNode.node_id = Node.number \
                AND Load.number = LoadNode.load_id \
                AND Node.component_id = ? ;"
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
    def load(self, node_load):
        """ """
        conn = create_connection(self._db_file)
        #
        if isinstance(node_load[-1], str):
            load_title = node_load.pop()
        else:
            load_title = None
        #
        # Load check
        with conn:
            node_id, load_id = get_load_specs(conn,
                                              node_name=self._node.name,
                                              load_name=self._name,
                                              component=self._component)      
        #
        #
        point_load = get_nodal_load(node_load)
        
        project = (load_id, None, node_id,
                   load_title, 'global', 'load',
                   *point_load,
                   None, None, None, None, None, None)
        #
        sql = 'INSERT INTO LoadNode(load_id, element_id, node_id, \
                                        title, system, type, \
                                        fx, fy, fz, mx, my, mz,\
                                        x, y, z, rx, ry, rz)\
                                    VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
        #
        with conn:         
            cur = conn.cursor()
            cur.execute(sql, project)        
        #1 / 0
    #
    def _get_loadX(self, conn, node_name:int|str):
        """ """
        load_data = get_load_data(conn, self._name, load_type='basic')
        load_id = load_data[0]
        #
        node = check_nodes(conn, node_name)
        try:
            node_id = node[0]
        except TypeError:
            raise IOError(f"Node {node_name} not found")
        # Node load
        cur = conn.cursor()
        cur.execute("SELECT * FROM LoadNode \
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
    #def __call__(self, node_id: str|int, load_id: str|int):
    #    self._node_id = node_id
    #    self._name = load_id
    #    return self    
    #
    #
    #@property
    #def load(self):
    #    """
    #    """
    #    1 / 0
    #    try:
    #        point_id = self._node_id
    #        return self._load[point_id]
    #    except :
    #        raise IndexError
    #
    #@load.setter
    #def load(self, node_load: list):
    #    """
    #    Point Load
    #    """
    #    node_name = self._node_id
    #    #if isinstance(values, dict):
    #    #    self._load[node_name] = values
    #    #elif isinstance(values[0], list):
    #    #    for value in values:
    #    #        self._load[node_name] = value
    #    #else:
    #    #    self._load[node_name] = values
    #    node_load.insert(0, 'load')
    #    self.__setitem__(node_name, node_load)
    #
#    
#
# ---------------------------------
#
class NodeLoadGlobalSQL(NodeMainSQL):
    __slots__ = ['_node', '_db_file', '_component'] 
    
    def __init__(self, component: int, db_file: str) -> None: # 
        """
        """
        self._component = component
        self._db_file = db_file
    #
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
            node = get_node(conn, node_name)
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
    def df(self):
        """ """
        conn = create_connection(self._db_file)
        df = get_nodedf(conn, load_name='*',
                        node_name='*',
                        component=self._component, 
                        beam_flag = True)
        return df
        
    @df.setter
    def df(self, df):
        """ """
        1 / 0
        conn = create_connection(self._db_file)
        with conn:
            query = (self._component, )
            table = "SELECT Node.name, Node.number FROM Node \
                     WHERE Node.component_id = ?;"
            cur = conn.cursor()
            cur.execute(table, query)
            nodes = cur.fetchall()
            nodes = {item[0]:item[1] for item in nodes}
            #
            #
            query = (self._component, )
            table = "SELECT Element.name, Element.number FROM Element \
                     WHERE Element.component_id = ?;"
            cur = conn.cursor()
            cur.execute(table, query)
            elements = cur.fetchall()
            elements = {item[0]:item[1] for item in elements}            
            #
            #
            query = ('basic', self._component, )
            table = "SELECT Load.name, Load.number FROM Load \
                         WHERE Load.level = ? \
                         AND Load.component_id = ?;"
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
# Operations
# ---------------------------------
#
def get_load_node(node_load):
    """ """
    # Load preparation
    if isinstance(node_load, dict):
        1 / 0
        point_load = get_nodal_load(node_load)
        load_title = point_load[-1]
        point_load.pop()
    else:
        if isinstance(node_load[-1], str):
            load_title = node_load.pop()
        else:
            load_title = None
        load_type = node_load.pop(0)
        point_load = get_nodal_load(node_load)
    #
    return load_type, point_load, load_title
#
#
def get_load_specs(conn, node_name: str|int,
                   load_name: str|int, component: int):
    """ """
    # Node check
    node = get_node(conn, node_name,
                    component=component)
    try:
        node_id = node.number
    except TypeError:
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
    return node_id, load_id
#
# ---------------------------------
# General sql operations
# ---------------------------------
#
#
def get_basic_load_id(conn, load_name:int|str):
    """ """
    cur = conn.cursor()
    cur.execute("SELECT * FROM Load\
                 WHERE name = {:} \
                 AND type = 'basic'".format(load_name))
    loads = cur.fetchone()
    return loads[0]
#
#
#
#class NodeItemSQL(NodeItem):
#    __slots__ = ['_load', '_displacement', '_mass', '_node_id',
#                 '_bd_file']
#
#    def __init__(self, load_name: str, bd_file: str):
#        """
#        """
#        super().__init__()
#        self._bd_file =  bd_file
#        self._load = NodeLoadSQL(load_name, "load", self._bd_file)
#        self._displacement = NodeLoadSQL(load_name, "displacement", self._bd_file)
#        self._mass = NodeLoadSQL(load_name, "mass", self._bd_file)
# 
#
#
#
def push_node_load(conn, load_name: str|int,
                   node_name: int|str,  
                   node_load: list|tuple|dict,
                   component: int,
                   system: str = "global"):
    """
    """
    #
    load_type, point_load, load_title = get_load_node(node_load)
    #
    node_id, load_id = get_load_specs(conn, node_name, load_name,
                                              component=component)
    # Load type check
    #
    if re.match(r"\b(point|load|node)\b", load_type, re.IGNORECASE):
        project = (load_id, None, node_id,
                   load_title, system, 'load',
                   *point_load,
                   None, None, None, None, None, None)

    elif re.match(r"\b(mass)\b", load_type, re.IGNORECASE):
        project = (load_id, None, node_id,
                   load_title, system, 'mass',
                   None, None, None,
                   *point_load,
                   None, None, None, None, None, None)

    elif re.match(r"\b(disp(lacement)?)\b", load_type, re.IGNORECASE):
        project = (load_id, None, 'displacement',
                   load_title, system, load_type,
                   None, None, None, None, None, None,
                   *point_load)

    else:
        raise IOError(f'node load type {load_type} not recognized')
    #
    # print('-->')
    #
    sql = 'INSERT INTO LoadNode(load_id, element_id, node_id, \
                                    title, system, type, \
                                    fx, fy, fz, mx, my, mz,\
                                    x, y, z, rx, ry, rz)\
                                VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
    cur = conn.cursor()
    cur.execute(sql, project)
    #print('-->')
#
#
# ---------------------------------
#
def get_SQLdata(conn, load_name: int|str,
                node_name: int|str,
                component: int, 
                beam_flag: bool = False):
    """ """
    #
    table = ""
    items = []
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
    table += "AND Component.number = ? ;"
    items.extend([component])
    #
    # Node load
    #
    nodal_load = get_Node(conn, items, table)
    #nodal_load.extend(get_Node(conn, items, table))
    #
    if beam_flag:
        nodal_load.extend(get_Beam(conn, items, table))
    #
    return nodal_load
#
def get_Node(conn, items: list, utils: str):
    """ """
    table = "SELECT Load.name AS load_name, \
                    Component.name AS component_name, \
                    Load.title AS load_title, \
                    Node.name AS node_name, \
                    LoadNode.element_id as element_name, \
                    LoadNode.* \
            FROM Load, Node, LoadNode, Component \
            WHERE LoadNode.load_id = Load.number \
            AND LoadNode.node_id = Node.number \
            AND LoadNode.element_id is NULL "
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
def get_Beam(conn, items: list, utils: str):
    """ """
    table = "SELECT Load.name AS load_name, \
                    Component.name AS component_name, \
                    Load.title AS load_title, \
                    Node.name AS node_name, \
                    Element.name as element_name, \
                    LoadNode.* \
            FROM Load, Node, Element, LoadNode, Component \
            WHERE LoadNode.load_id = Load.number \
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
#
#
def get_node_load(conn, load_name: int|str,
                  node_name: int|str,
                  beam_flag: bool = False):
    """ """
    #conn = create_connection(self._db_file)
    #load_data = get_load_data(conn, self._name, load_type='basic')
    #load_id = load_data[0
    
    table = "SELECT Node.name, LoadNode.* \
            FROM Load, Node, LoadNode \
            WHERE LoadNode.load_id = Load.number \
            AND LoadNode.node_id = Node.number "
            #AND Load.name = ?"
            #AND Node.name = ? ;"
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
        #if load_name:
        #    items = (load_name, node_name, )
        #else:
        #    items = (node_name, )
        items.extend([node_name])
        table += "AND Node.name = ? "    
    #
    # beam
    if beam_flag:
        #items.extend([''])
        table += "AND LoadNode.element_id is NULL"         

    #
    # Node load
    with conn:
        cur = conn.cursor()
        cur.execute(table, tuple(items))
        rows = cur.fetchall()
    #
    node_load = []
    for row in rows:
        if row[7] in ['mass']:
            1 / 0
        else:
            data = [*row[8:14],                      # [fx, fy, fz, mx, my, mz]
                    row[0], row[5], load_name,       # [name, title, load_name]
                    0, 0, 'load']                    # [system, load_complex, load_type]
            #
            # [fx, fy, fz, mx, my, mz, name, title, load_name, system, load_complex, load_type]
            node_load.append(PointNode._make(data))
    #
    return node_load         
#
#
def get_nodedf(conn, load_name: int|str,
               node_name: int|str,
               component: int, 
               beam_flag: bool):
    """nodes in dataframe format"""
    # FIXME : component
    ndata = get_SQLdata(conn, load_name=load_name,
                        node_name=node_name,
                        component=component, 
                        beam_flag = beam_flag)
    
    #
    cols = ['load_name', 'component_name', 
            'load_title', 'node_name', 'element_name',
            'number', 
            'load_id', 'element_id', 'node_id',
            'load_comment', 'load_system','load_type',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
            'x', 'y', 'z', 'rx', 'ry', 'rz',
            'psi', 'B', 'Tw']
    #
    # dataframe
    db = DBframework()
    df = db.DataFrame(data=ndata, columns=cols)
    df['load_level'] = 'basic'
    #df['component'] = component
    #
    df = df[['load_name', 'component_name', 
             'load_title', 'load_level',
             'load_id', 'load_system', 'load_comment',
             'element_name', 'node_name', 'load_type',
             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
             'x', 'y', 'z', 'rx', 'ry', 'rz',
             'psi', 'B', 'Tw']]
    return df