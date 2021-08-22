# Copyright (c) 2009-2021 fem2ufo

# 
# Python stdlib imports
from dataclasses import dataclass
from array import array
from collections import Counter
from collections.abc import Mapping
from math import dist
from typing import NamedTuple, Tuple, List, Iterator, Iterable, Union, Dict
from itertools import chain


# package imports
from steelpy.f2uModel.mesh.sqlite.nodes import get_node
from steelpy.f2uModel.material.matsql import get_materialSQL
from steelpy.f2uModel.sections.main import get_sectionSQL
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table
from steelpy.trave3D.preprocessor.assemble import (beam_stiffness, beam_Ks,
                                                   trans_3d_beam, Rmatrix)

#
#
@dataclass
class BeamElement:
    """ """
    __slots__ = ['name', 'db_file', 'type']
    
    def __init__(self, element_name:int, db_file:str) -> None:
        """
        """    
        self.name = element_name
        self.db_file = db_file
        self.type: str = "beam"
    #
    @property
    def number(self) -> int:
        """ """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name)
        return data[1]
    
    @number.setter
    def number(self, value:int) -> None:
        """"""
        1/0
        conn = create_connection(self.db_file)
        item = "number"
        with conn:
            update_element_item(conn, self.name, item, value)
    #
    @property
    def connectivity(self) -> List:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            connodes = get_connectivity(conn, self.name)
        return connodes

    @connectivity.setter
    def connectivity(self, nodes:List[int]) -> List:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            #push_connectivity(conn, self.name, nodes)
            update_connectivity(conn, self.name, nodes)
        #self._connectivity[self.index] = nodes
    #
    @property
    def material(self) -> List:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name)
        return data[4]

    @material.setter
    def material(self, material_name: str) -> None:
        """
        """
        conn = create_connection(self.db_file)
        item = "material"
        with conn:
            update_element_item(conn, self.name, item, material_name)

    #
    @property
    def section(self) -> List:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name)
        return data[5]

    @section.setter
    def section(self, section_name: str) -> None:
        """
        """
        conn = create_connection(self.db_file)
        item = "section"
        with conn:
            update_element_item(conn, self.name, item, self.name)

    #
    @property
    def beta(self):
        """beta angle roll"""
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name)
        return data[3]
    
    @beta.setter
    def beta(self, value):
        """beta angle roll"""
        conn = create_connection(self.db_file)
        item = "roll_angle"
        with conn:
            update_element_item(conn, self.name, item, self.name)
    #
    #
    def __str__(self) -> str:
        """ """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name)
        #title =  data[-1]
        if (title := data[-1]) == "NULL":
            title = ""        
        #
        return "{:8d} {:8d} {:8d} {:>12s} {:>12s} {: 6.4f} {:>6.3f} {:>12s}\n"\
               .format(self.name, *self.connectivity,
                       self.material, self.section, self.beta,
                       self.length, title)
    #
    #
    @property
    def DoF(self) -> List[ int ]:
        """
        """
        conn = create_connection(self.db_file)
        dof = [ ]
        for node_name in self.connectivity:
            node = get_node(conn, node_name=node_name)
            number = node[0] - 1
            dof.append(number) #  * 6
        return dof
    #
    @property
    def length(self) -> float:
        """
        """
        conn = create_connection(self.db_file)
        nodes = self.connectivity
        node1 = get_node(conn, node_name=nodes[0])
        node2 = get_node(conn, node_name=nodes[1])
        # _dx = _node1.x - _node2.x
        # _dy = _node1.y - _node2.y
        # _dz = _node1.z - _node2.z
        # dist2 = (_dx * _dx + _dy * _dy + _dz * _dz)**0.50
        return dist(node1[3:6], node2[3:6])
    #
    @property
    def unit_vector(self) -> List[ float ]:
        """
        """
        # TODO: get_node should be aligned with inmemmory
        conn = create_connection(self.db_file)
        node1 = get_node(conn, node_name=self.connectivity[0])
        node2 = get_node(conn, node_name=self.connectivity[1])
        dx = node2[3] - node1[3]
        dy = node2[4] - node1[4]
        dz = node2[5] - node1[5]
        # direction cosines
        L = dist(node1[3:6], node2[3:6])
        l = dx / L
        m = dy / L
        n = dz / L
        return [l, m, n]
    #
    @property
    def Kmatrix(self):
        """ """
        #conn = create_connection(self.db_file)
        material, section, beta = self._K_data()
        #section = get_sectionSQL(conn, self.section)
        #material = get_materialSQL(conn, self.material)
        # solve K matrix
        R = Rmatrix(*self.unit_vector, beta)
        # R = Rmatrix(*self.direction_cosines, self.beta)
        # K = beam_stiffness(self.length,
        #                   section.area,
        #                   section.J,
        #                   section.Iy,
        #                   section.Iz,
        #                   material.E,
        #                   material.G)
        K = beam_Ks(self.length,
                    section.area, section.J,
                    section.Iy, section.Iz,
                    material.E, material.G,
                    section.area, section.area)
        return trans_3d_beam(K, R)
    #
    def _K_data(self):
        """ """
        conn = create_connection(self.db_file)
        cur = conn.cursor()
        cur.execute ("SELECT * FROM tb_Elements\
                    WHERE tb_Elements.name = {:};".format(self.name))
        row = cur.fetchone()
        #
        #connodes = get_connectivity(conn, self.name)
        #data = [*row[4:], connodes]
        material = get_materialSQL(conn, row[4]) 
        section = get_sectionSQL(conn, row[5])
        beta = row[6]
        conn.close()
        return material, section, beta
    #
    @property
    def R(self):
        """
        Rotation matrix
        """
        if self.type in ['beam', 'truss']:
            return Rmatrix(*self.unit_vector, self.beta)
        else:
            raise IOError("no yet included")
    #
#
#
#
class ElementSQL(Mapping):
    __slots__ = ['db_file', '_labels']
    
    def __init__(self, db_file:str,
                 db_system:str="sqlite") -> None:
        """
        """
        self.db_file = db_file
        self._labels: array = array('I', [])
        # create node table
        self._create_table()
    #
    def __setitem__(self, element_number: int, parameters: List) -> None:
        """
        parameters = ['beam', node1, node2, material, section, roll_angle]
        """
        try:
            self._labels.index(element_number)
            raise Exception('element {:} already exist'.format(element_number))
        except ValueError:
            # default
            self._labels.append(element_number)
            # push to SQL
            conn = create_connection(self.db_file)
            with conn:
                self.push_element(conn, element_number, parameters)
                conn.commit()
    #
    def __getitem__(self, element_number: int):
        """ """
        try:
            self._labels.index(element_number)
            return BeamElement(element_number, self.db_file)
        except ValueError:
            raise IndexError(' ** element {:} does not exist'.format(element_number))    
    #
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
    def push_element(self, conn, element_number, parameters):
        """ """
        cur = conn.cursor()
        cur.execute("SELECT tb_Materials.name, tb_Materials.number FROM tb_Materials;")
        materials = cur.fetchall()
        materials = {item[0]:item[1] for item in materials}
        #
        #cur = conn.cursor()
        cur.execute("SELECT tb_Sections.name, tb_Sections.number FROM tb_Sections;")
        sections = cur.fetchall()
        sections = {item[0]:item[1] for item in sections}
        # connectivity
        push_connectivity(conn, element_number, parameters[1:3])
        #
        #try:
        roll_angle = parameters[5]
        #except IndexError:
        #    roll_angle = 0.0
        #print('-->')
        if (title := parameters[6]) == "NULL":
            title = None
        #
        project = (element_number, title, 
                   parameters[0],
                   materials[parameters[3]],
                   sections[parameters[4]],
                   roll_angle)
        #
        sql = 'INSERT INTO tb_Elements(name, title, type, material, section,\
                                       roll_angle)\
                                       VALUES(?,?,?,?,?,?)'
        #cur = conn.cursor()
        cur.execute(sql, project)
    #
    def _create_table(self) -> None:
        """ """
        _table_elements = "CREATE TABLE IF NOT EXISTS tb_Elements(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            name INTEGER NOT NULL,\
                            title TEXT,\
                            type TEXT NOT NULL,\
                            material INTEGER NOT NULL REFERENCES tb_Materials(number),\
                            section INTEGER NOT NULL REFERENCES tb_Sections(number),\
                            roll_angle DECIMAL);"
        #
        _table_connectivity = "CREATE TABLE IF NOT EXISTS tb_Connectivity(\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                element_name INTEGER NOT NULL REFERENCES tb_Elements(name),\
                                node_name INTEGER REFERENCES tb_Nodes(name),\
                                node_end INTEGER NOT NULL);"
        #
        _table_univectors = "CREATE TABLE IF NOT EXISTS tb_DirectionCosines(\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                element_name INTEGER NOT NULL REFERENCES tb_Elements(name),\
                                type TEXT NOT NULL);"
        #
        _table_offset = "CREATE TABLE IF NOT EXISTS tb_Eccentricities(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            element_name INTEGER NOT NULL REFERENCES tb_Elements(name),\
                            node_name INTEGER REFERENCES tb_Nodes(name),\
                            node_end INTEGER NOT NULL,\
                            system TEXT NOT NULL,\
                            x DECIMAL,\
                            y DECIMAL,\
                            z DECIMAL);"
        #
        conn = create_connection(self.db_file)
        create_table(conn, _table_elements)
        create_table(conn, _table_connectivity)
        create_table(conn, _table_offset)
        create_table(conn, _table_univectors)
    #
    #def iter_elements(self, arraysize=1000):
    #    """
    #    """
    #    conn = create_connection(self.db_file)
    #    cur = conn.cursor()
    #    # TODO: check if direction cosines given
    #    cur.execute("SELECT tb_Elements.name, tb_Elements.number, tb_Elements.type,\
    #                tb_Elements.roll_angle, tb_Elements.material, tb_Elements.section\
    #                FROM tb_Elements;" )
    #    #
    #    try:
    #        while True:
    #            elements = cur.fetchmany(arraysize)
    #            if not elements:
    #                break
    #            for element in elements:
    #                #cur.execute("SELECT tb_Connectivity.node_end, tb_Connectivity.node_name\
    #                #            FROM tb_Connectivity\
    #                #            WHERE tb_Connectivity.element_name = {:};".format(element[0]))
    #                #row = cur.fetchall()
    #                #connodes = [x for _, x in sorted(row)]
    #                connodes = get_connectivity(conn, element[0])
    #                data = [*element[0:6], connodes, self.db_file]
    #                yield BeamElement(data)
    #    except Exception as e:
    #        print(e)
    #    finally:
    #        conn.close()
    #
    @property
    def get_connectivities(self):
        """ """
        conn = create_connection(self.db_file)
        cur = conn.cursor()
        cur.execute( "SELECT tb_Elements.name FROM tb_Elements;")
        elements = cur.fetchall()
        connodes = []
        for element in elements:
            #cur.execute("SELECT tb_Connectivity.node_end, tb_Connectivity.node_name\
            #            FROM tb_Connectivity\
            #            WHERE tb_Connectivity.element_name = {:};".format(member[0]))
            #row = cur.fetchall()
            #connodes.append([x for _,x in sorted(row)])
            connodes.append(get_connectivity(conn, element[0]))
        conn.close()
        return connodes
    #
    #
    def get_number(self, start:int=0)-> Iterable[int]:
        """
        """
        try:
            n = max(self._labels)
        except ValueError:
            n = start
        #
        while True:
            n += 1
            yield n    
    #
    #
    def update_item(self, element_number:int, item:str, value:Union[float,int]):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            update_element_item(conn, element_number, item, value)
            #conn.commit()
    #
    @property
    def get_free_nodes(self):
        """
        find nodes not sharing elements
        """
        connectivities = self.get_connectivities
        #connectivities = [conn for conn in connectivities.values()]
        #column
        flat = list(chain.from_iterable(connectivities))
        return [k for k, v in Counter(flat).items() if v == 1]
#
#
def get_connectivity(conn, element_name):
    """ """
    cur = conn.cursor()
    cur.execute("SELECT tb_Connectivity.node_end, tb_Connectivity.node_name\
                FROM tb_Connectivity\
                WHERE tb_Connectivity.element_name = {:};".format(element_name))
    connodes = cur.fetchall()
    return [x for _, x in sorted(connodes)]
    #return connodes
#
def push_connectivity(conn, element_name, connectivity):
    """
    """
    cur = conn.cursor()
    for x, node in enumerate(connectivity):
        project = (element_name, node, x+1)
        sql = 'INSERT INTO  tb_Connectivity(element_name,\
                                            node_name, node_end)\
                                            VALUES(?,?,?)'
        cur.execute(sql, project)
    #return cur.lastrowid
#
def update_connectivity(conn, element_name, connectivity):
    """
    """
    cur = conn.cursor()
    for x, node in enumerate(connectivity):
        project = (node, element_name, x+1)
        sql = 'UPDATE tb_Connectivity SET node_name = ? \
               WHERE element_name = ?\
               AND node_end = ?'
        cur.execute(sql, project)
    #return cur.lastrowid
#
#
def update_element_item(conn, name, item, value):
    """ """
    project = (value, name)
    sql = 'UPDATE tb_Elements SET {:} = ? WHERE name = ?'.format(item)
    cur = conn.cursor()
    cur.execute(sql, project)
#
#
def get_element_data(conn, element_name):
    """ """
    cur = conn.cursor()
    cur.execute ("SELECT tb_Elements.name, tb_Elements.number, tb_Elements.type,\
                tb_Elements.roll_angle, tb_Materials.name, tb_Sections.name, tb_Elements.title\
                FROM tb_Elements, tb_Materials, tb_Sections\
                WHERE tb_Elements.name = {:} \
                AND tb_Elements.material = tb_Materials.number \
                AND tb_Elements.section = tb_Sections.number;".format(element_name))
    row = cur.fetchone()
    #
    connodes = get_connectivity(conn, element_name)
    data = [*row[:6], connodes, row[-1]]
    #conn.close ()
    return data
#
#