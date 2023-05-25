#
# Copyright (c) 2009-2023 fem2ufo
#

# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
#from itertools import chain
#import math
#from typing import NamedTuple
from operator import sub, add
from operator import itemgetter
#import os.path
from itertools import groupby
from math import dist

# package imports
from steelpy.f2uModel.mesh.process.process_sql import create_connection, check_nodes
#from steelpy.beam.main import BasicCalcs
from .nodes import get_node, NodeSQL
from steelpy.sections.sqlite.main import get_section, SectionSQL
from steelpy.material.sqlite.isotropic import  get_materialSQL, MaterialSQL
from ..process.elements.beam import BeamBasic, BeamItemBasic
#
from ..process.Kmatrix.stiffness import beam_stiffness
# (beam_Klocal, trans_3d_beam, Rmatrix, Rmatrix_new, trans3Dbeam)
#


class BeamSQL(BeamBasic):
    __slots__ = ['_labels', '_number',  '_title', 'db_file']
    
    def __init__(self, db_file:str) -> None:
        """
        beam elements
        """
        super().__init__()
        self.db_file = db_file
    #
    #
    def __setitem__(self, beam_name: int|str, parameters: list) -> None:
        """
        parameters = ['beam', node1, node2, material, section, roll_angle]
        """
        try:
            self._labels.index(beam_name)
            raise Exception('element {:} already exist'.format(beam_name))
        except ValueError:
            # default
            self._labels.append(beam_name)
            # push to SQL
            conn = create_connection(self.db_file)
            with conn:
                self.push_beam(conn, beam_name, parameters)
                #conn.commit()
    #
    def __getitem__(self, beam_name: int|str):
        """ """
        try:
            self._labels.index(beam_name)
            return BeamItemSQL(beam_name, self.db_file)
        except ValueError:
            raise IndexError(' ** element {:} does not exist'.format(beam_name))    
    #
    #
    def push_beam(self, conn, beam_name: int|str, parameters):
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
        #
        #try:
        roll_angle = parameters[4]
        #except IndexError:
        #    roll_angle = 0.0
        #print('-->')
        #if (title := parameters[5]) == "NULL":
        #    title = None
        title = None
        #
        project = (beam_name, title, 'beam', 
                   materials[str(parameters[2])],
                   sections[str(parameters[3])],
                   roll_angle)
        #
        sql = 'INSERT INTO tb_Elements(name, title, type, material_number, section_number,\
                                       roll_angle)\
                                       VALUES(?,?,?,?,?,?)'
        #cur = conn.cursor()
        cur.execute(sql, project)
        #
        # connectivity
        beam_number = cur.lastrowid
        push_connectivity(conn, beam_number, parameters[:2])
    #
    @property
    def get_connectivities(self):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            nodes = get_nodes_connec(conn)
        return nodes
    #
#
#
def get_nodes_connec(conn):
    """ """
    cur = conn.cursor()
    cur.execute("SELECT tb_Elements.number, tb_Connectivity.node_end, tb_Nodes.name \
                FROM tb_Connectivity, tb_Elements, tb_Nodes \
                WHERE tb_Elements.number = tb_Connectivity.element_number \
                AND tb_Nodes.number = tb_Connectivity.node_number;")
    connodes = cur.fetchall()
    connodes.sort(key = itemgetter(0))
    #
    groups = groupby(connodes, itemgetter(0))
    nodes = [[item[1:] for item in data]
             for (key, data) in groups]
    #
    nodes = [[x for _, x in sorted(items)]
             for items in nodes]
    return nodes
#
#
def get_beam_data(conn, element_name):
    """ """
    cur = conn.cursor()
    cur.execute ("SELECT tb_Elements.name, tb_Elements.number, tb_Elements.type,\
                tb_Elements.roll_angle, tb_Materials.name, tb_Sections.name, tb_Elements.title\
                FROM tb_Elements, tb_Materials, tb_Sections\
                WHERE tb_Elements.name = {:} \
                AND tb_Elements.material_number = tb_Materials.number \
                AND tb_Elements.section_number = tb_Sections.number;".format(element_name))
    row = cur.fetchone()
    #
    #element_number = row[1]
    connodes = get_connectivity(conn, element_name)
    data = [*row[:6], connodes, row[-1]]
    #conn.close ()
    return data
#
#
def push_connectivity(conn, element_number: int, connectivity: list):
    """
    """
    #nconn = [check_nodes(conn, node_name)[0]
    #        for node_name in connectivity]
    cur = conn.cursor()
    for x, item in enumerate(connectivity):
        node_number = check_nodes(conn, item)[0]
        project = (element_number, node_number, x+1)
        sql = 'INSERT INTO  tb_Connectivity(element_number,\
                                            node_number, node_end)\
                                            VALUES(?,?,?)'
        cur.execute(sql, project)
    #return cur.lastrowid
#
def get_connectivity(conn, element_name: int):
    """ """
    cur = conn.cursor()
    cur.execute("SELECT tb_Connectivity.node_end, tb_Nodes.name \
                FROM tb_Connectivity, tb_Nodes, tb_Elements \
                WHERE tb_Nodes.number = tb_Connectivity.node_number\
                AND tb_Elements.number = tb_Connectivity.element_number\
                AND tb_Elements.name = {:};".format(element_name))
    connodes = cur.fetchall()
    connodes = [x for _, x in sorted(connodes)]
    return connodes
#
def update_connectivity(conn, element_number: int, connectivity: list):
    """
    """
    #1 / 0
    cur = conn.cursor()
    for x, node in enumerate(connectivity):
        project = (node, element_number, x+1)
        sql = 'UPDATE tb_Connectivity SET node_number = ? \
               WHERE element_number = ?\
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
#
#
@dataclass
class BeamItemSQL(BeamItemBasic):
    """ """
    #__slots__ = ['name', 'db_file']
    
    def __init__(self, element_name:int, db_file:str) -> None:
        """
        """
        super().__init__(element_name)
        #
        self.db_file = db_file
        self.type: str = "beam"
        #self._f2u_sections = SectionSQL(db_file)
        #self._f2u_materials =  MaterialSQL(db_file)
        #self._f2u_nodes =  NodeSQL(db_file)
    #
    @property
    def number(self) -> int:
        """ """
        conn = create_connection(self.db_file)
        with conn:
            data = get_beam_data(conn, self.name)
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
    def connectivity(self) -> list:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            connodes = get_connectivity(conn, self.name)
        return connodes

    @connectivity.setter
    def connectivity(self, nodes:list[int]) -> list:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            #push_connectivity(conn, self.name, nodes)
            update_connectivity(conn, self.name, nodes)
        #self._connectivity[self.index] = nodes
    #
    @property
    def nodes(self) -> list:
        """
        """
        _nodes = []
        conn = create_connection(self.db_file)
        for _node in self.connectivity:
            with conn:
                _nodes.append(get_node(conn, node_name=_node))
        return _nodes    
    #
    @property
    def material(self) -> list:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            data = get_beam_data(conn, self.name)
            mat = get_materialSQL(conn, data[4])
        #return self._f2u_materials[data[4]]
        return mat

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
    def section(self) -> list:
        """
        """
        #BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        conn = create_connection(self.db_file)
        with conn:
            data = get_beam_data(conn, self.name)
            sect =  get_section(conn, data[5])
        #return self._f2u_sections[data[5]]
        return sect

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
            data = get_beam_data(conn, self.name)
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
            data = get_beam_data(conn, self.name)
        #title =  data[-1]
        if (title := data[-1]) == None:
            title = ""
        #
        #return "{:8d} {:8d} {:8d} {:>12s} {:>12s} {: 6.4f} {:>6.3f} {:>12s}\n"\
        #       .format(self.name, *self.connectivity,
        #               self.material, self.section, self.beta,
        #               self.length, title)
        return "{:8d} {:8d} {:8d} {:>12s} {:>12s} {: 1.2e} {:>1.3e} {:>12s}"\
               .format(self.name, *self.connectivity,
                       self.material.name, self.section.name,
                       self.beta, self.L, title)
    #
    #
    @property
    def DoF(self) -> list[ int ]:
        """
        """
        conn = create_connection(self.db_file)
        dof = [ ]
        for node_name in self.connectivity:
            node = get_node(conn, node_name=node_name)
            #number = node[0] - 1
            number = node.index
            dof.append(number) #  * 6
        return dof
    #
    @property
    def L(self) -> float:
        """
        """
        nodes = self.connectivity
        conn = create_connection(self.db_file)
        with conn:
            node1 = get_node(conn, node_name=nodes[0])
            node2 = get_node(conn, node_name=nodes[1])
        return dist(node1[:3], node2[:3])
        #return dist(node1[3:6], node2[3:6])
    #
    #
    @property
    def unit_vector(self) -> list[ float ]:
        """
        """
        nodes = self.connectivity
        conn = create_connection(self.db_file)
        with conn:
            node1 = get_node(conn, node_name=nodes[0])
            node2 = get_node(conn, node_name=nodes[1])
        # direction cosines
        L =  dist(node1[:3], node2[:3])
        uv = list(map(sub, node2[:3], node1[:3]))
        return [item / L for item in uv]       
    #
#
#