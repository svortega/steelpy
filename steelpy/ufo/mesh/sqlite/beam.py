#
# Copyright (c) 2009 steelpy
#
#
# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections import Counter
#from collections.abc import Mapping
from dataclasses import dataclass
#from itertools import chain
#import math
from typing import NamedTuple
#from operator import sub, add
#from operator import itemgetter
#import os.path
#from itertools import groupby
from math import dist
#
#
# package imports
from steelpy.sections.sqlite.utils import ShapeGeometrySQL #get_section 
from steelpy.material.sqlite.isotropic import  get_materialSQL
from steelpy.ufo.mesh.sqlite.nodes import get_node
from steelpy.ufo.mesh.sqlite.utils import (push_connectivity,
                                           get_element_data ,
                                           get_connectivity,
                                           update_connectivity,
                                           get_node_coord,
                                           get_unitvector, 
                                           update_element_item)
from steelpy.ufo.mesh.process.elements.beam import BeamBasic, BeamItemBasic
from steelpy.ufo.mesh.process.elements.bstiffness import Tmatrix, unitvec_0
from steelpy.utils.sqlite.utils import create_connection
#
import numpy as np
#
#
class BeamSQL(BeamBasic):
    __slots__ = ['_title', 'db_file', '_plane', '_component']
    
    def __init__(self, db_file:str,
                 component: int, 
                 plane: NamedTuple) -> None:
        """
        beam elements
        """
        super().__init__()
        self.db_file = db_file
        self._plane = plane
        self._component = component
    #
    @property
    def _labels(self):
        """ """
        query = ('beam', self._component)
        table = "SELECT name FROM Element \
                 WHERE type = ? AND component_id = ? "
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]    
    #
    def __setitem__(self, beam_name: int|str,
                    parameters: list) -> None:
        """
        beam_name : beam id
        parameters = [node1, node2, material, section, roll_angle]
        """
        try:
            self._labels.index(beam_name)
            raise Exception(f'element {beam_name} already exist')
        except ValueError:
            # default
            # push to SQL
            conn = create_connection(self.db_file)
            with conn:
                self._push_data(conn, beam_name, parameters)
    #
    def __getitem__(self, beam_name: int|str):
        """ """
        try:
            self._labels.index(beam_name)
            return BeamItemSQL(beam_name=beam_name,
                               plane=self._plane,
                               component=self._component,
                               db_file=self.db_file)
        except ValueError:
            raise IndexError(f' ** element {beam_name} not valid')    
    #
    # ---------------------------------------
    #
    def _push_data(self, conn, beam_name: int|str, parameters):
        """ """
        materials = self._pull_material(conn)
        #
        sections = self._pull_section(conn)
        #
        try:
            roll_angle = parameters[4]
        except IndexError:
            roll_angle = 0.0
        #
        #if (title := parameters[5]) == "NULL":
        #    title = None
        #try:
        #    title = parameters[5]
        #except IndexError:
        #    title = None
        #
        query = (beam_name, self._component, 'beam', 
                 materials[parameters[2]],
                 sections[parameters[3]],
                 roll_angle)
        #
        table = 'INSERT INTO Element(name, component_id, type, \
                                        material_id, section_id,\
                                        roll_angle)\
                                VALUES(?,?,?,?,?,?) ;'
        cur = conn.cursor()
        cur.execute(table, query)
        #
        # connectivity
        beam_number = cur.lastrowid
        node_id = push_connectivity(conn, beam_number, parameters[:2],
                                    component=self._component)
        #
        # Unit Vector
        coord = get_node_coord(conn, node_id)
        uvec = unitvec_0(nodei=coord[0], nodej=coord[1])
        self._push_unitvec(conn, element_id=beam_number,
                           unitvac=uvec)
        #
        #print('-->')
    #
    def _push_unitvec(self, conn, element_id: int, unitvac: list):
        """push unitvector to sql"""
        query = [(element_id, *item, x + 1, )
                 for x, item in enumerate(unitvac)]
        table = 'INSERT INTO ElementDirectionCosine( \
                             element_id, x, y, z, axis) \
                 VALUES(?,?,?,?,?)'
        #
        cur = conn.cursor()
        cur.executemany(table, query)
        #print('-->')
    #    
    #
    def _pull_material(self, conn):
        """ """
        query = (self._component, )
        table = "SELECT name, number \
                 FROM Material \
                 WHERE component_id = ?;"
        #
        cur = conn.cursor()
        cur.execute(table, query)
        materials = cur.fetchall()
        #
        materials = {item[0]:item[1] for item in materials}
        return materials
    #
    def _pull_section(self, conn):
        """ """
        query = (self._component, )
        table = "SELECT name, number \
                 FROM Section \
                 WHERE component_id = ?;"
        cur = conn.cursor()
        cur.execute(table, query)
        sections = cur.fetchall()
        #
        sections = {item[0]:item[1] for item in sections}
        return sections
    #
    #
    #
    def get_connectivities(self):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            nodes = get_nodes_connec(conn)
        return nodes
    #
    #
    # ---------------------------------------
    #
    @property
    def df(self):
        """ """
        print('--->')
        1 / 0
    
    @df.setter
    def df(self, df):
        """ """
        query = (self._component, )
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            table = "SELECT name, number FROM Material \
                     WHERE component_id = ? ;"
            #
            cur.execute(table, query)
            materials = cur.fetchall()
            materials = {item[0]:item[1] for item in materials}
            #
            #
            table = "SELECT .name, number FROM Section \
                     WHERE component_id = ? ;"
            #cur = conn.cursor()
            cur.execute(table, query)
            sections = cur.fetchall()
            sections = {item[0]:item[1] for item in sections}
        #
        df['component_id'] = self._component
        df['material_id'] = df['material_name'].apply(lambda x: materials[x])
        df['section_id'] = df['section_name'].apply(lambda x: sections[x])
        #
        if not df.columns.isin(['title']).any():
            df['title'] = None
            
        #
        # Element
        #
        mheader = ['name', 'component_id', 'type',
                  'material_id', 'section_id',
                  'roll_angle', 'title']
        members = df[mheader]
        #
        with conn:
            members.to_sql('Element', conn,
                           index_label=mheader, 
                           if_exists='append', index=False)
        #
        # ElementConnectivity
        #
        with conn:
            query = ('beam', self._component)
            table = "SELECT name, number FROM Element \
                     WHERE type = ? WHERE component_id = ? ;"
            #
            cur = conn.cursor()
            cur.execute(table, query)
            elements = cur.fetchall()
            elements = {item[0]:item[1] for item in elements}
            #
            #
            query = (self._component)
            table = "SELECT name, number FROM Node \
                     WHERE component_id = ?;"
            #
            cur = conn.cursor()
            cur.execute(table, query)
            nodes = cur.fetchall()
            nodes = {item[0]:item[1] for item in nodes}
        #
        nheader = ['element_id', 'node_id', 'node_end']
        df['element_id'] = df['name'].apply(lambda x: elements[x])
        #
        df['node_id'] = df['node_1'].apply(lambda x: elements[x])
        df['node_end'] = int(1)
        nodeconn = df[nheader]
        with conn:        
            nodeconn.to_sql('ElementConnectivity', conn,
                            index_label=nheader, 
                            if_exists='append', index=False)
            
        #
        df['node_id'] = df['node_2'].apply(lambda x: elements[x])
        df['node_end'] = int(2)
        nodeconn = df[nheader]
        with conn:        
            nodeconn.to_sql('ElementConnectivity', conn,
                            index_label=nheader, 
                            if_exists='append', index=False)        
        #
        #print('--->')
        #self._labels.extend(df['name'].tolist())
        #self._type.extend(df['type'].tolist())
        #1 / 0
#
#
@dataclass
class BeamItemSQL(BeamItemBasic):
    """ """
    __slots__ = ['name', '_releases', 'type',
                 'db_file', '_plane', '_component']
    
    def __init__(self, beam_name:int, plane: NamedTuple,
                 component: int, 
                 db_file:str) -> None:
        """
        """
        super().__init__(beam_name)
        self._plane = plane
        self.db_file = db_file
        self._component = component
    #
    @property
    def number(self) -> int:
        """ """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name,
                                    element_type='beam', 
                                    component=self._component)
        return data[1]
    
    #@number.setter
    #def number(self, value:int) -> None:
    #    """"""
    #    1/0
    #    conn = create_connection(self.db_file)
    #    item = "number"
    #    with conn:
    #        update_element_item(conn, self.name, item, value)
    #
    @property
    def connectivity(self) -> list:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            connodes = get_connectivity(conn, self.name,
                                        component=self._component)
        return connodes

    @connectivity.setter
    def connectivity(self, nodes:list[int]):
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
        nodes = []
        conn = create_connection(self.db_file)
        for node in self.connectivity:
            with conn:
                nodes.append(get_node(conn, node_name=node,
                                       component=self._component))
        return nodes
    #
    @property
    def material(self) -> list:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name,
                                    element_type='beam', 
                                    component=self._component)
            
            mat = get_materialSQL(conn, data[4],
                                  component=self._component)
        #
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
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name,
                                    element_type='beam', 
                                    component=self._component)
            
            #sect =  self._pull_section_section(conn, data[5],
            #                                   component=self._component)
            query = (data[5], self._component, )
            #table = f'SELECT Section.* FROM Section \
            #            WHERE Section.name = ? \
            #            AND component_id = ?'
            #
            table = f"SELECT Section.*, SectionGeometry.* \
                      FROM Section, SectionGeometry, Component \
                      WHERE Section.name = ? \
                      AND Component.number = ? \
                      AND Section.number = SectionGeometry.section_id ;"              
            #
            cur = conn.cursor()
            cur.execute(table, query)
            row = cur.fetchone()             
        #
        geometry = [*row[1:3], *row[11:]]
        #
        sect =  ShapeGeometrySQL(number=row[0], 
                                 name=row[1],
                                 geometry=geometry, 
                                 db_file=self.db_file)
        return sect

    @section.setter
    def section(self, section_name: str) -> None:
        """
        """
        conn = create_connection(self.db_file)
        item = "section"
        with conn:
            update_element_item(conn, self.name, item, section_name)
    #
    def _pull_section(self):
        """ get section """
        
    #
    #
    @property
    def beta(self):
        """beta angle roll"""
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name,
                                    element_type='beam', 
                                    component=self._component)
        return data[3]
    
    @beta.setter
    def beta(self, roll_angle:float):
        """beta angle roll"""
        conn = create_connection(self.db_file)
        item = "roll_angle"
        with conn:
            update_element_item(conn, self.name, item, roll_angle)
    #
    #
    def __str__(self) -> str:
        """ """
        beam_name = self.name
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, beam_name,
                                    element_type='beam', 
                                    component=self._component)
        #title =  data[-1]
        if (title := data[-1]) == None:
            title = ""
        #
        #return "{:8d} {:8d} {:8d} {:>12s} {:>12s} {: 6.4f} {:>6.3f} {:>12s}\n"\
        #       .format(self.name, *self.connectivity,
        #               self.material, self.section, self.beta,
        #               self.length, title)
        node1, node2 = self.connectivity
        return "{:>8s} {:>8s} {:>8s} {:>12s} {:>12s} {: 1.2e} {:>1.3e} {:}\n"\
               .format(str(beam_name), str(node1), str(node2),
                       str(self.material.name), str(self.section.name),
                       self.beta, self.L, title)
    #
    # ------------------------------------------------
    #
    @property
    def DoF(self) -> list[ int ]:
        """
        """
        conn = create_connection(self.db_file)
        dof = [ ]
        for node_name in self.connectivity:
            node = get_node(conn, node_name=node_name,
                            component=self._component)
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
            node1 = get_node(conn, node_name=nodes[0],
                             component=self._component)
            node2 = get_node(conn, node_name=nodes[1],
                             component=self._component)
        return dist(node1[:3], node2[:3])
        #return dist(node1[3:6], node2[3:6])
    #
    #
    @property
    def unit_vector(self) -> list[ float ]:
        """
        """
        #nodes = self.connectivity
        #conn = create_connection(self.db_file)
        #with conn:
        #    node1 = get_node(conn, node_name=nodes[0],
        #                     component=self._component)
        #    node2 = get_node(conn, node_name=nodes[1],
        #                     component=self._component)
        # direction cosines
        #L =  dist(node1[:3], node2[:3])
        #uv = list(map(sub, node2[:3], node1[:3]))
        #return [item / L for item in uv]
        #
        conn = create_connection(self.db_file)
        with conn:        
            unitvec = get_unitvector(conn, beam_id=self.name)
        #
        return unitvec
    #
    #
    def T3D(self):
        """ """
        #nodei, nodej = self.nodes
        #return Rmatrix(*self.unit_vector, self.beta)
        #return Rmatrix2(nodei, nodej, L=self.L)
        unitvec = np.array(self.unit_vector)
        return Tmatrix(dirCos=unitvec)
        #return Tr
#
#

