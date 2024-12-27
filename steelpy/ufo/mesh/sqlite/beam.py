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
#from typing import NamedTuple
from operator import sub #, add
#from operator import itemgetter
#import os.path
#from itertools import groupby
from math import dist
#
#
# package imports
from steelpy.sections.sqlite.utils import ShapeGeometrySQL
from steelpy.sections.utils.operations import ShapeGeometry
from steelpy.material.sqlite.isotropic import  get_materialSQL
from steelpy.ufo.mesh.sqlite.node import pull_node
from steelpy.ufo.mesh.sqlite.utils import (push_connectivity,
                                           get_element_data ,
                                           get_connectivity,
                                           get_connectivities,
                                           update_connectivity,
                                           get_node_coord,
                                           get_unitvector, 
                                           update_element_item,
                                           check_element,
                                           get_elements)

from steelpy.ufo.mesh.process.brotation import unitvec_0, unitvec_1
from steelpy.utils.sqlite.utils import create_connection
from steelpy.ufo.utils.beam import BeamBasic, BeamItemBasic
from steelpy.ufo.utils.element import get_beam_df
#
#
class BeamSQL(BeamBasic):
    __slots__ = ['_title', 'db_file', '_mesh_id']
    
    def __init__(self, db_file:str,
                 mesh_id: int,
                 name:str|int) -> None:
        """
        beam elements
        """
        super().__init__(name=name)
        self.db_file = db_file
        self._mesh_id = mesh_id
    #
    @property
    def _labels(self):
        """ """
        query = ('beam', self._mesh_id)
        table = "SELECT name FROM Element \
                 WHERE type = ? AND mesh_id = ? "
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
        #try:
        self._labels.index(beam_name)
        return BeamItemSQL(beam_name=beam_name,
                           #plane=self._plane,
                           mesh_id=self._mesh_id,
                           db_file=self.db_file)
        #except ValueError:
        #    raise IndexError(f' ** element {beam_name} not valid')    
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
        query = (beam_name, self._mesh_id, 'beam', 
                 materials[parameters[2]],
                 sections[parameters[3]],
                 roll_angle)
        #
        table = 'INSERT INTO Element(name, mesh_id, type, \
                                     material_id, section_id,\
                                     roll_angle)\
                                VALUES(?,?,?,?,?,?) ;'
        cur = conn.cursor()
        cur.execute(table, query)
        beam_number = cur.lastrowid
        #
        # connectivity
        node_id = push_connectivity(conn, beam_number, parameters[:2],
                                    mesh_id=self._mesh_id)
        #
        # Unit Vector
        coord = get_node_coord(conn, node_id)
        uvec = unitvec_0(nodei=coord[0], nodej=coord[1],
                         beta=roll_angle)
        
        #uvec2 = unitvec_1(nodei=coord[0], nodej=coord[1],
        #                  beta=roll_angle)
        
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
        query = (self._mesh_id, )
        table = "SELECT name, number \
                 FROM Material \
                 WHERE mesh_id = ?;"
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
        query = (self._mesh_id, )
        table = "SELECT name, number \
                 FROM Section \
                 WHERE mesh_id = ?;"
        cur = conn.cursor()
        cur.execute(table, query)
        sections = cur.fetchall()
        #
        sections = {item[0]:item[1] for item in sections}
        return sections
    #
    #
    #
    def get_connectivities(self)-> list[int]:
        """
        Return [element_id, node1, node2]
        """
        conn = create_connection(self.db_file)
        with conn:
            nodes = get_connectivities(conn, self._mesh_id)
        return nodes
    #
    #
    # ---------------------------------------
    #
    @property
    def df(self):
        """ """
        #print('elements df out')
        conn = create_connection(self.db_file)
        with conn:
            data = get_elements(conn,
                                mesh_id=self._mesh_id,
                                element_type='beam')
        
        header = ['name', 'number', 'type', 'material', 'section',
                  'node_1', 'node_2', 'node_3', 'node_4',
                  'roll_angle', 'title']
        return data[header]
    
    @df.setter
    def df(self, df):
        """ """
        # clean df
        df = get_beam_df(df)
        #
        query = (self._mesh_id, )
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            table = "SELECT name, number FROM Material \
                     WHERE mesh_id = ? ;"
            #
            cur.execute(table, query)
            materials = cur.fetchall()
            materials = {item[0]:item[1] for item in materials}
            #
            #
            table = "SELECT name, number FROM Section \
                     WHERE mesh_id = ? ;"
            #cur = conn.cursor()
            cur.execute(table, query)
            sections = cur.fetchall()
            sections = {item[0]:item[1] for item in sections}
        #
        df['mesh_id'] = self._mesh_id
        df['material_id'] = df['material'].apply(lambda x: materials[x])
        df['section_id'] = df['section'].apply(lambda x: sections[x])
        #
        if not df.columns.isin(['title']).any():
            df['title'] = None
            
        #
        # Element
        #
        mheader = ['name', 'mesh_id', 'type',
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
            query = ('beam', self._mesh_id)
            table = "SELECT name, number FROM Element \
                     WHERE type = ? AND mesh_id = ? ;"
            #
            cur = conn.cursor()
            cur.execute(table, query)
            elements = cur.fetchall()
            elements = {item[0]:item[1] for item in elements}
            #
            #
            query = (self._mesh_id, )
            table = "SELECT name, number FROM Node \
                     WHERE mesh_id = ? ;"
            #
            #cur = conn.cursor()
            cur.execute(table, query)
            nodes = cur.fetchall()
            nodes = {item[0]:item[1] for item in nodes}
        #
        nheader = ['element_id', 'node_id', 'node_end']
        df['element_id'] = df['name'].apply(lambda x: elements[x])
        #
        df['node_id'] = df['node1'].apply(lambda x: nodes[x])
        df['node_end'] = int(1)
        nodeconn = df[nheader]
        with conn:        
            nodeconn.to_sql('ElementConnectivity', conn,
                            index_label=nheader, 
                            if_exists='append', index=False)
            
        #
        df['node_id'] = df['node2'].apply(lambda x: nodes[x])
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
                 'db_file', '_mesh_id']
    
    def __init__(self, beam_name:int,
                 mesh_id: int, 
                 db_file:str) -> None:
        """
        """
        conn = create_connection(db_file)
        with conn:        
            beam = check_element(conn,
                                 element_name=beam_name,
                                 mesh_id=mesh_id)
        #
        super().__init__(beam_name)
        self.db_file = db_file
        self._mesh_id = mesh_id
        #print('--')
    #
    def __str__(self) -> str:
        """ """
        beam_name = self.name
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, beam_name,
                                    element_type='beam', 
                                    mesh_id=self._mesh_id)
        #title =  data[-1]
        if (title := data[-1]) == None:
            title = ""
        #
        node1, node2 = self.connectivity
        return "{:>8s} {:>8s} {:>8s} {:>12s} {:>12s} {: 1.2e} {:>1.3e} {:}\n"\
               .format(str(beam_name), str(node1), str(node2),
                       str(self.material.name), str(self.section.name),
                       self.beta, self.L, title)
    #    
    # ------------------------------------------------
    #
    @property
    def connectivity(self) -> list:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            connodes = get_connectivity(conn, self.name,
                                        mesh_id=self._mesh_id)
        return connodes

    @connectivity.setter
    def connectivity(self, nodes:list[int]):
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name,
                                    element_type='beam',
                                    mesh_id=self._mesh_id)
            nnodes = []
            for node in nodes:
                item = pull_node(conn, node_name=node,
                                 mesh_id=self._mesh_id)
                nnodes.append(item.number)
            update_connectivity(conn, data[1], nnodes)
        #self._connectivity[self.index] = nodes
    #
    #
    @property
    def material(self) -> list:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name,
                                    element_type='beam', 
                                    mesh_id=self._mesh_id)
            
            mat = get_materialSQL(conn, data[4],
                                  mesh_id=self._mesh_id)
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
                                    mesh_id=self._mesh_id)
            query = (data[5], self._mesh_id, )
            table = f"SELECT Section.*, SectionGeometry.* \
                      FROM Section, SectionGeometry, Mesh \
                      WHERE Section.name = ? \
                      AND Mesh.number = ? \
                      AND Section.number = SectionGeometry.section_id ;"              
            #
            cur = conn.cursor()
            cur.execute(table, query)
            row = cur.fetchone()
        #
        geometry = [*row[1:3], *row[11:]]
        shape = ShapeGeometry(section_type=geometry[2],
                              geometry=geometry)
        sect =  ShapeGeometrySQL(number=row[0], 
                                 name=row[1],
                                 geometry=shape,
                                 material=self.material, 
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
    #
    @property
    def beta(self):
        """beta angle roll"""
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name,
                                    element_type='beam', 
                                    mesh_id=self._mesh_id)
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
    # ------------------------------------------------
    #
    @property
    def number(self) -> int:
        """ """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn, self.name,
                                    element_type='beam', 
                                    mesh_id=self._mesh_id)
        return data[1]
    #
    @property
    def nodes(self) -> list:
        """
        """
        nodes = []
        conn = create_connection(self.db_file)
        with conn:
            for node in self.connectivity:
                nodes.append(pull_node(conn, node_name=node,
                                       mesh_id=self._mesh_id))
        return nodes
    #    
    @property
    def DoF(self) -> list[ int ]:
        """
        """
        conn = create_connection(self.db_file)
        dof = [ ]
        for node_name in self.connectivity:
            node = pull_node(conn, node_name=node_name,
                             mesh_id=self._mesh_id)
            #number = node[0] - 1
            #number = node.index
            dof.append(node.index) #  * 6
        return dof
    #
    @property
    def L(self) -> float:
        """
        """
        nodes = self.connectivity
        conn = create_connection(self.db_file)
        with conn:
            node1 = pull_node(conn, node_name=nodes[0],
                             mesh_id=self._mesh_id)
            node2 = pull_node(conn, node_name=nodes[1],
                             mesh_id=self._mesh_id)
        return dist(node1[:3], node2[:3])
        #return dist(node1[3:6], node2[3:6])
    #
    @property
    def unit_vector(self) -> list[ float ]:
        """
        """
        conn = create_connection(self.db_file)
        with conn:        
            unitvec = get_unitvector(conn,
                                     beam_name=self.name,
                                     mesh_id=self._mesh_id)
        return unitvec
    #
    @property
    def dircosines(self) -> list[ float ]:
        """
        """
        conn = create_connection(self.db_file)
        nodes = self.connectivity
        with conn:
            node1 = pull_node(conn, node_name=nodes[0],
                             mesh_id=self._mesh_id)
            node2 = pull_node(conn, node_name=nodes[1],
                             mesh_id=self._mesh_id)
        # direction cosines
        L = dist(node1[:3], node2[:3])
        uv = list(map(sub, node2[:3], node1[:3]))
        #return np.array([item / L for item in uv])
        return [item / L for item in uv]
    #
    #def T3D(self):
    #    """ """
    #    #nodei, nodej = self.nodes
    #    #return Rmatrix(*self.unit_vector, self.beta)
    #    #return Rmatrix2(nodei, nodej, L=self.L)
    #    unitvec = np.array(self.unit_vector)
    #    return Tmatrix(dirCos=unitvec)
    #    #return Tr
    #
    # ------------------------------------------------
    #    
#
#

