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
from operator import sub, add
#from operator import itemgetter
#import os.path
#from itertools import groupby
from math import dist, sqrt
#
#
# package imports
from steelpy.sections.sqlite.utils import ShapeGeometrySQL
from steelpy.sections.utils.operations import ShapeGeometry
from steelpy.material.sqlite.isotropic import  get_materialSQL
from steelpy.ufo.mesh.sqlite.utils import pull_node
from steelpy.ufo.mesh.sqlite.utils import (push_connectivity,
                                           get_element_data ,
                                           get_connectivity,
                                           get_connectivities,
                                           update_connectivity,
                                           get_node_coord,
                                           get_unitvector, 
                                           update_element_item,
                                           check_element,
                                           get_elements,
                                           push_dircos,
                                           get_roll_angle)

from steelpy.ufo.mesh.process.brotation import unitvec_0, unitvec_1
from steelpy.utils.sqlite.utils import create_connection
from steelpy.ufo.utils.beam import BeamBasic, BeamItemBasic, BeamItemDeformed
from steelpy.ufo.utils.element import get_beam_df
#from steelpy.utils.math.rotation import R_from_drot
#
import numpy as np
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
        #
        # Unit Vector
        coord = get_node_coord(conn, node_name=parameters[:2])
        uvec = unitvec_0(nodei=coord[0], nodej=coord[1],
                         beta=roll_angle)
        
        #uvec2 = unitvec_1(nodei=coord[0], nodej=coord[1],
        #                  beta=roll_angle)
        
        dircos_id = push_dircos(conn, unitvec=uvec,
                                roll_angle=roll_angle, 
                                mesh_id=self._mesh_id)        
        #
        # TODO: 
        eccentricity = [None, None]
        #
        query = (beam_name, self._mesh_id, 'beam', 
                 materials[parameters[2]],
                 sections[parameters[3]],
                 dircos_id)
        #
        table = 'INSERT INTO Element(name, mesh_id, type, \
                                     material_id, section_id,\
                                     dircosine_id)\
                                VALUES(?,?,?,?,?,?) ;'
        cur = conn.cursor()
        cur.execute(table, query)
        beam_id = cur.lastrowid
        #
        # connectivity
        push_connectivity(conn, beam_id,
                          connectivity=parameters[:2],
                          eccentricity=eccentricity, 
                          mesh_id=self._mesh_id)
        #
        #print('-->')
    #
    #def _push_unitvec(self, conn, element_id: int, unitvac: list):
    #    """push unitvector to sql"""
    #    query = [(element_id, *item, x + 1, )
    #             for x, item in enumerate(unitvac)]
    #    table = 'INSERT INTO ElementDirectionCosine( \
    #                         element_id, x, y, z, axis) \
    #             VALUES(?,?,?,?,?)'
    #    #
    #    cur = conn.cursor()
    #    cur.executemany(table, query)
    #    #print('-->')
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
        if (title := data.comment) == None:
            title = ""
        #
        node1, node2 = self.connectivity
        return "{:>8s} {:>8s} {:>8s} {:>12s} {:>12s} {: 1.2e} {:>1.3e} {:}\n"\
               .format(str(beam_name), str(node1), str(node2),
                       str(data.material.name), str(data.section.name),
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
            update_connectivity(conn, data.number, nnodes)
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
            
            #mat = get_materialSQL(conn, data.material,
            #                      mesh_id=self._mesh_id)
        return data.material

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
            #query = (data.section, self._mesh_id, )
            #table = f"SELECT Section.*, SectionGeometry.* \
            #          FROM Section, SectionGeometry, Mesh \
            #          WHERE Section.name = ? \
            #          AND Mesh.number = ? \
            #          AND Section.number = SectionGeometry.section_id ;"
            #
            #cur = conn.cursor()
            #cur.execute(table, query)
            #row = cur.fetchone()
        #
        #geometry = [*row[1:3], *row[11:]]
        #shape = ShapeGeometry(section_type=geometry[2],
        #                      geometry=geometry)
        #shape = data.section
        #sect =  ShapeGeometrySQL(number=shape.number,
        #                         name=shape.name,
        #                         geometry=shape,
        #                         material=data.material,
        #                         db_file=self.db_file)
        return data.section

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
            data = get_roll_angle(conn, self.name,
                                  mesh_id=self._mesh_id)
        return data
    
    @beta.setter
    def beta(self, roll_angle:float):
        """beta angle roll"""
        1 / 0
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
        return data.number
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
    def transform_unit(self, e1:list, e2:list|None=None, e3:list|None=None,
                       warnings:bool=False):
        '''
        Establish transformation matrix from e1 and temporary e2 or e3 vectors.
    
        Arguments
        -----------
        e1 : unit vector describing element longitudinal direction
        e2 : temporary unit vector describing a chosen vector that's perpendicular to the longitudinal direction (approximate y-direction)
        e3 : temporary unit vector describing a chosen vector that's perpendicular to the longitudinal direction (approximate z-direction)
            if both e2 and e3 are different from None, e2 is used (e3 disregarded)
    
        Returns:
            T : transformation matrix
        '''
    
        e1 = np.array(e1).flatten()
    
        if (e2 is not None) and (e3 is not None):
            e3 = None
    
        if e2 is not None:
            e2 = np.array(e2).flatten()
            e3_tmp = np.cross(e1, e2)  # Direction of the third unit vector
            if np.all(e3 == 0):
                if warnings:
                    print('Warning: e1 and e2 identical. Check orientations carefully!')
    
                e2 = e2 + 0.1
                e3 = np.cross(e1, e2)
    
            e2 = np.cross(e3_tmp, e1)  # Direction of the second unit vector
    
            e1 = e1 / np.linalg.norm(e1)  # Normalize the direction vectors to become unit vectors
            e2 = e2 / np.linalg.norm(e2)
            e3 = np.cross(e1, e2)
    
        elif e3 is not None:
            e3 = np.array(e3).flatten()
            e2_tmp = np.cross(e3, e1)  # Direction of the third unit vector
            e3 = np.cross(e1, e2_tmp)
            e1 = e1 / np.linalg.norm(e1)  # Normalize the direction vectors to become unit vectors
            e3 = e3 / np.linalg.norm(e3)
            e2 = np.cross(e3, e1)
    
        if e2 is None and e3 is None:
            raise ValueError('Specify either e2 or e3')
    
        T = np.vstack([e1, e2, e3])
    
        return T
    
    
    #    
    #
    # ------------------------------------------------
    #
    #
    def deformed(self, Fb: list, R: list,
                 L, uv, plane2D: bool):
        """ """
        #1 / 0
        #dim = 3
        #if plane2D:
        #    #dim = 2
        #    coord1_du = [*du[:2], 0.0]  # x,y,z
        #    coord2_du = [*du[3:5], 0.0] # x,y,z
        #    #
        #    rot1 = [0.0, 0.0, du[2]] # rx,ry,rz
        #    rot2 = [0.0, 0.0, du[5]]           
        #else:
        #    #dim = 3
        #    coord1_du = du[:3]
        #    coord2_du = du[6:9]
        #    #
        #    rot1 = du[3:6] # rx,ry,rz
        #    rot2 = du[9:]
        #
        #nodes = self.connectivity
        #conn = create_connection(self.db_file)
        #with conn:
        #    node1 = pull_node(conn, node_name=nodes[0],
        #                      mesh_id=self._mesh_id)
        #    node2 = pull_node(conn, node_name=nodes[1],
        #                      mesh_id=self._mesh_id)
        #
        #coord1 = list(map(add, node1[:3], coord1_du))
        #coord2 = list(map(add, node2[:3], coord2_du))
        #L = dist(coord1, coord2)
        #e1 = list(map(add, coord1, coord2))
        #e1 = [item / L for item in e1]
        #
        # -------------------------------------------------
        #
        #Tb = self.T3D()
        #T0 = Tb[:dim, :dim]
        #t0 = T0[2, :]
        #
        # -------------------------------------------------
        #
        #rot1, rot2 = Fb['rotation']
        #R1, R2 = Fb['R']
        #
        # update base vectors of element from rotations of nodes
        #R1 = R_from_drot(rot1) @ R[0]
        #R1 = R1 @ R1
        #R2 = R_from_drot(rot2) @ R[1]
        #R2 = R2 @ R2
        #
        #
        #e1 = Fb['e1']
        #L = Fb['L']
        #L = sqrt(sum([item ** 2 for item in e1]))
        #
        #e3_temp = R1 @ t0 + R2 @ t0
        #e2 = np.cross(e3_temp, e1) / np.linalg.norm(np.cross(e3_temp, e1))
        #uv = self.transform_unit(e1, e2=e2)
        #uv = Fb['uv']
        #
        # -------------------------------------------------
        # Return deformed beam item class
        beam_deformed = BeamItemDeformed(beam_name=self.name,
                                         material=self.material,
                                         section=self.section,
                                         #du=Fb, #['Fb_local']
                                         L=L, DoF=self.DoF, 
                                         unit_vector=uv,
                                         R = R)
        return beam_deformed
    #
    #
    def deformed2(self, Un: list):
        """ """
        # -------------------------------------------------
        # Return deformed beam item class
        beam_deformed = BeamItemDeformed(beam_name=self.name,
                                         material=self.material,
                                         section=self.section, 
                                         nodes=self.nodes, 
                                         roll_angle=self.beta)
        #
        beam_deformed.get_rot(du=Un)
        #
        #1 / 0
        return beam_deformed
#
#
#

