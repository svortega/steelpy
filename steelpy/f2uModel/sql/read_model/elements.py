#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
import sqlite3 as sqlite3
from math import dist
from typing import NamedTuple, Union, List, Tuple, Dict
from collections.abc import Mapping
import itertools
#

# package imports
from steelpy.trave3D.preprocessor.assemble import (beam_stiffness, beam_Ks,
                                                   trans_3d_beam, Rmatrix)
from steelpy.f2uModel.sql.operation.process_sql import create_connection


class BeamElement(NamedTuple):
    """ Cartesian coordinate system"""
    #__slots__ = ()
    name:Union[str, int]
    number:int
    type:str
    beta:float
    material: int
    section:int    
    connectivity: List[int]
    #direction_cosines: List[float]
    #length:float
    #DoF:List[int]
    
    #
    def DoF(self, nodes:Tuple) -> List[int]:
        """
        """
        dof = []
        for node_name in self.connectivity:
            dof.append((nodes[node_name].number) * 6)
        return dof
    #
    def length(self, nodes:Tuple) -> float:
        """
        """
        _node1 = nodes[self.connectivity[0]]
        _node2 = nodes[self.connectivity[-1]]
        #_dx = _node1.x - _node2.x
        #_dy = _node1.y - _node2.y
        #_dz = _node1.z - _node2.z
        #dist2 = (_dx * _dx + _dy * _dy + _dz * _dz)**0.50
        return dist(_node1[:3], _node2[:3])
    #
    def unit_vector(self, nodes) -> List[float]:
        """
        """
        _node1 = nodes[ self.connectivity[ 0 ] ]
        _node2 = nodes[ self.connectivity[-1 ] ]
        dx = _node2.x - _node1.x
        dy = _node2.y - _node1.y
        dz = _node2.z - _node1.z
        # direction cosines
        L = dist(_node1[:3], _node2[:3])
        l = dx / L
        m = dy / L
        n = dz / L
        return [l, m, n]
    #
    def Kmatrix(self, nodes:Dict, materials:Dict, sections:Dict):
        """ """
        section = sections[self.section]
        material = materials[self.material]
        # solve K matrix
        R = Rmatrix(*self.unit_vector(nodes), self.beta)
        #R = Rmatrix(*self.direction_cosines, self.beta)
        #K = beam_stiffness(self.length(nodes),
        #                   section.area, 
        #                   section.J,
        #                   section.Iy,
        #                   section.Iz,
        #                   material.E,
        #                   material.G)        
        K = beam_Ks(self.length(nodes),
                    section.area, section.J, 
                    section.Iy, section.Iz,
                    material.E, material.G,
                    section.area, section.area)
        return trans_3d_beam(K, R)
#
#
def get_elements(bd_file, component_name, arraysize=1000):
    """
    """
    conn = create_connection(bd_file)
    cur = conn.cursor()
    # TODO: check if direction cosines given
    #cur.execute("SELECT tb_elements.name, tb_elements.number, tb_elements.type,\
    #            tb_elements.roll_angle, tb_elements.material, tb_elements.section,\
    #            tb_connectivity.node_1, tb_connectivity.node_2,\
    #            tb_direction_cosines.C11, tb_direction_cosines.C22, tb_direction_cosines.C33\
    #            FROM tb_elements, tb_connectivity, tb_direction_cosines\
    #            WHERE  tb_elements.number = tb_direction_cosines.number\
    #            AND tb_elements.number = tb_connectivity.number;")
    #
    cur.execute("SELECT tb_Elements.name, tb_Elements.number, tb_Elements.type,\
                tb_Elements.roll_angle, tb_Elements.material, tb_Elements.section,\
                tb_Connectivity.node_1, tb_Connectivity.node_2\
                FROM tb_Elements, tb_Connectivity\
                WHERE tb_Elements.number = tb_Connectivity.number;")
    #
    try:
        while True:
            results = cur.fetchmany(arraysize)
            if not results:
                break
            for row in results:
                data = [*row[0:6], [*row[6:8]]]
                yield BeamElement._make(data)
    except sqlite3.Error as e:
        print(e)
    finally:
        conn.close()
    #
    #
    #while True:
    #    row = cur.fetchone()[0]
    #    data = [*row[0:6], [*row[6:8]]]
    #    yield BeamElement._make(data)
    #
    #rows = cur.fetchall()
    #elements = {}
    #for row in rows:
    #    #print(row)
    #    #data = [*row[0:6], [*row[6:8]], [*row[8:]]] # direction cosines
    #    data = [*row[0:6], [*row[6:8]]]
    #    elements[row[0]] = BeamElement._make(data)
    ##conn.close()
    ###print("--->")
    #return elements
#
#
def get_element(bd_file, item_name):
    """
    """
    conn = create_connection(bd_file)
    cur = conn.cursor()
    cur.execute("SELECT tb_Elements.name, tb_Elements.number, tb_Elements.type,\
                tb_Elements.roll_angle, tb_Elements.material, tb_Elements.section,\
                tb_Connectivity.node_1, tb_Connectivity.node_2\
                FROM tb_Elements, tb_Connectivity\
                WHERE tb_Elements.number = tb_Connectivity.number\
                AND tb_Elements.name = {:};".format(item_name))
    row = cur.fetchone()
    data = [*row[0:6], [*row[6:8]]]
    conn.close()
    return BeamElement._make(data)
    #item = cur.execute( "SELECT target FROM stringList WHERE target='specificString';" ).fetchall()
    #return isPresent == None    
    
#
#
def get_connectivities(conn):
    """
    """
    connectivities = {}
    cur = conn.cursor()
    cur = conn.execute("SELECT tb_Elements.name, tb_Elements.number,\
                       tb_Connectivity.node_1, tb_Connectivity.node_2\
                       FROM tb_Elements, tb_Connectivity\
                       WHERE tb_Elements.number = tb_Connectivity.number;")
    #cur.execute("SELECT * FROM tb_Connectivity;")
    rows = cur.fetchall()
    for row in rows:
        #print(row)
        connectivities[row[0]] = [row[2], row[3]]
    return connectivities
#
#
def ResultIter(cursor, arraysize=1000):
    'An iterator that uses fetchmany to keep memory usage down'
    while True:
        results = cursor.fetchmany(arraysize)
        if not results:
            break
        for result in results:
            yield result
#
#
class ElementSQL(Mapping):
    
    __slots__ = ['bd_file', 'arraysize', '_labels', '_conn']
    
    def __init__(self, bd_file) -> None:
        """
        """
        self.bd_file = bd_file
        self.arraysize=1000
        #
        self._conn = create_connection(self.bd_file)
        #conn.close()      
    
    def __setitem__(self, item_name):
        """
        """
        pass
    
    def __getitem__(self, item_name: int) -> Tuple:
        """
        node_name : node number
        """
        #conn = create_connection(self.bd_file)
        cur = self._conn.cursor()
        cur.execute("SELECT tb_Elements.name as ElemName, tb_Elements.number as ElemNumber,\
                    tb_Elements.type, tb_Elements.roll_angle,\
                    tb_Materials.name as MatName, tb_Sections.name as SecName,\
                    tb_Connectivity.node_1, tb_Connectivity.node_2\
                    FROM tb_Elements, tb_Connectivity, tb_Materials, tb_Sections\
                    WHERE tb_Elements.material = tb_Materials.number\
                    AND tb_Elements.section = tb_Sections.number\
                    AND tb_Elements.connectivity = tb_Connectivity.number\
                    AND tb_Elements.name = {:};".format(item_name))
        row = cur.fetchone()
        data = [*row[0:6], [*row[6:]]]
        #conn.close()
        return BeamElement._make(data)
    #
    def __iter__(self):
        """ """
        self._get_elements()
        return iter(self._labels)
    #
    def __len__(self) -> float:
        self._get_elements()
        return len(self._labels)
    #
    def __contains__(self, value) -> bool:
        return value in self._labels    
    #
    @property
    def iter_elements(self):
        """ """
        #conn = create_connection(self.bd_file)
        cur = self._conn.cursor()
        # TODO: check if direction cosines given
        #
        cur.execute("SELECT tb_Elements.name as ElemName, tb_Elements.number as ElemNumber,\
                    tb_Elements.type, tb_Elements.roll_angle,\
                    tb_Materials.name as MatName, tb_Sections.name as SecName,\
                    tb_Connectivity.node_1, tb_Connectivity.node_2\
                    FROM tb_Elements, tb_Connectivity, tb_Materials, tb_Sections\
                    WHERE tb_Elements.material = tb_Materials.number\
                    AND tb_Elements.section = tb_Sections.number\
                    AND tb_Elements.connectivity = tb_Connectivity.number;")
        #
        #try:
        while True:
            results = cur.fetchmany(self.arraysize)
            if not results:
                break
            for row in results:
                data = [*row[0:6], [*row[6:]]]
                yield BeamElement._make(data)
        #except sqlite3.Error as e:
        #    print(e)
        #finally:
        #    conn.close()
    #
    def _get_elements(self):
        """
        """
        cur = self._conn.cursor()
        cur.execute("SELECT tb_Elements.name FROM tb_Elements;")
        rows = cur.fetchall()
        self._labels = list(itertools.chain(*rows))