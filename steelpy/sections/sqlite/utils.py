#
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
#from typing import NamedTuple
#
# package imports
#
from steelpy.sections.utils.shape.main import SectionMain
from steelpy.utils.sqlite.utils import create_connection #, create_table
#from steelpy.sections.utils.operations import get_sect_prop_df
from steelpy.sections.utils.shape.utils import ShapeProperty
from steelpy.sections.utils.shape.stress import ShapeStressBasic
#from steelpy.utils.dataframe.main import DBframework
#
#
#-------------------------------------------------
#
class SectionMainSQL(SectionMain):
    """ """
    __slots__ = ['_tubular', '_solid', '_ibeam', '_box',
                 '_channel', '_tee', '_angle', '_default',
                 'db_file', '_mesh_id']
    
    def __init__(self, mesh_id:int, db_file: str):
        """
        """
        super().__init__()
        #
        self.db_file = db_file
        self._mesh_id = mesh_id        
    #
    # -----------------------------------------------
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._mesh_id, )
        table = "SELECT name from Section \
                 WHERE mesh_id = ? ;"
        #
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]
    #
    @property
    def _type(self):
        """ """
        query = (self._mesh_id, )
        table = "SELECT SectionGeometry.type \
                 FROM Section, SectionGeometry, Mesh \
                 WHERE Mesh.number = ? \
                 AND Section.number = SectionGeometry.section_id ; "
        conn = create_connection(self.db_file)
        #
        with conn:           
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]
    #
    @property
    def _title(self):
        """ """
        query = (self._mesh_id, )
        table = "SELECT title from Section \
                 WHERE mesh_id =? ; "
        #
        conn = create_connection(self.db_file)
        with conn:           
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]    
    #
    # ------------------------------------------------
    # SQL operations
    #
    #
    def push_section(self, data):
        """ """
        name = data[0]
        # Section ID
        #section = (name, self._mesh_id, *data[11:])       
        conn = create_connection(self.db_file)
        with conn:
            number = self._push_section(conn, section_id=name,
                                        properties=data[11:])
            # Geometry
            geometry = (number, *data[1:11])
            self._push_geometry(conn, geometry)
        #
        # Property
        item = self.__getitem__(name)
        properties =  item.properties(poisson=0)
        with conn:
            self._push_property(conn, number, properties)
        #
        return number    
    #    
    #
    def _pull_item(self, conn, name, item):
        """ """
        query = (name, self._mesh_id, )
        table = f'SELECT {item} FROM Section \
                 WHERE number = ? AND mesh_id = ?;'
        cur = conn.cursor()
        cur.execute(table, query)
        record = cur.fetchone()
        return record[0]
    #    
    #
    def _pull_section(self, conn, section_name:int|str):
        """
        """
        query = (section_name, self._mesh_id)
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
        #print("--->")
        #out = [*row[1:3], *row[6:]]
        #return [*row[1:3], *row[6:]]
        return [*row[1:3], *row[11:], *row[3:9]]
    #
    #
    def _push_section(self, conn, section_id: int,
                      properties :tuple|list) -> int:
        """
        name, mesh_id, title
        """
        query = (section_id, self._mesh_id, *properties)
        table = 'INSERT INTO  Section(name, mesh_id,\
                                        SA_inplane, SA_outplane,\
                                        shear_stress, build, compactness,\
                                        title)\
                    VALUES(?,?,?,?,?,?,?,?)'
        #
        cur = conn.cursor()
        cur.execute(table, query)
        return cur.lastrowid
    #    
    def _push_property(self, conn, section_id: int,
                       properties :tuple|list):
        """ """
        query = (section_id, *properties[:15])
        table = 'INSERT INTO  SectionProperty(section_id, \
                                               area, Zc, Yc,\
                                               Iy, Zey, Zpy, ry,\
                                               Iz, Zez, Zpz, rz,\
                                               J, Cw, \
                                               alpha_sy, alpha_sz)\
                                VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
        #
        #conn = create_connection(self.db_file)
        #with conn:
        #self._push_property_table(conn, number, properties)
        cur = conn.cursor()
        cur.execute(table, query)    
    #
    def _push_geometry(self, conn, geometry):
        """ """
        query = geometry
        table = 'INSERT INTO  SectionGeometry(section_id, type,\
                                            diameter, wall_thickness,\
                                            height, web_thickness,\
                                            top_flange_width, top_flange_thickness,\
                                            bottom_flange_width, bottom_flange_thickness,\
                                            fillet_radius)\
                            VALUES(?,?,?,?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(table, query)
    #
    #
    def _update_section(self, conn, section_id: int,
                        item: str, value: float):
        """ """
        query = (value, section_id, self._mesh_id, )
        table = f'UPDATE Section SET {item} = ? \
                  WHERE number = ? \
                  AND mesh_id = ?;'
        cur = conn.cursor()
        cur.execute(table, query)
    #
    #-------------------------------------------------
    #
    def _properties(self):
        """ """
        query = (self.number, )
        table = f'SELECT * FROM SectionProperty \
                 WHERE section_id = ?;'
        #
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            row = cur.fetchone()
        #
        return ShapeProperty(*row[2:])
    #
    def get_section(self, section_name):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            row = self._pull_section(conn, section_name)
        return row    
    #
#
#
#
# -----------------------------------------------
#
#
def get_sections(conn, mesh_id: int) -> list:
    """
    [name, number, mesh_id,
     SA_inplane, SA_outplane,
     shear_stress, build, compactness, title]
    """
    query = (mesh_id, )
    table = "SELECT * from Section \
             WHERE mesh_id = ?;"
    cur = conn.cursor()
    cur.execute(table, query)
    rows = cur.fetchall()
    rows = [[item[1], item[0], *item[3:]]
            for item in rows]    
    return rows
#
def get_properties(conn, mesh_id: int,
                    section_name: int|str|None = None) -> list:
    """
    [section_id, number,
    area, Zc, Yc,
    Iy, Zey, Zpy, ry,
    Iz, Zez, Zpz, rz,
    J, Cw, 
    alpha_sy, alpha_sz
    """
    query = (mesh_id, )
    table = "SELECT Section.name, SectionProperty.*\
             from Section, SectionProperty \
             WHERE Section.mesh_id  = ? \
             AND Section.number = SectionProperty.section_id ;"
    
    if section_name:
        query = (mesh_id, section_name,)
        table.join('AND Section.name = ?')
    
    cur = conn.cursor()
    cur.execute(table, query)
    rows = cur.fetchall()
    rows = [[item[0], item[1], *item[3:]]
            for item in rows]
    return rows
#
def get_geometry(conn, mesh_id: int,
                 section_name: int|str|None = None) -> list:
    """ """
    query = (mesh_id, )
    table = f"SELECT Section.name, SectionGeometry.* \
              FROM Section, SectionGeometry \
              WHERE Section.mesh_id  = ? \
              AND Section.number = SectionGeometry.section_id ;"
    if section_name:
        query = (mesh_id, section_name, )
        table.join('AND Section.name = ?')
    #
    cur = conn.cursor()
    cur.execute(table, query)
    rows = cur.fetchall()
    rows = [[item[0], item[1], *item[3:]]
            for item in rows]
    return rows
    
#
def get_section_df(conn, mesh_id: int) -> list:
    """
    [section_id, type, diameter, wall_thickness,
    height, web_thickness,
    top_flange_width, top_flange_thickness,
    bottom_flange_width, bottom_flange_thickness,
    fillet_radius,
    SA_inplane, SA_outplane,
    shear_stress, build, compactness, title,
    area, Zc, Yc,
    Iy, Zey, Zpy, ry,
    Iz, Zez, Zpz, rz,
    J, Cw, 
    alpha_sy, alpha_sz]
    """
    sect = get_sections(conn, mesh_id)
    sect = {item[0]: item[2:] for item in sect}
    #
    prop = get_properties(conn, mesh_id)
    prop = {item[0]: item[2:] for item in prop}
    #
    geom = get_geometry(conn, mesh_id)
    #
    out = [[item[0], *item[2:], *sect[item[0]], *prop[item[0]]]
           for item in geom]
    #
    return out
#
#
# def get_properties(self):
#    """
#    """
#    summary = {}
#    for key, item in self._sections.items():
#        summary[key] = item._get_properties()
#    return summary
#
#
def get_section2(conn, section_name):
    """ """
    with conn:
        row = _get_section(conn, section_name)
    return row[1:]
#
def _get_section(conn, section_id:int,
                 mesh_id: int):
    """
    """
    query = (section_id, mesh_id)
    table = "SELECT * from Section \
             WHERE number = ? AND mesh_id = ?;"
    cur = conn.cursor()
    cur.execute(table, query)
    row = cur.fetchone()
    #print("--->")
    return row
#
#
@dataclass
class ShapeGeometrySQL(ShapeStressBasic):
    name: str | int
    number: int
    geometry: tuple
    material: tuple
    db_file: str
    #
    def _properties(self, poisson: float):
        """ """
        query = (self.number, )
        table = f'SELECT * FROM SectionProperty \
                 WHERE section_id = ?;'
        #
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            row = cur.fetchone()
        #
        #
        alpha_sy, alpha_sz = self.geometry.alpha_s(poisson=poisson)
        #
        return ShapeProperty(*row[2:15], alpha_sy, alpha_sz)
    #
    #@property
    #def section(self):
    #    """ get section """
    #    #print('--')
    #    return ShapeGeometry(shape_type=self.geometry[2],
    #                         geometry=self.geometry)
    #
    @property
    def _stress(self):
        """ """
        #stress =  self.section._stress
        return self.geometry._stress
#
#
def get_section(conn, section_name: int|str,
              mesh_id: int):
    """ """
    #
    query = (section_name, mesh_id, )
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
    #
    shape = get_shape()
#    
# 
def get_sectionXX(conn, section_name: int|str,
                mesh_id: int):
    """ """
    1 / 0
    query = (section_name, mesh_id)
    table = f"SELECT Section.*, SectionGeometry.* \
              FROM Section, SectionGeometry, Mesh \
              WHERE Section.name = ? \
              AND Mesh.number = ? \
              AND Section.number = SectionGeometry.section_id ;"    
    #
    #cur = conn.cursor()
    #try:
    #    query = (section_name, mesh_id)
    #    table = "SELECT * from Section \
    #             WHERE Section.name = ? \
    #             AND mesh_id = ?;"        
    #    cur.execute(table, query)
    #    row = cur.fetchone()
    #
    #except :
    #    query = (mesh_id)
    #    table = "SELECT * from Section \
    #             WHERE mesh_id = ?;"        
    #    cur.execute(table, query)
    #    rows = cur.fetchall()
    #    for item in rows:
    #        if item[1] ==  section_name:
    #            row = item
    #            break
    #    #1 / 0
    #
    cur = conn.cursor()
    cur.execute(table, query)
    row = cur.fetchone()    
    #
    #geometry = [*row[1:3], *row[6:]]
    geometry = [*row[1:3], *row[11:], *row[3:9]]
    #
    #row = row[1:]
    shape_type = row[11]
    #geometry = row[7:]
    shape = ShapeGeometrySQL(number=row[0],
                             name=row[1],
                             mesh_id=row[3])
    #shape = get_shape(shape_type, geometry)
    return shape
#
#