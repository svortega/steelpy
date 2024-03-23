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
from steelpy.sections.utils.shape.utils import ShapeProperty, get_sect_prop_df
#from steelpy.sections.utils.shape.main import ShapeGeometry #, get_shape
from steelpy.sections.utils.shape.stress import ShapeStressBasic
from steelpy.utils.dataframe.main import DBframework
#
#
#-------------------------------------------------
#
class SectionMainSQL(SectionMain):
    """ """
    __slots__ = ['db_file']
    
    def __init__(self, component:int, db_file: str):
        """
        """
        self.db_file = db_file
        self._component = component
        self._properties = None        
    #
    # -----------------------------------------------
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._component, )
        table = "SELECT name from Section \
                 WHERE component_id = ? ;"
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
        query = (self._component, )
        table = "SELECT SectionGeometry.type \
                 FROM Section, SectionGeometry, Component \
                 WHERE Component.number = ? \
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
        query = (self._component, )
        table = "SELECT title from Section \
                 WHERE component_id =? ; "
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
        #section = (name, self._component, *data[11:])       
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
        query = (name, self._component, )
        table = f'SELECT {item} FROM Section \
                 WHERE number = ? AND component_id = ?;'
        cur = conn.cursor()
        cur.execute(table, query)
        record = cur.fetchone()
        return record[0]
    #    
    #
    def _pull_section(self, conn, section_name:int|str):
        """
        """
        query = (section_name, self._component)
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
        #print("--->")
        #out = [*row[1:3], *row[6:]]
        #return [*row[1:3], *row[6:]]
        return [*row[1:3], *row[11:], *row[3:9]]
    #
    #
    def _push_section(self, conn, section_id: int,
                      properties :tuple|list) -> int:
        """
        name, component_id, title
        """
        query = (section_id, self._component, *properties)
        table = 'INSERT INTO  Section(name, component_id,\
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
        query = (section_id, *properties)
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
        query = (value, section_id, self._component, )
        table = f'UPDATE Section SET {item} = ? \
                  WHERE number = ? \
                  AND component_id = ?;'
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
# -----------------------------------------------
#
class SectionItemSQLXX(SectionMainSQL):
    """ """
    __slots__ = ['db_file', '_component', '_properties']

    def __init__(self, component:int, db_file :str):
        """
        Shear Stress: MAXIMUM / AVERAGE
        section = (name, title, type,
                   diameter, wall_thickess,
                   height, web_thickness,
                   top_flange_width, top_flange_thickness,
                   bottom_flange_width, bottom_flange_thickness,
                   SA_inplane, SA_outplane,
                   shear_stress, build, compactness,)
        """
        super().__init__(db_file)
        #self.db_file = db_file
        self._component = component
        self._properties = None
    #
    # -----------------------------------------------
    #
    def update_itemX(self, item :str, value :float):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            self._update_item(conn, self.number, item, value)
    #
    #
    def get_itemX(self, item :str):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            value = self._get_item(conn, self.number, item)
        return value
    #
    #
    def push_propertyX(self, name):
        """ """
        index = self._labels.index(name)
        number = self._number[index]
        #try:
        properties = self.properties
        self._push_property(number, properties)
        #except:
        #    pass
    #
    def _push_property(self, section_id: int, properties :tuple):
        """ """
        query = (section_id, *properties)
        table = 'INSERT INTO  SectionProperty(section_id, \
                                               area, Zc, Yc,\
                                               Iy, Zey, Zpy, ry,\
                                               Iz, Zez, Zpz, rz,\
                                               J, Cw)\
                                VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
        #
        conn = create_connection(self.db_file)
        with conn:
            #self._push_property_table(conn, number, properties)
            cur = conn.cursor()
            cur.execute(table, query)            
    #
    #def _push_property_table(self, conn, section_id :int, properties :list[float]):
    #    """ """
    #    query = (section_id, *properties)
    #    table = 'INSERT INTO  SectionProperty(component_id, \
    #                                           area, Zc, Yc,\
    #                                           Iy, Zey, Zpy, ry,\
    #                                           Iz, Zez, Zpz, rz,\
    #                                           J, Cw)\
    #                            VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
    #    cur = conn.cursor()
    #    cur.execute(table, query)
    #
    #
    #
    def push_section(self, data):
        """ """
        name = data[0]
        # Section ID
        section = (name, self._component, *data[11:])       
        conn = create_connection(self.db_file)
        with conn:
            number = self._push_section(conn, section)
            # Geometry
            geometry = (number, *data[1:11])
            self._push_geometry(conn, geometry)
        #
        # Property
        item = self.__getitem__(name)
        properties =  item.properties()
        with conn:
            self._push_property(number, properties)
        #
        return number
    #
    #
    def get_sectionX(self, section_name):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            row = self._pull_section(conn, section_name)
        return row[1:]
    #
    #
    # -----------------------------------------------
    #
    @property
    def df(self):
        """ """
        1 / 0
        db = DBframework()
        conn = create_connection(self.db_file)
        #with conn:
        df = db.read_sql_query("SELECT * FROM Section", conn)          
        return df 
    
    @df.setter
    def df(self, df):
        """ """
        1 / 0
        df = get_sect_prop_df(df)
        conn = create_connection(self.db_file)
        with conn:
            df.to_sql('Section', conn,
                         index_label=['name', 'type', 'title',
                                      'diameter', 'wall_thickness',
                                      'height', 'web_thickness',
                                      'top_flange_width', 'top_flange_thickness',
                                      'bottom_flange_width', 'bottom_flange_thickness',
                                      'fillet_radius',
                                      'SA_inplane', 'SA_outplane',
                                      'shear_stress', 'build', 'compactness'], 
                         if_exists='append', index=False)
        #
        #self._labels.extend(df['name'].tolist())
        # TODO : fix numbering
        #nitems = len(self._number) + 1
        #self._number.extend([item + nitems for item in df.index])
        #print('-->')
#
#
# -----------------------------------------------
#
def get_propertyX(conn, component_name):
    """
    """
    sections = {}
    #
    query = (component_name, )
    table = "SELECT Section.name, Section.number, SectionProperty.*\
                from Section, SectionProperty, Component \
                WHERE Component.number = ? \
                AND Section.number = SectionProperty.section_id;"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    rows = cur.fetchall()
    for row in rows:
        # data = [row[0], *row[4:]]
        # print(row)
        sections[row[0]] = ShapeProperty._make(row[3:])
    # conn.close()
    # print("--->")
    1 / 0
    return sections
#
#
# -----------------------------------------------
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
                 component: int):
    """
    """
    query = (section_id, component)
    table = "SELECT * from Section \
             WHERE number = ? AND component_id = ?;"
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
        alpha_sy = self.geometry.alpha_s(poisson=poisson)
        #
        return ShapeProperty(*row[2:15], alpha_sy, alpha_sy)
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
              component: int):
    """ """
    #
    query = (section_name, component, )
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
    shape = get_shape()
#    
# 
def get_sectionXX(conn, section_name: int|str,
                component: int):
    """ """
    1 / 0
    query = (section_name, component)
    table = f"SELECT Section.*, SectionGeometry.* \
              FROM Section, SectionGeometry, Component \
              WHERE Section.name = ? \
              AND Component.number = ? \
              AND Section.number = SectionGeometry.section_id ;"    
    #
    #cur = conn.cursor()
    #try:
    #    query = (section_name, component)
    #    table = "SELECT * from Section \
    #             WHERE Section.name = ? \
    #             AND component_id = ?;"        
    #    cur.execute(table, query)
    #    row = cur.fetchone()
    #
    #except :
    #    query = (component)
    #    table = "SELECT * from Section \
    #             WHERE component_id = ?;"        
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
                             component=row[3])
    #shape = get_shape(shape_type, geometry)
    return shape
#
#