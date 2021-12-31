# 
# Copyright (c) 2019-2021 steelpy
#

# Python stdlib imports
#from dataclasses import dataclass
#import math
from typing import List
#

# package imports
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table



#
class SectionSQLite:
    """ """
    __slots__ = ['db_file']
    
    def __init__(self, db_file:str, section:tuple):
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
        self.db_file = db_file
        self._properties = None
        # create table
        conn = create_connection(self.db_file)
        with conn:        
            self.number = self._push_section_table(conn, section)     
    #
    #
    def update_item(self, item:str, value:float):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            self._update_item(conn, self.number, item, value)
            #conn.commit()
    #
    def _update_item(self, conn, name, item, value):
        """ """
        project = (value, name)
        sql = 'UPDATE tb_Sections SET {:} = ? WHERE number = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    #
    def get_item(self, item:str):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            value = self._get_item(conn, self.number, item)
        return value
    #
    def _get_item(self, conn, name, item):
        """ """
        project = (name,)
        sql = 'SELECT {:} FROM tb_Sections WHERE number = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
        record = cur.fetchone()
        return record[0]
    #
    #
    def push_property(self):
        """ """
        try:
            properties = self.properties
            self._push_property(properties)
        except:
            pass
    #
    def _push_property(self, properties:tuple):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            self._push_property_table(conn, self.number, properties)        
    #
    def _push_property_table(self, conn, section_number:int, properties:List[float]):
        """ """
        project = (section_number, *properties)
        sql = 'INSERT INTO  tb_SecProperties(number, area, Zc, Yc,\
                                            Iy, Zey, Zpy, ry,\
                                            Iz, Zez, Zpz, rz,\
                                            J, Cw)\
                                            VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
    #
    def _push_section_table(self, conn, section) -> int:
        """
        """
        project = section
        sql = 'INSERT INTO  tb_Sections(name, title, type,\
                                        diameter, wall_thickess,\
                                        height, web_thickness,\
                                        top_flange_width, top_flange_thickness,\
                                        bottom_flange_width, bottom_flange_thickness,\
                                        SA_inplane, SA_outplane,\
                                        shear_stress, build, compactness)\
                                        VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(sql, project)
        return cur.lastrowid
    #
    #
    #
#