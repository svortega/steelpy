#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#import os
#import zipfile
#import sqlite3 as sqlite3
#from sqlite3 import Error
#from datetime import datetime
#

# package imports
from steelpy.f2uModel.sql.operation.process_sql import create_table


#
def populate_section_table(conn, _section):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    if 'tubular' in _section.type.lower():
        project = (_section.name, _section.type,
                   _section.diameter.value, 
                   _section.thickness.value,
                   'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL',)
    if 'rectangular' in _section.type.lower():
        project = ( _section.name, _section.type,
                   'NULL', 'NULL',
                   _section.height.value, 'NULL', 
                   _section.width.value, 'NULL', 
                   _section.width.value, 'NULL',)
    else:
        project = (_section.name, _section.type,
                   'NULL', 'NULL',
                   _section.height.value, _section.web_thickness.value,
                   _section.top_flange_width.value, _section.top_flange_thickness.value,
                   _section.bottom_flange_width.value, _section.bottom_flange_thickness.value)        
    #
    #
    sql = 'INSERT INTO  tb_Sections(name, type, diameter, wall_thickess,\
                                    height, web_thickness,\
                                    top_flange_width, top_flange_thickness,\
                                    bottom_flange_width, bottom_flange_thickness)\
                                    VALUES(?,?,?,?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    return cur.lastrowid
#
#
def populate_property_table(conn, number, section):
    """ """
    project = (number, *section.properties)
    
    sql = 'INSERT INTO  tb_SecProperties(number, area, Zc, Yc,\
                                        Iy, Zey, Zpy, ry,\
                                        Iz, Zez, Zpz, rz,\
                                        J, Cw)\
                                        VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
#
#
#
table_sections =  "CREATE TABLE IF NOT EXISTS tb_Sections (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name INTEGER NOT NULL,\
                    type TEXT NOT NULL,\
                    diameter DECIMAL,\
                    wall_thickess DECIMAL,\
                    height DECIMAL,\
                    web_thickness DECIMAL,\
                    top_flange_width DECIMAL,\
                    top_flange_thickness DECIMAL,\
                    bottom_flange_width DECIMAL,\
                    bottom_flange_thickness DECIMAL);"
#
table_properties = "CREATE TABLE IF NOT EXISTS tb_SecProperties (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    area DECIMAL,\
                    Zc DECIMAL,\
                    Yc DECIMAL,\
                    Iy DECIMAL,\
                    Zey DECIMAL,\
                    Zpy DECIMAL,\
                    ry DECIMAL,\
                    Iz DECIMAL,\
                    Zez DECIMAL,\
                    Zpz DECIMAL,\
                    rz DECIMAL,\
                    J DECIMAL,\
                    Cw DECIMAL);"
#
def populate_sections(conn, sections):
    """
    """
    # sections
    create_table(conn, table_sections)
    create_table(conn, table_properties)
    for key, section in sections.items():
        number = populate_section_table(conn, section)
        populate_property_table(conn, number, section)
    conn.commit()
    #
