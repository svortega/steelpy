#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#import os
#import zipfile
#import sqlite3 as sqlite3
from sqlite3 import Error
#from datetime import datetime
#

# package imports
from steelpy.f2uModel.sql.operation.process_sql import create_table

#
def populate_material_table(conn, material):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    project = (material.name, material.type)
    
    sql = 'INSERT INTO  tb_Materials(name, type) VALUES(?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    return cur.lastrowid
#
def populate_linear_material_table(conn, material, number):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    project = (number,
               material.E.convert('pascal').value,
               material.Fy.convert('pascal').value, 
               material.Fu.convert('pascal').value, 
               material.G.convert('pascal').value, 
               material.density.convert('kilogram/metre^3').value, 
               material.poisson, material.alpha.value)
    
    sql = 'INSERT INTO  tb_MatElastoPlastic(material_number,\
                                            E, Fy, Fu, G, poisson, density , alpha)\
                                            VALUES(?,?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
#
def populate_curve_material_table(conn, material, number):
    """
    """
    _table_name = "tb_MatNonLinear"
    _column_disp = 'displacement_' + str(number)
    _column_force = 'force_' + str(number)
    #
    cur = conn.cursor()
    #
    try:
        sql = "ALTER TABLE {:} ADD COLUMN {:} DECIMAL".format(_table_name, _column_disp)
        cur.execute(sql)
        #
        sql = "ALTER TABLE {:} ADD COLUMN {:} DECIMAL".format(_table_name, _column_force)
        cur.execute(sql)
        #
        for x, _spring in enumerate(material._spring):
            cur.execute('UPDATE {:} SET {:}={:}, {:}={:} WHERE rowid = {:}'
                         .format(_table_name,
                                 _column_disp, _spring[0],
                                 _column_force, _spring[1], x+1))
    except Error:
        sql = "CREATE TABLE IF NOT EXISTS tb_MatNonLinear({:} DECIMAL,\
               {:} DECIMAL);".format(_column_disp, _column_force)
        cur.execute(sql)
        #
        sql = 'INSERT INTO {:} VALUES(?,?)'.format(_table_name, _column_disp, _column_force)
        cur.executemany(sql, material._spring)
    #
    #
    #
    #
    #cur = conn.cursor()
    #cur.execute(sql, project)
#
#
_table_materials = "CREATE TABLE IF NOT EXISTS tb_Materials (\
                    number INTEGER PRIMARY KEY,\
                    name INTEGER NOT NULL,\
                    type TEXT NOT NULL);"
#
_table_linear_material = "CREATE TABLE IF NOT EXISTS tb_MatElastoPlastic(\
                            number INTEGER PRIMARY KEY,\
                            material_number INTEGER NOT NULL,\
                            E DECIMAL NOT NULL,\
                            Fy DECIMAL NOT NULL,\
                            Fu DECIMAL NOT NULL,\
                            G DECIMAL NOT NULL,\
                            poisson DECIMAL NOT NULL,\
                            density DECIMAL NOT NULL, \
                            alpha DECIMAL NOT NULL);"
#
_table_curve_material = "CREATE TABLE IF NOT EXISTS tb_MatNonLinear(\
                         number INTEGER NOT NULL PRIMARY KEY,\
                         material_number INTEGER NOT NULL);"
#
def populate_materials(conn, materials):
    """
    """
    # materials
    create_table(conn, _table_materials)
    create_table(conn, _table_linear_material)
    create_table(conn, _table_curve_material)
    for key, _material in materials.items():
        number = populate_material_table(conn, _material)
        if _material.type == 'elastic' :
            populate_linear_material_table(conn, _material, number)
        elif _material.type == 'curve':
            populate_curve_material_table(conn, _material, number)
        else:
            print('mat to be fixed')
    conn.commit()
    #print('-->')
#