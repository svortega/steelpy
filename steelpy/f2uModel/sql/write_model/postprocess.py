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




def populate_element_stress_table(conn, force):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    project = (force)
    if load_type != "basic":
       project = (force) 
    
    sql = 'INSERT INTO  tb_PostElemStressIndex (bl_number, lc_number, element_name\
                                            node_1,\
                                            node_2, node_3, node_4,\
                                            node_4)\
                                           VALUES(?,?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    #return cur.lastrowid
#
#
def populate_stress_table(conn, stress_type, load_name, load):
    """
    """
    project = (load.fx, load.fy, load.fz, load.mx, load.my, load.mz)
    #
    sql = 'INSERT INTO {:} (hotspot_1, hotspot_2, hotspot_3,\
                            hotspot_4, hotspot_5, hotspot_6,\
                            hotspot_7, hotspot_8, hotspot_9,\
                            hotspot_10, hotspot_11, hotspot_12,\
                            hotspot_13, hotspot_14, hotspot_15,\
                            hotspot_16, hotspot_17, hotspot_18,\
                            hotspot_19, hotspot_20, hotspot_21,\
                            hotspot_22, hotspot_23, hotspot_24)\
                            VALUES(?,?,?,?,?,?,\
                                   ?,?,?,?,?,?,\
                                   ?,?,?,?,?,?,\
                                   ?,?,?,?,?,?)'.format(stress_type)
    #
    cur = conn.cursor()
    cur.execute(sql, project)
#
#
#
_table_member_stress_index = "CREATE TABLE IF NOT EXISTS tb_ResElemStressIndex (\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                system TEXT,\
                                point_1 DECIMAL,\
                                point_2 DECIMAL,\
                                point_3 DECIMAL,\
                                point_4 DECIMAL,\
                                point_5 DECIMAL,\
                                point_6 DECIMAL,\
                                point_7 DECIMAL,\
                                point_8 DECIMAL,\
                                point_9 DECIMAL,\
                                point_10 DECIMAL,\
                                point_12 DECIMAL);"
#
#
#
hot_spots = "hotspot_1, hotspot_2, hotspot_3,\
            hotspot_4, hotspot_5, hotspot_6,\
            hotspot_7, hotspot_8, hotspot_9,\
            hotspot_10, hotspot_11, hotspot_12,\
            hotspot_13, hotspot_14, hotspot_15,\
            hotspot_16, hotspot_17, hotspot_18,\
            hotspot_19, hotspot_20, hotspot_21,\
            hotspot_22, hotspot_23, hotspot_24"
#
_table_sigma_x = "CREATE TABLE IF NOT EXISTS tb_PostStressSigmaX (\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                {:});".format(hot_spots)
#
_table_sigma_y = "CREATE TABLE IF NOT EXISTS tb_PostStressSigmaY (\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                {:});".format(hot_spots)
#
_table_sigma_z = "CREATE TABLE IF NOT EXISTS tb_PostStressSigmaZ (\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                {:});".format(hot_spots)
#
_table_tau_x = "CREATE TABLE IF NOT EXISTS tb_PostStressTauX (\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                {:});".format(hot_spots)
#
_table_tau_y = "CREATE TABLE IF NOT EXISTS tb_PostStressTauY (\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                {:});".format(hot_spots)
#
_table_tau_z = "CREATE TABLE IF NOT EXISTS tb_PostStressTauZ (\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                {:});".format(hot_spots)
#
hot_spot_items = [_table_sigma_x, _table_sigma_y, _table_sigma_z,
                  _table_tau_x, _table_tau_y, _table_tau_z]
#
#
#
def populate_stress(conn, load_cases):
    """
    """
    #
    #
    create_table(conn, _table_member_stress_index)    
    for item in hot_spot_items:
        create_table(conn, item)
    #
    #
    #for load_type, load in load_cases.items():
    #    for item in load.values():
    #        for x in range(len(item.result[0])):
    #            #for memb in item:
    #            populate_Element_table(conn, load_type, item.name, 
    #                                   item.result[0][x][0], index=x+1)
    #            #force
    print('-->')