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
# Memeber Force
#
def populate_Element_table(conn, load_type, load_number, 
                           element_name, system, index):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    project = (load_number, "NULL", element_name, system, *index)
    if load_type != "basic":
        project = ("NULL", load_number, element_name, system,  *index)
    
    sql = 'INSERT INTO  tb_ResElements (bl_number, lc_number,\
                                        element_name, system,\
                                        point_1,\
                                        point_2, point_3, point_4,\
                                        point_5)\
                                        VALUES(?,?,?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    #return cur.lastrowid
#
def populate_element_force_table(conn, load):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    project = load
    sql = 'INSERT INTO tb_ResElemForce (fx, fy, fz, mx, my, mz)\
                                        VALUES(?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    return cur.lastrowid
#
def populate_element_deflec_table(conn, load_name, deflection):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    #
    project = (deflection)
    #
    sql = 'INSERT INTO tb_ResElemForce (fx, fy, fz, mx, my, mz)\
                                      VALUES(?,?,?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    return cur.lastrowid
#

#
#
#
_table_beam_results = "CREATE TABLE IF NOT EXISTS tb_ResElements (\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                bl_number INTEGER REFERENCES tb_LoadBasic (number),\
                                lc_number INTEGER REFERENCES tb_LoadCombination (number),\
                                element_name INTEGER NOT NULL REFERENCES tb_Elements (name),\
                                system TEXT NOT NULL,\
                                point_1 INTEGER REFERENCES tb_ResElemForce (number),\
                                point_2 INTEGER REFERENCES tb_ResElemForce (number),\
                                point_3 INTEGER REFERENCES tb_ResElemForce (number),\
                                point_4 INTEGER REFERENCES tb_ResElemForce (number),\
                                point_5 INTEGER REFERENCES tb_ResElemForce (number));"
#
#
_table_beam_forces = "CREATE TABLE IF NOT EXISTS tb_ResElemForce (\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                fx DECIMAL NOT NULL,\
                                fy DECIMAL NOT NULL,\
                                fz DECIMAL NOT NULL,\
                                mx DECIMAL NOT NULL,\
                                my DECIMAL NOT NULL,\
                                mz DECIMAL NOT NULL);"
#
#
_table_beam_deflection = "CREATE TABLE IF NOT EXISTS tb_ResElemDeflection (\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                x DECIMAL,\
                                y DECIMAL,\
                                z DECIMAL,\
                                rx DECIMAL,\
                                ry DECIMAL,\
                                rz DECIMAL);"
#

#
def populate_force(conn, mforces):
    """
    """
    #
    #cur = conn.cursor()
    #cur.execute("SELECT tb_Elements.name, tb_Elements.number FROM tb_Elements;")
    #elements = cur.fetchall()
    #elements = {item[0]:item[1] for item in elements}    
    #
    create_table(conn, _table_beam_results)
    create_table(conn, _table_beam_forces)
    create_table(conn, _table_beam_deflection)
    #
    #for load_type, load in load_cases.items():
    for key, forces in mforces.items():
        for mname, memb in forces.items.items():
            system = "local"
            end_force = memb[system][:6]
            end1 = populate_element_force_table(conn, end_force) # local
            end_force = memb[system][6:]
            end2 = populate_element_force_table(conn, end_force) # local
            #for memb in item:
            index = [end1, "NULL", "NULL", "NULL",  end2]
            populate_Element_table(conn, forces.load_type, forces.name, 
                                   mname, system, index)
            #
            system = "global"
            end_force = memb[system][:6]
            end1 = populate_element_force_table(conn, end_force) # local
            end_force =  memb[system][6:]
            end2 = populate_element_force_table(conn, end_force) # local
            #for memb in item:
            index = [end1, "NULL", "NULL", "NULL",  end2]
            populate_Element_table(conn, forces.load_type, forces.name, 
                                   mname, system, index)            
    #print('-->')
    conn.commit()
#
#
# Node Results
#
#
# Nodal Reactions
#
def populate_reacc_table(conn, load_type, load_number, node_name, reacc):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    project = (load_number, "NULL", node_name, *disp)
    if load_type != "basic":
        project = ("NULL", load_number, node_name, *disp)
    #
    sql = 'INSERT INTO tb_ResNodeReactions(bl_number, lc_number, node_name,\
                                                Fx, Fy, Fz, Mx, My, Mz)\
                                                VALUES(?,?,?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
#
_table_reacctions = "CREATE TABLE IF NOT EXISTS tb_ResNodeReactions (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        bl_number INTEGER REFERENCES tb_LoadBasic (number),\
                        lc_number INTEGER REFERENCES tb_LoadCombination (number),\
                        node_name INTEGER NOT NULL REFERENCES tb_Nodes (name),\
                        Fx DECIMAL,\
                        Fy DECIMAL,\
                        Fz DECIMAL,\
                        Mx DECIMAL,\
                        My DECIMAL,\
                        Mz DECIMAL);"
#
def populate_reacctions(conn, nodes, node_reac):
    """
    """
    #TODO: fix this
    #create_table(conn, _table_reacctions)
    for key, disp in node_reac.items():
        for node in nodes.values():
            item = disp.displacement[node.number]
            populate_reacc_table(conn, disp.load_type, disp.name, 
                                 node.name, item)
    #print("-->")
    conn.commit()
#
#
# Nodal displacement
#
def populate_disp_table(conn, load_type, load_number, node_name, disp):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    project = (load_number, "NULL", node_name, *disp)
    if load_type != "basic":
        project = ("NULL", load_number, node_name, *disp)
    #
    sql = 'INSERT INTO tb_ResNodeDisp(bl_number, lc_number, node_name,\
                                    x, y, z, rx, ry, rz)\
                                    VALUES(?,?,?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
#
_table_disp = "CREATE TABLE IF NOT EXISTS tb_ResNodeDisp (\
                number INTEGER PRIMARY KEY NOT NULL,\
                bl_number INTEGER REFERENCES tb_LoadBasic (number),\
                lc_number INTEGER REFERENCES tb_LoadCombination (number),\
                node_name INTEGER NOT NULL REFERENCES tb_Nodes (name),\
                x DECIMAL,\
                y DECIMAL,\
                z DECIMAL,\
                rx DECIMAL,\
                ry DECIMAL,\
                rz DECIMAL);"
#
def populate_displacements(conn, nodes, node_disp):
    """
    """
    create_table(conn, _table_reacctions)
    create_table(conn, _table_disp)
    for key, disp in node_disp.items():
        for node in nodes.values():
            item = disp.items[node.number]
            populate_disp_table(conn, disp.load_type, disp.name, 
                                node.name, item)
    #print("-->")
    conn.commit()
#
#
#