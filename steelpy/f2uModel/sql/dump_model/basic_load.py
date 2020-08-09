#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#

# package imports
from steelpy.f2uModel.sql.operation.process_sql import create_table


#
def populate_master_load(conn, load_number, load, load_type):
    """
    """
    project = (load_number, load.name, load_type)
    sql = 'INSERT INTO tb_LoadBasic(number, name, type) \
           VALUES(?,?,?)'
    cur = conn.cursor()
    cur.execute(sql, project)
    #return cur.lastrowid
#
def populate_node_load_table(conn, load_type, load_number, load_name, node):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    #
    system = "global"
    if node.system != 0:
       system = "local" 
    #
    project = (load_number,
               node.name, system, node.load_complex,
               node.fx, node.fy, node.fz, node.mx, node.my, node.mz,
               'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL')
    #
    sql = 'INSERT INTO tb_LoadNode(load_number,\
                                    node_name, system, load_complex, \
                                    fx, fy, fz, mx, my, mz,\
                                    fxi, fyi, fzi, mxi, myi, mzi)\
                                VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
#
def populate_beam_line_load_table(conn, load_type, load_number, load_name, memb):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    system = "global"
    if memb.system != 0:
       system = "local" 
    #
    project = (load_number, 
               memb.name, system, memb.load_complex,
               memb.L1, memb.qx1, memb.qy1, memb.qz1,
               'NULL', 'NULL', 'NULL',
               memb.L2, memb.qx2, memb.qy2, memb.qz2,
               'NULL', 'NULL', 'NULL')
    #
    sql = 'INSERT INTO tb_LoadBeamLine(load_number,\
                                         element_name, system, load_complex,\
                                         Lnode1, qx1, qy1, qz1, qx1i, qy1i, qz1i,\
                                         Lnode2, qx2, qy2, qz2, qx2i, qy2i, qz2i)\
                                         VALUES(?,?,?,?,\
                                                ?,?,?,?,?,?,?,\
                                                ?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
#
def populate_beam_point_load_table(conn, load_type, load_number, load_name, memb):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    system = "global"
    if memb.system != 0:
       system = "local" 
    #    
    project = (load_number,
               memb.name, system, memb.load_complex, 
               memb.distance, 
               memb.fx, memb.fy, memb.fz,
               memb.mx, memb.my, memb.mz,
               'NULL', 'NULL', 'NULL',
               'NULL', 'NULL', 'NULL')
    #
    sql = 'INSERT INTO tb_LoadBeamPoint(load_number,\
                                          element_name, system, load_complex,\
                                          Lnode1, fx, fy, fz, mx, my, mz,\
                                          fxi, fyi, fzi, mxi, myi, mzi)\
                                         VALUES(?,?,?,?,\
                                                ?,?,?,?,?,?,?,\
                                                ?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
#
#
#
_table_load_index = "CREATE TABLE IF NOT EXISTS tb_LoadBasic (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name TEXT NOT NULL,\
                    type TEXT NOT NULL);"
#
_table_node_load = "CREATE TABLE IF NOT EXISTS tb_LoadNode(\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    load_number INTEGER NOT NULL REFERENCES tb_LoadBasic (number),\
                    node_name INTEGER NOT NULL REFERENCES tb_Nodes (name),\
                    system TEXT,\
                    load_complex INTEGER NOT NULL,\
                    fx DECIMAL,\
                    fy DECIMAL,\
                    fz DECIMAL,\
                    mx DECIMAL,\
                    my DECIMAL,\
                    mz DECIMAL,\
                    fxi DECIMAL,\
                    fyi DECIMAL,\
                    fzi DECIMAL,\
                    mxi DECIMAL,\
                    myi DECIMAL,\
                    mzi DECIMAL);"

#
_table_element_line_load = "CREATE TABLE IF NOT EXISTS tb_LoadBeamLine(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            load_number INTEGER NOT NULL REFERENCES tb_LoadBasic (number),\
                            element_name INTEGER NOT NULL REFERENCES tb_Elements (name),\
                            system TEXT,\
                            load_complex INTEGER NOT NULL,\
                            Lnode1 DECIMAL,\
                            qx1 DECIMAL,\
                            qy1 DECIMAL,\
                            qz1 DECIMAL,\
                            qx1i DECIMAL,\
                            qy1i DECIMAL,\
                            qz1i DECIMAL,\
                            Lnode2 DECIMAL,\
                            qx2 DECIMAL,\
                            qy2 DECIMAL,\
                            qz2 DECIMAL,\
                            qx2i DECIMAL,\
                            qy2i DECIMAL,\
                            qz2i DECIMAL);"
#
_table_element_point_load = "CREATE TABLE IF NOT EXISTS tb_LoadBeamPoint(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            load_number INTEGER NOT NULL REFERENCES tb_LoadBasic (number),\
                            element_name INTEGER NOT NULL REFERENCES tb_Elements (name),\
                            system TEXT,\
                            load_complex INTEGER NOT NULL,\
                            Lnode1 DECIMAL,\
                            fx DECIMAL,\
                            fy DECIMAL,\
                            fz DECIMAL,\
                            mx DECIMAL,\
                            my DECIMAL,\
                            mz DECIMAL,\
                            fxi DECIMAL,\
                            fyi DECIMAL,\
                            fzi DECIMAL,\
                            mxi DECIMAL,\
                            myi DECIMAL,\
                            mzi DECIMAL);"
#


def populate_basic_loading(conn, basic_load):
    """
    """
    #
    #load_data = loading.data
    #
    load_type = "basic"
    #
    # master load
    create_table(conn, _table_load_index)
    create_table(conn, _table_node_load)
    create_table(conn, _table_element_line_load)
    create_table(conn, _table_element_point_load)
    #
    for load_name, load_basic in basic_load.items():
        populate_master_load(conn, load_name, load_basic, load_type)
    #
    # node
    #create_table(conn, _table_node_load)
    #for load_name, load_basic in basic_load.items():
        for basic in load_basic.point_node.values():
            for node in basic:
                populate_node_load_table(conn, load_type, 
                                         load_name, 
                                         load_basic.name, node)
    #
    # beam element
    #create_table(conn, _table_element_line_load)
    #create_table(conn, _table_element_point_load)
    #for load_name, load_basic in basic_load.items():
        for basic in load_basic.udl_beam.values():
            for element in basic:
                populate_beam_line_load_table(conn, load_type, 
                                              load_name, 
                                              load_basic.name, element)
        for basic in load_basic.point_beam.values():
            for element in basic:
                populate_beam_point_load_table(conn, load_type, 
                                               load_name, 
                                               load_basic.name, element)
    #
    conn.commit()
    #print('-->')