#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#

# package imports
from steelpy.f2uModel.sql.operation.process_sql import create_table

#
#
def populate_master_load(conn, load, load_type):
    """
    """
    project = (load.number, load.name, "combination", load.level)
    sql = 'INSERT INTO tb_LoadCombination (number, name, type, level) \
           VALUES(?,?,?,?)'
    cur = conn.cursor()
    cur.execute(sql, project)
    #return cur.lastrowid
#
def populate_lc_lv1_tb(conn, comb_number,
                       bl_number, factor):
    """
    """
    project = (comb_number, bl_number, "NULL", factor)
    
    sql = 'INSERT INTO tb_LoadCombIndex (\
           comb_number, bl_number, lc_number, factor)\
           VALUES(?,?,?,?)'
    cur = conn.cursor()
    cur.execute(sql, project)
#
#
def populate_lc_lv2_tb(conn, comb_number,
                       lc_number, factor):
    """
    """
    project = (comb_number, "NULL",  lc_number, factor)
    
    sql = 'INSERT INTO tb_LoadCombIndex (\
           comb_number, bl_number, lc_number, factor)\
           VALUES(?,?,?.?)'
    cur = conn.cursor()
    cur.execute(sql, project)
#
#
_table_load_index = "CREATE TABLE IF NOT EXISTS tb_LoadCombination (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name TEXT  NOT NULL,\
                    type TEXT  NOT NULL,\
                    level INTEGER NOT NULL);"
#
_table_comb_load = "CREATE TABLE IF NOT EXISTS tb_LoadCombIndex (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    comb_number INTEGER NOT NULL,\
                    bl_number INTEGER,\
                    lc_number INTEGER,\
                    factor DECIMAL NOT NULL);"
#
#
#
def populate_load_combination(conn, combined_load):
    """
    """
    #
    #load_data = loading.data
    #
    #
    # master load
    create_table(conn, _table_load_index)
    create_table(conn, _table_comb_load)
    for lc_name, load_comb in combined_load.items():
        level = load_comb.level
        populate_master_load(conn, load_comb, level)  
        for key, item in load_comb.basic_load.items():
            populate_lc_lv1_tb(conn, lc_name, key, item)
    #
    conn.commit()
    #print('-->')    