#
# Copyright (c) 2009 steelpy
# 
# Python stdlib imports

#
# package imports
#

#
#
def get_load_data(conn, load_name:int|str, load_level: str,
                  component: int):
    """ """
    query = (load_name, load_level, component)
    table = f"SELECT * FROM Load \
                WHERE name = ? \
                AND level = ? \
                AND mesh_id = ? ;"
    
    cur = conn.cursor()
    cur.execute (table, query)
    loads = cur.fetchone()
    return loads
#
#
#
def push_basic(conn, load_name: str|int,
               component: int, 
               design_load: str,
               #title: str|None = None,
               step_type: str|None = None):
    """ """
    load = pull_load(conn, load_name,
                     load_level='basic',
                     component=component)
    try:
        load_id = load[0]
    except IndexError:
        raise IOError(f'load {load_name} not found')
    
    query = (load_id, step_type, design_load)
    table = 'INSERT INTO LoadBasic(\
                load_id, step_type, design_load)\
                VALUES(?,?,?)'
    with conn:
        cur = conn.cursor()
        cur.execute(table, query)
        idx = cur.lastrowid
        #
        #if not idx:
        #    query = (load_id, load_type, )
        #    table = 'SELECT number\
        #             FROM LoadBasic \
        #             WHERE load_id = ? \
        #             AND load_type = ? '
        #    cur.execute(table, query)
        #    idx = cur.fetchone()
        #    idx = idx[0]
    #if not idx:
    #    idx = pull_basic(conn, load_id, load_type)
    return idx
#
def pull_basic(conn, load_name: str|int):
    """ """
    query = (load_name,)
    table = 'SELECT LoadBasic.* \
             FROM Load, LoadBasic \
             WHERE LoadBasic.load_id = Load.number \
             AND Load.name = ?'
    with conn:
        cur = conn.cursor()
        cur.execute(table, query)
        idx = cur.fetchone()
        #idx = idx[0]
    return idx
#
def pull_load(conn, load_name: str|int,
              load_level: str,
              component: int,
              item: str = '*'):
    """ """
    query = (load_name, load_level, component)
    table = f"SELECT {item} \
              FROM Load \
              WHERE name = ?\
              AND level = ? \
              AND mesh_id = ? ;"    
    with conn:
        cur = conn.cursor()
        cur.execute(table, query)
        idx = cur.fetchone()
        #idx = idx[0]
    return idx
#
#
#
def get_load_basics(conn, component: int):
    """ """
    query = (component, )
    #with conn:
    cur = conn.cursor()
    #
    table = "SELECT Node.name, Node.number FROM Node \
            WHERE Node.mesh_id = ? ;"
    cur.execute(table, query)
    nodes = cur.fetchall()
    nodes = {item[0]:item[1] for item in nodes}
    #
    table = "SELECT Element.name, Element.number FROM Element \
             WHERE Element.mesh_id = ? ;"
    cur.execute(table, query)
    elements = cur.fetchall()
    elements = {item[0]:item[1] for item in elements}            
    #
    query = ('basic', component, )
    table = "SELECT Load.name, Load.number FROM Load \
             WHERE Load.level = ? \
             AND Load.mesh_id = ? ;"
    cur.execute(table, query)
    basic = cur.fetchall()
    basic = {item[0]:item[1] for item in basic}
    #
    return nodes, elements, basic
#
#
def pull_basic_load_id(conn, load_name:int|str):
    """ """
    cur = conn.cursor()
    cur.execute("SELECT * FROM Load\
                 WHERE name = {:} \
                 AND type = 'basic'".format(load_name))
    loads = cur.fetchone()
    return loads[0]
#