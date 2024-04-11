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
                WHERE name = ? AND level = ? AND component_id = ? ;"
    
    cur = conn.cursor()
    #if isinstance(load_name, str):
    #    cur.execute(f"SELECT * FROM Load \
    #                  WHERE name = '{load_name}' AND level = '{load_level}' ;")
    #else:
    #    cur.execute(f"SELECT * FROM Load \
    #                  WHERE name = {load_name} AND level = '{load_level}' ;")
    #
    cur.execute (table, query)
    loads = cur.fetchone()
    return loads
#
#
