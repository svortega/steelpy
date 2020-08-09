#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#
#

# package imports
from steelpy.f2uModel.sql.operation.process_sql import create_table

#
def populate_node_table(conn, node):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    if node.system == 'cylindrical':
        project = (node.name, node.number, node.system,  
                   'NULL', 'NULL', node.z,  
                   node.r, node.theta, 'NULL')
    elif node.system == 'spherical':
        project = (node.name, node.number, node.system,  
                   'NULL', 'NULL', 'NULL',  
                   node.r, node.theta, node.phi)
    else:
        project = (node.name, node.number, node.system,  
                   node.x.value, node.y.value, node.z.value,  
                  'NULL', 'NULL', 'NULL')
    #
    sql = 'INSERT INTO tb_Nodes (name, number, type,\
                                x, y, z, r,\
                                theta, phi)\
                                VALUES(?,?,?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    #return cur.lastrowid
#
def populate_boundary_table(conn, fixity):
    """
    """
    project = (fixity.name, *fixity[:6])
    
    sql = 'INSERT INTO tb_Boundaries(node_name,\
                                     x, y, z, rx, ry, rz)\
                                     VALUES(?,?,?,?,?,?,?)'
    cur = conn.cursor()
    cur.execute(sql, project)
    #return cur.lastrowid
#
#
#
_table_nodes = "CREATE TABLE IF NOT EXISTS tb_Nodes (\
                name INTEGER PRIMARY KEY NOT NULL,\
                number INTEGER NOT NULL,\
                type TEXT NOT NULL,\
                x DECIMAL,\
                y DECIMAL,\
                z DECIMAL,\
                r DECIMAL,\
                theta DECIMAL,\
                phi DECIMAL);"
#
_table_boundary = "CREATE TABLE IF NOT EXISTS tb_Boundaries(\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    node_name INTEGER NOT NULL REFERENCES tb_Nodes (name),\
                    x DECIMAL, y DECIMAL, z DECIMAL,\
                    rx DECIMAL, ry DECIMAL, rz DECIMAL);"
#
#
def get_tb_nodes(_table_name):
    """
    """
    _table = "CREATE TABLE IF NOT EXISTS {:} (\
                number NOT NULL,\
                x, y, z\
                r, theta, phi);".format(_table_name)
    
    return _table
#
#
def populate_nodes(conn, nodes):
    """
    """
    # nodes
    create_table(conn, _table_nodes)
    for _name, _node in nodes.items():
        populate_node_table(conn, _node)
    conn.commit()
    # boundaries
    #create_table(conn, _table_boundary)
    #for _boundary in mesh.boundaries:
    #    populate_boundary_table(conn, _boundary)
    #
#
#
def populate_boundaries(conn, boundaries):
    """
    """
    create_table(conn, _table_boundary)
    for key, boundary in boundaries.node.items():
        populate_boundary_table(conn, boundary)
    conn.commit()
    #print("--")
#