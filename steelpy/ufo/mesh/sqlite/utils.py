#
# Copyright (c) 2009 steelpy
#
#
# Python stdlib imports
from __future__ import annotations
#from array import array
from operator import itemgetter
#
# package imports
from steelpy.utils.dataframe.main import DBframework
#
#
#
# --------------------------------------------
# Nodes
#
#
def check_nodes(conn, node_name: int|str, component: int):
    """check if node exist"""
    query = (node_name, component, )
    table = f"SELECT * FROM Node \
                WHERE name = ? \
                AND mesh_id= ? ;"
    
    cur = conn.cursor()
    #if isinstance(node_name, str):
    #    cur.execute(table, query)
    #else:
    #    cur.execute(f"SELECT * FROM Node \
    #                  WHERE name = {node_name} ;")
    #
    cur.execute(table, query)
    node = cur.fetchone()
    
    if not node:
        raise IOError(f'Node {node_name} not found')    
    return node    
#
def pull_node_id(conn, name: int|str, mesh_name: str|int):
    """ """
    query = (name, mesh_name, )
    table = f"SELECT Node.number \
                FROM Node, Mesh \
                WHERE Node.name = ? \
                AND Mesh.name = ? \
                AND Node.mesh_id = Mesh.number ;"

    cur = conn.cursor()
    cur.execute(table, query)
    node = cur.fetchone()
    if not node:
        raise IOError(f'Node {name} not found')
    return node[0]
#
def pull_node_mesh(conn, mesh_name: str|int):
    """ """
    query = (mesh_name, )
    table = f"SELECT Node.number, Node.name \
                FROM Node, Mesh \
                WHERE Mesh.name = ? \
                AND Node.mesh_id = Mesh.number ;"

    cur = conn.cursor()
    cur.execute(table, query)
    node = cur.fetchall()
    if not node:
        raise IOError(f'Mesh {mesh_name} not found')
    return node
#
def get_nodes_connec(conn, component: int):
    """ """
    query = ('beam', component)
    table = 'SELECT Element.number, ElementConnectivity.node_end, Node.name \
                FROM ElementConnectivity, Element, Node, Mesh \
                WHERE Element.number = ElementConnectivity.element_id \
                AND Node.number = ElementConnectivity.node_id \
                AND Element.type = ? \
                AND Mesh.number = ? ;'
    #
    cur = conn.cursor()
    cur.execute(table, query)
    connodes = cur.fetchall()
    connodes.sort(key = itemgetter(0))
    #
    1 / 0
    groups = groupby(connodes, itemgetter(0))
    nodes = [[item[1:] for item in data]
             for (key, data) in groups]
    #
    nodes = [[x for _, x in sorted(items)]
             for items in nodes]
    return nodes
#
def get_node_coord(conn, node_id: int):
    """get node coordinates"""
    query = (*node_id, )
    table = f"SELECT x,y,z FROM NodeCoordinate \
              WHERE node_id "
    table += " IN ({:})".format(",".join(["?"] * len(node_id)))
    #
    cur = conn.cursor()
    cur.execute(table, query)
    coord = cur.fetchall()
    return coord
#
# --------------------------------------------
# Element
#
#
def check_element(conn, element_name: int|str, component: int):
    """ """
    query = (element_name, component, )
    table = f"SELECT * FROM Element \
                       WHERE name = ? \
                       AND mesh_id = ? ;"
    cur = conn.cursor()
    cur.execute(table, query)
    row = cur.fetchone()
    if not row:
        raise IOError(f'Element {element_name} not found')    
    return row
#
def pull_element_id(conn, name: int|str,
                    mesh_name: int|str):
    """ """
    query = (name, mesh_name, )
    table = f"SELECT Element.number \
              FROM Element, Mesh \
              WHERE Element.name = ? \
              AND Mesh.name = ? \
              AND Element.mesh_id = Mesh.number ;"
    cur = conn.cursor()
    cur.execute(table, query)
    row = cur.fetchone()
    if not row:
        raise IOError(f'Element {name} not found')
    return row[0]
#
def pull_element_mesh(conn, mesh_name: int|str):
    """ """
    query = (mesh_name, )
    table = f"SELECT Element.number, Element.name\
              FROM Element, Mesh\
              WHERE Element.mesh_id = Mesh.number\
              AND Mesh.name = ?;"
    #
    #
    cur = conn.cursor()
    cur.execute(table, query)
    row = cur.fetchall()
    #
    if not row:
        raise IOError(f'Mesh {mesh_name} not found')
    return row    
#
def get_elements(conn, component: int,
                 element_type: str|int|None = None):
    """ """
    query = (component, )
    table = "SELECT Element.name, Element.number, Element.type,\
                    Material.name, Section.name, \
                    Element.roll_angle, Element.title\
            FROM Element, Material, Section, Mesh\
            WHERE Element.material_id = Material.number \
            AND Element.section_id = Section.number \
            AND Mesh.number = ? "
    #
    if element_type:
        table += 'AND Element.type = ?'
        #query.extend([element_type])
        query = (component, element_type, )
    table += ';'
    #
    cur = conn.cursor()
    cur.execute(table, query)
    rows = cur.fetchall()
    #
    db = DBframework()
    header = ['name', 'number', 'type',
              'material', 'section', 'roll_angle', 'title']
    membdf = db.DataFrame(data=rows, columns=header)
    membdf.set_index('name', inplace=True)
    #
    connodes = get_connectivities(conn, component=component)
    conndf = db.DataFrame(data=connodes, columns=['name', 'nodes', 'end'])
    conndf = conndf.pivot(index='name', columns='end', values='nodes')
    #conndf.reset_index(inplace=True)
    #conndf.set_index('name')
    #conndf.rename(columns={1: 'node_1', 2: 'node_2'}, inplace=True)
    #
    membdf = membdf.join(conndf)
    membdf.reset_index(inplace=True)
    membdf.rename(columns={1: 'node_1', 2: 'node_2'}, inplace=True)
    membdf[['node_3', 'node_4']] = None
    #
    return membdf
#
#
def get_element_data(conn, element_name: str|int,
                     element_type: str, 
                     component: int):
    """ """
    query = (element_name, element_type, component, )
    table = "SELECT Element.name, Element.number, Element.type,\
                    Element.roll_angle, Material.name, \
                    Section.name, Element.title \
            FROM Element, Material, Section, Mesh \
            WHERE Element.name = ? \
            AND Element.type = ? \
            AND Mesh.number = ?  \
            AND Element.material_id = Material.number \
            AND Element.section_id = Section.number \
            AND Element.mesh_id =  Mesh.number ;"
    #
    cur = conn.cursor()
    cur.execute (table, query)
    row = cur.fetchone()
    #
    #element_id = row[1]
    connodes = get_connectivity(conn, element_name,
                                component=component)
    data = [*row[:6], connodes, row[-1]]
    #conn.close ()
    return data
#
def update_element_item(conn, name: str|int, item: str,
                        value, component: int):
    """ """
    query = (value, name, component, )
    table = f'UPDATE Element SET {item} = ? \
             WHERE name = ? \
             AND mesh_id = ? ;'
    cur = conn.cursor()
    cur.execute(table, query)
#
#
def push_connectivity(conn, element_id: int,
                      connectivity: list, component: int):
    """
    """
    items = []
    cur = conn.cursor()
    for x, item in enumerate(connectivity):
        node_id = check_nodes(conn, item,
                              component=component)
        items.append(node_id[0])
        query = (element_id, node_id[0], x+1)
        table = 'INSERT INTO  ElementConnectivity(element_id,\
                                            node_id, node_end)\
                                VALUES(?,?,?)'
        cur.execute(table, query)
    #return cur.lastrowid
    return items
#
def get_connectivity(conn, element_name: int|str,
                     component: int):
    """ """
    query = (element_name, component)
    table = "SELECT ElementConnectivity.node_end, Node.name \
                FROM ElementConnectivity, Node, Element, Mesh \
            WHERE Node.number = ElementConnectivity.node_id\
            AND Element.number = ElementConnectivity.element_id\
            AND Element.mesh_id =  Mesh.number \
            AND Element.name = ? \
            AND Mesh.number = ? ;"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    connodes = cur.fetchall()
    connodes = [x for _, x in sorted(connodes)]
    return connodes
#
def update_connectivity(conn, element_id: int,
                        connectivity: list):
    """
    """
    cur = conn.cursor()
    for x, node in enumerate(connectivity):
        query = (node, element_id, x+1)
        table = ('UPDATE ElementConnectivity SET node_id = ? \
                    WHERE element_id = ? \
                    AND node_end = ?;')
        cur.execute(table, query)
    #return cur.lastrowid
#
def get_connectivities(conn, component: int):
    """
    Return [element_id, node_id, node_end]
    """
    query = (component, )
    table = "SELECT Element.name, Node.name, ElementConnectivity.node_end \
                FROM Element, Node, ElementConnectivity \
             WHERE Element.number = ElementConnectivity.element_id \
             AND Node.number = ElementConnectivity.node_id \
             AND Element.mesh_id = ?;"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    connodes = cur.fetchall()
    #xxx = [x for i, j, x in sorted(connodes)]
    #memb = defaultdict(list)
    #for item in connodes:
    #    memb[item[0]].append()
    #return [x for _, x in sorted(connodes)]
    return connodes
#
def get_unitvector(conn, beam_name: int,
                   component: int):
    """get direction cosines"""
    query = (beam_name, component, )
    table = 'SELECT ElementDirectionCosine.x, \
                    ElementDirectionCosine.y, \
                    ElementDirectionCosine.z \
             FROM Element, ElementDirectionCosine \
             WHERE Element.number = ElementDirectionCosine.element_id \
             AND element.name = ? \
             AND Element.mesh_id = ? \
             ORDER BY ElementDirectionCosine.axis ASC;'
    cur = conn.cursor()
    cur.execute (table, query)
    items = cur.fetchall()
    return items
#
# --------------------------------------------
# Mesh
#
def pull_mesh_id(conn, name:str|int, component: int):
    """ """
    query = (name, component, )
    table = 'SELECT number FROM Mesh \
             WHERE name = ? AND component_id = ?'
    cur = conn.cursor()
    cur.execute (table, query)
    items = cur.fetchone()
    return items[0]
#
#
def pull_results_mesh(conn, mesh_name: int|str):
    """ """
    query = (mesh_name, )
    table = f"SELECT Result.number, Result.name \
              FROM Result, Mesh \
              WHERE Mesh.name = ? \
              AND Result.mesh_id = Mesh.number ;"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    row = cur.fetchall()
    return row
#
#
# --------------------------------------------
# Load
#
def pull_load_id(conn, name: int|str,
                 mesh_name: int|str):
    """ """
    query = (name, mesh_name, )
    table = f"SELECT Load.number \
              FROM Load, Mesh \
              WHERE Load.name = ? \
              AND Mesh.name = ? \
              AND Load.mesh_id = Mesh.number ;"
    cur = conn.cursor()
    cur.execute(table, query)
    row = cur.fetchone()
    if not row:
        raise IOError(f'Load {name} not found')
    return row[0]
#
def pull_load_mesh(conn, mesh_name: int|str):
    """ """
    query = (mesh_name, )
    table = f"SELECT Load.number, Load.name \
              FROM Load, Mesh \
              WHERE Mesh.name = ? \
              AND Load.mesh_id = Mesh.number ;"
    cur = conn.cursor()
    cur.execute(table, query)
    row = cur.fetchall()
    if not row:
        raise IOError(f'Mesh {mesh_name} not found')
    return row