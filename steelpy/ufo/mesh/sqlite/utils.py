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
                AND component_id= ? ;"
    
    cur = conn.cursor()
    #if isinstance(node_name, str):
    #    cur.execute(table, query)
    #else:
    #    cur.execute(f"SELECT * FROM Node \
    #                  WHERE name = {node_name} ;")
    #
    cur.execute(table, query)
    node = cur.fetchone()
    return node    
#
#
def get_nodes_connec(conn, component: int):
    """ """
    query = ('beam', component)
    table = 'SELECT Element.number, ElementConnectivity.node_end, Node.name \
                FROM ElementConnectivity, Element, Node, Component \
                WHERE Element.number = ElementConnectivity.element_id \
                AND Node.number = ElementConnectivity.node_id \
                AND Element.type = ? \
                AND Component.number = ? ;'
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
                       AND component_id = ? ;"
    cur = conn.cursor()
    #if isinstance(element_name, str):
    #    cur.execute (f"SELECT * FROM Element \
    #                   WHERE name = '{element_name}';")
    #else:
    #    cur.execute (f"SELECT * FROM Element\
    #                   WHERE name = {element_name};")
    cur.execute(table, query)
    row = cur.fetchone()
    return row
#
#
def get_elements(conn, component: int):
    """ """
    query = (component, )
    table = "SELECT Element.name, Element.type,\
                    Material.name, Section.name, \
                    Element.roll_angle, Element.title\
            FROM Element, Material, Section, Component\
            WHERE Element.material_id = Material.number \
            AND Element.section_id = Section.number \
            AND Component.number = ? ;"
    #
    cur = conn.cursor()
    cur.execute (table, query)
    rows = cur.fetchall()
    #
    db = DBframework()
    header = ['name', 'type', 'material', 'section', 'roll_angle', 'title']
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
            FROM Element, Material, Section, Component \
            WHERE Element.name = ? \
            AND Element.type = ? \
            AND Component.number = ?  \
            AND Element.material_id = Material.number \
            AND Element.section_id = Section.number \
            AND Element.component_id =  Component.number ;"
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
             AND component._number = ? ;'
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
                FROM ElementConnectivity, Node, Element, Component \
            WHERE Node.number = ElementConnectivity.node_id\
            AND Element.number = ElementConnectivity.element_id\
            AND Element.component_id =  Component.number \
            AND Element.name = ? \
            AND Component.number = ? ;"
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
        table = 'UPDATE ElementConnectivity SET node_id = ? \
                    WHERE element_id = ?\
                    AND node_end = ?;'
        cur.execute(table, query)
    #return cur.lastrowid
#
def get_connectivities(conn, component: int):
    """ """
    query = (component, )
    table = "SELECT Element.name, Node.name, ElementConnectivity.node_end \
                FROM Element, Node, ElementConnectivity \
                WHERE Element.number = ElementConnectivity.element_id \
                AND Node.number = ElementConnectivity.node_id \
                AND Element.component_id = ?;"
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
             AND Element.component_id = ? \
             ORDER BY ElementDirectionCosine.axis ASC;'
    cur = conn.cursor()
    cur.execute (table, query)
    items = cur.fetchall()
    return items
#
# --------------------------------------------
#
#
#
#
#
# --------------------------------------------
#
#
#