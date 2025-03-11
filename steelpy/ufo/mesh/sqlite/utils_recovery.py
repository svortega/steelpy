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
from steelpy.ufo.utils.node import NodePoint
from steelpy.ufo.utils.boundary import BoundaryItem
from steelpy.utils.dataframe.main import DBframework
#
import numpy as np
#
# --------------------------------------------
# Nodes
#
#
def check_nodes(conn, node_name: int|str, mesh_id: int):
    """check if node exist"""
    query = (node_name, mesh_id, )
    table = f"SELECT * FROM Node \
                WHERE name = ? \
                AND mesh_id= ? ;"
    
    cur = conn.cursor()
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
def get_nodes_connec(conn, mesh_id: int):
    """ """
    query = ('beam', mesh_id)
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
def get_node_coord(conn, node_name: list):
    """get node coordinates"""
    query = (*node_name, )
    table = f"SELECT NodeCoordinate.x, NodeCoordinate.y, NodeCoordinate.z \
              FROM Node, NodeCoordinate \
              WHERE NodeCoordinate.node_id = Node.number \
              AND  Node.name "
    table += " IN ({:})".format(",".join(["?"] * len(node_name)))
    #
    cur = conn.cursor()
    cur.execute(table, query)
    coord = cur.fetchall()
    return coord
#
#
#
def pull_node(conn, node_name:int|str, mesh_id: int, item:str='*'):
    """ """
    data = pull_node_item(conn, node_name, mesh_id, item)
    boundary = pull_node_boundary(conn,
                                  node_name=node_name,
                                  mesh_id=mesh_id)
    #
    node = NodePoint(*data, boundary=boundary)
    return node.system()
#
def pull_node_item(conn, node_name:int|str, mesh_id: int, item:str='*'):
    """ """
    project = (node_name, mesh_id)
    table = f'SELECT NodeCoordinate.{item}, \
                     Node.title, Node.idx \
            FROM Node, NodeCoordinate \
            WHERE Node.number =  NodeCoordinate.node_id\
            AND Node.name = ? \
            AND Node.mesh_id = ?'
    cur = conn.cursor()
    cur.execute(table, project)
    record = cur.fetchone()
    return [*project, *record[1:]]
#
def pull_node_boundary(conn, node_name: int|str,
                       mesh_id: int, item:str="*"):
    """
    """
    #boundary_type = 'restrain'
    data = pull_boundary(conn, mesh_id,
                         node_name,
                         #boundary_type,
                         item=item)
    
    try:
        data = data[0]
        boundary = BoundaryItem(*data[6:12],
                                number=data[5],
                                name=data[2],
                                node=node_name,
                                boundary_type=data[3])
    except IndexError:
        boundary = None
    #
    return boundary
#
def pull_boundary(conn, mesh_id: int,
                  node_name: int|str|None = None,
                  boundary_type: str|None = None, 
                  item:str = "*"):
    """pull all boundary data"""
    #
    project = [mesh_id]
    #
    table = f'SELECT Node.idx, Node.name, \
             Boundary.name, Boundary.type, \
             BoundaryFixity.{item} \
             FROM Node, BoundaryFixity, Boundary, Mesh \
             WHERE Node.boundary_id = Boundary.number \
             AND BoundaryFixity.boundary_id = Boundary.number \
             AND Node.mesh_id = Boundary.mesh_id \
             AND Mesh.number = ?'
    #
    if node_name:
        table += 'AND Node.name = ?'
        project.extend([node_name])
    #
    if boundary_type:
        table += 'AND Boundary.type = ?'
        project.extend([boundary_type])    
    #
    table += ';'
    #
    cur = conn.cursor()
    cur.execute(table, tuple(project))
    data = cur.fetchall()
    return data
#
#
# --------------------------------------------
# Element
#
#
def check_element(conn, element_name: int|str, mesh_id: int):
    """ """
    if isinstance(element_name, (int, np.integer)):
        element_name = int(element_name)

    query = (element_name, mesh_id, )
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
def get_elements(conn, mesh_id: int,
                 element_type: str|int|None = None):
    """ """
    query = (mesh_id, )
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
        query = (mesh_id, element_type, )
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
    connodes = get_connectivities(conn, mesh_id=mesh_id)
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
                     mesh_id: int):
    """ """
    query = (element_name, element_type, mesh_id, )
    table = "SELECT Element.name, Element.number, Element.type,\
                    Material.name, \
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
    1 / 0
    #
    #element_id = row[1]
    connodes = get_connectivity(conn, element_name,
                                mesh_id=mesh_id)
    data = [*row[:5], connodes, row[-1]]
    #conn.close ()
    return data
#
def update_element_item(conn, name: str|int, item: str,
                        value, mesh_id: int):
    """ """
    query = (value, name, mesh_id, )
    table = f'UPDATE Element SET {item} = ? \
             WHERE name = ? \
             AND mesh_id = ? ;'
    cur = conn.cursor()
    cur.execute(table, query)
#
#
def push_connectivity(conn, element_id: int,
                      connectivity: list,
                      eccentricity: list, 
                      mesh_id: int):
    """
    """
    items = []
    cur = conn.cursor()
    for x, item in enumerate(connectivity):
        node_id = check_nodes(conn, item,
                              mesh_id=mesh_id)
        #
        items.append(node_id[0])
        query = (element_id, node_id[0], eccentricity[x], x+1)
        table = 'INSERT INTO  ElementConnectivity(element_id,\
                                                  node_id, \
                                                  eccentricity_id,\
                                                  node_end)\
                                VALUES(?,?,?,?)'
        cur.execute(table, query)
    #return cur.lastrowid
    #return items
#
def get_connectivity(conn, element_name: int|str,
                     mesh_id: int):
    """ """
    query = (element_name, mesh_id)
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
def get_connectivities(conn, mesh_id: int):
    """
    Return [element_id, node_id, node_end]
    """
    query = (mesh_id, )
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
#
def push_dircos(conn, unitvec: list,
                mesh_id: int,
                roll_angle: float|None) -> int:
    """
    directional cosines


    """
    query = (mesh_id, roll_angle)
    table = ('INSERT INTO DirCosine(mesh_id, roll_angle) \
              VALUES(?,?)' )
    #
    cur = conn.cursor()
    cur.execute (table, query)
    dircos_id = cur.lastrowid    
    #
    push_unitvector(conn, dircos_id, unitvec)
    #row = cur.lastrowid
    return dircos_id
#
def push_unitvector(conn, dircos_id: int, 
                    unitvec: list):
    """ """
    query = [(dircos_id, *item, x + 1, )
             for x, item in enumerate(unitvec)]
    
    table = 'INSERT INTO DirCosineUnitVector( \
                        dircosine_id, x, y, z, axis) \
             VALUES(?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.executemany(table, query)
#
#
def get_unitvector(conn, beam_name: int,
                   mesh_id: int):
    """get direction cosines"""
    query = (beam_name, mesh_id, )
    table = 'SELECT DirCosineUnitVector.x, \
                    DirCosineUnitVector.y, \
                    DirCosineUnitVector.z \
             FROM Element, DirCosine, DirCosineUnitVector \
             WHERE Element.dircosine_id = DirCosine.number \
             AND DirCosine.number = DirCosineUnitVector.dircosine_id \
             AND element.name = ? \
             AND Element.mesh_id = ? \
             ORDER BY DirCosineUnitVector.axis ASC;'
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
