#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
from typing import Tuple #NamedTuple, 
#

# package imports
from steelpy.f2uModel.mesh.node import CoordCylindrical, CoordSpherical, CoordCartesian
from steelpy.f2uModel.mesh.boundary import BoundaryItem

def get_coordinate_system(data:Tuple):
    """
    """
    if 'cylindrical' in data[2].lower():
        return CoordCylindrical(r=data[6], theta=data[7], z=data[5],
                                number=data[1], name=data[0])
    elif 'spherical' in data[2].lower():
        return CoordSpherical(r=data[6], theta=data[7], phi=data[8],
                              number=data[1], name=data[0])
    else:
        return CoordCartesian(x=data[3], y=data[4], z=data[5],
                              number=data[1], name=data[0])
#
#
def get_nodes(conn, component_name):
    """
    """
    nodes = {}
    cur = conn.cursor()
    cur.execute("SELECT * FROM tb_Nodes;")
    rows = cur.fetchall()
    for row in rows:
        #print(row)
        nodes[row[0]] = get_coordinate_system(row)
    #conn.close()
    #print("--->")
    return nodes
#
#
def get_boundaries(conn, component_name):
    """
    """
    boundary = {}
    cur = conn.cursor()
    cur.execute("SELECT * FROM tb_Boundaries;")
    rows = cur.fetchall()
    for row in rows:
        boundary[row[1]] = BoundaryItem(x=row[2], y=row[3], z=row[4],
                                        rx=row[5], ry=row[6], rz=row[7],
                                        number=row[0], name=row[1])
    return boundary