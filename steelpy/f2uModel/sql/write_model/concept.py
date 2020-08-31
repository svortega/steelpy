#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#

# package imports
from steelpy.f2uModel.sql.operation.process_sql import create_table

#
_table_concepts = "CREATE TABLE IF NOT EXISTS tb_Concepts (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name TEXT NOT NULL,\
                    type TEXT NOT NULL,\
                    point_conn INTEGER NOT NULL REFERENCES tb_Points(number),\
                    step_number INTEGER NOT NULL);" # REFERENCES tb_Steps(number)
#
_table_points = "CREATE TABLE IF NOT EXISTS tb_PointConn (\
                number INTEGER PRIMARY KEY NOT NULL,\
                point_1 INTEGER NOT NULL REFERENCES tb_Nodes(name), \
                point_2 INTEGER NOT NULL REFERENCES tb_Nodes(name), \
                point_3 INTEGER REFERENCES tb_Nodes(name), \
                point_4 INTEGER REFERENCES tb_Nodes(name));"
#
_table_steps = "CREATE TABLE IF NOT EXISTS tb_Steps (\
                number INTEGER PRIMARY KEY NOT NULL,\
                concept INTEGER NOT NULL REFERENCES tb_Nodes(name), \
                step_number INTEGER NOT NULL REFERENCES tb_Nodes(name), \
                element INTEGER REFERENCES tb_Nodes(name), \
                material INTEGER REFERENCES tb_Material(name),\
                section INTEGER REFERENCES tb_Sections(name));"
#
#
def populate_points_table(conn, points):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    project = (points[0].name, points[1].name, 'NULL', 'NULL')
    
    sql = 'INSERT INTO  tb_PointConn(point_1, point_2,\
                                    point_3, point_4)\
                                VALUES(?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    return cur.lastrowid
#
def populate_steps_table(conn, concept, step_number, 
                         element, material, section):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    number = step_number + 1
    project = (concept, number, element, material, section)    
    
    sql = 'INSERT INTO  tb_Steps(concept, step_number,\
                                element, material, section)\
                                VALUES(?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    #return cur.lastrowid
#
def populate_concept_table(conn, concept, points):
    """
    Create a new project into the projects table
    :param conn:
    :param project:

    :return: project id
    """
    project = (concept.name, concept.type, 
               points, len(concept.step))
    #
    sql = 'INSERT INTO tb_Concepts(name, type, point_conn, step_number)\
                                   VALUES(?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    #return cur.lastrowid
#
#
def populate_concepts(conn, concepts):
    """
    """
    #
    cur = conn.cursor()
    cur.execute("SELECT tb_Materials.name, tb_Materials.number FROM tb_Materials;")
    material = cur.fetchall()
    material = {item[0]:item[1] for item in material}
    #
    cur = conn.cursor()
    cur.execute("SELECT tb_Sections.name, tb_Sections.number FROM tb_Sections;")
    sections = cur.fetchall()
    sections = {item[0]:item[1] for item in sections}
    #
    create_table(conn, _table_concepts)
    create_table(conn, _table_points)
    create_table(conn, _table_steps)
    #
    for key, item in concepts.items():
        # connectivity
        point_no = populate_points_table(conn, item.connectivity)
        material_no = material[item.material.name]
        section_no = sections[item.section.name]        
        #
        if item.type in ["beam", "truss"]:
            populate_concept_table(conn, item, point_no)
        #
        for number, step in enumerate(item.step):
            populate_steps_table(conn, key, number, step._mesh, 
                                 material_no, section_no)
        #
    #print('---?')
#
#