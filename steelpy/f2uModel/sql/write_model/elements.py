#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#

# package imports
from steelpy.f2uModel.sql.operation.process_sql import create_table

#
def populate_beam_element_table(conn, element, material,
                                section, univector, 
                                coonectivity, offset):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    #univector = element.number
    #try:
    #    univector = element.direction_cosines.number
    #except ValueError:
    #    univector = 'NULL'
    #    
    #try:
    #    offset = element.offset_index + 1
    #except ValueError:
    #    offset = 'NULL'
    #
    #
    project = (element.name, element.number, element.type,
               material, section, 
               coonectivity, element.beta,
               univector, offset)
    #
    sql = 'INSERT INTO tb_Elements(name, number, type, material, section,\
                                   connectivity, roll_angle,\
                                   direction_cosines, eccentricities)\
                                   VALUES(?,?,?,?,?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    #return cur.lastrowid
#
def populate_coonectivity_table(conn, connectivity):
    """
    Create a new project into the projects table
    :param conn:
    :param project:
    
    :return: project id
    """
    project = (connectivity[0],
               connectivity[1], 'NULL', 'NULL')    
    
    sql = 'INSERT INTO  tb_Connectivity(node_1, node_2,\
                                        node_3, node_4)\
                                        VALUES(?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    return cur.lastrowid
#
def populate_univector_table(conn, univectors):
    """
    """
    #project = (component, number, 1, 
    #           univectors[0][0], univectors[0][1], univectors[0][2],
    #           univectors[1][0], univectors[1][1], univectors[1][2],
    #           univectors[2][0], univectors[2][1], univectors[2][2])    
    # TODO: need fix
    #if univectors.type == 3:
    #    project = (component, univectors.number,
    #               univectors.type, 
    #               univectors.x, 'NULL', 'NULL',
    #               'NULL', univectors.y, 'NULL',
    #               'NULL', 'NULL', univectors.z)
    #else:
    #    project = (univectors.number, univectors.type, *univectors[:6])
    #
    project = (1, 
               univectors[0], 'NULL', 'NULL',
               'NULL', univectors[1], 'NULL',
               'NULL', 'NULL', univectors[2])    
    
    sql = 'INSERT INTO tb_DirectionCosines(type, \
                                            C11, C12, C13, C21, C22, C23, C31, C32, C33)\
                                            VALUES(?,?,?,?,?,?,?,?,?,?)'
    cur = conn.cursor()
    cur.execute(sql, project)
    return cur.lastrowid
#
def populate_offset_table(conn, offset):
    """
    """
    project = (offset.system, *offset[:3])
    
    sql = 'INSERT INTO tb_Eccentricities(system, x, y, z)\
                                         VALUES(?,?,?,?)'
    cur = conn.cursor()
    cur.execute(sql, project)
    return cur.lastrowid
#
def populate_offset_index_table(conn, element):
    """
    """
    #_index = element.offset_index
    _offsets = element.eccentricities
    _node_number =  len(_offsets)
    #project = [_index + 1, _node_number]
    project = [_node_number]
    
    for x in range(4):
        try:
            project.append(_offsets[x].number)
        except IndexError:
            project.append('NULL')
    
    sql = 'INSERT INTO  tb_EccIndex(ecc_number,\
                                    ecc_1, ecc_2, \
                                    ecc_3, ecc_4) VALUES(?,?,?,?,?)'
    #
    cur = conn.cursor()
    cur.execute(sql, project)
    return cur.lastrowid
#
#
#
#
_table_beam_elements = "CREATE TABLE IF NOT EXISTS tb_Elements (\
                        name INTEGER PRIMARY KEY NOT NULL,\
                        number INTEGER NOT NULL,\
                        type TEXT NOT NULL,\
                        material INTEGER NOT NULL REFERENCES tb_Materials (number),\
                        section INTEGER NOT NULL REFERENCES tb_Sections (number),\
                        connectivity INTEGER NOT NULL REFERENCES tb_Connectivity (number),\
                        roll_angle DECIMAL NOT NULL,\
                        direction_cosines INTEGER REFERENCES tb_DirectionCosines (number),\
                        eccentricities INTEGER REFERENCES tb_Eccentricities (number));" 
#
_table_connectivity = "CREATE TABLE IF NOT EXISTS tb_Connectivity (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        node_1 DECIMAL NOT NULL REFERENCES tb_Nodes (name), \
                        node_2 DECIMAL NOT NULL REFERENCES tb_Nodes (name), \
                        node_3 DECIMAL REFERENCES tb_Nodes (name), \
                        node_4 DECIMAL REFERENCES tb_Nodes (name));"
#
_table_univectors = "CREATE TABLE IF NOT EXISTS tb_DirectionCosines (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    type TEXT NOT NULL,\
                    C11 DECIMAL, C12 DECIMAL, C13 DECIMAL,\
                    C21 DECIMAL, C22 DECIMAL, C23 DECIMAL,\
                    C31 DECIMAL, C32 DECIMAL, C33 DECIMAL);"
#
_table_offset = "CREATE TABLE IF NOT EXISTS tb_Eccentricities (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        system TEXT NOT NULL,\
                        x DECIMAL, \
                        y DECIMAL, \
                        z DECIMAL);"
#
_table_offset_index = "CREATE TABLE IF NOT EXISTS tb_EccIndex (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        ecc_number INTEGER NOT NULL,\
                        ecc_1 DECIMAL NOT NULL REFERENCES tb_Eccentricities (number), \
                        ecc_2 DECIMAL NOT NULL REFERENCES tb_Eccentricities (number), \
                        ecc_3 DECIMAL REFERENCES tb_Eccentricities (number), \
                        ecc_4 DECIMAL REFERENCES tb_Eccentricities (number));"
#
#
#
def populate_elements(conn, elements):
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
    # Elements
    create_table(conn, _table_beam_elements)
    create_table(conn, _table_connectivity)
    create_table(conn, _table_offset)
    create_table(conn, _table_offset_index)
    create_table(conn, _table_univectors)
    #
    for _name, _element in elements.items():
        # TODO: offsets
        offset_number = 'NULL'
        # univectors
        univec_no = 'NULL'
        #univec_no = populate_univector_table(conn, _element.number, 
        #                                     _element.unit_vector)
        # connectivity
        coonec_no = populate_coonectivity_table(conn, _element.connectivity)
        #
        material_no = material[_element.material]
        section_no = sections[_element.section]
        #
        if _element.type in ["beam", "truss"]:
            populate_beam_element_table(conn, _element, material_no, 
                                        section_no, univec_no, 
                                        coonec_no, offset_number)
        #
        #
        #
        # offsets
        #create_table(conn, _table_offset)
        #for _offset in mesh.eccentricities.values():
        #    populate_offset_table(conn, _offset, component)
        #
        # offset index
        #create_table(conn, _table_offset_index)
        #for _element in mesh.elements:
        #    try:
        #        _element.offset_index
        #        populate_offset_index_table(conn, _element)
        #    except ValueError:
        #        continue
    #
    conn.commit()
    #
#
