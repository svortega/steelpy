#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#from typing import Tuple #NamedTuple, 
#

# package imports
from steelpy.f2uModel.sections.process.io_sections import PropertyOut

def get_sections(conn, component_name):
    """
    """
    sections = {}
    cur = conn.cursor()
    cur.execute("SELECT tb_Sections.name, tb_Sections.number, tb_SecProperties.*\
                from tb_Sections, tb_SecProperties\
                WHERE  tb_SecProperties.number = tb_Sections.number;")
    rows = cur.fetchall()
    for row in rows:
        #data = [row[0], *row[4:]]
        #print(row)
        sections[row[0]] = PropertyOut._make(row[3:])
    #conn.close()
    #print("--->")
    return sections