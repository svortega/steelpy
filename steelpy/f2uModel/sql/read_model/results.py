#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
#from typing import Union, NamedTuple, List
#

# package imports
from steelpy.trave3D.postprocessor.operations import Results

#
#
#
def get_node_displacements(conn):
    """ """
    #
    cur = conn.cursor()
    cur.execute("SELECT Count(*) FROM tb_Nodes")
    n_nodes = cur.fetchone()[0]
    #
    # Basic Load
    #
    cur.execute("SELECT name, number FROM tb_LoadBasic")
    bloads = cur.fetchall()
    #
    basic = {}
    for bl in bloads:
        cur.execute("SELECT tb_LoadBasic.name as LoadName,\
                    tb_Nodes.number as NodeNumber,\
                    tb_ResNodeDisp.*\
                    FROM tb_LoadBasic, tb_Nodes, tb_ResNodeDisp\
                    where tb_Nodes.name = tb_ResNodeDisp.node_name\
                    and tb_LoadBasic.number = tb_ResNodeDisp.bl_number\
                    and tb_LoadBasic.number = {:};".format(bl[1]))
        #
        rows = cur.fetchall()
        node_res = [0] * n_nodes
        for row in rows:
            node_res[row[1]] = row[6:]
        
        basic[bl[1]] = Results(name=bl[1], title=bl[0],
                               load_type="basic", items=node_res)
    #
    #
    # Combinations
    #
    cur.execute("SELECT name, number FROM tb_LoadCombination")
    bloads = cur.fetchall()    
    #
    comb = {}
    for bl in bloads:
        cur.execute("SELECT tb_LoadCombination.name as LoadName,\
                    tb_Nodes.number as NodeNumber,\
                    tb_ResNodeDisp.*\
                    FROM tb_LoadCombination, tb_Nodes, tb_ResNodeDisp\
                    where tb_Nodes.name = tb_ResNodeDisp.node_name\
                    and tb_LoadCombination.number = tb_ResNodeDisp.lc_number\
                    and tb_LoadCombination.number = {:};".format(bl[1]))
        rows = cur.fetchall()
        node_res = [0] * n_nodes
        for row in rows:
            node_res[row[1]] = row[6:]
        #
        comb[bl[1]] = Results(name=bl[1], title=bl[0],
                              load_type="combination", items=node_res)
    #
    displacement = {"basic":basic, "combination":comb}
    return displacement
#
#
#
def get_element_forces(conn):
    """ """
    cur = conn.cursor()
    cur.execute("SELECT name, number FROM tb_LoadBasic")
    loads = cur.fetchall()
    #sections = {item[0]:item[1] for item in sections} 
    cur = conn.cursor()
    #
    memf_basic = {}
    for item in loads:
        basic = {} 
        member_load = get_forces_basic(cur, system="local", load_number=item[1])
        basic = {"local": member_load}
        #
        member_load = get_forces_basic(cur, system="global", load_number=item[1])
        basic.update({"global": member_load})
        #print("--")
        memf_basic[item[1]] = Results(name=item[1], title=item[0],
                                      load_type="basic", items = basic)        
    #
    cur = conn.cursor()
    cur.execute("SELECT name, number FROM tb_LoadCombination")
    loads = cur.fetchall()
    cur = conn.cursor()
    #
    memf_comb = {}
    for item in loads:
        comb = {}
        member_load = get_forces_comb(cur, system="local", load_number=item[1])
        comb = {"local": member_load}
        #
        member_load = get_forces_comb(cur, system="global", load_number=item[1])
        comb.update({"global": member_load})
        #print("--")        
        #
        #
        # memf_basic[bload] = Results()
        memf_comb[item[1]] = Results(name=item[1], title=item[0],
                                     load_type="basic", items = comb)
    #
    #member_load[element.name] = {"global":ngforce, "local":nlforce}
    #print("--")
    return {"basic":memf_basic, "combination":memf_comb}
#
#
def iter_items(rows):
    """ """
    member_load = {}
    for row in rows:
        try:
            member_load[row[0]].extend(list(row[2:]))
        except KeyError:
            member_load[row[0]] = list(row[2:])
    return member_load
    
#
def get_forces_basic(cur, system, load_number):
    """ """
    #member_load = {}
    cur.execute("SELECT tb_Elements.name as ElemName,\
                tb_ResElemForce.*\
                FROM tb_ResElements, tb_ResElemForce, tb_Elements\
                WHERE (tb_ResElements.point_1 = tb_ResElemForce.number\
                       OR tb_ResElements.point_5 = tb_ResElemForce.number)\
                AND tb_ResElements.element_name = tb_Elements.name\
                AND tb_ResElements.system = '{:}'\
                AND tb_ResElements.bl_number = {:};"
                .format(system, load_number))
    rows = cur.fetchall()
    #
    #for row in rows:
    #    try:
    #        member_load[row[0]].extend(list(row[2:]))
    #    except KeyError:
    #        member_load[row[0]] = list(row[2:])
    return iter_items(rows)
#
def get_forces_comb(cur, system, load_number):
    """ """
    #member_load = {}
    cur.execute("SELECT tb_Elements.name as ElemName,\
                tb_ResElemForce.*\
                FROM tb_ResElements, tb_ResElemForce, tb_Elements\
                WHERE (tb_ResElements.point_1 = tb_ResElemForce.number\
                       OR tb_ResElements.point_5 = tb_ResElemForce.number)\
                AND tb_ResElements.element_name = tb_Elements.name\
                AND tb_ResElements.system = '{:}'\
                AND tb_ResElements.lc_number = {:};"
                .format(system, load_number))
    rows = cur.fetchall()
    return iter_items(rows)
#