#
# Copyright (c) 2009-2020 fem2ufo
#

# Python stdlib imports
from typing import Dict # NamedTuple, 
from dataclasses import dataclass
#

# package imports
from steelpy.f2uModel.load.node import PointNode
from steelpy.f2uModel.load.element import PointBeam, LineBeam


@dataclass
class LoadTypes:
    """ """
    __slots__ = ['_nodal_load', '_nodal_mass', '_nodal_displacement',
                 '_beam_line', '_beam_point','name', 'number']
    
    #_nodal_load:Dict
    #_nodal_mass:Dict
    #_nodal_displacement:Dict
    #_beam_line:Dict
    #_beam_point:Dict
    name:str
    number:int


def get_basic_load(conn, component_name):
    """
    """
    cur = conn.cursor()
    cur.execute("SELECT name, number FROM tb_LoadBasic")
    loads = cur.fetchall()
    #
    basic_load = {}
    cur = conn.cursor()
    for item in loads:
        basic_load[item[1]] = LoadTypes(name=item[0], number=item[1])
        # nodal load
        cur.execute("SELECT tb_LoadNode.*\
                     FROM tb_LoadBasic, tb_LoadNode, tb_Nodes\
                     WHERE tb_LoadBasic.number = tb_LoadNode.load_number\
                     AND tb_Nodes.name = tb_LoadNode.node_name\
                     AND tb_LoadBasic.number = {:};".format(item[1])) 
        rows = cur.fetchall()
        nodal_load = {}
        for row in rows:
            #print(row)
            data = [*row[5:11], row[2], row[2], *row[3:5]]
            try:
                nodal_load[row[2]].append(PointNode._make(data))
            except KeyError:
                nodal_load[row[2]] = [PointNode._make(data)]
        basic_load[item[1]]._nodal_load = nodal_load
        #
        # beam load
        cur.execute("SELECT tb_LoadBeamPoint.*\
                     FROM tb_LoadBasic, tb_LoadBeamPoint, tb_Elements\
                     WHERE tb_LoadBasic.number = tb_LoadBeamPoint.load_number\
                     AND tb_Elements.name = tb_LoadBeamPoint.element_name\
                     AND tb_LoadBasic.number = {:};".format(item[1])) 
        rows = cur.fetchall()
        beam_point = {}
        for row in rows:
            #print(row)
            data = [*row[6:12], row[5], row[2], row[2], *row[3:5]]
            try:
                beam_point[row[2]].append(PointBeam._make(data))
            except KeyError:
                beam_point[row[2]] = [PointBeam._make(data)]
        basic_load[item[1]]._beam_point = beam_point
        #
        # beam concentrated load
        cur.execute("SELECT tb_LoadBeamLine.*\
                    FROM tb_LoadBasic, tb_LoadBeamLine, tb_Elements\
                    WHERE tb_LoadBasic.number = tb_LoadBeamLine.load_number\
                    AND tb_Elements.name = tb_LoadBeamLine.element_name\
                    AND tb_LoadBasic.number = {:};".format(item[1])) 
        rows = cur.fetchall()
        beam_line = {}
        for row in rows:
            print(row)
            data = [*row[6:9], *row[13:16], row[5], row[12],
                    row[2], row[2], *row[3:5]]
            try:
                beam_line[row[2]].append(LineBeam._make(data))
            except KeyError:
                beam_line[row[2]] = [LineBeam._make(data)]
        basic_load[item[1]]._beam_line = beam_line        
    #
    #print("-->")
    return basic_load