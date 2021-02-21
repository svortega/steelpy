#
# Copyright (c) 2009-2021 fem2ufo
# 

# Python stdlib imports
from array import array
from collections.abc import Mapping
#from dataclasses import dataclass
from typing import NamedTuple, Tuple, List, Iterator, Dict, Iterable, ClassVar, Union
#import re

# package imports
# package imports
from steelpy.trave3D.preprocessor.assemble import Rmatrix
from steelpy.trave3D.processor.operations import zeros, trns_3Dv, zeros_vector
from steelpy.process.math.vector import Vector
#
#
# ---------------------------------
#
#
class BasicLoad(NamedTuple):
    """ Basic load transfer"""
    name: int
    number:int
    title: str
    nodal_load: list
    member_load: dict
#
def get_basic_load(basic_load, elements, nodes, 
                   materials, sections):
    """
    """
    member_nload = {}
    nnp = len(nodes)
    for load_name, lcase in basic_load.items():
        nodal_load = zeros_vector(nnp, 6)
        # nodal load
        items = lcase.point_node.get_group_nodes()
        for lnode in items:
            index = nodes[lnode.number].index
            nodal_load[index] += Vector(lnode[:6])
        # beam line load
        for key, item in lcase.line_beam.items():
            element = elements[key]
            material = materials[element.material]
            section = sections[element.section].properties
            #
            end_nodes = element.connectivity
            n_index0 = nodes[end_nodes[0]].index
            n_index1 = nodes[end_nodes[1]].index            
            for line_load in item:
                res = line_load.node_equivalent(element, material, section)
                # global nodal load
                update_node_force(nodal_load, n_index0, n_index1, res[0])
                # local beam nodal load
                update_member_force(member_nload, key, res[1])
        # beam point load
        for key, item in lcase.point_beam.items():
            element = elements[key]
            material = materials[element.material]
            section = sections[element.section].properties
            #
            end_nodes = element.connectivity
            n_index0 = nodes[end_nodes[0]].index
            n_index1 = nodes[end_nodes[1]].index            
            for point_load in item:
                res = point_load.node_equivalent(element, material, section)
                # global nodal load
                update_node_force(nodal_load, n_index0, n_index1, res[0])
                # local beam nodal load
                update_member_force(member_nload, key, res[1])
        #
        yield BasicLoad(load_name, lcase.number,  lcase.title,
                        nodal_load, member_nload)
#
#
def update_node_force(nodal_load:List, index0:int, 
                      index1:int, gnload:List):
    """ """
    nodal_load[index0] += Vector(gnload[:6])
    nodal_load[index1] += Vector(gnload[6:])
#
def update_member_force(member_nload:Dict, name:int, nload:List):
    """ """
    try:
        member_nload[name] += nload
    except KeyError:
        member_nload[name] = nload
#
#
def beam_eqX(eq_lnloads, nodal_load, elements, nodes):
    """ """
    #nmb = len(elements)
    member_nload = {} # zeros(nmb, 12)
    for item in eq_lnloads:
        ename = item[-1]
        element = elements[ename]
        #ename = element.name # element.index
        end_nodes = element.connectivity #element.connectivity
        #
        # global nodal load
        univec = element.unit_vector
        R = Rmatrix(*univec, element.beta)
        gnload = trns_3Dv(item[:12], R)
        #
        # axial
        gnload[0] *= -1
        gnload[6] *= -1
        #
        # rotation X
        #gnload[3] *= -1
        #gnload[9] *= -1
        #
        # rotation Y
        gnload[4] *= -1
        gnload[10] *= -1
        #
        # rotation Z
        gnload[5] *= -1
        gnload[11] *= -1
        #
        # node 1
        n_index = nodes[end_nodes[0]].index
        nodal_load[n_index] += Vector(gnload[:6])
        # node 2
        n_index = nodes[end_nodes[1]].index
        nodal_load[n_index] += Vector(gnload[6:])
        #
        # axial
        item[0] *= -1
        item[6] *= -1
        #
        # rotation X
        #item[3] *= -1
        #item[9] *= -1
        #
        # rotation Y
        item[4] *= -1
        item[10] *= -1
        #
        item[5] *= -1
        item[11] *= -1
        #
        try:
            #member_nload[ename] = [member_nload[ename][x] + item[x]
            #                       if x == 0 or x == 6
            #                       else member_nload[ename][x] - item[x]
            #                       for x in range(12)]
            #member_nload[ename] += Vector([item[x] if x == 0 or x == 6
            #                               else -item[x] for x in range(12)])
            member_nload[ename] += Vector(item[:12])
        except KeyError:
            #member_nload[ename] = Vector([item[x] if x == 4 or x == 10
            #                              else - item[x] for x in range(12)])
            member_nload[ename] = Vector(item[:12])
            #member_nload[ename] = Vector([_item*-1 for _item in item[:12]])
    return member_nload
#
def nodal_loadX(nodes, load_case):
    """
    """
    nnp = len(nodes)
    #nodal_load = zeros(nnp, 6)
    nodal_load = zeros_vector(nnp, 6)
    for load_name, lcase in load_case.items():
        items = lcase._nodal_load.get_items_sum
        for node in items:
            index = nodes._labels.index(node.name)
            nodal_load[index] += Vector(node[:6])
    return nodal_load
#
def beam_loadX(elements, nodes, load_case,
              nodal_load):
    """
    equivalent local nodal load
    """
    eq_lnloads = []
    for lcase in load_case.values():
        eq_lnloads.extend(lcase._beam_udl.get_nodal_load(elements))
        eq_lnloads.extend(lcase._beam_point.get_nodal_load(elements))
    member_nload = beam_udl(eq_lnloads, nodal_load,
                            elements, nodes)
    return member_nload
#
def beam_udlX(eq_lnloads, nodal_load,
             elements, nodes):
    """ """
    nmb = len(elements)
    member_nload = zeros(nmb, 12)
    for item in eq_lnloads:
        element = elements[item[-1]]
        index = element.index
        end_nodes = element.connectivity
        #
        # update member nodal load based on local udl
        #member_nload[index] = [member_nload[index][x] - item[x]
        #                       for x in range(12)]
        #
        member_nload[index] = [member_nload[index][x] + item[x]
                               if x == 0 or x == 6
                               else member_nload[index][x] - item[x]
                               for x in range(12)]
        #
        # global nodal load
        #gload = [-it for it in item[:12]]
        gnload = trns_3Dv(item[:12], element.R)
        gnload[0] *= -1
        gnload[6] *= -1
        #
        gnload[4] *= -1
        gnload[10] *= -1
        # node 1
        n_index = nodes._labels.index(end_nodes[0])
        nodal_load[n_index] += Vector(gnload[:6])
        #nodal_load[n_index] = [nodal_load[n_index][x] + _nload
        #                       for x, _nload in enumerate(gnload[:6])]
        # node 2
        n_index = nodes._labels.index(end_nodes[1])
        nodal_load[n_index] += Vector(gnload[6:])
        #nodal_load[n_index] = [nodal_load[n_index][x] + _nload
        #                       for x, _nload in enumerate(gnload[6:])]
    return member_nload
#
#
#