# 
# Copyright (c) 2009-2020 fem2ufo
#
# Python stdlib imports
from __future__ import annotations
from typing import NamedTuple
import pickle

# package imports
from steelpy.trave3D.preprocessor.assemble import Rmatrix
from steelpy.trave3D.processor.operations import zeros, trns_3Dv, zeros_vector
from steelpy.process.math.vector import Vector
#
# ---------------------
#
def process_load(basic_load:dict, jbc:list, 
                 nodes, elements):
    """
    """
    print("** Processing Load")
    #
    node_load = nodal_load(nodes, basic_load)
    #
    member_load = beam_load(elements, nodes, 
                            basic_load, node_load)
    #
    loads = [node_load[i][j] 
             for i in range(len(jbc)) for j in range(6) 
             if jbc[i][j] != 0]    
    return loads, member_load
#
def nodal_load(nodes, load_case):
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
def beam_load(elements, nodes, load_case, 
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
def beam_udl(eq_lnloads, nodal_load, 
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
# ---------------------
#
class BasicLoad(NamedTuple):
    """ Basic load transfer"""
    name: str
    title: str
    nodal_load: list
    member_load: dict
#
def get_basic_load(elements:dict, nodes:dict, materials:dict, 
                   sections:dict, basic_load:dict):
    """
    """   
    #
    #file = open( "mesh.f2u", "rb" )
    #nodes = pickle.load(file )
    #boundaries = pickle.load(file )
    #elements = pickle.load(file )
    #materials = pickle.load(file )
    #sections = pickle.load(file )
    ##free_nodes = pickle.load(file )
    #file.close()
    #
    nnp = len(nodes)
    #elements = pickle.load(open( "memb.p", "rb" ))
    #
    for load_name, lcase in basic_load.items():
        nodal_load = zeros_vector(nnp, 6)
        items = lcase._nodal_load.get_items_sum
        for lnode in items:
            index = nodes[lnode.name].number
            nodal_load[index] += Vector(lnode[:6])
        #
        eq_lnloads = lcase._beam_line.get_nodal_load(elements, nodes, materials, sections)
        eq_lnloads.extend(lcase._beam_point.get_nodal_load(elements, nodes, materials, sections))
        member_nload = beam_eq(eq_lnloads, #member_nload,
                               nodal_load, elements, nodes)
        #
        #nloads = [nodal_load[i][j] 
        #          for i in range(len(jbc)) for j in range(6) 
        #          if jbc[i][j] != 0]
        yield BasicLoad(load_name, lcase.name, nodal_load, member_nload)
    #
    #return loads, member_nload
#
#
#
def beam_eq(eq_lnloads, # member_nload,
            nodal_load, elements, nodes):
    """ """
    #nmb = len(elements)
    member_nload = {} # zeros(nmb, 12)
    for item in eq_lnloads:
        ename = item[-1]
        element = elements[ename]
        #ename = element.name # element.index
        end_nodes = element.connectivity #element.connectivity
        #
        #
        # global nodal load 
        univec = element.unit_vector(nodes)
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
        n_index = nodes[end_nodes[0]].number
        nodal_load[n_index] += Vector(gnload[:6])            
        # node 2
        n_index = nodes[end_nodes[1]].number
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
#