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
#
class BasicLoadBasic(Mapping):
    
    def __init__(self):
        """
        """
        self._labels: array = array("I", [])
        self._title: List[str] = []
        self._number: array = array("I", [])        
        self.gravity = 9.80665 # m/s^2
    #
    def __len__(self) -> int:
        return len(self._labels)
    #
    def __iter__(self):
        """
        """
        return iter(self._labels)
    #
    @property
    def g(self):
        """"""
        return self.gravity
    #
    def __str__(self, units:str="si") -> str:
        """ """
        unit_lenght = " m"
        unit_force = "  N"
        unit_bm = "N*m"
        unit_fl = "N/m"
        #
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += f"{35*' '}BASIC LOAD\n"
        output += "\n"
        output += f"--- Nodal Load\n"
        output += f"Node Number{14*' '}fx [{unit_force}] fy [{unit_force}] fz [{unit_force}] Coord System\n"
        output += f"Load Comment{13*' '}mx [{unit_bm}] my [{unit_bm}] mz [{unit_bm}] Complex\n"
        output += "\n"
        output += f"--- Beam Line Load\n"
        output += f"Element Number{4*' '}L1[{unit_lenght}] qx1[{unit_fl}] qy1[{unit_fl}] qz1[{unit_fl}] Coord System\n"
        output += f"Load Comment{6*' '}L2[{unit_lenght}] qx2[{unit_fl}] qy2[{unit_fl}] qz2[{unit_fl}] Complex\n"
        output += "\n"
        output += f"--- Beam Point Load\n"
        output += f"Element Number{4*' '}L1[{unit_lenght}] fx [{unit_force}] fy [{unit_force}] fz [{unit_force}] Coord System\n"
        output += f"Load Comment{13*' '}mx [{unit_bm}] my [{unit_bm}] mz [{unit_bm}] Complex\n"
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for key in self._labels:
            item = self.__getitem__(key)
            output += f"Load Name : {item.name:12.0f}  Number : {item.number:8.0f}  Title : {item.title}\n"
            # node load
            if item.node:
                output += "--- Nodal Load \n"
                for node_name, points in item.node.items():
                    for point in points:
                        output += point.__str__()
            # beam line
            if item.line_beam:
                output += "--- Beam Line Load\n"
                for beam_name, lines in item.line_beam.items():
                    for line in lines:
                        output += line.__str__()
            # beam point
            if item.point_beam:
                output += "--- Beam Point Load\n"
                for beam_name, points in item.point_beam.items():
                    for point in points:
                        output += point.__str__()
            output += "\n"
        #print('---')
        return output
    #
    #
    def get_basic_load(self, load_name:str, elements, nodes,
                       boundaries, materials, sections):
        """
        """
        member_nload = {}
        nnp = len(nodes)
        bnodes = boundaries.node
        #
        lcase = self.__getitem__(load_name)
        nodal_load = zeros_vector(nnp, 6)
        # nodal load
        items = lcase.point_node.get_group_nodes()
        for lnode in items:
            index = nodes[lnode.number].index
            nodal_load[index] += Vector(lnode[:6])
        # beam line load
        for key, item in lcase.line_beam.items():
            element = elements[key]
            end_nodes = element.connectivity
            n_index0 = nodes[end_nodes[0]].index
            n_index1 = nodes[end_nodes[1]].index
            #
            beam = element.beam()
            for line_load in item:
                bload = line_load.local_system(element.R)
                beam.load.line = [bload[1], bload[2],
                                  bload[7], bload[8],
                                  line_load.L0, line_load.L1]
                #
                #res = line_load.node_equivalent(element, material,
                #                                section, bnodes)
                # global nodal load
                #update_node_force(nodal_load, n_index0, n_index1, res[0])
                # local beam nodal load
                #update_member_force(member_nload, key, res[1])
            #
            boundary = []
            for x, node in enumerate(end_nodes):
                fixity = bnodes[node]
                if fixity:
                    boundary.append(fixity[:6])
                else:
                    boundary.append((1,1,1,1,1,1))
            beam.supports(boundary)
            # get total reactions
            reactions = beam.reactions
            # convert local to global
            gnload, bnload = update_line_load(element, reactions)
            # global nodal load
            update_node_load(nodal_load, n_index0, n_index1, gnload)
            # local beam nodal load
            update_member_load(member_nload, key, bnload, beam)
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
        return BasicLoad(load_name, lcase.number,  lcase.title,
                        nodal_load, member_nload)
        #print('---')
    #
    def get_member_force(self, load_name:str, elements, nodes,
                         reaction, displacement):
        """
        reaction : node load in global system
        ndisp : node displacement in global system
        """
        lcase = self.__getitem__(load_name)
        # beam line load
        beam_res = {}
        for key, item in lcase.line_beam.items():
            element = elements[key]
            end_nodes = element.connectivity
            n_index0 = nodes[end_nodes[0]].index
            n_index1 = nodes[end_nodes[1]].index
            #
            nload = [*reaction[n_index0], *reaction[n_index1]]
            nload = trns_3Dv ( nload, element.R )
            #
            ndisp = displacement[n_index0] + displacement[n_index1]
            ndisp = trns_3Dv ( ndisp, element.R )
            #
            # [FV,FM,Fw,Ftheta]
            bres1 = [nload[1], nload[5], ndisp[1], ndisp[5]]
            bres2 = [nload[2], nload[4], ndisp[2], ndisp[4]]
            #
            beam = element.beam ()
            for line_load in item:
                bload = line_load.local_system(element.R)
                beam.load.line = [bload[1], bload[2],
                                  bload[7], bload[8],
                                  line_load.L0, line_load.L1]
            #
            #beam.reactions = [[bres1, bres2],[[0,0,0,0], [0,0,0,0]]]
            beam_res[key] = beam.response([bres1, bres2])
            #resp
            #1/0
        #print ( '---' )
        return beam_res
#
#
#
def update_node_load(nodal_load:List, index0:int,
                      index1:int, gnload:List):
    """ """
    nodal_load[index0] += Vector(gnload[:6])
    nodal_load[index1] += Vector(gnload[6:])
#
def update_member_load(member_nload:Dict, name:int, nload:List, beam):
    """ """
    try:
        member_nload[name].end_force += nload
    except KeyError:
        member_nload[name] = beam
        member_nload[name].end_force = nload
#
#
def update_line_load(element, res):
    """ """
    lnload = [0] * 12
    # lnload[ 0 ] # axial
    lnload[ 1 ] = res[0][0][0] # y
    lnload[ 5 ] = res[0][0][1] # mz
    lnload[ 2 ] = res[1][0][0] # z
    lnload[ 4 ] = res[1][0][1] # my
    # End 2
    lnload[ 7 ]  = -1*res[0][1][0]  # y
    lnload[ 11 ] = -1*res[0][1][1]  # mz
    lnload[ 8 ]  = -1*res[1][1][0]  # z
    lnload[ 10 ] = -1*res[1][1][1]  # my
    #
    gnload = trns_3Dv(lnload, element.R)
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
    #gnload = Vector(gnload)
    bnload = beam_eq(lnload)
    return gnload, bnload
#
def beam_eq(eq_lnloads):
    """ """
    item = eq_lnloads
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
    return Vector(item)