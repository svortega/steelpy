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
    def get_basic_load(self, elements, nodes, 
                       materials, sections):
        """
        """
        member_nload = {}
        nnp = len(nodes)
        #basic = [item for item in self.__iter__()]
        for load_name in self.__iter__():
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
        #print('---') 
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
#