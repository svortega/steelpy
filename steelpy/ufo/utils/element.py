#
# Copyright (c) 2009 steelpy
# 
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
from collections.abc import Mapping
#from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
#from typing import NamedTuple
import re
#from math import isclose
#import os
#
# package imports
import steelpy.utils.io_module.text as common
from steelpy.utils.dataframe.main import DBframework
#
#

class ElementMain(Mapping):
    """ """
    __slots__ = ['_name']

    def __init__(self, name:int|str) -> None:
        """
        Manages f2u elements
        """
        self._name = name
    #
    def __len__(self) -> int:
        return len(self._labels)

    def __iter__(self):
        return iter(self._labels)

    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    # ---------------------------------
    #
    def __str__(self, units:str="si") -> str:
        """ """
        lenght = ' m'
        space = " "
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += f"{33*space}ELEMENTS\n"
        output += "\n"
        output += (f"Element     Node1    Node2 {4*space} Material  Section {4*space}")
        output += (f"Beta {3*space} Lenght {2*space} Title")
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for key in self._labels:
            element = self.__getitem__(key)
            output += element.__str__()
            #output += "\n"
        return output
    #
    # ---------------------------------
    #
    #@property
    def _beam(self, values:None|list=None,
              df=None):
        """
        """
        beam_items = []
        if values:
            1 / 0
            if isinstance(values, list):
                for item in values:
                    beam_name = item[0]
                    try:
                        self._labels.index(beam_name)
                        raise Exception('element {:} already exist'.format(beam_name))
                    except ValueError:
                        #element_type = 'beam'
                        # default
                        #self._labels.append(element_name)
                        #self._type.append(element_type)
                        self._beams[beam_name] = item[1:]
                        beam_items.append(beam_name)
        #
        # dataframe input
        try:
            columns = list(df.columns)
            1 / 0
            self._beams.df = df
        except AttributeError:
            pass
        #
        return beam_items
    #
    # ---------------------------------
    #
    def max_bandwidth(self,  jbc):
        """
        calculate max bandwidth
        ------------------------  
        npi : connectivity end 1
        npj : connectivity end 2
        jbc : nodes freedom
        nel: number of elements
        if we
        npi ,npj, jbc, nel
        """
        #TODO : plates 
        ibndm4 = [0]
        for key in self.self._labels:
            element = self.__getitem__(key)
            nodes = element.connectivity
            #
            bc1 = jbc.loc[nodes[0]].tolist()
            bc2 = jbc.loc[nodes[1]].tolist()
            ieqn = bc1 + bc2
            try:
                ibndm4.append(max([abs(ieqn1 - ieqn2)
                                   for x, ieqn1 in enumerate(ieqn) if ieqn1 > 0
                                   for ieqn2 in ieqn[x+1:] if ieqn2 > 0]))
            except ValueError:
                continue
        #
        return max(ibndm4) + 1
    #
    #
    def get_number(self, start: int = 1):
        """
        """
        try:
            n = max(self._labels) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1    
    #
    # ---------------------------------
    #
    @property
    def df(self):
        """nodes in dataframe format"""
        # TODO : additional elements tobe added
        bdf = self._beams.df
        return bdf

    @df.setter
    def df(self, df):
        """nodes in dataframe format"""
        #df = get_element_df(df)
        columns = list(df.columns)
        for key in columns:
            if re.match(r"\b((element(s)?(_|-|\s*)?)?type)\b", key, re.IGNORECASE):
                df['type'] = df[key].apply(lambda x: find_element_type(x))
                break
        #
        group = df.groupby("type")
        #
        try:
            group = group.get_group("beam")
            self._beams.df = group
        except KeyError:
            raise IOError('Element df not valid')
#
#
#
def get_element(parameters:tuple|list|dict)->list:
    """ """
    # [element_type, material]
    element_type = None 
    if isinstance(parameters, dict):
        matkeys = list(parameters.keys())
        for key in matkeys:
            if re.match(r"\b((element(s)?(_|-|\s*)?)?type)\b", key, re.IGNORECASE):
                element_type = find_element_type(parameters.pop(key)) # parameters[key])
                #parameters.pop(key)
                break   
    elif isinstance(parameters, (list, tuple)):
        parameters = list(parameters)
        element_type = parameters.pop(0)
        element_type = find_element_type(element_type)
    else:
        raise IOError(f' element input not valid')
    #
    if  element_type:
        match element_type:
            case "beam":
                data = get_beam(parameters)
            case "shell":
                raise NotImplementedError
    else:
        raise IOError(f'Element type not given')
    
    return [element_type, *data]
#
#
def find_element_type(word_in:str) -> str:
    """
    """
    key = {"beam": r"\b(beam|viga)\b",
           "shell": r"\b(shell(s)?|plate(s)?)\b"}
    #try:
    match = common.find_keyword(word_in, key)
    #except IOError:
    #    match = get_node(word_in)
    return match
#
def find_element_item(word_in:str) -> str:
    """ """
    key = {"name": r"\b(id|name|element(s)?|member(s)?)\b",
           "material": r"\b(material(s)?((_|-|\s*)?name)?)\b",
           "section": r"\b(section(s)?((_|-|\s*)?name)?)\b",
           "title": r"\b(title)\b",
           "type": r"\b((element(s)?(_|-|\s*)?)?type)\b"}
    try:
        match = common.find_keyword(word_in, key)
    except IOError:
        match = get_node(word_in)
    return match
#
def find_beam_item(word_in:str) -> str:
    """ """
    key = {"name": r"\b(id|name|element(s)?|member(s)?)\b",
           "material": r"\b(material(s)?((_|-|\s*)?name)?)\b",
           "section": r"\b(section(s)?((_|-|\s*)?name)?)\b",
           "roll_angle": r"\b(alpha|roll(_|-|\s*)?angle)\b",
           #"title": r"\b(title)\b",
           "type": r"\b((element(s)?(_|-|\s*)?)?type)\b"}
    try:
        match = common.find_keyword(word_in, key)
    except IOError:
        try:
            match = get_node(word_in)
        except IOError:
            match = get_coord(word_in)
    return match
#
def get_beam(parameters:tuple|list|dict) -> list:
    """ [node1, node2, material, section, roll_angle,
         title, concept_idx] """
    if isinstance(parameters, dict):
        #keys = list(parameters.keys())
        blist = get_beam_dict(parameters)
    elif isinstance(parameters, (list, tuple)):
        blist = get_beam_list(parameters)
    #
    if not all(blist[:4]):
        raise IOError('beam basic input missing')
    #
    return blist
#
def get_beam_dict(parameters:dict) -> list:
    """ [node1, node2, material, section, roll_angle,
         title, concept_idx] """
    blist = [None] * 7
    blist[4] = 0.0 # roll angle
    if isinstance(parameters, dict):
        #keys = list(parameters.keys())
        for key, item in parameters.items():
            if re.match(r"\b((nod(e|o)|end|coord(inate)?)(_|-|\s*)?1)\b", key, re.IGNORECASE):
                blist[0] = item
            elif re.match(r"\b((nod(e|o)|end|coord(inate)?)(_|-|\s*)?2)\b", key, re.IGNORECASE):
                blist[1] = item
            elif re.match(r"\b(material(s)?((_|-|\s*)?name)?)\b", key, re.IGNORECASE):
                blist[2] = item
            elif re.match(r"\b(section(s)?((_|-|\s*)?name)?)\b", key, re.IGNORECASE):
                blist[3] = item
            elif re.match(r"\b(alpha|roll(_|-|\s*)?angle)\b", key, re.IGNORECASE):
                blist[4] = item            
            elif re.match(r"\b(title)\b", key, re.IGNORECASE):
                blist[5] = item
    #
    return blist
#
def get_beam_list(parameters:tuple|list,
                  steps:int = 7) -> list:
    """ [node1, node2, material, section, roll_angle,
         title, concept_idx] """
    blist = [None] * steps
    blist[4] = 0.0 # roll angle
    # basic data
    for x in range(steps):
        try:
            blist[x] = parameters[x]
        except IndexError:
            pass
    #
    #for x in range(4, len(parameters)):
    #    
    #
    #1 / 0
    return blist
#
def get_node(word_in:str) -> str:
    """ """
    key = {"node1": r"\b((nod(e|o)|end|point)(_|-|\s*)?1)\b",
           "node2": r"\b((nod(e|o)|end|point)(_|-|\s*)?2)\b", 
           "node3": r"\b((nod(e|o)|end|point)(_|-|\s*)?3)\b", 
           "node4": r"\b((nod(e|o)|end|point)(_|-|\s*)?4)\b"}
    match = common.find_keyword(word_in, key)
    return match      
#
def get_coord(word_in:str) -> str:
    """ """
    key = {"x1": r"\b((coord(inate)?)?(_|-|\s*)?x(_|-|\s*)?(0|1|start))\b",
           "y1": r"\b((coord(inate)?)?(_|-|\s*)?y(_|-|\s*)?(0|1|start))\b",
           "z1": r"\b((coord(inate)?)?(_|-|\s*)?z(_|-|\s*)?(0|1|start))\b",
           #
           "x2": r"\b((coord(inate)?)?(_|-|\s*)?x(_|-|\s*)?(2|end))\b",
           "y2": r"\b((coord(inate)?)?(_|-|\s*)?y(_|-|\s*)?(2|end))\b",
           "z2": r"\b((coord(inate)?)?(_|-|\s*)?z(_|-|\s*)?(2|end))\b"}    
    match = common.find_keyword(word_in, key)
    return match
#
def get_beam_df(df: DBframework.DataFrame):
    """ [node1, node2, material, section, roll_angle,
         title, concept_idx] """
    columns = list(df.columns)
    header = {item:find_beam_item(item) for item in columns}
    df.rename(columns=header, inplace=True)
    #
    if "roll_angle" not in columns:
        df["roll_angle"] = 0.0
    #
    if "title" not in columns:
        df["title"] = None
    #
    if "concept_idx" not in columns:
        df["concept_idx"] = None
    #
    return df
#
#
#@functools.lru_cache(maxsize=2048)
def find_element_data(word_in: str) -> str:
    """
    Identify beam data from user
    """
    key: dict = {"number": r"\b(number|mesh)\b",
                  "name": r"\b(name|label)\b",
                  "connectivity": r"\b(connectivity|node(s)?|joint(s)?)\b",
                  "material": r"\b(material)\b",
                  "section": r"\b(geometry|section|shape)\b",
                  "hinges": r"\b(hinge(s)?)\b",
                  "group": r"\b(group|set)(s)?\b",
                  # "boundary": r"\b(boundar(y|ies))\b",
                  "type": r"type"}

    match = common.find_keyword(word_in, key)
    return match
#