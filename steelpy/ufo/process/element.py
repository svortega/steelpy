#
# Copyright (c) 2009 steelpy
# 
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
#from collections.abc import Mapping
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
           #"Tee section": r"\b(t(ee)?((_|-|\s*)?section)?)\b",
           #"Tubular section": r"\b(tub(ular)?|pipe|chs((_|-|\s*)?section)?)\b",
           "shell": r"\b(shell(s)?|plate(s)?)\b",}
    match = common.find_keyword(word_in, key)
    return match
#
def find_beam_item(word_in:str) -> str:
    """ """
    key = {"name": r"\b(id|name|element(s)?|member(s)?)\b",
           "node1": r"\b((nod(e|o)|end|coord(inate)?)(_|-|\s*)?1)\b",
           "node2": r"\b((nod(e|o)|end|coord(inate)?)(_|-|\s*)?2)\b", 
           "material": r"\b(material(s)?((_|-|\s*)?name)?)\b",
           "section": r"\b(section(s)?((_|-|\s*)?name)?)\b",
           "roll_angle": r"\b(alpha|roll(_|-|\s*)?angle)\b",
           "title": r"\b(title)\b",
           "type": r"\b((element(s)?(_|-|\s*)?)?type)\b"}
    match = common.find_keyword(word_in, key)
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
    