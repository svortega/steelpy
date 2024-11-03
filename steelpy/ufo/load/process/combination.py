#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from array import array
from collections.abc import Mapping
#from collections import defaultdict
#from typing import NamedTuple
import re

# package imports
import steelpy.utils.io_module.text as common
from steelpy.utils.dataframe.main import DBframework

#
class LoadCombinationBasic(Mapping):
    __slots__ = ['_labels', '_title', '_number',
                 '_index', '_basic', '_combination', '_metocean']
    
    def __init__(self):
        """
        """
        self._labels: list[str|int] = []
        self._title: list[str] = []
        self._number: array = array("I", [])
    #
    def __len__(self) -> int:
        return len(self._labels)
    #
    def __iter__(self):
        """
        """
        #comb =  list(dict.fromkeys(self._labels))
        #return iter(comb)
        return iter(self._labels)
    #
    def get_number(self, start:int=1):
        """
        """
        try:
            n = len(self._labels) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1    
    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += "{:}LOAD COMBINATIONS\n".format(30*" ")
        #output += "--- Basic Load \n"
        #output += f"Load Type [Basic/Combination]\n"
        output += f"Load Type   Name{' '*10} Factor\n"
        #output += "\n"
        #output += "\n"
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for key in self._labels:
            lcase = self.__getitem__(key)
            output += f"Load Name : {str(key):12s}  Number : {lcase.number:8.0f}  Title : {lcase.title}\n"
            try:
                #output += f"--- Basic\n"
                for basic_name, factor in lcase.basic.items():
                    output += f"Basic       {str(basic_name):12s} {factor: 1.4e}\n"
                    #name, basic
            except TypeError:
                pass
            # comb
            try:
                #output += f"--- Combination\n"
                for comb_name, factor in lcase.combination.items():
                    output += f"Combination {str(comb_name):12s} {factor: 1.4e}\n"
            except TypeError:
                continue
            #name, comb
            output += "\n"
        #print('---')
        return output
    #
    def to_basic(self):
        """ """
        # get combination of combination and convert them to basic loads
        comb = {}
        for key, item in self._combination.items():
            cbasic = {}
            # basic load
            for comb_name, factor in item._basic.items():
                cbasic[comb_name] = factor
            # comb of comb
            for comb_name, factor in item._combination.items():
                for precomb, prefactor in self._combination[comb_name]._basic.items():
                    try:
                        cbasic[precomb] += prefactor * factor
                    except KeyError:
                        cbasic[precomb] = prefactor * factor
            comb[key] = cbasic
        #
        # Combination formed by basic loads only
        dftemp = []
        for key, item in self._combination.items():
            #mesh_name = item.mesh_name
            for bl_name, factor in comb[key].items():
                dftemp.append([key, item.number, 'combination',
                               item.title, item.mesh_name,
                               bl_name, factor])
        #
        header = ['load_name', 'load_id','load_type',
                  'load_title', 'mesh_name', 
                  'basic_load', 'factor']
        db = DBframework()
        df_comb = db.DataFrame(data=dftemp, columns=header, index=None)
        return df_comb 
    #
    # ----------------------------
    #
    @property
    def df(self):
        """ """
        1 / 0
    
    @df.setter
    def df(self, df):
        """ """
        df = get_comb_df(df)
        group = df.groupby("name")
        for cname, items in group:
            for item in items.itertuples():
                #cname =  item.name
                try:
                    self.__getitem__(cname)
                except KeyError:
                    self.__setitem__(cname, item.title)
                #
                if item.basic:
                    self._combination[cname]._basic[item.basic] = item.factor
                
                if item.combination:
                    self._combination[cname]._combination[item.combination] = item.factor             
#
#
def find_comb_item(word_in:str) -> str:
    """ """
    key = {"name": r"\b((load)?(_|-|\s*)?(comb(ination)?)?(_|-|\s*)?name)\b",
           "basic": r"\b(basic(_|-|\s*)?(load)?)\b",
           "combination": r"\b((load)?(_|-|\s*)?comb(ination)?(_|-|\s*)?(load)?)\b",
           "factor": r"\b((load)?(_|-|\s*)?factor)\b",
           "title": r"\b(title|comment)\b",}
    match = common.find_keyword(word_in, key)
    return match    
#
#
def get_comb_df(df:DBframework.DataFrame):
    """ """
    try:
        columns = list(df.columns)
        header = {key: find_comb_item(key)
                  for key in columns}
        #
        df.rename(columns=header, inplace=True)
    except AttributeError:
        raise IOError('Combination df not valid')
    #
    columns = list(df.columns)
    if 'title' not in columns:
        df['title'] = df['name']
    #
    if 'basic' not in columns:
        df['basic'] = None
    #
    if 'combination' not in columns:
        df['combination'] = None
    #
    return df
#
def get_comb_list(values: list|tuple, steps: int = 5) -> list:
    """ [name, type, load_id, factor, title] """
    output = [None] * steps
    for x, item in enumerate(values):
        try:
            output[x] = item
        except IndexError:
            pass
    #
    return output
#
def get_comb_dict(values: dict, steps: int = 5) -> list:
    """ [name, type, load_id, factor, title] """
    output = [None] * steps
    
    for key, item in values.items():
        if re.match(r"\b((load)?(_|-|\s*)?(comb(ination)?)?(_|-|\s*)?name)\b", key, re.IGNORECASE):
            output[0] = item

        elif re.match(r"\b(basic(_|-|\s*)?(load)?)\b", key, re.IGNORECASE):
            output[1] = 'basic'
            output[2] = item
            
        elif re.match(r"\b((load)?(_|-|\s*)?comb(ination)?(_|-|\s*)?(load)?)\b", key, re.IGNORECASE):
            output[1] = 'combination'
            output[2] = item
            
        elif re.match(r"\b((load)?(_|-|\s*)?factor)\b", key, re.IGNORECASE):
            output[3] = item
            
        elif re.match(r"\b(title|comment)\b", key, re.IGNORECASE):
            output[4] = item
    
    return output
#
#