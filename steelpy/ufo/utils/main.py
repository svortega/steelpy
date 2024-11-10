# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
from collections.abc import Mapping
#from typing import NamedTuple
import re

#
# package imports
#from steelpy.ufo.process.element import find_element_item
from steelpy.utils.dataframe.main import DBframework


class ufoBasicModel:
    __slots__ = ['_component']
    
    def __init__(self, component:str|int) -> None:
        """ """
        self._component = component
    #   
    # --------------------
    # Common
    # --------------------
    #    
    def material(self, values: None | list | tuple | dict = None,
                  df = None):
        """
        """
        if values :
            if isinstance(values, dict):
                mname = values['name']
                if isinstance(mname, (list | tuple)):
                    db = DBframework()
                    self._materials.df = db.DataFrame(values)
                else:
                    mname = values.pop('name')
                    self._materials[mname] = values

            elif isinstance(values, (list, tuple)):
                if isinstance(values[0], (list, tuple)):
                    for item in values:
                        self._materials[item[0]] = item[1:]
                elif isinstance(values[0], dict):
                    for item in values:
                        self._materials['name'] = item
                else:
                    self._materials[values[0]] = values[1:]

            else:
                raise IOError('Material input data not valid')
        #
        try:
            #df.columns
            columns = list(df.columns)
            self._materials.df = df
        except AttributeError:
            pass
        #
        return self._materials
    #
    #
    def section(self, values: None|list|tuple|dict = None,
                 df=None):
        """
        """
        if values:
            if isinstance(values, dict):
                sname = values['name']
                if isinstance(sname, (list | tuple)):
                    db = DBframework()
                    newdf = db.DataFrame(values)
                    self._sections.df = newdf
                else:
                    self._sections[sname] = values
            elif isinstance(values, (list, tuple)):
                if isinstance(values[0], (list, tuple)):
                    for item in values:
                        self._sections[item[0]] = item[1:]
                elif isinstance(values[0], dict):
                    for item in values:
                        self._sections[item['name']] = item
                else:
                    self._sections[values[0]] = values[1:]
            else:
                raise IOError('Section input data not valid')
        #
        # dataframe input
        try:
            #df.columns
            columns = list(df.columns)
            self._sections.df = df
        except AttributeError:
            pass
        #
        return self._sections
    #
    #
    def group(self):
        """
        """
        return self._groups    
    #
    # --------------------
    # Mesh items
    # -------------------- 
    #
    #
    def element(self, values:None|list|tuple|dict=None,
                 df=None):
        """
        """
        if values:
            if isinstance(values, dict):
                #columns = list(values.keys())
                #columns = {item:find_element_item(item)
                #           for item in columns}
                ename = values['name']
                if isinstance(ename, (list | tuple)):
                    db = DBframework()
                    edf = db.DataFrame(values)
                    self._elements.df = edf
                else:
                    sname = values.pop('name')
                    self._elements[sname] = values                
            elif isinstance(values, (list, tuple)):
                if isinstance(values[0], (list,tuple)):
                    for value in values:
                        self._elements[value[0]] = value[1:]
                elif isinstance(values[0], dict):
                    for item in values:
                        self._elements[item['name']] = item                
                else:
                    self._elements[values[0]] = values[1:]
            else:
                raise IOError('Element input data not valid')
        #
        # dataframe input
        try:
            #df.columns
            columns = list(df.columns)
            self._elements.df = df
        except AttributeError:
            pass
        #
        return self._elements
    #
    #
    def boundary(self, values: None | list | tuple | dict = None,
                 df=None):
        """
        """
        if values:
            if isinstance(values, dict):
                bname = values['name']
                if isinstance(bname, (list | tuple)):
                    db = DBframework()
                    self._boundaries.df = db.DataFrame(values)
                else:
                    sname = values.pop('name')
                    self._boundaries[sname] = values   

            elif isinstance(values, (list, tuple)):
                if isinstance(values[0], (list, tuple)):
                    for item in values:
                        self._boundaries[item[0]] = item[1:]

                elif isinstance(values[0], dict):
                    for item in values:
                        self._boundaries[item['name']] = item
                else:
                    self._boundaries[values[0]] = values[1:]

            else:
                raise IOError('Boundary input data not valid')
        #
        # dataframe input
        try:
            df.columns
            self._boundaries.df = df
        except AttributeError:
            pass
            #
        return self._boundaries

    #
    #
    # --------------------
    # Load
    # -------------------- 
    #
    def load(self, values:None|list|tuple=None,
             df=None):
        """
        """
        #
        if isinstance(values, (list, tuple)):
            if isinstance(values[0], (list,tuple)):
                for item in values:
                    if re.match(r"\b(basic(\_)?(load)?)\b", item[0], re.IGNORECASE):
                        self._load.basic(item[1:])
                    elif re.match(r"\b(comb(ination)?(\_)?(load)?)\b", item[0], re.IGNORECASE):
                        self._load.combination(item[1:])
                    else:
                        raise IOError(f'load {item[0]}')
            else:
                if re.match(r"\b(basic(\_)?(load)?)\b", values[0], re.IGNORECASE):
                    self._load.basic(values[1:])
                elif re.match(r"\b(comb(ination)?(\_)?(load)?)\b", values[0], re.IGNORECASE):
                    self._load.combination(values[1:])
                else:
                    raise IOError(f'load {values[0]}')                
        #
        # dataframe input
        try:
            df.columns
            #self._boundaries.df(df)
        except AttributeError:
            pass
        #
        return self._load
#
#
#
class ModelClassBasic(Mapping):
    """ """
    __slots__ = ['_component']
    
    #
    def __init__(self, component: str|int):
        """
        """
        self._component = component
    #  
    #
    def __len__(self) -> int:
        return len(self._label)

    def __iter__(self):
        """
        """
        return iter(self._label)
    
    def __contains__(self, value) -> bool:
        return value in self._label
    #
    #
