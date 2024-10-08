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


class ufoBasicModel:
    #__slots__ = ['_name', '_title', ]
    
    #def __init__(self, name:str, title: str) -> None:
    #    """ """
    #    self._name = name
    #    self._title = title    
    
    # --------------------
    # Common
    # --------------------
    #    
    def material(self, values: None|list|dict=None,
                  df=None):
        """
        """
        if isinstance(values, (list, tuple)):
            #1/0
            if isinstance(values[0], (list, tuple)):
                for item in values:
                    self._materials[item[0]] = item[1:]
            else:
                self._materials[values[0]] = values[1:]
        #
        try:
            #df.columns
            columns = list(df.columns)
            header = {}
            for key in columns:
                if re.match(r"\b(id|name|material(s)?)\b", key, re.IGNORECASE):
                    header[key] = 'name'
                
                elif re.match(r"\b(type)\b", key, re.IGNORECASE):
                    header[key] = 'type'
                
                elif re.match(r"\b(fy|yield)\b", key, re.IGNORECASE):
                    header[key] = 'Fy'
                
                elif re.match(r"\b(fu)\b", key, re.IGNORECASE):
                    header[key] = 'Fu'
                
                elif re.match(r"\b(E(\-|\s*)?(modulus)?)\b", key, re.IGNORECASE):
                    header[key] = 'E'
                
                elif re.match(r"\b(poisson|nu)\b", key, re.IGNORECASE):
                    header[key] = 'poisson'
                
                elif re.match(r"\b(density)\b", key, re.IGNORECASE):
                    header[key] = 'density'
                
                elif re.match(r"\b(thermal(\-|\s*)?(x)?|alpha)\b", key, re.IGNORECASE):
                    header[key] = 'alpha'                
            #
            df['type'] = df['type'].apply(lambda x: 'elastic'
                                          if re.match(r"\b(elastic|linear(\_elastic)?)\b", x, re.IGNORECASE)
                                          else x)            
            #
            mat = df[header.keys()].copy()
            mat.rename(columns=header, inplace=True)
            try:
                mat['Fu'].replace(to_replace=[''], value=[mat['Fy'] / 0.75], inplace=True)
            except KeyError:
                mat['Fu'] = mat['Fy'] / 0.75
            #
            self._materials.df = mat
        except AttributeError:
            pass
        #
        return self._materials
    #
    #
    def section(self, values: None|list|dict=None,
                 df=None):
        """
        """
        if isinstance(values, (list, tuple)):
            #1/0
            if isinstance(values[0], (list, tuple)):
                for item in values:
                    self._sections[item[0]] = item[1:]
            else:
                self._sections[values[0]] = values[1:]
        #
        # dataframe input
        try:
            #df.columns
            columns = list(df.columns)
            header = {}
            for key in columns:
                if re.match(r"\b(id|name|section(s)?)\b", key, re.IGNORECASE):
                    header[key] = 'name'
                
                elif re.match(r"\b(type)\b", key, re.IGNORECASE):
                    header[key] = 'type'
                
                elif re.match(r"\b(title)\b", key, re.IGNORECASE):
                    header[key] = 'title'
                #
                # tubular
                #
                elif re.match(r"\b(d(iamet(er|re))?)\b", key, re.IGNORECASE):
                    header[key] = 'diameter'
                
                elif re.match(r"\b((wall(\_|\s*)?)?t(hickness)?(\_|\-|\s*)?(w)?)\b", key, re.IGNORECASE):
                    header[key] = 'wall_thickness'
                #
                # PG/box/tee/channel
                #
                elif re.match(r"\b(height)\b", key, re.IGNORECASE):
                    header[key] = 'height'
                
                elif re.match(r"\b(web_thickness)\b", key, re.IGNORECASE):
                    header[key] = 'web_thickness'
                
                elif re.match(r"\b((top(\_|\s*)?flange(\_|\s*)?)?width)\b", key, re.IGNORECASE):
                    header[key] = 'top_flange_width'
                
                elif re.match(r"\b(top(\_|\s*)?flange(\_|\s*)?thickness)\b", key, re.IGNORECASE):
                    header[key] = 'top_flange_thickness'
                
                elif re.match(r"\b(bottom(\_|\s*)?flange(\_|\s*)?width)\b", key, re.IGNORECASE):
                    header[key] = 'bottom_flange_width'
                
                elif re.match(r"\b(bottom(\_|\s*)?flange(\_|\s*)?thickness)\b", key, re.IGNORECASE):
                    header[key] = 'bottom_flange_thickness'
                
                elif re.match(r"\b(fillet(\_|\s*)?radius)\b", key, re.IGNORECASE):
                    header[key] = 'fillet_radius'
                #
                # ops
                #
                elif re.match(r"\b(SA(\_|\s*)?inplane)\b", key, re.IGNORECASE):
                    header[key] = 'SA_inplane'
                
                elif re.match(r"\b(SA(\_|\s*)?outplane)\b", key, re.IGNORECASE):
                    header[key] = 'SA_outplane'
                
                elif re.match(r"\b(shear(\_|\s*)?stress)\b", key, re.IGNORECASE):
                    header[key] = 'shear_stress'
                
                elif re.match(r"\b(build)\b", key, re.IGNORECASE):
                    header[key] = 'build'
                
                elif re.match(r"\b(compactness)\b", key, re.IGNORECASE):
                    header[key] = 'compactness'                
            #
            #
            sect = df[header.keys()].copy()
            sect.rename(columns=header, inplace=True)#
            #
            if not 'title' in sect:
                sect['title'] = None
            #
            self._sections.df = sect
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
    def element(self, values:None|list|tuple=None,
                 df=None):
        """
        """
        if isinstance(values, (list, tuple)):
            if isinstance(values[0], (list,tuple)):
                for value in values:
                    self._elements[value[0]] = value[1:]
            else:
                self._elements[values[0]] = values[1:]
        #
        #
        # dataframe input
        try:
            #df.columns
            columns = list(df.columns)
            header = {}
            for key in columns:
                if re.match(r"\b(id|name|element(s)?|member(s)?)\b", key, re.IGNORECASE):
                    header[key] = 'name'
                
                elif re.match(r"\b(type)\b", key, re.IGNORECASE):
                    header[key] = 'type'
                
                elif re.match(r"\b(title)\b", key, re.IGNORECASE):
                    header[key] = 'title'
                #
                #
                elif re.match(r"\b(material(s)?(\_|\s*)?(name|id)?)\b", key, re.IGNORECASE):
                    header[key] = 'material_name'
                
                elif re.match(r"\b(section(s)?(\_|\s*)?(name|id)?)\b", key, re.IGNORECASE):
                    header[key] = 'section_name'
                #
                #
                elif re.match(r"\b((node|point|end)?(\_|\s*)?1)\b", key, re.IGNORECASE):
                    header[key] = 'node_1'
                
                elif re.match(r"\b((node|point|end)?(\_|\s*)?2)\b", key, re.IGNORECASE):
                    header[key] = 'node_2'
                
                elif re.match(r"\b((node|point|end)?(\_|\s*)?3)\b", key, re.IGNORECASE):
                    header[key] = 'node_3'
                
                elif re.match(r"\b((node|point|end)?(\_|\s*)?4)\b", key, re.IGNORECASE):
                    header[key] = 'node_4'
                #
                #
                elif re.match(r"\b(alpha|roll(\_|\s*)?angle)\b", key, re.IGNORECASE):
                    header[key] = 'roll_angle'               
            #
            #
            members = df[header.keys()].copy()
            members.rename(columns=header, inplace=True) 
            self._elements.df = members
        except AttributeError:
            pass
        return self._elements
    #
    #
    def boundary(self, values: None | list | tuple | dict = None,
                 df=None):
        """
        """
        # 1 / 0
        if values:
            if isinstance(values, dict):
                1 / 0

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
    __slots__ = ['_labels']
    
    #
    def __init__(self):
        """
        """
        self._labels:list = []
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
