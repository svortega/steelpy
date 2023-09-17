# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass

#
# package imports


class BasicModel:
    
    # --------------------
    # Common
    # --------------------
    #    
    def materials(self, values: None|list|dict=None,
                  df=None):
        """
        """
        if isinstance(values, list):
            #1/0
            if isinstance(values[0], list):
                for item in values:
                    self._materials[item[0]] = item[1:]
            else:
                self._materials[values[0]] = values[1:]
        #
        try:
            df.columns
            self._materials.df = df
        except AttributeError:
            pass
        #
        return self._materials
    #
    #
    def sections(self, values: None|list|dict=None,
                 df=None):
        """
        """
        if isinstance(values, list):
            #1/0
            if isinstance(values[0], list):
                for item in values:
                    self._sections[item[0]] = item[1:]
            else:
                self._sections[values[0]] = values[1:]
        #
        # dataframe input
        try:
            df.columns
            self._sections.df = df
        except AttributeError:
            pass
        #
        return self._sections
    #
    #
    def groups(self):
        """
        """
        return self._groups    
    #
    # --------------------
    # Mesh items
    # -------------------- 
    #
    #
    def elements(self, values:None|list|tuple=None,
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
            df.columns   
            self._elements.df(df)
        except AttributeError:
            pass
        return self._elements    
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
        #self._load = Load(plane=self._plane,
        #                  mesh_type=self.data_type,
        #                  db_file=self.db_file)
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
