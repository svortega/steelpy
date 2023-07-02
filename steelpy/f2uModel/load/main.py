# 
# Copyright (c) 2009-2023 fem2ufo
#

# Python stdlib imports
from __future__ import annotations
import re
# from typing import NamedTuple, Dict, List, Iterable, Union

# package imports
# steelpy.f2uModel.load
#from .inmemory.main import LoadingInmemory
from .inmemory.basic_load import BasicLoad
from .inmemory.combination import LoadCombination
from .inmemory.timehistory import TimeHistory
#
#from .sqlite.main import LoadingSQL
from .sqlite.basic_load import BasicLoadSQL
from .sqlite.combination import LoadCombSQL

#
#
class Load:
    """
    """
    __slots__ = ['_labels', '_number', '_df_nodal',
                 '_basic', '_combination']

    def __init__(self, nodes, elements,
                 mesh_type: str,
                 db_file: str | None = None):
        """
        """
        self._number = []
        self._labels = []
        if mesh_type != "inmemory":
            #self._load = LoadingSQL(db_file=db_file,
            #                        db_system=mesh_type)
            self._basic = BasicLoadSQL(db_file)
            self._combination = LoadCombSQL(db_file)
        else:
            #self._load = LoadingInmemory(nodes=nodes,
            #                             elements=elements)
            self._basic = BasicLoad(nodes=nodes, elements=elements)
            self._combination = LoadCombination(basic_load=self._basic)

    #
    # @property
    # def basic(self):
    #    """
    #    """
    #    return self._load._basic

    # @basic.setter
    def basic(self, values: None|list|tuple = None,
              df = None):
        """
        """
        if isinstance(values, (list, tuple)):
            number = len(self._basic) #+ 1
            if isinstance(values[0], (list,tuple)):
                for item in values:
                    name = item[0]
                    try:
                        index = self._labels.index(name)
                        number = self._number[index]
                    except ValueError:
                        number += 1
                        self._basic[number] = name
                        self._number.append(number)
                        self._labels.append(name)
                    #
                    if re.match(r"\b(node(s)?|point(s)?)\b", item[1], re.IGNORECASE):
                        self._basic[number].node(item[2:])
                    
                    elif re.match(r"\b(beam(s)?)\b", item[1], re.IGNORECASE):
                        self._basic[number].beam(item[2:])
                    
                    else:
                        raise IOError(f"Basic load type {item[1]} not available")
            else:
                name = values[0]
                try:
                    index = self._labels.index(name)
                    number = self._number[index]
                except ValueError:
                    number += 1
                    self._basic[number] = name
                    self._number.append(number)
                    self._labels.append(name)
                #
                if re.match(r"\b(node(s)?|point(s)?)\b", values[1], re.IGNORECASE):
                    self._basic[number].node(values[2:])
                
                elif re.match(r"\b(beam(s)?)\b", values[1], re.IGNORECASE):
                    self._basic[number].beam(values[2:])                
        #
        # dataframe input
        try:
            df.columns
            #self._sections.df(df)
        except AttributeError:
            pass
        #
        return self._basic

    #
    # @property
    def combination(self, values: None | list = None,
                    df = None):
        """
        """
        if isinstance(values, list):
            number = len(self._combination)
            for item in values:
                name = item[0]
                try:
                    index = self._labels.index(name)
                    number = self._number[index]
                except ValueError:
                    number += 1
                    self._combination[number] = name
                    self._number.append(number)
                    self._labels.append(name)
                #
                if re.match(r"\b(basic(_)?(load)?)\b", item[1], re.IGNORECASE):
                    self._combination[number]._basic[item[2]] = item[3]
                
                elif re.match(r"\b(comb(ination)?(_)?(load)?)\b", item[1], re.IGNORECASE):
                    self._combination[number]._combination[item[2]] = item[3]
                
                else:
                    raise IOError(f"Combination load type {item[1]} not available")
        #
        #
        # dataframe input
        try:
            df.columns
            #self._sections.df(df)
        except AttributeError:
            pass
        #
        return self._combination

    #
    def time_history(self):
        """
        """
        return self.th

    #
    def mass(self):
        """
        :return:
        """
        return self._mass

    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += self._basic.__str__()
        output += self._combination.__str__()
        return output
#
#
