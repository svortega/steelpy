# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
import re
from typing import NamedTuple

# package imports
from steelpy.f2uModel.mesh.sqlite.nodes import NodeSQL
from steelpy.f2uModel.mesh.sqlite.elements import ElementsSQL
# steelpy.f2uModel.load
#from .inmemory.main import LoadingInmemory
from .concept.main import BasicLoadConcept, LoadCombConcept
#from .concept.combination import LoadCombConcept
#from .inmemory.timehistory import TimeHistory
from steelpy.f2uModel.load.concept.wave_load import MetoceanLoadIM
#
#from .sqlite.main import LoadingSQL
from .sqlite.main import BasicLoadSQL, LoadCombSQL
from steelpy.f2uModel.load.sqlite.wave_load import MetoceanLoadSQL
#
from steelpy.f2uModel.plot.main import PlotLoad

#
#
class MeshLoad:
    """
    """
    __slots__ = ['_df_nodal', '_plane',
                 '_hydro', '_basic', '_combination',
                 '_elements', '_nodes', '_boundaries']

    def __init__(self, #nodes, elements, boundaries, 
                 plane: NamedTuple, 
                 mesh_type: str,
                 db_file: str | None = None):
        """
        """
        #
        #self._number = []
        #self._labels = []        
        # mesh components
        #self._nodes = nodes
        #self._elements = elements
        #self._boundaries = boundaries
        self._plane = plane
        #
        # mesh 
        self._nodes = NodeSQL(db_system=mesh_type,
                              plane=self._plane,
                              db_file=db_file)
        #
        self._elements = ElementsSQL(plane=self._plane,
                                     db_system=mesh_type,
                                     db_file=db_file)
        #
        #
        self._basic = BasicLoadSQL(db_file, plane=self._plane)
        self._combination = LoadCombSQL(db_file, plane=self._plane)
        #
        self._hydro = MetoceanLoadSQL(db_file=db_file, plane=self._plane)
    #
    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += self._basic.__str__()
        output += self._combination.__str__()
        return output
    #
    # ----------------------------
    #
    @property
    def plane(self) -> NamedTuple:
        """
        """
        return self._plane
        
    @plane.setter
    def plane(self, plane: NamedTuple) -> None:
        """
        """
        self._plane = plane
        self._basic._plane = self._plane
        self._hydro._plane = self._plane
        self._combination._plane = self._plane

    #
    # -----------------------------------------------
    #
    def basicXX(self, values: None|list|tuple = None,
              df = None):
        """
        """
        if isinstance(values, (list, tuple)):
            1 / 0
            #number = len(self._basic) #+ 1
            #if isinstance(values[0], (list,tuple)):
            #    for item in values:
            #        name = item[0]
            #        try:
            #            self._basic[name]
            #        except IOError:
            #            self._basic[name] = name
            #        #
            #        if re.match(r"\b(node(s)?|point(s)?)\b", item[1], re.IGNORECASE):
            #            self._basic[name].node(item[2:])
            #        
            #        elif re.match(r"\b(beam(s)?)\b", item[1], re.IGNORECASE):
            #            self._basic[name].beam(item[2:])
            #        
            #        else:
            #            raise IOError(f"Basic load type {item[1]} not available")
            #else:
            #    name = values[0]
            #    try:
            #        self._basic[name]
            #    except IOError:
            #        self._basic[name] = name
            #    #
            #    if re.match(r"\b(node(s)?|point(s)?)\b", values[1], re.IGNORECASE):
            #        self._basic[name].node(values[2:])
            #    
            #    elif re.match(r"\b(beam(s)?)\b", values[1], re.IGNORECASE):
            #        self._basic[name].beam(values[2:])
        
        #
        # dataframe input
        try:
            columns = list(df.columns)
            1 / 0
            #self._sections.df(df)
        except AttributeError:
            pass            
        #
        # name input
        #if name:
        #    try:
        #        return self._basic[name]
        #    except IOError:
        #        number = len(self._basic) + 1
        #        #self._number.append(number)
        #        #self._labels.append(name)                
        #        if not title:
        #            title = f"{name}_{number}"
        #        self._basic[name] = title
        #        return self._basic[name]
        #
        return self._case
    #
    #
    def basic(self, values:tuple|list|None=None,
              df=None):
        """ """
        if isinstance(values, (list, tuple)):
            #number = len(self._basic) #+ 1
            if isinstance(values[0], (list,tuple)):
                for item in values:
                    name = item[0]
                    try:
                        self._basic[name]
                    except IOError:
                        self._basic[name] = name
                    #
                    if re.match(r"\b(node(s)?|point(s)?)\b", item[1], re.IGNORECASE):
                        self._basic[name].node(item[2:])
                    
                    elif re.match(r"\b(beam(s)?)\b", item[1], re.IGNORECASE):
                        self._basic[name].beam(item[2:])
                    
                    else:
                        raise IOError(f"Basic load type {item[1]} not available")
            else:
                name = values[0]
                try:
                    self._basic[name]
                except IOError:
                    self._basic[name] = name
                #
                if re.match(r"\b(node(s)?|point(s)?)\b", values[1], re.IGNORECASE):
                    self._basic[name].node(values[2:])
                
                elif re.match(r"\b(beam(s)?)\b", values[1], re.IGNORECASE):
                    self._basic[name].beam(values[2:])
        
        #
        # dataframe input
        try:
            columns = list(df.columns)
            1 / 0
            #self._sections.df(df)
        except AttributeError:
            pass         
        #
        return self._basic 
    #    
    #
    # ----------------------------
    #
    #
    def metocean(self, condition:dict|list|None=None,
                 df=None):
        """
        design_load : max_BS
        criterion : select design load based on local (member) or global system
        """
        #
        #if isinstance(condition, dict):
        if condition:
            #cases = []
            for key, item in condition.items():
                self._hydro[key] = item
                
            #1 / 0
            #if isinstance(values[0], (list, tuple)):
            #    for value in values:            
            #        self._hydro[value[0]] = value[1:]
            #else:
            #    self._hydro[values[0]] = [*values[1:], 'local']
        #elif isinstance(condition, (list, tuple)):
        #    1 / 0
        
        #
        # dataframe input
        try:
            columns = list(df.columns)
            header = {}
            1 / 0
        except AttributeError:
            pass        
        #
        return self._hydro
    #
    # ----------------------------
    #
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
    # ----------------------------
    #
    def time_history(self):
        """
        """
        return self.th

    #
    # ----------------------------
    #
    #def mass(self):
    #    """
    #    :return:
    #    """
    #    return self._mass
    #
    #
    #
    # --------------------
    # Plotting
    # --------------------
    #
    #@property
    def plot(self, figsize:tuple = (10, 10)):
        """ """
        return PlotLoad(cls=self, figsize=figsize)
#
#
#
class ConceptLoad:
    
    __slots__ = ["_basic", "_combination", '_hydro', '_mass',
                 '_nodes', '_elements', '_boundaries']
    
    def __init__(self, points, elements, boundaries) -> None:
        """
        """
        self._nodes = points
        self._elements = elements
        self._boundaries = boundaries
        #
        self._basic = BasicLoadConcept(points=self._nodes,
                                       elements=self._elements)
        self._hydro = MetoceanLoadIM(elements=self._elements)
        #self.th = TimeHistory()
        self._combination = LoadCombConcept(basic_load=self._basic)
        self._mass = LoadCombConcept(basic_load=self._basic)
    #
    #
    def basic(self):
        """
        """
        return self._basic
    #
    #
    def combination(self):
        """
        """
        return self._combination
    #
    #
    def metocean(self, condition:dict|list|None=None,
                 df=None):
        """
        design_load : max_BS
        criterion : select design load based on local (member) or global system
        """
        #
        #if isinstance(condition, dict):
        if condition:
            #cases = []
            for key, item in condition.items():
                self._hydro[key] = item
                
            #1 / 0
            #if isinstance(values[0], (list, tuple)):
            #    for value in values:            
            #        self._hydro[value[0]] = value[1:]
            #else:
            #    self._hydro[values[0]] = [*values[1:], 'local']
        #elif isinstance(condition, (list, tuple)):
        #    1 / 0
        
        #
        # dataframe input
        try:
            columns = list(df.columns)
            header = {}
            1 / 0
        except AttributeError:
            pass        
        #
        return self._hydro
    #
    #
    #def time_history(self):
    #    """
    #    """
    #    return self.th
    #
    def mass(self):
        """
        """
        return self._mass
    #
        #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += self._basic.__str__()
        output += self._combination.__str__()
        return output
    #
    # --------------------
    # Plotting
    # --------------------
    #
    #@property
    def plot(self, figsize:tuple = (10, 10)):
        """ """
        return PlotLoad(cls=self, figsize=figsize)    
#
