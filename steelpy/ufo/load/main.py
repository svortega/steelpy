# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
import re
from typing import NamedTuple

# package imports
from steelpy.ufo.mesh.sqlite.nodes import NodeSQL
from steelpy.ufo.mesh.sqlite.elements import ElementsSQL
from steelpy.ufo.load.sqlite.wave_load import MetoceanLoadSQL
from steelpy.ufo.load.sqlite.main import BasicLoadSQL, LoadCombSQL
#
from steelpy.ufo.load.concept.main import BasicLoadConcept, LoadCombConcept
from steelpy.ufo.load.concept.wave_load import MetoceanLoadIM
#from steelpy.ufo.load.concept.combination import LoadCombConcept
#from steelpy.ufo.load.inmemory.timehistory import TimeHistory
from steelpy.ufo.plot.main import PlotLoad
# 
#
from steelpy.utils.sqlite.utils import create_connection, create_table
#
#
#
class MeshLoad:
    """
    """
    __slots__ = ['_df_nodal', '_component', 
                 '_hydro', '_basic', '_combination',
                 '_elements', '_nodes', '_boundaries',
                 '_db_file']

    def __init__(self, 
                 mesh_type: str,
                 component: str|int, 
                 db_file: str | None = None):
        """
        """
        self._component = component
        self._db_file = db_file
        #
        # mesh 
        self._nodes = NodeSQL(db_system=mesh_type,
                              component=component, 
                              db_file=db_file)
        #
        self._elements = ElementsSQL(db_system=mesh_type,
                                     component=component,
                                     db_file=db_file)
        #
        #
        self._basic = BasicLoadSQL(db_file,
                                   component=component)
        #
        self._combination = LoadCombSQL(db_file,
                                        component=component)
        #
        self._hydro = MetoceanLoadSQL(db_file=db_file,
                                      component=component)
        #
        conn = create_connection(self._db_file)
        with conn: 
            self._new_table(conn)        
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
    def _new_table(self, conn):
        """ """
        table = "CREATE TABLE IF NOT EXISTS Load(\
                number INTEGER PRIMARY KEY NOT NULL,\
                name NOT NULL,\
                mesh_id INTEGER NOT NULL REFERENCES Mesh(number), \
                level TEXT NOT NULL,\
                title TEXT,\
                input_type TEXT, \
                input_file TEXT);"
        create_table(conn, table)
    #
    #
    # ----------------------------
    #
    #@property
    #def plane(self) -> NamedTuple:
    #    """
    #    """
    #    return self._plane
    #    
    #@plane.setter
    #def plane(self, plane: NamedTuple) -> None:
    #    """
    #    """
    #    self._plane = plane
    #    self._basic._plane = self._plane
    #    self._hydro._plane = self._plane
    #    self._combination._plane = self._plane
    #
    #
    # -----------------------------------------------
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
                 '_point', '_elements', '_boundaries', '_component']
    
    def __init__(self, points, elements, boundaries, component: int) -> None:
        """
        """
        self._point = points
        self._elements = elements
        self._boundaries = boundaries
        self._component = component
        #
        self._basic = BasicLoadConcept(points=self._point,
                                       elements=self._elements,
                                       component=component)
        #
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
