# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
import re
#from typing import NamedTuple

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
from steelpy.ufo.load.process.combination import (get_comb_dict,
                                                  get_comb_list)

from steelpy.ufo.load.process.load_case import (get_bload_list,
                                                get_bload_dict)

from steelpy.ufo.plot.main import PlotLoad
# 
#
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework
#
#
#
class MasterLoad:
    """ """
    __slots__ = ['_component']
    
    def __init__(self, component: str|int):
        """
        """
        self._component = component    
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
    def _comb_setup(self, values: list):
        """" values : [name, type, load_id, factor, title] """
        cname = values[0]
        ltype = values[1]
        lname = values[2]
        factor = float(values[3])
        try:
            title = values[4]
        except IndexError:
            title = cname        
        #
        try:
            self._combination[cname]
        except KeyError:
            self._combination[cname] = title
        #
        if re.match(r"\b(basic(_|-|\s*)?(load)?)\b", ltype, re.IGNORECASE):
            self._combination[cname]._basic[lname] = factor
        
        elif re.match(r"\b((load)?(_|-|\s*)?comb(ination)?(_|-|\s*)?(load)?)\b", ltype, re.IGNORECASE):
            self._combination[cname]._combination[lname] = factor
        
        else:
            raise IOError(f"Combination load type {ltype} not available")        
    #
    def combination(self, values: None|list|tuple|dict = None,
                    df = None):
        """
        """
        if values:
            if isinstance(values, dict):
                comb = get_comb_dict(values)
                sname = comb[2] # lname
                if isinstance(sname, (list | tuple)):
                    db = DBframework()
                    newdf = db.DataFrame(values)
                    self._combination.df = newdf
                else:
                    self._comb_setup(values=comb)
            elif isinstance(values, (list, tuple)):
                for value in values:
                    if isinstance(value, (list, tuple)):
                        comb = get_comb_list(value)
                    elif isinstance(value, dict):
                        comb = get_comb_dict(value)
                    else:
                        raise IOError(f'load format not valid')
                    #
                    self._comb_setup(values=comb)
        #
        # dataframe input
        try:
            columns = list(df.columns)
            self._combination.df = df
        except AttributeError:
            pass
        #
        return self._combination
    #
    # -----------------------------------------------
    #
    def _basic_setup(self, values: list):
        """ [name, item, item_id, load_type, values] """
        name = values[0]
        item = values[1]
        load = values[2:]
        #
        try:
            self._basic[name]
        except IOError:
            self._basic[name] = name
        #
        if re.match(r"\b(node(s)?|point(s)?)\b", item, re.IGNORECASE):
            self._basic[name].node(load)
        
        elif re.match(r"\b(beam(s)?)\b", item, re.IGNORECASE):
            self._basic[name].beam(load)
        
        else:
            raise IOError(f"Basic load type {item[1]} not available")        
    #
    def basic(self, values:None|list|tuple|dict = None,
              df=None):
        """ """
        if values:
            if isinstance(values, dict):
                bload = get_bload_dict(values)
                itype = bload[2] # item_type: node/beam
                if isinstance(itype, (list | tuple)):
                    db = DBframework()
                    newdf = db.DataFrame(values)
                    self._basic.df = newdf
                else:
                    self._basic_setup(values=bload)
            
            elif isinstance(values, (list, tuple)):
                for value in values:
                    if isinstance(value, (list, tuple)):
                        bload = get_bload_list(value)
                    elif isinstance(value, dict):
                        bload = get_bload_dict(value)
                    else:
                        raise IOError(f'load format not valid')
                    #
                    self._basic_setup(values=bload)
        
        #
        # dataframe input
        try:
            columns = list(df.columns)
            self._basic.df = df
        except AttributeError:
            pass         
        #
        return self._basic 
    #    
    #
    # ----------------------------
    #
    #
    def metocean(self, condition:None|list|tuple|dict = None,
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
    #def time_history(self):
    #    """
    #    """
    #    return self.th
    #
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
    # ----------------------------
    # Plotting
    # ----------------------------
    #
    #@property
    def plot(self, figsize:tuple = (10, 10)):
        """ """
        return PlotLoad(cls=self, figsize=figsize)
#
#
#
class MeshLoad(MasterLoad):
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
        super().__init__(component)
        #self._component = component
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
    #def __str__(self) -> str:
    #    """ """
    #    output = "\n"
    #    output += self._basic.__str__()
    #    output += self._combination.__str__()
    #    return output
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
#
#
class ConceptLoad(MasterLoad):
    
    __slots__ = ["_basic", "_combination", '_hydro', '_mass',
                 '_point', '_elements', '_boundaries', '_component']
    
    def __init__(self, points, elements, boundaries,
                 component: int|str) -> None:
        """
        """
        super().__init__(component)
        #
        self._point = points
        self._elements = elements
        self._boundaries = boundaries
        #self._component = component
        #
        self._basic = BasicLoadConcept(points=self._point,
                                       elements=self._elements,
                                       component=component)
        #
        self._hydro = MetoceanLoadIM(elements=self._elements)
        #self.th = TimeHistory()
        self._combination = LoadCombConcept(basic_load=self._basic,
                                            component=component)
        self._mass = LoadCombConcept(basic_load=self._basic,
                                     component=component)
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
    #def __str__(self) -> str:
    #    """ """
    #    output = "\n"
    #    output += self._basic.__str__()
    #    output += self._combination.__str__()
    #    return output
    #
    # --------------------
    # Plotting
    # --------------------
    #
    #@property
    #def plot(self, figsize:tuple = (10, 10)):
    #    """ """
    #    return PlotLoad(cls=self, figsize=figsize)    
#
