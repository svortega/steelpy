#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
#from collections import defaultdict
#from dataclasses import dataclass
#from operator import sub, add
from typing import NamedTuple
import re
#
# package imports
#from steelpy.ufo.load.process.utils import find_load_type
from steelpy.utils.dataframe.main import DBframework
#
#
#
# ---------------------------------
#
#
class BasicLoad(NamedTuple):
    """ Basic load transfer"""
    name: int
    number: int
    title: str
    nodal_load: list
    member_load: dict
#
#
class BasicLoadMain(Mapping):
    
    def __init__(self):
        """
        """
        pass
    #
    # -----------------------------------------------
    #    
    def __len__(self) -> int:
        return len(self._labels)

    def __iter__(self):
        """
        """
        return iter(self._labels)
    
    def __contains__(self, value):
        return value in self._labels
    #
    # -----------------------------------------------
    #
    def get_number(self, start:int=1):
        """
        """
        try:
            n = max(self._number) + 1
        except ValueError:
            n = start
        #
        while True:
            yield n
            n += 1
    #
    #
    #  
#
#
class LoadCaseBasic(BasicLoadMain):
    __slots__ = ['_labels', '_title','_number', 'gravity']
    
    def __init__(self):
        """
        """
        super().__init__()
        self.gravity = 9.80665  # m/s^2
    #
    # -----------------------------------------------
    #
    def __str__(self, units: str = "si") -> str:
        """ """
        unit_lenght = " m"
        unit_force = "  N"
        unit_bm = "N*m"
        unit_fl = "N/m"
        #
        output = "\n"
        output += "{:}\n".format(80 * "_")
        output += "\n"
        output += f"{35 * ' '}BASIC LOAD\n"
        output += "\n"
        output += f"--- Node\n"
        output += f"Node Name{16 * ' '}fx [{unit_force}] fy [{unit_force}] fz [{unit_force}] System Complex\n"
        output += f"Load{21* ' '}mx [{unit_bm}] my [{unit_bm}] mz [{unit_bm}] Comment\n"
        output += "\n"
        output += f"--- Beam \n"
        output += f"Element Name{6 * ' '}L1[{unit_lenght}] qx1[{unit_fl}] qy1[{unit_fl}] qz1[{unit_fl}] System Complex\n"
        output += f"Line Load{9 * ' '}L2[{unit_lenght}] qx2[{unit_fl}] qy2[{unit_fl}] qz2[{unit_fl}] Comment\n"
        output += "\n"
        output += f"--- Beam \n"
        output += f"Element Name{6 * ' '}L1[{unit_lenght}] fx [{unit_force}] fy [{unit_force}] fz [{unit_force}] System Complex\n"
        output += f"Point Load{15 * ' '}mx [{unit_bm}] my [{unit_bm}] mz [{unit_bm}] Comment\n"
        output += "\n"
        output += f"--- Gravity/Selfweight\n"
        output += "\n"
        output += "{:}\n".format(80 * ".")
        output += "\n"
        loadname = list(dict.fromkeys(self._labels))
        for key in loadname:
            lcase = self.__getitem__(key)
            output += lcase.__str__()
            output += "\n"
        # print('---')
        return output
    #   
    #
    # -----------------------------------------------
    #
    @property
    def g(self):
        """"""
        return self.gravity

    #
    #
    # -----------------------------------------------
    #
    def _set_load(self, values, steps: int):
        """ """
        columns = list(values.keys())
        for key in columns:
            if re.match(r"\b(load(s)?((_|-|\s*)?id|name)?)\b", key, re.IGNORECASE):
                values['load'] = check_column(values[key], steps=steps)
                for x, load_name in enumerate(values['load']):
                    load_title = f"{load_name}_{x + 1}"
                    try:
                        self.__setitem__(load_name, load_title)
                    except Warning:
                        pass
                break
        #
        #db = DBframework()
        #dfnew = db.DataFrame(data=values)        
        return values
    #
    #
    def node(self, values:tuple|list|None=None,
             df=None):
        """ """
        if values:
            # Input data for specific basic node load
            if isinstance(values, dict):
                nodeid = values['node']
                if isinstance(nodeid, (list, tuple)):
                    nitems = len(nodeid)
                    values = self._set_load(values, steps=nitems)
                    db = DBframework()
                    dfnew = db.DataFrame(data=values)                    
                    #
                    basic = dfnew.groupby(['load'])
                    for x, load_name in enumerate(basic.groups):
                        load = basic.get_group((load_name,)) #.copy()
                        bload = self.__getitem__(load_name)
                        bload._node.df = load                    
                else:
                    load_name = values['load']
                    try:
                        self.__setitem__(load_name, load_name)
                    except Warning:
                        pass
                    #
                    bload = self.__getitem__(load_name)
                    bload._node[nodeid] = values
            
            elif isinstance(values, (list, tuple)):
                if isinstance(values[0], (list, tuple, dict)):
                    for item in values:
                        if isinstance(item, dict):
                            load_name = item['load']
                            nodeid = item['node']
                            load = item
                        elif isinstance(item, (list, tuple)):
                            load_name = item[0]
                            nodeid = item[1]
                            load = item[2:]
                            #nodeid = item.pop(1)
                            #load_name = item.pop(0)
                        #
                        try:
                            self.__setitem__(load_name, load_name)
                        except Warning:
                            pass
                        #
                        bload = self.__getitem__(load_name)
                        bload._node[nodeid] = load
                else:
                    load_name = values[0]
                    try:
                        self.__setitem__(load_name, load_name)
                    except Warning:
                        pass
                    #
                    bload = self.__getitem__(load_name)
                    nodeid = values[1]
                    bload._node[nodeid] = values[2:]
        #
        #
        # dataframe input
        try:
            columns = list(df.columns)
            nodeid = df['node']
            bitems = len(nodeid)
            values = self._set_load(df, steps=bitems)
            #
            # Set basic load if doesn't exist
            basic = values.groupby(['load'])
            for x, load_name in enumerate(basic.groups):
                bload = self.__getitem__(load_name)
                nload = basic.get_group((load_name,))
                bload._node.df = nload
            # End loop, return none cos no specific load_name
            #1 / 0
            return
        except AttributeError:
            pass
        #
        #
        return self._nodes
    #
    #
    def beam(self, values:tuple|list|dict|None=None,
             df=None):
        """ """
        if values:
            if isinstance(values, dict):
                beamid = values['beam']
                if isinstance(beamid, (list, tuple)):
                    bitems = len(beamid)
                    values = self._set_load(values, steps=bitems)
                    db = DBframework()
                    dfnew = db.DataFrame(data=values)                    
                    #
                    basic = dfnew.groupby(['load'])
                    for x, load_name in enumerate(basic.groups):
                        load = basic.get_group((load_name,)) #.copy()
                        bload = self.__getitem__(load_name)
                        bload._beam.df = load
                else:
                    load_name = values['load']
                    try:
                        self.__setitem__(load_name, load_name)
                    except Warning:
                        pass
                    bload = self.__getitem__(load_name)                    
                    bload._beam[beamid] = values
            elif isinstance(values, (list, tuple)):
                if isinstance(values[0], (list, tuple, dict)):
                    for item in values:
                        if isinstance(item, dict):
                            load_name = item['load']
                            beamid = item['beam']
                            load =  item
                        elif isinstance(item, (list, tuple)):
                            load_name = item[0]
                            beamid = item[1]
                            load =  item[2:]
                            #beamid = item.pop(1)
                            #load_name = item.pop(0)
                        # check if basic load name exist
                        try:
                            self.__setitem__(load_name, load_name)
                        except Warning:
                            pass
                        # push beam load data
                        bload = self.__getitem__(load_name)
                        bload._beam[beamid] = load
                else:
                    load_name = values[0]
                    try:
                        self.__setitem__(load_name, load_name)
                    except Warning:
                        pass
                    bload = self.__getitem__(load_name)
                    beamid = values[1]
                    bload._beam[beamid] = values[2:]        
        #
        # dataframe input
        try:
            columns = list(df.columns)
            beamid = df['beam']
            bitems = len(beamid)
            values = self._set_load(df, steps=bitems)
            # Set basic load if doesn't exist
            basic = values.groupby(['load'])
            for x, load_name in enumerate(basic.groups):
                load = basic.get_group((load_name,)) #.copy()
                bload = self.__getitem__(load_name)
                bload._beam.df = load            
            #
            return 
        except AttributeError:
            pass
        #
        #print('beam load')
        return self._beams
#  
#
def check_column(col: list|tuple|str|int, steps: int) -> list:
    """ """
    if isinstance(col, (list, tuple)):
        if len(col) != steps:
            raise IOError(f'items {len(col)} <> {steps}')
        return col
    else:
        return [col for _ in range(steps)]
# 
#
#
class LoadTypeBasic(BasicLoadMain):
    __slots__ = ['name', 'title', 'number', 
                 '_node', '_beam', '_selfweight']
    
    def __init__(self):
        """
        """
        super().__init__()
    #
    def __setitem__(self, load_name: str | int,
                    properties: list[float]) -> None:
        """
        """
        try:
            self._labels.index(load_name)
            self.name = load_name
            self.title = properties[0]
            self.number = properties[1]
        except ValueError:
            raise Exception(f'Load {load_name} already exist')
    
    
    def __getitem__(self, load_name: int|str):
        """
        node_name : node number
        """
        try:
            index = self._labels.index(load_name)
            #load_type = self._type[index]
            return self
        except ValueError:
            raise KeyError(f'   *** Load {load_name} does not exist')
    #    
    #
    # -----------------------------------------------
    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += f"Load Name : {str(self.name):12s}  Number : {self.number:8.0f}  Title : {self.title}\n"
        # node load
        if self._node:
            output += f"--- Node\n"
            output += self._node.__str__()
        # beam line
        if self._beam:
            output += f"--- Beam \n"
            output += self._beam.__str__()
        #
        if self._selfweight:
            output += f"--- Gravity/Selfweight\n"
            output += self._selfweight.__str__()
        #
        #output += "\n"
        # print('---')
        return output
    #
    #
    # -----------------------------------------------
    #
    #@property
    def gravity(self, values:list|None=None):
        """
        The self weight form allows you to specify multipliers to 
        acceleration due to gravity (g) in the X, Y, and Z axes. 
        If switched on, the default self weight acts in the Y axis 
        with a magnitude and sign of -1."""
        #
        if isinstance(values, list):
            if isinstance(values[0], list):
                for value in values:
                    self._selfweight[value[0]] = value[1:]
            else:
                self._selfweight[values[0]] = values[1:]
        return self._selfweight
    
    def selfweight(self, values:list|None=None):
        """
        The self weight form allows you to specify multipliers to 
        acceleration due to gravity (g) in the X, Y, and Z axes. 
        If switched on, the default self weight acts in the Y axis 
        with a magnitude and sign of -1."""
        #
        return self.gravity(values) 
    #
    #
    # -----------------------------------------------
    #
    #
    def node(self, values:tuple|list|None=None,
             df=None):
        """ Nodal load"""
        if values:
            # Input data for specific basic node load
            if isinstance(values, dict):
                nodeid = values['node']
                if isinstance(nodeid, (list, tuple)):
                    db = DBframework()
                    dfnew = db.DataFrame(data=values)
                    dfnew['load'] = self.name
                    self._node.df = dfnew
                else:
                    self._node[nodeid] = values
            elif isinstance(values, (list, tuple, dict)):
                if isinstance(values[0], (list, tuple, dict)):
                    for item in values:
                        if isinstance(item, dict):
                            nodeid = item['node']
                            load = item
                        elif isinstance(item, (list, tuple)):
                            nodeid = item[0]
                            load = item[1:]
                        #
                        self._node[nodeid] = load
                else:
                    self._node[values[0]] = values[1:]
        #
        # dataframe input
        try:
            columns = list(df.columns)
            df['load'] = self.name
            self._node.df = df            
        except AttributeError:
            pass
        #
        return self._node
    #
    #
    def beam(self, values:tuple|list|dict|None=None,
             df=None):
        """ beam loading """
        if values:
            if isinstance(values, dict):
                beamid = values['beam']
                if isinstance(beamid, (list, tuple)):
                    db = DBframework()
                    dfnew = db.DataFrame(data=values)
                    dfnew['load'] = self.name
                    self._beam.df = dfnew
                else:
                    self._beam[beamid] = values
            elif isinstance(values, (list, tuple)):
                if isinstance(values[0], (list, tuple, dict)):
                    for item in values:
                        if isinstance(item, dict):
                            beamid = item['beam']
                            load =  item
                        elif isinstance(item, (list, tuple)):
                            beamid = item[0]
                            load =  item[1:]
                        #
                        self._beam[beamid] = load
                else:
                    self._beam[values[0]] = values[1:]
        # dataframe input
        try:
            columns = list(df.columns)
            df['load'] = self.name
            self._beam.df = df
            return 
        except AttributeError:
            pass
        #
        return self._beam
    #    
    #
    # -----------------------------------------------
    #
    #

#
#