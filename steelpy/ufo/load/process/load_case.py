#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
from collections.abc import Mapping
#from collections import defaultdict
#from dataclasses import dataclass
#from operator import sub, add
#import time
from typing import NamedTuple
import re
#
# package imports
import steelpy.utils.io_module.text as common
#from steelpy.ufo.load.process.utils import find_load_type
from steelpy.utils.dataframe.main import DBframework
#
from steelpy.ufo.load.process.beam.utils import find_BeamLoad_item
from steelpy.ufo.load.process.node.utils import find_NodeLoad_item
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
class BasicLoadRoot(Mapping):
    __slots__ = ['name', 'gravity']
    
    def __init__(self, name:str|int):
        """
        """
        self.name = name
        self.gravity = 9.80665  # m/s^2
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
    # -----------------------------------------------
    #
    @property
    def g(self):
        """"""
        return self.gravity

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
    # -----------------------------------------------
    #    
    #   
#
#
class BasicLoadCase(BasicLoadRoot):
    __slots__ = ['name','gravity']
    
    def __init__(self, name:int|str):
        """
        """
        super().__init__(name)
        #self.gravity = 9.80665  # m/s^2
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
                columns = list(values.keys())
                header = {item: find_NodeLoad_item(item)
                          for item in columns}
                values ={header[key]: item
                         for key, item in values.items()}
                #
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
                        load_title = values['title']
                    except KeyError:
                        load_title = values['load']
                    #
                    try:
                        self.__setitem__(load_name, load_title)
                    except Warning:
                        pass
                    #
                    bload = self.__getitem__(load_name)
                    bload._node[nodeid] = values
            
            elif isinstance(values, (list, tuple)):
                if isinstance(values[0], (list, tuple, dict)):
                    for item in values:
                        if isinstance(item, dict):
                            header = {item: find_NodeLoad_item(item)
                                      for item in item}
                            update = {header[key]: item
                                      for key, item in item.items()}
                            load_name = update['load']
                            try:
                                load_title = update['title']
                            except KeyError:
                                load_title = update['load']
                            node_id = update['node']
                            load = update
                        elif isinstance(item, (list, tuple)):
                            load_name = item[0]
                            if isinstance(item[-1], str):
                                load_title = item[-1]
                            else:
                                load_title = load_name
                            node_id = item[1]
                            load = item[2:]
                        #
                        try:
                            self.__setitem__(load_name, load_title)
                        except Warning:
                            pass
                        #
                        bload = self.__getitem__(load_name)
                        bload._node[node_id] = load
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
            header = {key: find_NodeLoad_item(key)
                      for key in columns}
            df.rename(columns=header, inplace=True)
            #nodeid = df['node']
            #bitems = len(nodeid)
            #values = self._set_load(df, steps=bitems)
            #
            columns = list(df.columns)
            #if not 'title' in columns:
            #    df['title'] = df['load']
            #
            # Set basic load if doesn't exist
            basic = df.groupby(['load'])
            for x, load_name in enumerate(basic.groups):
                # TODO : title
                load_title = f"{load_name}_{x + 1}"
                try:
                    self.__setitem__(load_name, load_title)
                except Warning:
                    pass
                #
                bload = self.__getitem__(load_name)
                nload = basic.get_group((load_name,))
                bload._node.df = nload
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
                columns = list(values.keys())
                header = {item: find_BeamLoad_item(item)
                          for item in columns}
                values ={header[key]: item
                         for key, item in values.items()}
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
                        load_title = values['title']
                    except KeyError:
                        load_title = values['load']
                    try:
                        self.__setitem__(load_name, load_title)
                    except Warning:
                        pass
                    bload = self.__getitem__(load_name)                    
                    bload._beam[beamid] = values
            elif isinstance(values, (list, tuple)):
                if isinstance(values[0], (list, tuple, dict)):
                    for item in values:
                        if isinstance(item, dict):
                            header = {item: find_BeamLoad_item(item)
                                        for item in item}
                            update = {header[key]: item
                                      for key, item in item.items()}
                            load_name = update['load']
                            try:
                                load_title = update['title']
                            except KeyError:
                                load_title = update['load']
                            beam_id = update['beam']
                            load =  update
                        elif isinstance(item, (list, tuple)):
                            # fix load title
                            load_name = item[0]
                            if isinstance(item[-1], str):
                                load_title = item[-1]
                            else:
                                load_title = load_name
                            beam_id = item[1]
                            load =  item[2:]
                        # check if basic load name exist
                        try:
                            self.__setitem__(load_name, load_title)
                        except Warning:
                            pass
                        # push beam load data
                        bload = self.__getitem__(load_name)
                        bload._beam[beam_id] = load
                else:
                    load_name = values[0]
                    if isinstance(values[-1], str):
                        load_title = values[-1]
                    else:
                        load_title = load_name
                    try:
                        self.__setitem__(load_name, load_title)
                    except Warning:
                        pass
                    bload = self.__getitem__(load_name)
                    beamid = values[1]
                    bload._beam[beamid] = values[2:]
        #
        # dataframe input
        try:
            columns = list(df.columns)
            header = {key: find_bload_item(key)
                      for key in columns}
            df.rename(columns=header, inplace=True)            
            #beamid = df['beam']
            #bitems = len(beamid)
            #values = self._set_load(df, steps=bitems)
            # Set basic load if doesn't exist
            basic = df.groupby(['load'])
            for x, load_name in enumerate(basic.groups):
                # TODO : title
                load_title = f"{load_name}_{x + 1}"
                try:
                    self.__setitem__(load_name, load_title)
                except Warning:
                    pass
                #
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
    # -----------------------------------------------
    #
    #
    #def functionX(self, steps:int,
    #             Pa:float=0.0, factor:float=1):
    #    """utils element load"""
    #    #
    #    dftemp = self._beams.load_function(steps=steps,
    #                                       Pa=Pa, factor=factor)
    #    #
    #    # Axial   [FP, blank, blank, Fu]
    #    # torsion [T, B, Psi, Phi, Tw]
    #    # Bending [V, M, theta, w]
    #    #
    #    # [Fx, Fy, Fz, Mx, My, Mz]
    #    # [V, M, w, theta]
    #    header = ['load_name', 'mesh_name',
    #              'load_comment', 'load_type',
    #              'load_level', 'load_system',
    #              'element_name', 'length',
    #              'axial', 'torsion', 'VM_inplane', 'VM_outplane']
    #    #
    #    #          'FP', 'blank1', 'blank2', 'Fu',
    #    #          'T', 'B', 'Psi', 'Phi', 'Tw',
    #    #          'Vy', 'Mz', 'theta_y', 'w_y',
    #    #          'Vz', 'My', 'theta_z', 'w_z']
    #    df = DBframework()
    #    dfload = df.DataFrame(data=dftemp, columns=header, index=None)
    #    #return load_func
    #    return dfload
    #
    #
    # -----------------------------------------------
    #
    def Pn(self):
        """
        Returns:
            Global nodal force vector
        """
        P = self._nodes.force
        return P
    #
    # -----------------------------------------------
    #
    @property
    def df(self):
        """ """
        1 / 0
    
    @df.setter
    def df(self, df):
        """ """
        columns = list(df.columns)
        header = {key: find_bload_item(key)
                  for key in columns}
        df.rename(columns=header, inplace=True)
        #
        # input node
        try:
            grpnode = df.groupby('node')
            for key, item in grpnode:
                print(key)
                self.node(df=item)
        except KeyError:
            pass        
        #
        # input beam
        try:
            grpbeam = df.groupby('beam')
            for key, item in grpbeam:
                print(key)
                self.beam(df=item)
        except KeyError:
            pass
# 
#  
#
#
#
# ---------------------------------
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
# ---------------------------------
#
def find_bload_item(word_in:str) -> str:
    """ """
    key = {"load": r"\b((basic)?(_|-|\s*)?(load)?(_|-|\s*)?name)\b",
           "node": r"\b(node(s)?|point(s)?)\b",
           "beam": r"\b(beam(s)?)\b",
           #"factor": r"\b((load)?(_|-|\s*)?factor)\b",
           "title": r"\b(title|comment)\b",}
    try:
        match = common.find_keyword(word_in, key)
        return match
    except IOError:
        return word_in
#
#
def get_bload_list(values: list|tuple, steps: int = 3) -> list:
    """ [name, item, item_id, values] """
    output = [None] * steps
    for x, item in enumerate(values):
        try:
            output[x] = item
        except IndexError:
            pass
    #
    output.extend(values[steps:])
    return output
#
#
def get_bload_dict(values: dict, steps: int = 3) -> list:
    """ [name, item, item_id, values] """
    output = [None] * steps
    valtemp = values.copy()
    for key, item in values.items():
        if re.match(r"\b((basic)?(_|-|\s*)?(load)?(_|-|\s*)?name)\b", key, re.IGNORECASE):
            output[0] = item
            del valtemp[key]

        elif re.match(r"\b(node(s)?|point(s)?)\b", key, re.IGNORECASE):
            output[1] = 'node'
            output[2] = item
            del valtemp[key]
            #valtemp
            
        elif re.match(r"\b(beam(s)?)\b", key, re.IGNORECASE):
            output[1] = 'beam'
            output[2] = item
            del valtemp[key]
            #valtemp            
            
        #elif re.match(r"\b((load)?(_|-|\s*)?factor)\b", key, re.IGNORECASE):
        #    output[3] = item
        #    
        #elif re.match(r"\b(title|comment)\b", key, re.IGNORECASE):
        #    output[4] = item
    #
    for key, items in valtemp.items():
        output.extend([key, *items])
    return output
#
#