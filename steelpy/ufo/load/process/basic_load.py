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
# TODO: wrape pandas
import pandas as pd
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
    # -----------------------------------------------
    # Loading Operations
    # -----------------------------------------------
    #
    def process(self, elements, steps:int):
        """process element load"""
        #
        beams = elements.beam()
        #
        #beamload = self.beam()
        #dftemp = beamload.beam_function(beams=beams, steps=steps)
        #
        dftemp = []
        ##load_func = {}
        for name in self._labels:
            lcase = self.__getitem__(name)
            #lcase =  self._basic[name]
        #    #
        #    #bfunc = {}
        #    #bfunc = []
        #    # -------------------------------
        #    # beam line load (line & point)
        #    # -------------------------------
        #    #
            beamload = lcase.beam()
            for key, items in beamload.items():
                try:
                    1/(len(items.line) + len(items.point))
                    dftemp.extend(items.beam_function(beams=beams, steps=steps))
                    #bfunc.update(items.beam_function(beams=beams,
                    #                                 steps=steps))
                except ZeroDivisionError:
                    continue
        #    #
        #    # -------------------------------
        #    # wave line load
        #    # -------------------------------
        #    #
        #    waveload = lcase.wave()
        #    #bfunc.update(waveload.beam_load(steps=steps))
        #    dftemp.extend(waveload.beam_load(steps=steps))
        #    #wload = waveload.beam_load(steps=steps)
        #    #if wload:
        #    #    bfunc.update(wload)
        #    #
        #    # -------------------------------
        #    #
        #    #load_func[name] = bfunc
        #    #
        #    # -------------------------------
        #    # plate load
        #    # -------------------------------
        #    #
        #
        #
        #print('---')
        #
        # Axial   [FP, blank, blank, Fu]
        # torsion [T, B, Psi, Phi, Tw]
        # Bending [V, M, theta, w]        
        #
        # [Fx, Fy, Fz, Mx, My, Mz]
        # [V, M, w, theta]
        header = ['load_name', 'component_name',
                  'load_comment', 'load_type', 
                  'load_level', 'load_system',
                  'element_name', 'node_end',
                  'axial', 'torsion', 'VM_inplane', 'VM_outplane']
        #
        #          'FP', 'blank1', 'blank2', 'Fu',
        #          'T', 'B', 'Psi', 'Phi', 'Tw',
        #          'Vy', 'Mz', 'theta_y', 'w_y',
        #          'Vz', 'My', 'theta_z', 'w_z']
        dfload = pd.DataFrame(data=dftemp, columns=header, index=None)
        #return load_func
        return dfload
    #
    #
    def FER(self, elements):
        """Convert element load to global fixed end reactions"""
        #
        beams = elements.beam()
        load_name = list(self.keys())
        for lname in load_name:
            lcase = self.__getitem__(lname)
            #lcase =  self._basic[lname]
            # -------------------------------
            # beam line load (line & point)
            # -------------------------------
            lcase._beam.fer(beams=beams)
            # -------------------------------
            # plates
            # -------------------------------
            #
            #
            # -------------------------------
            # wave
            # -------------------------------
            # FIXME: broken 
            #lcase._wave.fer()
            #
        #
        print('---> element2node')
    #
    #
    def Fn(self):
        """
        Global matrix consisting of summation of force & displacement 
        """
        dfnodal = self._nodes.df
        dfnodal = dfnodal.reindex(columns=['load_name', 'component_name', 
                                           'load_id', 'load_level',
                                          'load_title', 'load_comment', 'load_system',
                                          'element_name', 'node_name',
                                          'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                                          'x', 'y', 'z', 'rx', 'ry', 'rz'])
        # Select FER's force & displacement
        # [Fx, Fy, Fz, Mx, My, Mz] and [x, y, z, rx, ry, rz]
        columns = [*self._plane.hforce, *self._plane.hdisp]
        # Sum total 
        dfnodal = (dfnodal.groupby(['load_name', 'component_name',
                                    'load_id', 'load_level',
                                    'load_title','load_system',
                                    'node_name'])
                   [columns].sum())
        dfnodal.reset_index(inplace=True)
        return dfnodal
        
    #
    #def node_df(self): 
    #    """
    #    Return beams' end node load in dataframe form
    #    """
    #    #
    #    dftemp = self._nodes.df
        #
        #1 / 0
        #load_name = list(self.keys())
        #
        #def update_df(dfval, dfin, number: int, title: str):
        #    """ """
        #    # update dataframe
        #    dfin['load_id'] = number
        #    dfin['load_title'] = title
        #    #
        #    try:
        #        dfval = pd.concat([dfval, dfin], ignore_index=True, sort=False)
        #    except UnboundLocalError:
        #        dfval = dfin
        #    return dfval
        #
        #dftemp = None
        #for lname in load_name:
        #    lcase = self.__getitem__(lname)
        #    # -------------------------------
        #    # nodal load
        #    # 
        #    df2 = lcase._node.df
        #    if not df2.empty:
        #        dftemp = update_df(dftemp, df2, lcase.number, lcase.title)
        #    #
        #    # -------------------------------
        #    # beam line load (line & point)
        #    #
        #    #df2 = lcase._beam._node_eq.df
        #    #if not df2.empty:
        #    #    dftemp = update_df(dftemp, df2, lcase.number, lcase.title)
        #    #
        #    # -------------------------------
        #    # wave load (line & point)  
        #    #
        #    #df2 = lcase._wave._node_eq.df
        #    #if not df2.empty:
        #    #    dftemp = update_df(dftemp, df2, lcase.number, lcase.title)
        #
        #dftemp = dftemp.reindex(columns=['load_name', 'load_id', 'load_level',
        #                                 'load_title', 'load_comment', 'load_system',
        #                                 'element_name', 'node_name',
        #                                 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'])
        #print('-->??')
        #1 / 0
        #return dftemp
    #    
    #
    def wave_process(self):
        """ Wave los postprocessing"""
        # FIXME: wave broken
        1 / 0
        print('# --> wave2beamload')
        load_name = list(self.keys())
        for lname in load_name:
            lcase = self.__getitem__(lname)
            #try:
            lcase._wave.process(lname)
            #except (KeyError, AttributeError, IndexError):
            #    continue
            #
            #except AttributeError:
            #    break            
        print('# --> End')
    #    
#
#
class LoadCaseBasic(BasicLoadMain):
    #__slots__ = ['_labels', '_title','_number', 'gravity']
    
    def __init__(self):
        """
        """
        #self._labels: list = []
        #self._title: list[str] = []
        #self._number: array = array("I", [])
        self.gravity = 9.80665  # m/s^2
        #1 / 0
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
    #def __delitem__(self, load_name: str|int):
    #    """
    #    """
    #    load_name = str(load_name)
    #    del self._basic[load_name]    
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
    #
    def node(self, values:tuple|list|None=None,
             df=None):
        """ """
        # Input data for specific basic node load 
        if isinstance(values, (list, tuple)):
            #1 / 0
            if isinstance(values[0], (list, tuple)):
                for item in values:
                    load_name = item[0]
                    try:
                        self.__setitem__(load_name, load_name)
                    except Warning:
                        pass
                    #
                    bload = self.__getitem__(load_name)
                    nodeid = item[1]
                    bload._node[nodeid] = item[2:]
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
            header = {}
            for key in columns:
                if re.match(r"\b(id|name|load(s)?)\b", key, re.IGNORECASE):
                    header[key] = 'name'
                    
                elif re.match(r"\b(type)\b", key, re.IGNORECASE):
                    header[key] = 'type'
                
                elif re.match(r"\b(title|comment)\b", key, re.IGNORECASE):
                    header[key] = 'title'
                #
                elif re.match(r"\b(node(s)?(\_)?(name|id)?)\b", key, re.IGNORECASE):
                    header[key] = 'node'                 
                #
                # Load 
                header = read_point(key, header)
                #
                # Displacement/mass
                header = read_dispmass(key, header)
            #
            node_load = df[header.keys()].copy()
            node_load.rename(columns=header, inplace=True)            
            #
            # Set basic load if doesn't exist
            basic = node_load.groupby(['name'])
            for x, load_name in enumerate(basic.groups):
                #load_title = f"{load_name}_{x + 1}"
                try:
                    self.__setitem__(load_name, load_name)
                except Warning:
                    pass
                #
                bload = self.__getitem__(load_name)
                nload = basic.get_group(load_name)
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
    def beam(self, values:tuple|list|None=None,
             df=None):
        """ """
        if isinstance(values, (list, tuple)):
            #1 / 0
            if isinstance(values[0], (list, tuple)):
                for item in values:
                    load_name = item[0]
                    try:
                        self.__setitem__(load_name, load_name)
                    except Warning:
                        pass
                    #
                    bload = self.__getitem__(load_name)
                    beamid = item[1]
                    bload._beam[beamid] = item[2:]
            else:
                load_name = values[0]
                try:
                    self.__setitem__(load_name, load_name)
                except Warning:
                    pass
                #
                bload = self.__getitem__(load_name)
                beamid = values[1]
                bload._beam[beamid] = values[2:]
        #
        # dataframe input
        try:
            columns = list(df.columns)
            header = {}
            for key in columns:
                if re.match(r"\b(id|name|load(s)?)\b", key, re.IGNORECASE):
                    header[key] = 'name'
                    
                elif re.match(r"\b(type)\b", key, re.IGNORECASE):
                    header[key] = 'type'
                
                elif re.match(r"\b(title|comment)\b", key, re.IGNORECASE):
                    header[key] = 'title'
                #
                elif re.match(r"\b(beam(s)?(\_)?(name|id)?)\b", key, re.IGNORECASE):
                    header[key] = 'beam'
                #
                # Load
                #
                elif re.match(r"\b(a|L(0|p))\b", key, re.IGNORECASE):
                    header[key] = 'L0'
                
                elif re.match(r"\b(b|L1)\b", key, re.IGNORECASE):
                    header[key] = 'L1'
                #
                # Point
                #
                header = read_point(key, header)
                #
                # Line
                #
                header = read_line(key, header)
            #
            beam_load = df[header.keys()].copy()
            beam_load.rename(columns=header, inplace=True)
            #
            #
            # Set basic load if doesn't exist
            basic = beam_load.groupby(['name'])
            for x, load_name in enumerate(basic.groups):
                load_title = f"{load_name}_{x + 1}"
                try:
                    self.__setitem__(load_name, load_title)
                except Warning:
                    pass
                #
                bload = self.__getitem__(load_name)
                loadx = basic.get_group(load_name).copy()
                bload._beam.df = loadx
                #
            #
            return 
        except AttributeError:
            pass
        #
        #print('beam load')
        #1 / 0
        return self._beams
    #    
    #

#
#
class LoadTypeBasic(BasicLoadMain):
    __slots__ = ['name', 'title', 'number', 
                 '_node', '_beam', '_selfweight']
    #def __init__(self, name: str | int, number: int, title: str):
    #    """
    #    """
    #    self.name = name
    #    self.number = number
    #    self.title = title
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
    def node(self, values:tuple|list|None=None,
             df=None):
        """ Nodal load"""
        # Input data for specific basic node load 
        if isinstance(values, (list, tuple)):
            try:
                nodal = self._node
            except AttributeError:
                raise IOError('Basic Load name is required')
            #
            if isinstance(values[0], (list, tuple)):
                for value in values:
                    #load = [self.name] + value[1:]
                    nodal[value[0]] = value[1:]
            else:
                #load = [self.name] + values[1:]
                nodal[values[0]] = values[1:]
        #
        #
        # dataframe input
        try:
            columns = list(df.columns)
            1 / 0
        except AttributeError:
            pass
        #
        #
        return self._node
    #
    #
    def beam(self, values:tuple|list|None=None,
             df=None):
        """ beam loading """
        if isinstance(values, (list, tuple)):
            if isinstance(values[0], (list, tuple)):
                for value in values:
                    #try:
                    #    self._f2u_beams[value[0]]
                    #except KeyError:
                    #    raise IOError(f"beam {value[0]} not found")
                    #
                    self._beam[value[0]] = value[1:]
            else:
                #try:
                #    self._f2u_beams[values[0]]
                #except KeyError:
                #    raise IOError(f"beam {values[0]} not found")
                #
                self._beam[values[0]] = values[1:]
        # dataframe input
        try:
            columns = list(df.columns)
            header = {}
            for key in columns:
                if re.match(r"\b(id|name|load(s)?)\b", key, re.IGNORECASE):
                    header[key] = 'name'
                    
                elif re.match(r"\b(type)\b", key, re.IGNORECASE):
                    header[key] = 'type'
                
                elif re.match(r"\b(title|comment)\b", key, re.IGNORECASE):
                    header[key] = 'title'
                #
                elif re.match(r"\b(beam(s)?(\_)?(name|id)?)\b", key, re.IGNORECASE):
                    header[key] = 'beam'
                #
                # Load
                #
                elif re.match(r"\b(a|L(0|p))\b", key, re.IGNORECASE):
                    header[key] = 'L0'
                
                elif re.match(r"\b(b|L1)\b", key, re.IGNORECASE):
                    header[key] = 'L1'
                #
                # Point
                #
                header = read_point(key, header)
                #
                # Line
                #
                header = read_line(key, header)
            #
            beam_load = df[header.keys()].copy()
            beam_load.rename(columns=header, inplace=True)
            #
            #
            # Set basic load if doesn't exist
            basic = beam_load.groupby(['name'])
            for x, load_name in enumerate(basic.groups):
                load_title = f"{load_name}_{x + 1}"
                try:
                    self.__setitem__(load_name, load_title)
                except Warning:
                    pass
                #
                bload = self.__getitem__(load_name)
                loadx = basic.get_group(load_name).copy()
                bload._beam.df = loadx
                #
            #
            return 
        except AttributeError:
            pass
        #
        #print('beam load')
        #1 / 0
        return self._beam
    #    
    #
    # -----------------------------------------------
    #

#
#def read_main(key, header):
#    """ """
#    pass
#
def read_point(key, header):
    """ """
    if re.match(r"\b(fx)\b", key, re.IGNORECASE):
        header[key] = 'fx'
    
    elif re.match(r"\b(fy)\b", key, re.IGNORECASE):
        header[key] = 'fy'
    
    elif re.match(r"\b(fz)\b", key, re.IGNORECASE):
        header[key] = 'fz'
    
    elif re.match(r"\b(mx)\b", key, re.IGNORECASE):
        header[key] = 'mx'
    
    elif re.match(r"\b(my)\b", key, re.IGNORECASE):
        header[key] = 'my'
    
    elif re.match(r"\b(mz)\b", key, re.IGNORECASE):
        header[key] = 'mz'
    #
    return header
#
def read_dispmass(key, header):
    """ """
    if re.match(r"\b(x)\b", key, re.IGNORECASE):
        header[key] = 'x'
    
    elif re.match(r"\b(y)\b", key, re.IGNORECASE):
        header[key] = 'y'
    
    elif re.match(r"\b(z)\b", key, re.IGNORECASE):
        header[key] = 'z'
    
    elif re.match(r"\b(rx)\b", key, re.IGNORECASE):
        header[key] = 'rx'
    
    elif re.match(r"\b(ry)\b", key, re.IGNORECASE):
        header[key] = 'ry'
    
    elif re.match(r"\b(rz)\b", key, re.IGNORECASE):
        header[key] = 'rz'
    #
    return header
#
def read_line(key, header):
    """ """
    if re.match(r"\b(qx0)\b", key, re.IGNORECASE):
        header[key] = 'qx0'
    
    elif re.match(r"\b(qy0)\b", key, re.IGNORECASE):
        header[key] = 'qy0'
    
    elif re.match(r"\b(qz0)\b", key, re.IGNORECASE):
        header[key] = 'qz0'
    
    elif re.match(r"\b(qx1)\b", key, re.IGNORECASE):
        header[key] = 'qx1'
    
    elif re.match(r"\b(qy1)\b", key, re.IGNORECASE):
        header[key] = 'qy1'
    
    elif re.match(r"\b(qz1)\b", key, re.IGNORECASE):
        header[key] = 'qz1'
    #
    return header
#
#