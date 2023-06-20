#
# Copyright (c) 2009-2023 fem2ufo
# 

# Python stdlib imports
from __future__ import annotations
from array import array
from collections.abc import Mapping
#from collections import defaultdict
from dataclasses import dataclass
#from operator import sub, add
from typing import NamedTuple
# import re
#
# package imports
import pandas as pd
# steelpy
#from ....process.math.operations import zeros, trns_3Dv #, zeros_vector , linspace
#from ....process.math.vector import Vector
#from steelpy.formulas.beam.main import BasicCalcs
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
class BasicLoadBasic(Mapping):

    def __init__(self):
        """
        """
        self._labels: list = []
        self._title: list[str] = []
        self._number: array = array("I", [])
        self.gravity = 9.80665  # m/s^2

    #
    def __len__(self) -> int:
        return len(self._labels)

    #
    def __iter__(self):
        """
        """
        return iter(self._labels)

    #
    @property
    def g(self):
        """"""
        return self.gravity

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
        #output += f"Element Number{4 * ' '}L1[{unit_lenght}] fx [{unit_force}] fy [{unit_force}] fz [{unit_force}] System Complex\n"
        output += "\n"
        output += "{:}\n".format(80 * ".")
        output += "\n"
        loadname = list(dict.fromkeys(self._labels))
        for key in loadname:
            lcase = self.__getitem__(key)
            output += f"Load Name : {str(lcase.name):12s}  Number : {lcase.number:8.0f}  Title : {lcase.title}\n"
            # node load
            if lcase._node:
                output += f"--- Node\n"
                output += lcase._node.__str__()
            # beam line
            if lcase._beam:
                output += f"--- Beam \n"
                output += lcase._beam.__str__()
            #
            if lcase._selfweight:
                output += f"--- Gravity/Selfweight\n"
                output += lcase._selfweight.__str__()
            #
            output += "\n"
        # print('---')
        return output
    #
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
    # ---------------------
    # Loading Operations
    # ---------------------
    #
    def process(self, elements, steps:int):
        """process element load"""
        #
        load_func = {}
        #dftemp = []
        for name in self._labels:
            lcase = self.__getitem__(name)
            #
            bfunc = {}
            #bfunc = []
            # -------------------------------
            # beam line load (line & point)
            # -------------------------------
            #
            beamload = lcase.beam()
            for key, items in beamload.items():
                try:
                    1/(len(items.line) + len(items.point))
                    bfunc.update(items.beam_function(steps=steps))
                except ZeroDivisionError:
                    continue
            #
            # -------------------------------
            # wave line load
            # -------------------------------
            #
            waveload = lcase.wave()
            waveload.beam_load()
            #for key, items in waveload.items():
                #try:
               #1/len(items.line())
                #bfunc.update(items.beam_function(steps=steps))
                #except ZeroDivisionError:
                #continue

            #
            # -------------------------------
            #
            load_func[name] = bfunc
            #
            # -------------------------------
            # plate load
            # -------------------------------
            #         
        #print('---')
        # [Fx, Fy, Fz, Mx, My, Mz]
        # [V, M, w, theta]
        #header = ['load_name', 'load_type', 'load_title', 'system',
        #          'element_name', 'node_end',
        #          'Fx_V', 'Fx_M', 'Fx_w', 'Fx_theta',
        #          'Fy_V', 'Fy_M', 'Fy_w', 'Fy_theta',
        #          'Fz_V', 'Fz_M', 'Fz_w', 'Fz_theta',
        #          'Mx_V', 'Mx_M', 'Mx_w', 'Mx_theta',
        #          'My_V', 'My_M', 'My_w', 'My_theta',
        #          'Mz_V', 'Mz_M', 'Mz_w', 'Mz_theta']
        #dfload = pd.DataFrame(data=dftemp, columns=header, index=None)
        return load_func
    #
    #
    def FER(self):
        """Convert element load to global fixed end reactions"""
        load_name = list(self.keys())
        for lname in load_name:
            lcase = self.__getitem__(lname)
            # -------------------------------
            # beam line load (line & point)
            # -------------------------------
            lcase._beam.fer()
            # -------------------------------
            # plates
            # -------------------------------
            #
            #
            # -------------------------------
            # wave
            # -------------------------------
            #
            lcase._wave.fer()
            #
        #
        print('---> element2node')
    #
    #
    def node_df(self): 
        """
        Return beams' end node load in dataframe form
        """
        #1 / 0
        load_name = list(self.keys())
        #
        def update_df(dfval, dfin, number: int, title: str):
            """ """
            # update dataframe
            dfin['load_number'] = number
            dfin['load_title'] = title
            #
            try:
                dfval = pd.concat([dfval, dfin], ignore_index=True, sort=False)
            except UnboundLocalError:
                dfval = dfin
            return dfval
        #
        dftemp = None
        for lname in load_name:
            lcase = self.__getitem__(lname)
            # -------------------------------
            # nodal load
            # 
            df2 = lcase._node.load.df
            if not df2.empty:
                dftemp = update_df(dftemp, df2, lcase.number, lcase.title)
            #
            # -------------------------------
            # beam line load (line & point)
            #
            df2 = lcase._beam._node_eq.df
            if not df2.empty:
                dftemp = update_df(dftemp, df2, lcase.number, lcase.title)
            #
            # -------------------------------
            # wave load (line & point)  
            #
            df2 = lcase._wave._node_eq.df
            if not df2.empty:
                dftemp = update_df(dftemp, df2, lcase.number, lcase.title)
        #
        dftemp = dftemp.reindex(columns=['load_name', 'load_number', 'load_type',
                                         'load_title', 'load_comment', 'load_system',
                                         'element_name', 'node_name',
                                         'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'])
        #print('-->??')
        return dftemp
    #    
    #
    def wave_process(self):
        """ Wave los postprocessing"""
        print('# --> wave2beamload')
        load_name = list(self.keys())
        for lname in load_name:
            lcase = self.__getitem__(lname)
            try:
                lcase._wave.process(lname)
            except KeyError:
                continue
        print('# --> End')
    #
#
#
class LoadTypeBasic:
        
    def __init__(self, name: str | int, number: int, title: str):
        """
        """
        self.name = name
        self.number = number
        self.title = title
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
    #
    #@property
    #def node(self):
    #    """ """
    #    return self._node

    #@node.setter
    def node(self, values:list|None=None):
        """ """
        #
        if isinstance(values, list):
            if isinstance(values[0], list):
                for value in values:
                    self._node[value[0]] = value[1:]
            else:
                self._node[values[0]] = values[1:]
        #
        return self._node

    #
    #@beam.setter
    def beam(self, values:list|None=None,
             df=None):
        """ """
        if isinstance(values, list):
            if isinstance(values[0], list):
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
        #
        # dataframe input
        try:
            df.columns
            #1/0
            self._beam.df(df)
        except AttributeError:
            pass
        #
        return self._beam
    #
    #@property
    #def wave(self):
    #    """ """
    #    return self._wave
    #
    #@wave.setter
    def wave(self, wave_load=None, design_load: str = 'max_BSOTM'):
        """
        design_load : max_BSOTM
        """
        if wave_load:
            self._wave[self.name] = [wave_load, design_load, self.title]
        #
        return self._wave
# 
#