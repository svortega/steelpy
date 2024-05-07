#
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from typing import NamedTuple
from dataclasses import dataclass
import re
#from collections.abc import Mapping


# package imports
from steelpy.metocean.wave.utils.main import WaveBasic
from steelpy.metocean.wave.regular.main import RegularWaves
from steelpy.utils.sqlite.utils import create_table # create_connection, 

#
#
@dataclass
class Wave(WaveBasic):
    __slots__ = ['_regular', '_iregular', '_spectrum',
                 '_rho_w', '_labels','_type', '_db_file',
                 '_criteria']
    
    def __init__(self, criteria: str, rho_w:float, db_file: str):
        """
        """
        super().__init__(db_file)
        #
        self._criteria = criteria
        self._rho_w: float = rho_w  # kg / m^3
        #self._db_file = db_file
        #
        self._labels: list[str|int] = []
        self._type: list = []
        #
        self._regular = RegularWaves(criteria=self._criteria, 
                                     db_file=self.db_file)
        #self._iregular_wave = WaveIrregular()
        # self._spectrum = Sprectrum()
        #
    #
    #
    def __setitem__(self, name: int|str,
                    wave_data: list|tuple|dict) -> None:
        """ wave input data"""
        try:
            self._labels.index(name)
            raise Exception(f'Wave {name} already exist')
        except ValueError:
            wave_type = wave_data[0]
            #properties = get_sect_properties(properties[1:])
            self._labels.append(name)
            self._type.append(wave_type)
            #conn = create_connection(self.db_file)
            #with conn:
            #    # name, type, title, criteria_id, 
            #    wave_data = (name, wave_type, None, self._criteria,)
            #    self._push_data(conn, wave_data)
        #
        #wave_type = wave_data[0]
        #self._input(name, wave_type)
        #
        if re.match(r"\b(regular(_wave)?)\b", wave_type, re.IGNORECASE):
            self._regular[name] = wave_data
        elif re.match(r"\b(iregular(_wave)?)\b", wave_type, re.IGNORECASE):
            #self._iregular[name] = wave_data
            raise NotImplementedError
        else:
            raise ValueError("wave type invalid")

    def __getitem__(self, name: int|str):
        """
        node_name : node number
        """
        try:
            index = self._labels.index(name)
            wave_type = self._type[index]
        except ValueError:
            raise KeyError(f'   *** Wave {name} does not exist')
        #
        if re.match(r"\b(regular(_wave)?)\b", wave_type, re.IGNORECASE):
            return self._regular[name]
        elif re.match(r"\b(iregular(_wave)?)\b", wave_type, re.IGNORECASE):
            # return self._iregular[name]
            raise NotImplementedError
        else:
            raise ValueError("wave type invalid")
    #
    #def _input(self, name:int|str, wave_data:str):
    #    """ """
    #    try:
    #        self._labels.index(name)
    #        raise Exception(f'Wave {name} already exist')
    #    except ValueError:
    #        wave_type = wave_data
    #        self._labels.append(name)
    #        self._type.append(wave_type)
    #
    #
    #@property
    #def spectrum(self):
    #    """
    #    """
    #    return self._spectrum
    #
    #@property
    #def irregular(self):
    #    """
    #    """
    #    return self._iregular_wave
    #
    #@property
    def regular(self, values:None|list=None,
                df=None):
        """
        """
        if values:
            print('-->')
            1/0
        else:
            try:
                df.columns
                self._regular.df(df)
            except AttributeError:
                pass
        return self._regular
    #
    #
    # -----------------------------------
    # SQL ops
    #
    def _create_table(self, conn) -> None:
        """ """
        # Wave main
        table = "CREATE TABLE IF NOT EXISTS Wave (\
                number INTEGER PRIMARY KEY NOT NULL,\
                name NOT NULL,\
                type TEXT NOT NULL, \
                title TEXT, \
                criteria_id INTEGER NOT NULL REFERENCES Criteria(number));"
        create_table(conn, table)
    #
    #
    #def _push_data(self, conn, wave_data):
    #    """
    #    Create a new project into the projects table
    #    """
    #    #project = (wave_data[0], wave_data[1])
    #    table = 'INSERT INTO Wave(name, type, tile, criteria_id) \
    #             VALUES(?,?,?,?,?)'
    #    # push
    #    cur = conn.cursor()
    #    cur.execute(table, wave_data)
    #    wave_id = cur.lastrowid
    #    print('--')
    #
    # -----------------------------------
    # Operations
    #
    def solve(self, surface_points:int = 36,
              depth_points:int = 100):
        """ """
        print('--> regular solve')
        self._regular.solve(surface_points=surface_points,
                            depth_points=depth_points)
        #1 / 0
    #
    #
    #
#
#