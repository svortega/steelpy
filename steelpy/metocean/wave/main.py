#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations

# Python stdlib imports
#from typing import NamedTuple
#from dataclasses import dataclass
import re
from collections.abc import Mapping
#from math import tau
#
from steelpy.metocean.wave.regular.main import RegularWaves
from steelpy.utils.sqlite.utils import create_connection, create_table

#
#
class WaveMain(Mapping):
    __slots__ = ['_regular', '_iregular', '_spectrum',
                 '_rho_w', '_labels','_type', '_db_file']
    
    def __init__(self, rho_w:float, db_file: str):
        """
        """
        self._db_file = db_file
        #
        self._labels: list[str|int] = []
        self._type: list = []
        self._rho_w: float = rho_w  # kg / m^3
        #
        self._regular = RegularWaves(db_file=self._db_file)
        #self._iregular_wave = WaveIrregular()
        # self._spectrum = Sprectrum()
        # create node table
        conn = create_connection(self._db_file)
        with conn:
            self._create_table(conn)        
    #
    #
    def __setitem__(self, name: int|str,
                    wave_data: list|tuple|dict) -> None:
        #try:
        #    self._labels.index(name)
        #    raise Exception(f'Wave {name} already exist')
        #except ValueError:
        #    wave_type = wave_data[0]
        #    #properties = get_sect_properties(properties[1:])
        #    self._labels.append(name)
        #    self._type.append(wave_type)
        wave_type = wave_data[0]
        self._input(name, wave_type)
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
    def _input(self, name:int|str, wave_data:str):
        """ """
        try:
            self._labels.index(name)
            raise Exception(f'Wave {name} already exist')
        except ValueError:
            wave_type = wave_data
            self._labels.append(name)
            self._type.append(wave_type)
    #
    def __len__(self):
        return len(self._labels)

    def __iter__(self):
        return iter(self._labels)

    def __contains__(self, value):
        return value in self._labels
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
        table = "CREATE TABLE IF NOT EXISTS tb_Wave (\
                number INTEGER PRIMARY KEY NOT NULL,\
                name NOT NULL,\
                type TEXT NOT NULL);"
        create_table(conn, table)
    #
    #
    # -----------------------------------
    # Operations
    #
    def solve(self, surface_points:int = 36,
              depth_points:int = 100):
        """ """
        #print('--> regular')
        self._regular.solve(surface_points=surface_points,
                            depth_points=depth_points)
        #1 / 0
    #
    #
    #
#
#