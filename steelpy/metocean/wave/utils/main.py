# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass

# package imports
from steelpy.utils.sqlite.main import ClassBasicSQL
#from steelpy.utils.sqlite.utils import create_connection, create_table



#
@dataclass
class WaveBasic(ClassBasicSQL):
    """ """
    __slots__ = ['db_file', '_component']
    #
    def __init__(self, component: str|int, db_file: str):
        """
        """
        super().__init__(db_file)
        self._component = component
    #
    #
    def _push_wave(self, conn, wave_data: list) -> int:
        """
        wave_data: [name, type, title]
        """
        project = (*wave_data, self._criteria)
        table = 'INSERT INTO Wave(name, type, title, criteria_id) \
                 VALUES(?,?,?,?)'
        #
        cur = conn.cursor()
        cur.execute(table, project)
        wave_id = cur.lastrowid
        return wave_id
    #
    def _pull_wave(self, conn, wave_name: str|int,
                   wave_type: str) -> tuple:
        """ """
        project = (wave_name, wave_type, self._criteria)
        table = "SELECT * FROM Wave \
                 WHERE name = ? \
                 AND type = ?\
                 AND criteria_id = ?"
        cur = conn.cursor()
        cur.execute(table, project)
        wave = cur.fetchone()
        return wave
#