# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass

# package imports
from steelpy.utils.sqlite.main import ClassBasicSQL
#from steelpy.process.units.units import Units
#from steelpy.utils.sqlite.utils import create_connection, create_table



#
@dataclass
class WaveBasic(ClassBasicSQL):
    """ """
    #__slots__ = ['db_file']
    #
    def __init__(self, db_file: str):
        """
        """
        super().__init__(db_file)
    #
    #
    def _push_wave(self, conn, wave_data):
        """ get wave data"""
        #
        project = (*wave_data, None, self._criteria)
        table = 'INSERT INTO Wave(name, type, title, criteria_id) \
                 VALUES(?,?,?,?)'
        #
        #push
        cur = conn.cursor()
        cur.execute(table, project)
        wave_id = cur.lastrowid
        #
        #project = (self._criteria, )
        #table = "SELECT * FROM Wave \
        #         WHERE name = ? AND type = 'regular'"
        #cur = conn.cursor()
        #cur.execute(table, project)
        #wave = cur.fetchone()
        #wave
        #1 / 0
        return wave_id
    #
    def _pull_wave(self, conn, wave_name):
        """ """
        project = (wave_name, self._criteria)
        table = "SELECT * FROM Wave \
                 WHERE name = ? AND criteria_id = ? \
                 AND type = 'regular'"
        cur = conn.cursor()
        cur.execute(table, project)
        wave = cur.fetchone()
        return wave
#