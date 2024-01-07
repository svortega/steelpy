# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
#from operator import itemgetter
import re

#
# package imports
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.metocean.hydrodynamic.utils.main import HydroItem, HydroBasic


class WaveKinFactor(HydroBasic):
    """
    """
    __slots__ = ['db_file']
    
    def __init__(self, db_file: str):
        """
        """
        super().__init__(db_file)
    #
    # ------------------
    #
    def __getitem__(self, name: str|int):
        """
        """
        #
        conn = create_connection(self.db_file)
        with conn:
            data = self._pull_data(conn, name=name)        
        #
        if not data:
            raise IOError(f'WKF {name} not found')
        #
        return WKFitem(name=name,
                      db_file=self.db_file)
    #
    #
    # ------------------
    # SQL ops
    # ------------------
    #
    def _create_table(self, conn) -> None:
        """ """
        # Main
        table = "CREATE TABLE IF NOT EXISTS WaveKinFactor (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        name NOT NULL,\
                        type TEXT NOT NULL, \
                        title TEXT);"
        create_table(conn, table)
        # Profile
        table = "CREATE TABLE IF NOT EXISTS WaveKinFactorProfile (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        wkf_id NOT NULL REFERENCES WaveKinFactor(number),\
                        elevation DECIMAL NOT NULL,\
                        factor DECIMAL NOT NULL);"
        create_table(conn, table)
    #
    def _push_data(self, conn, data):
        """
        Create a new project into the projects table
        """
        cur = conn.cursor()
        table = 'INSERT INTO WaveKinFactor(name, type, title) \
                 VALUES(?,?,?)'
        # push
        cur = conn.cursor()
        cur.execute(table, data)
    #
    def _pull_data(self, conn, name: str|int,
                  item: str = "*"):
        """ """
        #
        project = (name,)
        sql = 'SELECT {:} FROM WaveKinFactor WHERE name = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
        data = cur.fetchone()
        return data
    #

#
#
@dataclass
class WKFitem(HydroItem):
    __slots__ = ['name', '_db_file']
    #
    def __init__(self, name: int|str,
                 db_file: str):
        """ """
        super().__init__(name=name, db_file=db_file)
    #
    #
    def _push_profile(self, conn, profile_data):
        """ """
        mg_name = (self.name, )
        #
        table = f"UPDATE WaveKinFactor \
                 SET type = 'profile' \
                 WHERE name = ?"
        cur = conn.cursor()
        cur.execute(table, mg_name)
        #
        cur = conn.cursor()
        table = 'SELECT * FROM WaveKinFactor WHERE name = ?'
        cur.execute(table, mg_name)
        values = cur.fetchone()        
        #
        profile = tuple((values[0], *item, ) for item in profile_data)
        cur = conn.cursor()
        table = 'INSERT INTO WaveKinFactorProfile(wkf_id, \
                                              elevation, factor) \
                                              VALUES(?,?,?)'
        # push
        cur = conn.cursor()
        cur.executemany(table, profile)
    #