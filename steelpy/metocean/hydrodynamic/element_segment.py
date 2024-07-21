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
from steelpy.metocean.hydrodynamic.utils import HydroItem, HydroBasic
#
import numpy as np


class ElementSegmentation(HydroBasic):
    """
    Wave integration points on element
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
            raise IOError(f'Element Segmentation {name} not found')
        #
        return WIPitem(name=name,
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
        table = "CREATE TABLE IF NOT EXISTS ElementSegment (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        name NOT NULL,\
                        type TEXT NOT NULL, \
                        title TEXT);"
        create_table(conn, table)
        # Profile
        table = "CREATE TABLE IF NOT EXISTS ElementSegmentProfile (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        element_segment_id NOT NULL REFERENCES ElementSegment(number),\
                        elevation DECIMAL NOT NULL,\
                        segment_number INTEGER NOT NULL);"
        create_table(conn, table)
    #
    def _push_data(self, conn, data):
        """
        Create a new project into the projects table
        """
        cur = conn.cursor()
        table = 'INSERT INTO ElementSegment(name, type, title) \
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
        sql = 'SELECT {:} FROM ElementSegment WHERE name = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
        data = cur.fetchone()
        return data
    #

#
#
@dataclass
class WIPitem(HydroItem):
    __slots__ = ['name', '_db_file']
    #
    def __init__(self, name: int|str,
                 db_file: str):
        """ """
        super().__init__(name=name, db_file=db_file)
    #
    def nelev(self, beam, up: str):
        """"""
        n1, n2 = beam.nodes
        n1 = getattr(n1, up)
        n2 = getattr(n2, up)
        zmax = np.maximum(n1, n2)
        zmin = np.minimum(n1, n2)
        point = [zmax, zmin]
        #
        prof = self.profile
        elev = list(reversed([item[0] for item in prof]))
        value = list(reversed([item[1] for item in prof]))
        #
        nelev = np.interp(point, elev, value)
        nelev = int(nelev.max())
        #
        section = beam.section.geometry
        bmax = int(beam.L // section.Dh)
        #
        return max(nelev, bmax)
    #
    def _push_profile(self, conn, profile_data):
        """ """
        mg_name = (self.name, )
        #
        table = f"UPDATE ElementSegment \
                 SET type = 'profile' \
                 WHERE name = ?"
        cur = conn.cursor()
        cur.execute(table, mg_name)
        #
        cur = conn.cursor()
        table = 'SELECT * FROM ElementSegment WHERE name = ?'
        cur.execute(table, mg_name)
        values = cur.fetchone()        
        #
        profile = tuple((values[0], *item, ) for item in profile_data)
        cur = conn.cursor()
        table = 'INSERT INTO ElementSegmentProfile(element_segment_id, \
                                                    elevation, segment_number) \
                                                    VALUES(?,?,?)'
        # push
        cur = conn.cursor()
        cur.executemany(table, profile)
    #
    def _pull_item(self, conn):
        """ """
        item_name = (self.name, )
        table = 'SELECT * FROM ElementSegment WHERE name = ?'
        cur = conn.cursor()
        cur.execute(table, item_name)
        item = cur.fetchone()
        return item
    #
    def _pull_profile(self, conn):
        """get profile data"""
        wip = self._pull_item(conn)
        #
        if re.match(r"\b(profile)\b", wip[2], re.IGNORECASE):
            wip_name = (wip[0], )
            cur = conn.cursor()
            table = 'SELECT * FROM ElementSegmentProfile \
                     WHERE element_segment_id = ?'
            cur.execute(table, wip_name)
            wip_profile = cur.fetchall()
            wip_profile = [item[2:] for item in wip_profile]
        else:
            1 / 0
        #
        #print('-->')
        return wip_profile