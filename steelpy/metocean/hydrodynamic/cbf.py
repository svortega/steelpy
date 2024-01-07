# 
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
#from operator import itemgetter

#
# package imports
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.metocean.hydrodynamic.utils.main import HydroItem, HydroBasic


class CurrentBlockFactor(HydroBasic):
    """
    """
    __slots__ = ['db_file']
    
    def __init__(self, db_file: str):
        """
        """
        super().__init__(db_file)
    #
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
        return CBFitem(name=name,
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
        table = "CREATE TABLE IF NOT EXISTS CurrentBlockageFactor (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        name NOT NULL,\
                        title TEXT NOT NULL, \
                        type TEXT, \
                        factor DECIMAL);"
        create_table(conn, table)
        # Profile
        table = "CREATE TABLE IF NOT EXISTS CurrentBlockFactorProfile (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        wkf_id NOT NULL REFERENCES CurrentBlockageFactor(number),\
                        elevation DECIMAL NOT NULL,\
                        factor DECIMAL NOT NULL);"
        create_table(conn, table)
    #
    def _push_data(self, conn, mg_data):
        """
        Create a new project into the projects table
        """
        cur = conn.cursor()
        table = 'INSERT INTO CurrentBlockageFactor(name, title, type) \
                 VALUES(?,?,?)'
        # push
        cur = conn.cursor()
        cur.execute(table, mg_data)
    #
    def _pull_data(self, conn, name: str|int,
                  item: str = "*"):
        """ """
        #
        project = (name,)
        sql = 'SELECT {:} FROM CurrentBlockageFactor WHERE name = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
        data = cur.fetchone()
        return data
    #
    # ------------------
    #
    #
    def get_dataX(self, values):
        """ [title, density, constant/profile, thickness] """
        outval = [None, self._rho_w, None, None]
        #mgprofile = outval[2]
        #
        if isinstance(values, dict):
            1/0
        else:
            #if isinstance(values[-1], str):
            #    mgprofile = values.pop()
            #
            outval[0] = values[0]
            #
            try:
                outval[1] = values[1].convert('kilogram/metre^3').value
            except IndexError:
                pass
            #
            try:
                outval[3] = values[2].value
                #mgprofile = 'uniform'
                outval[2] = 'uniform'
            except IndexError:
                #outval[2] = 'profile'
                pass
        #
        #if re.match(r"\b(constant|uniform|mean)\b", mgprofile, re.IGNORECASE):
        #    outval[2] = 'uniform'
        #
        return outval
    #
#
#
#
@dataclass
class CBFitem(HydroItem):
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
        table = f"UPDATE CurrentBlockageFactor \
                 SET type = 'profile' \
                 WHERE name = ?"
        cur = conn.cursor()
        cur.execute(table, mg_name)
        #
        cur = conn.cursor()
        table = 'SELECT * FROM CurrentBlockageFactor WHERE name = ?'
        cur.execute(table, mg_name)
        values = cur.fetchone()        
        #
        profile = tuple((values[0], *item, ) for item in profile_data)
        cur = conn.cursor()
        table = 'INSERT INTO CurrentBlockFactorProfile(wkf_id, \
                                                elevation, factor) \
                                                VALUES(?,?,?)'
        # push
        cur = conn.cursor()
        cur.executemany(table, profile)
    #