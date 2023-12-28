# 
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
from collections.abc import Mapping
#from array import array
from dataclasses import dataclass
from operator import itemgetter
#from typing import NamedTuple
import re

# package imports
from steelpy.metocean.hydrodynamic.utils import HydroBasic, HydroItem
from steelpy.utils.sqlite.utils import create_connection, create_table
#from steelpy.utils.units.main import Number
import numpy as np
#
#
#
#
class Current(HydroBasic):
    """
    """
    __slots__ = ['_current', 'db_file']
    
    def __init__(self, db_file: str):
        """
        """
        super().__init__(db_file)
        #
        self._current:dict = {}
        #
        # create table
        #conn = create_connection(self.db_file)
        #with conn:
        #    self._create_table(conn)        
    #
    #
    def __getitem__(self, name):
        """
        """
        #
        conn = create_connection(self.db_file)
        with conn:
            data = self._pull_data(conn, current_name=name)
        #
        if not data:
            raise IOError(f'Current {name} not found')
        #
        return CurrentItem(name=name,
                           db_file=self.db_file)

    #
    #
    def get_data(self, values):
        """ [title, profile (linear/exponential/user), velocity_top, velocity_bottom] """
        outval = [None, None, None, None]
        #cprofile = outval[2]
        if isinstance(values, dict):
            1/0
        else:
            #if isinstance(values[-1], str):
            #    cprofile = values.pop()
            outval[0] = values[0]
            #
            try:
                outval[1] = values[1]
                if re.match(r"\b(user)\b", outval[1], re.IGNORECASE):
                    pass
                
                else:
                    try:
                        outval[2] = values[2].value
                        try:
                            outval[3] = values[3].value
                        #outval[idx] = item.value
                        except IndexError:
                            outval[3] = outval[2]
                            outval[1] = 'uniform' 
                    
                    except IndexError:
                        raise IOError('velocity_top missing')
            except IndexError:
                pass
        #
        #if re.match(r"\b(constant|uniform|mean)\b", cprofile, re.IGNORECASE):
        #    outval[2] = 'uniform'
        #else:
        #    outval[2] = 'exponential'
        return outval
    #
    #
    # ------------------
    # SQL ops
    # ------------------
    #
    def _create_table(self, conn) -> None:
        """ """
        # Main
        table = "CREATE TABLE IF NOT EXISTS tb_Current (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    type TEXT NOT NULL,\
                    title TEXT);"
        create_table(conn, table)
        # Profile
        table = "CREATE TABLE IF NOT EXISTS tb_CurrentProfile (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    current_number NOT NULL REFERENCES tb_Current(number),\
                    elevation DECIMAL NOT NULL,\
                    velocity DECIMAL NOT NULL);"
        create_table(conn, table)
    #
    #
    def _push_data(self, conn, current_data):
        """
        Create a new project into the projects table
        """
        cur = conn.cursor()
        #
        table = 'INSERT INTO tb_Current(name, type, \
                                        title) \
                                        VALUES(?,?,?)'
        # push
        cur = conn.cursor()
        cur.execute(table, current_data)
    #
    #
    def _pull_data(self, conn, current_name: str|int,
                  item: str = "*"):
        """ """
        #
        project = (current_name,)
        sql = 'SELECT {:} FROM tb_Current WHERE name = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
        data = cur.fetchone()
        return data    
    #
    # --------------------------------------
    #
    def df(self, df, columns:dict|None=None):
        """ """
        grpname = df.groupby(['name'])
        for key, item in grpname:
            #velmax = item["velocity"].max()
            profile = item[["elevation", "velocity"]].values.tolist()
            self.__setitem__(name=key[0], value=profile)
            #item["velocity"] /= value
            #item["velocity"] =  [vel.value / velmax.value for vel in item["velocity"]]
            #profile = item[["elevation", "velocity"]].values.tolist()
            #zlevel =  item["zlevel"].values.tolist()
            #self._current[key[0]].profile = profile
            
        #print('-->')
#
#
#
@dataclass
class CurrentItem(HydroItem):
    __slots__ = ['name', '_db_file',
                 '_depth', '_eta', 'zd']
    
    def __init__(self, name:int|str, db_file: str):
        """ """
        #self.name = name
        #self.tvelocity = vel_top
        #self.bvelocity = vel_bottom
        #self.cprofile = profile
        #self._profile: list = []
        self._depth: float = 0.0
        #
        #self._db_file = db_file
        super().__init__(name=name, db_file=db_file)
    #
    #@property
    #def profile(self):
    #    """ """
    #    if not self._profile:
    #        Vct = self.tvelocity
    #        #zd = (self.zd + self._depth) / max(self._depth + self._eta)
    #        zd = self.zd[::-1]
    #        val = (zd + self._depth) / max(self._depth + self._eta)
    #        #1/0
    #        val = np.power(np.abs(val), 1.0 / 7.0)
    #        profile = []
    #        for idx, item in enumerate(zd):
    #            profile.append([item, val[idx]])
    #        #self._profile = [self.zd, np.power(np.abs(Vct * zd), 1.0 / 7.0)]
    #        self._profile = profile
    #    1/0
    #    print('-->')
    #    return self._profile
    #
    #@profile.setter
    #def profile(self, value: list[list]):
    #    """ """
    #    prof = []
    #    for item in value:
    #        #try:
    #        elev = item[0].value
    #        #except AttributeError:
    #        #    elev = item[0]
    #        velocity = item[1].value
    #        prof.append([elev, velocity])
    #    #
    #    prof.sort(key=itemgetter(0), reverse=True)
    #    #
    #    conn = create_connection(self._db_file)
    #    with conn:        
    #        self._push_profile(conn, profile_data=prof)
    #    #
    #    #self._profile = prof
    #    #1 / 0
    #    #print('-->')
    #
    #
    def _pull_item(self, conn):
        """ """
        item_name = (self.name, )
        cur = conn.cursor()
        table = 'SELECT * FROM tb_Current WHERE name = ?'
        cur.execute(table, item_name)
        item = cur.fetchone() 
        return item    
    #
    #
    def _pull_profile(self, conn):
        """get profile data"""
        item = self._pull_item(conn)
        #
        if re.match(r"\b(profile)\b", item[2], re.IGNORECASE):
            mg_name = (item[0], )
            cur = conn.cursor()
            table = 'SELECT * FROM tb_CurrentProfile WHERE current_number = ?'
            cur.execute(table, mg_name)
            profile = cur.fetchall()
            profile = [item[2:] for item in profile]
        else:
            1 / 0
        #
        #print('-->')
        return profile    
    #
    def _push_profile(self, conn, profile_data):
        """ """
        #current_name = (self.name, )
        #
        #table = f"UPDATE tb_Current \
        #         SET type = 'profile' \
        #         WHERE name = ?"
        #cur = conn.cursor()
        #cur.execute(table, current_name)        
        #
        #cur = conn.cursor()
        #table = 'SELECT * FROM tb_Current WHERE name = ?'
        #cur.execute(table, current_name)
        #current = cur.fetchone()        
        #
        current = self._pull_item(conn)
        #
        # Converting velocity to SI units
        #try:
        #    profile_data = [[item[0], item[1].value] for item in profile_data]
        #except AttributeError:
        #    pass        
        #
        profile = tuple((current[0], *item, ) for item in profile_data)
        cur = conn.cursor()
        table = 'INSERT INTO tb_CurrentProfile(current_number, \
                                                elevation, velocity) \
                                                VALUES(?,?,?)'
        # push
        cur = conn.cursor()
        cur.executemany(table, profile)
        #print('--->')
    #
    def seastate(self, d:float, z:list, eta:list):
        """ """
        #try:
        #    1/ self._depth
        #except ZeroDivisionError:
        self._depth = d
        #
        self._eta = eta
        self.zd = z
        #print('-->')
    #
    def Vc2(self, Vct: float, d: float, Z):
        """Current Velocity"""
        vc =  np.zeros((Z.shape))
        vc[:, :, :] = Vct
        vcr = np.power(np.abs(Vct * (Z + d), 1.0/7.0))
        vidx = Z < 0
        vc[vidx] = vcr[vidx]
        return vc
    #
    def Vc3(self, d: float, z: list, eta: list):
        """Current Velocity"""
        Vct = self.tvelocity
        zd = (z + d) / d
        vc = np.zeros((eta.size, zd.size))
        vc[:, :] = Vct
        #vcr1 = np.abs(Vct *  d)
        #
        ze = np.zeros((eta.size, zd.size))
        ze[:, :] = z
        zebool = np.transpose(ze[:, :].T < eta)
        #
        vcr = np.zeros((eta.size, zd.size))
        vcr[:, :] = np.power(np.abs(Vct * zd), 1.0/7.0)
        #
        vc[zebool] = vcr[zebool]
        return vc
    #
    def Vc(self, d: float, eta: list, zdepth: list):
        """Current Velocity"""
        Vct = self.tvelocity
        elev = list(reversed([item[0] for item in self.profile]))
        value = list(reversed([item[1] * Vct for item in self.profile]))
        #
        vc = np.zeros((eta.size, zdepth.size))
        #vc[:, :] = Vct
        for i, item in enumerate(eta):
            for j, point in enumerate(zdepth):
                vc[i,j] = np.interp(point, elev, value)
                                    #right=self.tvelocity)
        return vc
    #
    #
    def get_profile(self, eta: list, zdepth: list):
        """Current Velocity"""
        #Vct = self.tvelocity
        profile = self.profile
        elev = list(reversed([item[0] for item in profile]))
        value = list(reversed([item[1] for item in profile]))
        #
        vc = np.zeros((eta.size, zdepth.size))
        #vc[:, :] = Vct
        for i, item in enumerate(eta):
            for j, point in enumerate(zdepth):
                vc[i,j] = np.interp(point, elev, value)
                                    #right=self.tvelocity)
        return vc    
    #
    #
#
#
#
class MeanCurrent:
    __slots__ = ['c_type', '_current']
    
    def __init__(self):
        """ """
        self.c_type = 1
        self._current = 0
    #
    @property
    def Euler(self):
        """ Eularian mean """
        self.c_type = 1
         
    @Euler.setter
    def Euler(self, value:Units):
        """ Eularian mean """
        self._current = value.convert('metre/second').value
        self.c_type = 1     
    #
    @property
    def Stokes(self):
        """ Stokes depth integrated mean """
        self.c_type = 2
    
    @Stokes.setter
    def Stokes(self, value:Units):
        """ Stokes depth integrated mean """
        self._current = value.convert('metre/second').value
        self.c_type = 2