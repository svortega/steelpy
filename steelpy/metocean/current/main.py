# 
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
#from collections.abc import Mapping
#from array import array
from dataclasses import dataclass
#from operator import itemgetter
from typing import NamedTuple
import re

# package imports
from steelpy.metocean.current.utils.main import CurrentBasic
from steelpy.metocean.hydrodynamic.utils import HydroItem
from steelpy.utils.sqlite.utils import create_connection, create_table
#
import numpy as np
#
#
#
#
class Current(CurrentBasic):
    """
    """
    __slots__ = ['db_file', '_criteria']
    
    def __init__(self, criteria: str|int, db_file: str):
        """
        """
        super().__init__(db_file)
        self._criteria = criteria
    #
    #
    def __getitem__(self, name):
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            data = self._pull_data(conn, name=name)
        #
        if not data:
            raise IOError(f'Current {name} not found')
        #
        return CurrentItem(name=name,
                           criteria= self._criteria, 
                           db_file=self.db_file)
    #
    #
    def get_dataX(self, values):
        """ [profile (linear/exponential/user), velocity_top, velocity_bottom] """
        outval = [None, None, None, None]
        #cprofile = outval[2]
        if isinstance(values, dict):
            1/0
        else:
            #if isinstance(values[-1], str):
            #    cprofile = values.pop()
            #outval[0] = values[0]
            #
            try:
                outval[0] = values[0]
                if re.match(r"\b(user)\b", outval[0], re.IGNORECASE):
                    pass
                
                else:
                    try:
                        outval[1] = values[1].value
                        try:
                            outval[2] = values[2].value
                        #outval[idx] = item.value
                        except IndexError:
                            outval[2] = outval[1]
                            outval[0] = 'uniform' 
                    
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
    # --------------------------------------
    # SQL ops
    #
    def _create_table(self, conn) -> None:
        """ """
        # Main
        table = "CREATE TABLE IF NOT EXISTS Current (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    type TEXT NOT NULL,\
                    title TEXT, \
                    criteria_id INTEGER NOT NULL REFERENCES Criteria(number));"
        create_table(conn, table)
        # Profile
        table = "CREATE TABLE IF NOT EXISTS CurrentProfile (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    current_id NOT NULL REFERENCES Current(number),\
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
        table = 'INSERT INTO Current(name, type,\
                                        title, criteria_id)\
                                        VALUES(?,?,?,?)'
        # push
        cur = conn.cursor()
        cur.execute(table, current_data)
    #
    #
    #def _pull_data(self, conn, current_name: str|int,
    #              item: str = "*"):
    #    """ """
    #    #
    #    project = (current_name,)
    #    sql = 'SELECT {:} FROM Current WHERE name = ?'.format(item)
    #    cur = conn.cursor()
    #    cur.execute(sql, project)
    #    data = cur.fetchone()
    #    return data
    #
    def _pull_data(self, conn, name):
        """ """
        project = (name, self._criteria)
        table = "SELECT * FROM Current \
                 WHERE name = ? AND criteria_id = ?"
        cur = conn.cursor()
        cur.execute(table, project)
        data = cur.fetchone()
        return data     
    #
    # --------------------------------------
    #
    def df(self, df, columns:dict|None=None):
        """ """
        grouptype = df.groupby(['profile'])
        #grpname = df.groupby(['name'])
        for ctype, dfdata in grouptype:
            if re.match(r"\b(constant|linear)\b", ctype[0], re.IGNORECASE):
                for item in dfdata.itertuples():
                    profile = [item.profile, item.elevation, item.velocity]
                    self.__setitem__(name=item.name, value=profile)
                #print('---')
                
            elif re.match(r"\b(user|profile)\b", ctype[0], re.IGNORECASE):
                grpname = df.groupby(['name'])
                for key, item in grpname:
                    #velmax = item["velocity"].max()
                    profile = item[['elevation', 'velocity']].values.tolist()
                    self.__setitem__(name=key[0], value=profile)
                    #item["velocity"] /= value
                    #item["velocity"] =  [vel.value / velmax.value for vel in item["velocity"]]
                    #profile = item[["elevation", "velocity"]].values.tolist()
                    #zlevel =  item["zlevel"].values.tolist()
                    #self._current[key[0]].profile = profile
            else:
                raise IOError(f' Profile {ctype[0]} not valid')
        #print('-->')
#
#
#
@dataclass
class CurrentItem(HydroItem):
    __slots__ = ['name', '_db_file', '_criteria', 
                 '_depth', '_eta', '_grid', '_cbf']
    
    def __init__(self, name:int|str, criteria: str|int,
                 db_file: str):
        """ """
        #self._depth: float = 0.0
        self._criteria = criteria
        super().__init__(name=name, db_file=db_file)
    #
    #
    # -------------------------------------------
    #
    #
    def _pull_profile(self, conn):
        """get profile data"""
        item = self._pull_current(conn)
        #
        #1 / 0
        if re.match(r"\b(profile)\b", item[2], re.IGNORECASE):
            mg_name = (item[0], )
            cur = conn.cursor()
            table = 'SELECT * FROM CurrentProfile WHERE current_id = ?'
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
        #table = f"UPDATE Current \
        #         SET type = 'profile' \
        #         WHERE name = ?"
        #cur = conn.cursor()
        #cur.execute(table, current_name)        
        #
        #cur = conn.cursor()
        #table = 'SELECT * FROM Current WHERE name = ?'
        #cur.execute(table, current_name)
        #current = cur.fetchone()        
        #
        current = self._pull_current(conn)
        #
        # Converting velocity to SI units
        #try:
        #    profile_data = [[item[0], item[1].value] for item in profile_data]
        #except AttributeError:
        #    pass        
        #
        profile = tuple((current[0], *item, ) for item in profile_data)
        cur = conn.cursor()
        table = 'INSERT INTO CurrentProfile(current_id, \
                                                elevation, velocity) \
                                                VALUES(?,?,?)'
        # push
        cur = conn.cursor()
        cur.executemany(table, profile)
        #print('--->')
    #
    #
    def _push_current(self, conn, data):
        """ get wave data"""
        #
        project = (*data, None, self._criteria)
        table = 'INSERT INTO Current(name, type, title, criteria_id) \
                 VALUES(?,?,?,?,?)'
        #
        #push
        cur = conn.cursor()
        cur.execute(table, project)
        number = cur.lastrowid
        #
        return number
    #
    def _pull_current(self, conn):
        """ """
        project = (self.name, self._criteria)
        table = "SELECT * FROM Current \
                 WHERE name = ? AND criteria_id = ?"
        cur = conn.cursor()
        cur.execute(table, project)
        data = cur.fetchone()
        return data    
    #
    #
    # -------------------------------------------
    #
    def seastate(self, grid:list, eta:list, cbf: list):
        """ """
        #try:
        #    1/ self._depth
        #except ZeroDivisionError:
        #self._depth = d
        #
        self._eta = eta
        self._grid = grid
        self._cbf = cbf
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
    def Vc4(self, d: float, eta: list, zdepth: list):
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
    def get_profile(self, eta: list, zdepth: list, cbf: list):
        """Current Velocity
        
        eta : surface elevation
        zdepth : 
        cbf : Currrent blockage Factor
        """
        #
        # TODO: must be better ways to do this
        profile = self.profile
        elev = list(reversed([item[0] for item in profile]))
        value = list(reversed([item[1] for item in profile]))
        #
        vc = np.zeros((eta.size, zdepth.size))
        #vc[:, :] = Vct
        for i, item in enumerate(eta):
            for j, point in enumerate(zdepth):
                if point > item:
                    vc[i,j] = 0.0
                else:
                    vc[i,j] = np.interp(point, elev, value) * cbf[j]
                                        #right=self.tvelocity)
        return vc    
    #
    #
    def Vc(self, uvector: tuple):
        """
        global system 
        """
        cbf = self._cbf.get_profile(self._grid)
        #cbf *= 0
        vn = self.get_profile(self._eta, self._grid, cbf)
        #
        Vn = uvector.x * vn
        Vnx = vn - uvector.x * Vn
        Vny = - uvector.y * Vn
        Vnz = - uvector.z * Vn
        return CurrVel(Vnx, Vny, Vnz, vn)
#
#
#
class CurrVel(NamedTuple):
    """

    Un : Current components of velocity x
    Vn : Current components of velocity y
    Wn : Current components of velocity z
    vn : Current velocity normal to the cylinder axis
    """
    Un: list
    Vn: list
    Wn: list
    vn : list
    #rho: float
    #
    #def fdn(self, D, cd, UX, vn, rho: float):
    #    """
    #    D  : Member diametre
    #    cd : Drag coefficient
    #    Ux : Instantaneus velocity resolved normal to the member
    #    Vn : Fluid velocity normal to the cylinder axis
    #    """
    #    # drag load per unit length
    #    Fdn = 0.5 * rho * cd * D * UX * vn
    #    return Fdn
    #
    #def FDn(self, Dt:float, Cd:float, rho: float):
    #    """
    #    Component of drag force per unit of cylinder length
    #
    #    Dt : Diametre tubular
    #    Cd : Drag coefficient
    #    
    #    Return:
    #    FDn [x,y,z]
    #    """
    #    FDnx = self.fdn(Dt, Cd, self.Un, self.vn, rho)
    #    FDny = self.fdn(Dt, Cd, self.Vn, self.vn, rho)
    #    FDnz = self.fdn(Dt, Cd, self.Wn, self.vn, rho)
    #    return FDnx, FDny, FDnz    
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