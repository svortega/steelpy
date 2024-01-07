# 
# Copyright (c) 2009 fem2ufo
#

# Python stdlib imports
from __future__ import annotations
#from array import array
from collections.abc import Mapping
import re
from typing import NamedTuple
from dataclasses import dataclass
#from operator import itemgetter


# package imports
from steelpy.metocean.hydrodynamic.utils.main import (BasicProperty, HydroBasic,
                                                      get_list, HydroItem)
from steelpy.utils.sqlite.utils import create_connection, create_table
import numpy as np

#
class Direction(NamedTuple):
    """
    Direction dependent coefficients
    """
    x:float
    y:float
    z:float    
#
class CdCm(NamedTuple):
    """
    Morison parameters CD & Cm
    """
    #Cdx:float
    #Cdy:float
    Cd:float
    #
    #Cmx:float
    #Cmy:float
    Cm:float
    #
    #number:int
    #sets:List
    #case:str
    #
    @property
    def drag(self):
        return Direction(self.Cdx, self.Cdy, self.Cdz)
    
    @property
    def mass(self):
        return Direction(self.Cmx, self.Cmy, self.Cmz)  
#
#
#
#
# ---------------------------------------------
#
class CdCmCoefficients(HydroBasic):
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
            raise IOError(f'CdCm {name} not found')
        #
        return CdCmitem(name=name,
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
        table = "CREATE TABLE IF NOT EXISTS CdCm (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        name NOT NULL,\
                        type TEXT NOT NULL, \
                        title TEXT);"
        create_table(conn, table)
        # Profile
        table = "CREATE TABLE IF NOT EXISTS CdCmProfile (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        cdcm_id NOT NULL REFERENCES CdCm(number),\
                        elevation DECIMAL NOT NULL,\
                        cd_coefficient DECIMAL NOT NULL, \
                        cm_coefficient DECIMAL NOT NULL);"
        create_table(conn, table)
    #
    def _push_data(self, conn, mg_data):
        """
        Create a new project into the projects table
        """
        cur = conn.cursor()
        table = 'INSERT INTO CdCm(name, type, title) \
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
        sql = 'SELECT {:} FROM CdCm WHERE name = ?'.format(item)
        cur = conn.cursor()
        cur.execute(sql, project)
        data = cur.fetchone()
        return data
    #
#
#
class CdCmCoefficientsXX(Mapping):
    """
    """
    __slots__ =  ['_cdcm']
    
    def __init__(self):
        """
        """
        self._cdcm: dict = {}
    
    def __getitem__(self, cdcm_name):
        """
        """
        return self._cdcm[cdcm_name]
    
    def __setitem__(self, name, values:list|tuple|dict|str) -> None:
        """
        rule
        specified
        diametre
        """
        cdcm_type = values.pop(0)
        #values = self._get_value(values)
        if re.match(r"\b((rule(\_)?)(api|iso))\b", cdcm_type, re.IGNORECASE):
            self._cdcm[name] = DiametreFunction()
            self._cdcm[cdcm_name].rule = cdcm_type
        elif re.match(r"\b(specified|coefficients)\b", cdcm_type, re.IGNORECASE):
            self._cdcm[name] = SpecifiedFunction()
            self._cdcm[name].set_cdcm(cdcm=values)
        elif re.match(r"\b(diamet(re|re))\b", cdcm_type, re.IGNORECASE):
            self._cdcm[name] = DiametreFunction()
        elif re.match(r"\b(depth(\_)?profile)\b", cdcm_type, re.IGNORECASE):
            self._cdcm[name] = DepthProfileFunction()
        else:
            raise IOError("CdCm type {:} not implemented".format(cdcm_type))         
    #
    def _get_value(self, value:list|tuple|dict|str):
        """ """
        if isinstance(value, (list, tuple)):
            # [name, cdx, cdy, cdz, cmx, cmy, cmz]
            value = get_list(value, steps)
        elif isinstance(value, dict):
            value = get_dic(value)
        else:
            raise Exception('   *** input format not recognized')
        return value

    #
    #
    def __len__(self) -> float:
        return len(self._cdcm)

    def __iter__(self) -> Iterator[Tuple]:
        return iter(self._cdcm)
    
    def __contains__(self, value) -> bool:
        return value in self._cdcm
    #
    #
    #
    def getCdCm(self, Z, HTs: float, condition:int):
        """ """
        cm = np.zeros((Z.shape))
        cm += 1.2
        cm[Z > HTs] = 1.6
        #cm[Z <= 2] = 1.2
        #
        cd = np.zeros((Z.shape))
        # switch condition
        if condition == 1:
            cd += 1.15
        else:
            #elif condition == 2:
            cd += 1.05
            cd[Z > HTs] = 0.65
            #cd[Z <= 2] = 1.05
        #
        return cd, cm
    #
#
#
@dataclass
class CdCmitem(HydroItem):
    __slots__ = ['name', '_db_file']
    #
    def __init__(self, name: int|str,
                 db_file: str):
        """ """
        super().__init__(name=name, db_file=db_file)
    #
    #
    def _pull_item(self, conn):
        """ """
        mg_name = (self.name, )
        cur = conn.cursor()
        table = 'SELECT * FROM CdCm WHERE name = ?'
        cur.execute(table, mg_name)
        values = cur.fetchone()
        return values
    #
    #
    def _push_profile(self, conn, profile_data):
        """ """
        #mg_name = (self.name, )
        #
        #table = f"UPDATE CdCm \
        #         SET type = 'profile' \
        #         WHERE name = ?"
        #cur = conn.cursor()
        #cur.execute(table, mg_name)
        #
        #cur = conn.cursor()
        #table = 'SELECT * FROM CdCm WHERE name = ?'
        #cur.execute(table, mg_name)
        #cdcm = cur.fetchone()
        #
        cdcm = self._pull_item(conn)
        #
        profile = tuple((cdcm[0], *item, ) for item in profile_data)
        cur = conn.cursor()
        table = 'INSERT INTO CdCmProfile(cdcm_id, \
                                            elevation, \
                                            cd_coefficient, \
                                            cm_coefficient) \
                                            VALUES(?,?,?,?)'
        # push
        cur = conn.cursor()
        cur.executemany(table, profile)
    #
    def _pull_profile(self, conn):
        """get profile data"""
        cdcm = self._pull_item(conn)
        #
        if re.match(r"\b(profile)\b", cdcm[2], re.IGNORECASE):
            item_name = (cdcm[0], )
            cur = conn.cursor()
            table = 'SELECT * FROM CdCmProfile WHERE cdcm_id = ?'
            cur.execute(table, item_name)
            profile = cur.fetchall()
            profile = [item[2:] for item in profile]
        else:
            1 / 0
        #
        #print('-->')
        return profile    
    #
    #
    def getCdCm(self, Z, HTs: float, condition:int):
        """ """
        #
        cdcm_profile = self.profile
        #
        cm = np.zeros((Z.shape))
        cm += 1.2
        cm[Z > HTs] = 1.6
        #cm[Z <= 2] = 1.2
        #
        cd = np.zeros((Z.shape))
        # switch condition
        if condition == 1:
            cd += 1.15
        else:
            #elif condition == 2:
            cd += 1.05
            cd[Z > HTs] = 0.65
            #cd[Z <= 2] = 1.05
        #
        return cd, cm
    #
    #
    def get_profile(self, Z):
        """ """
        prof = self.profile
        elev = list(reversed([item[0] for item in prof]))
        cdprof = list(reversed([item[1] for item in prof]))
        cmprof = list(reversed([item[1] for item in prof]))
        #
        cd = np.zeros((Z.shape))
        cm = np.zeros((Z.shape))
        for i, point in enumerate(Z):
            cd[i] = np.interp(point, elev, cdprof)
            cm[i] = np.interp(point, elev, cmprof)
        #
        return cd, cm
#
#
class DiametreFunction(BasicProperty):
    __slots__ = ['_type', '_diameter', '_rule']
    
    def __init__(self):
        """
        """
        BasicProperty.__init__(self)
        self._type = 'diameter_function'
        self._diameter:Dict = {}
    
    def diameter_function(self, diameter, smooth, rough):
        """
        """
        self._diameter[diameter] = [smooth, rough]
    #
    @property
    def rule(self):
        """
        """
        return self._rule
    
    @rule.setter
    def rule(self, value:str):
        """
        """
        self._rule = value
#
class SpecifiedFunction(BasicProperty):
    __slots__ = ['_type', '_coefficients']
    
    def __init__(self):
        """
        """
        BasicProperty.__init__(self)
        self._type = 'specified_function'
    
    #@property
    def set_cdcm(self, cdcm:list|dict):
        """
        """
        if type(cdcm) == list:
            self._coefficients = CdCm._make(cdcm)
        elif type(cdcm) == dict:
            self._coefficients = CdCm.Fixity(**cdcm)
        else:
            self._coefficients = CdCm(Cd, Cm)
#
class DepthProfileFunction(BasicProperty):
    __slots__ = ['_type', '_profile']
    
    def __init__(self):
        """
        """
        BasicProperty.__init__(self)
        self._type = 'depth_profile'
        self._profile:List = []
    #
    #
    def depth_profile(self, depth, cdcm):
        """
        """
        #_coefficient = MorrisonCoefficient.set_cdcm(self, cdcm)
        self._profile.extend([depth, cdcm])
   