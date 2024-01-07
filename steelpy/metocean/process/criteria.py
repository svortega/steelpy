# 
# Copyright (c) 2019 steelpy
# 
# Python stdlib imports
from __future__ import annotations
#from collections.abc import Mapping
from datetime import datetime as dt
from dataclasses import dataclass
from typing import NamedTuple
#import os
import re


# package imports
from steelpy.utils.sqlite.main import ClassBasicSQL
from steelpy.metocean.wave.main import Wave
from steelpy.metocean.wind.main import Wind
from steelpy.metocean.current.main import Current
from steelpy.metocean.process.conditions import HydroCondition
from steelpy.utils.sqlite.utils import create_connection, create_table


#
class HydroCriteria(ClassBasicSQL):
    """ f2u Concept model Class """
    __slots__ = ['_labels', '_item',
                 'db_file', '_properties', '_rho_w']
    #
    def __init__(self, rho_w, properties, db_file: str):
        """
        """
        super().__init__(db_file=db_file)
        #self._name = name
        #self.db_file = db_file
        self._labels = list = []
        #
        self._properties = properties
        self._rho_w = rho_w
    #
    def __setitem__(self, name: int|str, title: int|str) -> None:
        """
        """
        try:
            self._labels.index(name)
            raise Exception('Concept {:} already exist'.format(name))
        except ValueError:
            self._labels.append(name)
            #
            conn = create_connection(self.db_file)
            with conn:
                self._push_data(conn, criteria=name, title=title)
    #
    def __getitem__(self, name: int|str):
        """ """
        try:
            self._labels.index(name)
            #return self._item[name]
            conn = create_connection(self.db_file)
            with conn:            
                data = self._pull_data(conn, name)
            #
            #1 / 0
            return CriteriaItem(name=name,
                                number=data.number,
                                rho_w=data.rho, 
                                properties= self._properties,
                                db_file=self.db_file)
        except ValueError:
            raise IndexError(' ** Concept {:} does not exist'.format(name))         
    #
    # ------------------
    # SQL ops
    #
    def _create_table(self, conn) -> None:
        """ """
        table = "CREATE TABLE IF NOT EXISTS Criteria (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    type TEXT NOT NULL,\
                    units TEXT NOT NULL,\
                    rho_water DECIMAL NOT NULL,\
                    date TEXT, \
                    title TEXT);"
        #
        create_table(conn, table)
    #
    #
    def _push_data(self, conn, criteria: str,
                   title: str|None = None):
        """
        Create a new project into the projects table
        """
        #
        table = 'INSERT INTO Criteria(name, type, units,\
                                    rho_water, date, \
                                    title)\
                                    VALUES(?,?,?,?,?,?)'
        #
        time=dt.now().strftime('%Y-%m-%d')
        data = (criteria, 'metocean', 'si', self._rho_w, time, title)
        # push
        cur = conn.cursor()
        row = cur.execute(table, data)
        #print('-->')
        #return row.lastrowid
    #
    def _pull_data(self, conn, name: str|int,
                   item: str = "*"):
        """ """
        project = (name, )
        table = "SELECT * FROM Criteria \
                 WHERE name = ?"
        cur = conn.cursor()
        cur.execute(table, project)
        case = cur.fetchone()
        #print('-->')
        return CriteriaOut._make(case)
    #
    # ------------------
    # Operations
    #
    def solve(self, surface_points:int = 36,
              depth_points:int = 100):
        """ """
        for name in self._labels:
            item = self.__getitem__(name)
            item.solve(surface_points, depth_points)
        #print('---')
#
#
class CriteriaOut(NamedTuple):
    """ """
    number: int
    name: str | int
    criteria_type: str | None
    units: str
    rho: float
    date: str
    #criteria: str | int
    title: str | None
#
#
@dataclass
class CriteriaItem:
    """ f2u Concept model Class """
    __slots__ = ['_wave', '_current', '_wind',
                 '_condition', 
                 '_name', '_number']
    #
    def __init__(self, name:str, number: int, 
                 rho_w, properties, db_file: str):
        """
        """
        self._name = name
        self._number = number
        #
        self._wave = Wave(criteria=self._number, 
                          rho_w=rho_w,
                          db_file=db_file)
        
        self._wind = Wind(criteria=self._number,
                          db_file=db_file)
        
        self._current = Current(criteria=self._number,
                                db_file=db_file)
        #
        self._condition = HydroCondition(criteria=self._number,
                                         properties=properties,
                                         db_file=db_file)
    #
    # ------------------------------------------
    #
    def wave(self, values:None|list=None,
             df=None):
        """ """
        if values:
            print('-->')
            1/0
        else:
            try:
                columns = list(df.columns)
                grpwave = df.groupby('type')
                for wtype, item in grpwave:
                    if re.match(r"\b(regular(_wave)?)\b", wtype, re.IGNORECASE):
                        self._wave._regular.df(df=item)
                    elif re.match(r"\b(iregular(_wave)?)\b", wtype, re.IGNORECASE):
                        # self._wave._iregular.df(df=item)
                        raise NotImplementedError
                    else:
                        raise ValueError("wave type invalid")
            except AttributeError:
                pass
        #
        return self._wave
    #
    #@property
    def wind(self, values:None|list=None,
             df=None):
        """
        """
        if values:
            print('-->')
            1/0
        else:
            try:
                df.columns            
                self._wind.df(df)
            except AttributeError:
                pass
        #
        return self._wind
    #
    #@property
    def current(self, values:None|list=None,
                df=None):
        """
        """
        if values:
            print('-->')
            1/0
        else:
            try:
                columns = list(df.columns)
                self._current.df(df)
            except AttributeError:
                pass
        return self._current
    #
    #
    # ------------------------------------------
    #@property
    def condition(self, values:None|list=None,
                 df=None):
        """
        """
        if values:
            print('-->')
            1/0
        else:
            try:
                df.columns            
                self._condition.df(df)
            except AttributeError:
                pass        
        return self._condition
    #
    # ------------------------------------------
    #
    def solve(self, surface_points:int = 36,
              depth_points:int = 100):
        """ """
        # wave
        self._wave.solve(surface_points=surface_points,
                         depth_points=depth_points)
        #
        # current
        #
        # wind
        #print('-->')
        #1 / 0
    #    
#
#