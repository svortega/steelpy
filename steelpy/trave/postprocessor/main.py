#
# Copyright (c) 2009 steelpy
#
from __future__ import annotations
# Python stdlib imports
#from dataclasses import dataclass
#from typing import NamedTuple
#import pickle
#from typing import NamedTuple
from datetime import datetime as dt

#
# package imports
from steelpy.trave.postprocessor.operations import MainProcess
#from .output import  Node, Beam
#from steelpy.utils.dataframe.main import DBframework
from steelpy.utils.sqlite.main import ClassMainSQL
from steelpy.utils.sqlite.utils import create_table
#
from steelpy.trave.postprocessor.solution import UnSQL
#
#
# --------------------
# Results
# --------------------
#
#
class PostProcess(ClassMainSQL):
    __slots__ = ['_mesh', '_process', '_solution', 'db_file']
    
    def __init__(self, mesh, name:str|None,
                 sql_file:str|None) -> None:
        """
        """
        self._mesh = mesh
        super().__init__(name=name, sql_file=sql_file)
        #
        self._process = MainProcess(mesh=mesh,
                                    db_file=self.db_file)
        #
        self._solution =  UnSQL(load=self._mesh._load,
                                db_file=self.db_file)        
    #
    #def mesh(self, mesh):
    #    """ """
    #    self._mesh = mesh
    #    self._plane = self._mesh._plane
    #
    # --------------------
    # Results
    # --------------------
    #
    @property
    def Un(self):
        """ Nodal displacement solution"""
        return self._solution.df
    #
    @Un.setter
    def Un(self, df):
        """Nodal displacement solution"""
        self._solution.df = df
    #
    #
    # --------------------
    #
    #
    def results(self, beam_steps):
        """ """
        res = self._process.results(Un=self.Un,
                                    beam_steps=beam_steps)
        #print('-->')
        return res
    #
    # ------------------
    # SQL ops
    # ------------------
    #
    def _create_table(self, conn) -> None:
        """ """
        # conn = create_connection(self.db_file)
        table = "CREATE TABLE IF NOT EXISTS tb_Main (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    type TEXT,\
                    mesh_name TEXT,\
                    units TEXT NOT NULL,\
                    plane TEXT NOT NULL,\
                    date TEXT);"
        #
        create_table(conn, table)
        #
        #
        table = 'INSERT INTO tb_Main(name, type, mesh_name, \
                                     units, plane, date)\
                                     VALUES(?,?,?,?,?,?)'
        #
        time=dt.now().strftime('%Y-%m-%d')
        plane = '3D'
        if self._mesh._plane.plane2D:
            plane = '2D'
        data = (self._name, None, self._mesh._name, # name, type, mesh_name
                'si', plane, time)                  # units, plane, date
        # push
        cur = conn.cursor()
        cur.execute(table, data)