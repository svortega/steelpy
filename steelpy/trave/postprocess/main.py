#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass

import numpy as np

#from typing import NamedTuple
#import pickle
#from typing import NamedTuple
#from datetime import datetime as dt

#
# package imports
from steelpy.trave.postprocess.operations import MainProcess
from steelpy.utils.sqlite.utils import create_table, create_connection
from steelpy.utils.dataframe.main import DBframework
#
#
# --------------------
# Results
# --------------------
#
#
class PostProcess: # (ClassMainSQL)
    __slots__ = ['_mesh', '_process', '_Un',
                 'db_file', '_name']
    
    def __init__(self, mesh, name:str|None) -> None:
        """
        """
        self._mesh = mesh
        if not name:
            self._name = self._mesh._name
        #
        self.db_file = self._mesh.db_file
        
        #super().__init__(name=name, sql_file=sql_file)
        #
        #
        self._Un = UnSQL(mesh=self._mesh,
                         db_file=self.db_file)
        #
        self._process = MainProcess(name=self._name,
                                    mesh=mesh,
                                    db_file=self.db_file)
        #
        #self._solution =  UnSQL(load=self._mesh._load,
        #                        db_file=self.db_file)
        #
        
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
        return self._Un
    #
    #@Un.setter
    #def Un(self, df):
    #    """Nodal displacement solution"""
    #    self._Un = df
    #
    #
    # --------------------
    #
    #
    def results(self, beam_steps: int,
                Pdelta: bool):
        """ """
        res = self._process.results(Un=self.Un.df,
                                    beam_steps=beam_steps,
                                    Pdelta=Pdelta)
        #print('-->')
        return res
    #
    # ------------------
    # SQL ops
    # ------------------
    #
    #def _create_table(self, conn) -> None:
    #    """ """
    #    # conn = create_connection(self.db_file)
    #    table = "CREATE TABLE IF NOT EXISTS Component (\
    #                number INTEGER PRIMARY KEY NOT NULL,\
    #                name NOT NULL,\
    #                type TEXT,\
    #                mesh_name TEXT,\
    #                units TEXT NOT NULL,\
    #                plane TEXT NOT NULL,\
    #                date TEXT);"
    #    #
    #    create_table(conn, table)
    #    #
    #    #
    #    table = 'INSERT INTO Component(name, type, mesh_name, \
    #                                 units, plane, date)\
    #                                 VALUES(?,?,?,?,?,?)'
    #    #
    #    time=dt.now().strftime('%Y-%m-%d')
    #    plane = '3D'
    #    if self._mesh._plane.plane2D:
    #        plane = '2D'
    #    data = (self._name, None, self._mesh._name, # name, type, mesh_name
    #            'si', plane, time)                  # units, plane, date
    #    # push
    #    cur = conn.cursor()
    #    cur.execute(table, data)
#
#
#
@dataclass
class UnSQL:
    """ """
    def __init__(self, mesh, db_file) -> None:
        self.mesh = mesh
        self.db_file = db_file
        #
        conn = create_connection(self.db_file)
        with conn:
            self._create_table(conn)
    #
    # ---------------------------------
    #
    @property
    def df(self):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            #Un = get_Undf(conn)
            ndata = self._pull_data(conn)
        #
        # dataframe
        cols = ['number', 'load_name', 'component_name',
                'load_level', 'load_system', 
                'node_name',
                'x', 'y', 'z', 'rx', 'ry', 'rz']         
        db = DBframework()
        Un = db.DataFrame(data=ndata, columns=cols)
        #Un = Un.copy()
        Un.fillna(value=float(0.0), inplace=True)
        Un.drop(labels=['number'], axis=1, inplace=True)
        return Un.astype({'x': 'float64', 'y': 'float64', 'z': 'float64',
                          'rx': 'float64', 'ry': 'float64', 'rz': 'float64'})
    
    @df.setter
    def df(self, df):
        """ """
        conn = create_connection(self.db_file)
        #
        #1 / 0
        # TODO: here combine with load displacement
        #
        # TODO: check this works
        # combination
        #load_comb = self.mesh._load.combination()
        #df_comb = load_comb.to_basic()        
        # displacments
        #df = self._update_ndf(dfnode=df, dfcomb=df_comb)
        #
        #
        header = ['load_name', 'component_name',
                  'load_level', 'load_system',
                  'node_name',
                  'x', 'y', 'z', 'rx', 'ry', 'rz']
        try: # 3d plane
            nodeconn = df[header].copy()
        except KeyError: # 2d plane
            header = ['load_name', 'component_name', 
                      'load_level', 'load_system', 
                      'node_name','x', 'y', 'rz']
            nodeconn = df[header].copy()
            #nodeconn[['z', 'rx', 'ry']] = float(0)
        #
        #nodeconn.replace(to_replace=[''], value=[float(0)], inplace=True)
        #nodeconn['element_id'] = df['element_id']
        #nodeconn['title'] = df['title']
        #
        with conn:
            nodeconn.to_sql('ResultNodeU', conn,
                            index_label=header, 
                            if_exists='append', index=False)
        #
        #print('nodal disp results saved')
    #
    def _update_ndf(self, dfnode, dfcomb, 
                   values:list[str] = ['x', 'y', 'z', 'rx', 'ry', 'rz']):
        """
        Update node displacements to include lcomb
        """
        db = DBframework()
        # group basic load by name
        ndgrp = dfnode.groupby('load_name')
        # get combinations with basic loads 
        #
        combgrp = dfcomb.groupby('load_name')
        for key, combfactors in combgrp:
            for row in combfactors.itertuples():
                comb = ndgrp.get_group(row.basic_load).copy()
                comb.loc[:, values] *= row.factor
                comb['load_level'] = 'combination'
                comb['load_name'] = row.load_name
                comb['load_id'] = row.load_id
                comb['load_title'] = row.load_title
                #
                try:
                    dftemp = db.concat([dftemp, comb], ignore_index=True)
                except UnboundLocalError:
                    dftemp = comb
        #
        try:
            #check = dftemp.groupby(['node_name', 'c']).sum().reset_index()
            dftemp = dftemp.groupby(['load_name', 'load_id','load_level',
                                     'load_title', 'load_system','node_name'],
                                      as_index=False)[values].sum()
            #test
            dfnode = db.concat([dfnode, dftemp], ignore_index=True)
        except UnboundLocalError:
            pass
        #
        return dfnode #, memb_comb
    #
    # ---------------------------------
    # SQL ops
    #
    def _create_table(self, conn) -> None:
        """ """
        # Node displacement solution
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultNodeU (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
                        component_name NOT NULL, \
                        load_level TEXT NOT NULL,\
                        load_system TEXT NOT NULL,\
                        node_name NOT NULL,\
                        x DECIMAL,\
                        y DECIMAL,\
                        z DECIMAL,\
                        rx DECIMAL,\
                        ry DECIMAL,\
                        rz DECIMAL);"
        #
        create_table(conn, table_nodes)
    #
    #
    def _pull_data(self, conn):
        """ """
        #table = "SELECT number AS number, load_name AS load_name, \
        #                component_name AS component_name,  load_level AS load_level, \
        #                load_system AS load_system, node_name AS node_name, \
        #                IFNULL(x, 0.0) AS x, IFNULL(y, 0.0) AS y, IFNULL(z, 0.0) AS z,\
        #                IFNULL(rx, 0.0) AS rx, IFNULL(ry, 0.0) AS ry, IFNULL(rz, 0.0) AS rz \
        #        FROM ResultNodeU"
        #
        table = "SELECT ResultNodeU.* FROM ResultNodeU"
        #
        # Node load
        with conn:
            cur = conn.cursor()
            cur.execute(table)
            rows = cur.fetchall()
        return rows
    #    
#