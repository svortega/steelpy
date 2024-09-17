#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
import time
#from typing import NamedTuple
#import re

# package imports
# steelpy.f2uModel
from steelpy.ufo.load.process.actions import SelfWeight
from steelpy.ufo.load.process.basic_load import LoadCaseBasic, LoadTypeBasic 
from steelpy.ufo.load.sqlite.beam import BeamLoadItemSQL, BeamLoadGloabalSQL
from steelpy.ufo.load.sqlite.node import  NodeLoadItemSQL, NodeLoadGlobalSQL
#
# steelpy.f2uModel
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework
#
#
#
class BasicLoadSQL(LoadCaseBasic):
    __slots__ = ['db_file', '_component']

    #
    def __init__(self, db_file:str, #plane: NamedTuple,
                 component: int) -> None:
        """
        """
        super().__init__()
        #
        self.db_file = db_file
        #self._plane = plane
        self._component = component
        #
        self._nodes = NodeLoadGlobalSQL(component=self._component,
                                        db_file=self.db_file)
        #
        self._beams = BeamLoadGloabalSQL(component=self._component,
                                         #plane=self._plane, 
                                         db_file=self.db_file)
        #
        conn = create_connection(self.db_file)
        with conn: 
            self._new_table(conn)        
    #
    @property
    def _labels(self):
        """ """
        query = ('basic', self._component, )
        table = "SELECT Load.name \
                  FROM Load, LoadBasic \
                  WHERE Load.level = ? \
                  AND LoadBasic.load_id = Load.number\
                  AND Load.mesh_id = ? ;"
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        items = set([item[0] for item in items])
        return list(items)
    #
    #
    # -----------------------------------------------
    #
    def __setitem__(self, load_name:int|str, load_title:str) -> None:
        """
        """
        try:
            self._labels.index(load_name)
            raise Warning(f'    *** warning load name {load_name} already exist')
        except ValueError:
            conn = create_connection(self.db_file)
            with conn:
                self._push_load(conn, load_name, load_title)
    #           
    def __getitem__(self, load_name: str|int):
        """
        """
        try:
            index = self._labels.index(load_name)
            return LoadTypeSQL(load_name=load_name,
                               component=self._component, 
                               bd_file=self.db_file)
        except ValueError:
            raise IOError("load case not defined")

    #
    # -----------------------------------------------
    #
    def _push_load(self, conn, load_name:int|str, load_title:str):
        """ """
        query = (load_name,  self._component, "basic", load_title, 'ufo', None)
        table = 'INSERT INTO Load(name, mesh_id, level, title, \
                 input_type, input_file) \
                 VALUES(?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.execute(table, query)
        idx = cur.lastrowid
        #return cur.lastrowid
        query = (idx, None)
        table = 'INSERT INTO LoadBasic(load_id, step_type) \
                 VALUES(?,?)'
        cur = conn.cursor()
        cur.execute(table, query)        
    #
    #
    def _new_table(self, conn):
        """ """
        # -------------------------------------
        # Main
        table = "CREATE TABLE IF NOT EXISTS LoadBasic(\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    load_id INTEGER NOT NULL REFERENCES Load(number),\
                    step_type TEXT, \
                    design_load TEXT);"
        # UNIQUE(load_id, load_type, step)
        create_table(conn, table)
        # -------------------------------------
        # Node Load
        table = "CREATE TABLE IF NOT EXISTS LoadNode(\
                number INTEGER PRIMARY KEY NOT NULL,\
                basic_id INTEGER NOT NULL REFERENCES LoadBasic(number),\
                node_id INTEGER NOT NULL REFERENCES Node(number),\
                title TEXT,\
                system TEXT NOT NULL,\
                type TEXT NOT NULL,\
                fx DECIMAL,\
                fy DECIMAL,\
                fz DECIMAL,\
                mx DECIMAL,\
                my DECIMAL,\
                mz DECIMAL,\
                x DECIMAL,\
                y DECIMAL,\
                z DECIMAL,\
                rx DECIMAL,\
                ry DECIMAL,\
                rz DECIMAL,\
                psi DECIMAL,\
                B DECIMAL,\
                Tw DECIMAL,\
                step DECIMAL);"
        create_table(conn, table)
        # -------------------------------------
        # line 
        table = "CREATE TABLE IF NOT EXISTS LoadBeamLine(\
                number INTEGER PRIMARY KEY NOT NULL,\
                basic_id INTEGER NOT NULL REFERENCES LoadBasic(number),\
                element_id INTEGER NOT NULL REFERENCES Element(number),\
                title TEXT,\
                system INTEGER NOT NULL,\
                type TEXT NOT NULL,\
                L0 DECIMAL,\
                qx0 DECIMAL,\
                qy0 DECIMAL,\
                qz0 DECIMAL,\
                qt0 DECIMAL,\
                L1 DECIMAL,\
                qx1 DECIMAL,\
                qy1 DECIMAL,\
                qz1 DECIMAL,\
                qt1 DECIMAL,\
                step DECIMAL);"
        create_table(conn, table)
        # -------------------------------------
        # point
        table = "CREATE TABLE IF NOT EXISTS LoadBeamPoint(\
                number INTEGER PRIMARY KEY NOT NULL,\
                basic_id INTEGER NOT NULL REFERENCES LoadBasic(number),\
                element_id INTEGER NOT NULL REFERENCES Element(number),\
                title TEXT,\
                system INTEGER NOT NULL,\
                type TEXT NOT NULL,\
                L0 DECIMAL,\
                fx DECIMAL,\
                fy DECIMAL,\
                fz DECIMAL,\
                mx DECIMAL,\
                my DECIMAL,\
                mz DECIMAL,\
                x DECIMAL,\
                y DECIMAL,\
                z DECIMAL,\
                rx DECIMAL,\
                ry DECIMAL,\
                rz DECIMAL, \
                step DECIMAL);"
        create_table(conn, table)
        # -------------------------------------
        # FER
        table = "CREATE TABLE IF NOT EXISTS LoadBeamFER(\
                number INTEGER PRIMARY KEY NOT NULL,\
                basic_id INTEGER NOT NULL REFERENCES LoadBasic(number),\
                element_id INTEGER REFERENCES Element(number),\
                node_id INTEGER NOT NULL REFERENCES Node(number),\
                title TEXT,\
                system TEXT NOT NULL,\
                type TEXT NOT NULL,\
                fx DECIMAL,\
                fy DECIMAL,\
                fz DECIMAL,\
                mx DECIMAL,\
                my DECIMAL,\
                mz DECIMAL,\
                x DECIMAL,\
                y DECIMAL,\
                z DECIMAL,\
                rx DECIMAL,\
                ry DECIMAL,\
                rz DECIMAL,\
                Psi DECIMAL,\
                B DECIMAL,\
                Tw DECIMAL, \
                step DECIMAL);"
        create_table(conn, table)        
    #
    # -----------------------------------------------
    #
    # -----------------------------------------------
    # Loading Operations
    # -----------------------------------------------
    #
    def function(self, steps:int,
                 Pa:float=0.0, factor:float=1):
        """process element load"""
        #
        dftemp = self._beams.load_function(steps=steps,
                                           Pa=Pa, factor=factor)
        #
        # Axial   [FP, blank, blank, Fu]
        # torsion [T, B, Psi, Phi, Tw]
        # Bending [V, M, theta, w]        
        #
        # [Fx, Fy, Fz, Mx, My, Mz]
        # [V, M, w, theta]
        header = ['load_name', 'mesh_name',
                  'load_comment', 'load_type', 
                  'load_level', 'load_system',
                  'element_name', 'node_end',
                  'axial', 'torsion', 'VM_inplane', 'VM_outplane']
        #
        #          'FP', 'blank1', 'blank2', 'Fu',
        #          'T', 'B', 'Psi', 'Phi', 'Tw',
        #          'Vy', 'Mz', 'theta_y', 'w_y',
        #          'Vz', 'My', 'theta_z', 'w_z']
        df = DBframework()
        dfload = df.DataFrame(data=dftemp, columns=header, index=None)
        #return load_func
        return dfload
    #
    #
    def FER(self, elements):
        """
        Fixed End Reactions (FER) 
        Convert element load to global node loads
        """
        # TODO : loop slow
        #print('--- FER Operation')
        start_time = time.time()
        beams = elements.beam()
        load_name = list(self.keys())
        for lname in load_name:
            lcase = self.__getitem__(lname)
            #
            #
            # -------------------------------
            # beam line load (line & point)
            # -------------------------------
            #
            lcase._beam.fer(beams=beams)
            #
            # -------------------------------
            # plates
            # -------------------------------
            #
            #
            # -------------------------------
            # 
            # -------------------------------
        #
        uptime = time.time() - start_time
        print(f"** FER Operation: {uptime:1.4e} sec")
    #
    #
    def ENL(self):
        """Equivalent Nodal Loads """
        #
        # beam FER
        conn = create_connection(self.db_file)
        with conn:
            dfbeam = pull_ENL_df(conn, self._component)
        return dfbeam
    #
    #
    def Fn(self, plane):
        """
        Global matrix consisting of summation of force & displacement 
        """
        # FIXME: step needs to be sorted
        columns = [*plane.hforce, *plane.hdisp]
        headgrp = ['load_name', 'mesh_name',
                   'load_id', 'load_level',
                   'load_title','load_system',
                   'node_name', 'node_index']
        # beam FER
        conn = create_connection(self.db_file)
        with conn:
            dfbeam = pull_FER_df(conn, self._component)
        #
        dfbeam = dfbeam.groupby(headgrp)[columns].sum()
        #
        # Node load
        dfnodal = self._nodes.df
        dfnodal = dfnodal.groupby(headgrp)[columns].sum()
        #
        #dfnodal[columns] = dfnodal[columns].add(dfbeam[columns], fill_value=0)
        Fn_df = dfnodal.add(dfbeam, fill_value=0)
        Fn_df.reset_index(inplace=True)
        #
        return Fn_df.reindex(columns=['load_name', 'mesh_name', 
                                      'load_id', 'load_level',
                                      'load_title','load_system',
                                      #'load_title', 'load_comment', 'load_system',
                                      #'element_name', 'node_name', 'node_index',
                                      'node_name', 'node_index', 
                                      'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                                      'x', 'y', 'z', 'rx', 'ry', 'rz'])
#
#
def pull_ENL_df(conn, component: int):
    """ Equivalent Nodal Loads """
    df = pull_FER_data(conn, component)
    df = df[['load_name', 'mesh_name', 
             'load_title', 'load_level',
             'load_id', 'load_system', 'load_comment',
             'element_name',
             'node_name', 'node_index', 
             'load_type',
             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
             'Psi', 'B', 'Tw', 'step']]
    return df
#
#
def pull_FER_df(conn, component: int):
    """ """
    df = pull_FER_data(conn, component)
    #
    return df[['load_name', 'mesh_name', 
             'load_title', 'load_level',
             'load_id', 'load_system', 'load_comment',
             'element_name',
             'node_name', 'node_index', 
             'load_type',
             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
             'x', 'y', 'z', 'rx', 'ry', 'rz']]
#
#
def pull_FER_data(conn, component: int):
    """ """
    query = (component, )
    table = "SELECT Load.name AS load_name, \
                    Mesh.name AS mesh_name, \
                    Load.title AS load_title, \
                    Node.name AS node_name, \
                    Node.idx as node_index, \
                    Element.name as element_name, \
                    LoadBeamFER.* \
            FROM Load, Node, Element, LoadBasic, LoadBeamFER, Mesh \
            WHERE LoadBeamFER.basic_id = LoadBasic.number \
            AND LoadBasic.load_id = Load.number \
            AND LoadBeamFER.node_id = Node.number \
            AND LoadBeamFER.element_id =  Element.number \
            AND Load.mesh_id = Mesh.number \
            AND Mesh.number = ?;"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    rows = cur.fetchall()
    #
    #
    cols = ['load_name', 'mesh_name', 
            'load_title',
            'node_name', 'node_index', 
            'element_name',
            'number', 
            'load_id', 'element_id', 'node_id',
            'load_comment', 'load_system','load_type',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
            'x', 'y', 'z', 'rx', 'ry', 'rz',
            'Psi', 'B', 'Tw', 'step']
    #
    # dataframe
    db = DBframework()
    df = db.DataFrame(data=rows, columns=cols)
    df['load_level'] = 'basic'    
    #
    return df
#      
#
#
class LoadTypeSQL(LoadTypeBasic):
    """
    """
    __slots__ = ['_node', '_beam', '_selfweight', '_component', 
                 'name', 'number', 'title', '_db_file']

    def __init__(self, load_name: str|int,
                 component: int, 
                 bd_file:str):
        """
        """
        self.name = load_name
        self._db_file = bd_file
        self._component = component
        #
        self.number, self.title = self._load_spec(self.name)
        #
        self._node = NodeLoadItemSQL(load_name=load_name,
                                     component=component, 
                                     db_file=self._db_file)
        #
        self._beam = BeamLoadItemSQL(load_name=load_name,
                                     component=component, 
                                     db_file=self._db_file)
        #
        self._selfweight = SelfWeight()          
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._component, )
        table = 'SELECT Load.name \
                 FROM Load, LoadBasic\
                 WHERE LoadBasic.load_id = Load.number \
                 AND mesh_id = ?;'
        conn = create_connection(self._db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]    
    #
    def _load_spec(self, load_name: str|int):
        """ """
        query = (load_name, self._component)
        table = 'SELECT Load.number, Load.title\
                 FROM Load, LoadBasic \
                 WHERE Load.name = ? \
                 AND LoadBasic.load_id = Load.number \
                 AND mesh_id = ?;'
        #
        conn = create_connection(self._db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchone()
        #
        return items
    #
    @property
    def _component_name(self):
        """ component name """
        query = (self._component, )
        table = 'SELECT name \
                 FROM Component WHERE number = ?;'
        conn = create_connection(self._db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchone()
        return items[0]
    #  
#
#
# ---------------------------------
#
