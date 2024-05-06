#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
from typing import NamedTuple
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
    
    __slots__ = ['db_file', '_plane', '_component']
                # '_nodes', '_beams']
    #
    def __init__(self, db_file:str, plane: NamedTuple,
                 component: int) -> None:
        """
        """
        super().__init__()
        #
        self.db_file = db_file
        self._plane = plane
        self._component = component
        #
        self._nodes = NodeLoadGlobalSQL(component=self._component,
                                        db_file=self.db_file)
        #
        self._beams = BeamLoadGloabalSQL(component=self._component,
                                         plane=self._plane, 
                                         db_file=self.db_file)
        #
        conn = create_connection(self.db_file)
        with conn: 
            self._create_table(conn)        
    #
    @property
    def _labels(self):
        """ """
        query = ('basic', self._component, )
        table = "SELECT Load.name FROM Load \
                  WHERE level = ? AND component_id = ?;"
        conn = create_connection(self.db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchall()
        return [item[0] for item in items]
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
                               plane=self._plane,
                               component=self._component, 
                               bd_file=self.db_file)
        except ValueError:
            raise IOError("load case not defined")

    #
    # -----------------------------------------------
    #
    def _push_load(self, conn, load_name:int|str, load_title:str):
        """ """
        #load_name = str(load_name)
        project = (load_name,  self._component, "basic", load_title)
        table = 'INSERT INTO Load(name, component_id, level, title) \
                 VALUES(?,?,?,?)'
        cur = conn.cursor()
        cur.execute(table, project)
        #return cur.lastrowid
    #
    #
    def _create_table(self, conn):
        """ """
        table = "CREATE TABLE IF NOT EXISTS Load(\
                number INTEGER PRIMARY KEY NOT NULL,\
                name NOT NULL,\
                component_id INTEGER NOT NULL REFERENCES Component(number), \
                level TEXT NOT NULL,\
                title TEXT );"
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
        header = ['load_name', 'component_name',
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
        print('--- FER Operation')
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
    def Fn(self):
        """
        Global matrix consisting of summation of force & displacement 
        """
        columns = [*self._plane.hforce, *self._plane.hdisp]
        headgrp = ['load_name', 'component_name',
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
        return Fn_df.reindex(columns=['load_name', 'component_name', 
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
    #
    df = df[['load_name', 'component_name', 
             'load_title', 'load_level',
             'load_id', 'load_system', 'load_comment',
             'element_name',
             'node_name', 'node_index', 
             'load_type',
             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
             #'x', 'y', 'z', 'rx', 'ry', 'rz',
             'Psi', 'B', 'Tw']]
    return df
#
#
def pull_FER_df(conn, component: int):
    """ """
    df = pull_FER_data(conn, component)
    #
    return df[['load_name', 'component_name', 
             'load_title', 'load_level',
             'load_id', 'load_system', 'load_comment',
             'element_name',
             'node_name', 'node_index', 
             'load_type',
             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
             'x', 'y', 'z', 'rx', 'ry', 'rz']]
             #'Psi', 'B', 'Tw']]
#
#
def pull_FER_data(conn, component: int):
    """ """
    query = (component, )
    table = "SELECT Load.name AS load_name, \
                    Component.name AS component_name, \
                    Load.title AS load_title, \
                    Node.name AS node_name, \
                    Node.mesh_idx as node_index, \
                    Element.name as element_name, \
                    LoadBeamFER.* \
            FROM Load, Node, Element, LoadBeamFER, Component \
            WHERE LoadBeamFER.load_id = Load.number \
            AND LoadBeamFER.node_id = Node.number \
            AND LoadBeamFER.element_id =  Element.number \
            AND Load.component_id = Component.number \
            AND Component.number = ?;"
    #
    cur = conn.cursor()
    cur.execute(table, query)
    rows = cur.fetchall()
    #
    #
    cols = ['load_name', 'component_name', 
            'load_title',
            'node_name', 'node_index', 
            'element_name',
            'number', 
            'load_id', 'element_id', 'node_id',
            'load_comment', 'load_system','load_type',
            'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
            'x', 'y', 'z', 'rx', 'ry', 'rz',
            'Psi', 'B', 'Tw']
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
                 'name', 'number', 'title', '_db_file', '_plane']

    def __init__(self, load_name: str|int,
                 plane: NamedTuple,
                 component: int, 
                 bd_file:str):
        """
        """
        self.name = load_name
        self._db_file = bd_file
        self._plane = plane
        self._component = component
        #
        self.number, self.title = self._load_spec(self.name)
        #
        self._node = NodeLoadItemSQL(load_name=load_name,
                                     component=component, 
                                     db_file=self._db_file)
        #
        self._beam = BeamLoadItemSQL(load_name=load_name,
                                     plane=self._plane,
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
                 FROM Load WHERE component_id = ?;'
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
                 FROM Load \
                 WHERE Load.name = ? \
                 AND component_id = ?;'
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
