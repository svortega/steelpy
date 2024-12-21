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
from steelpy.ufo.load.process.load_case import BasicLoadCase, BasicLoadType 
from steelpy.ufo.load.sqlite.beam import BeamLoadItemSQL, BeamLoadGloabalSQL
from steelpy.ufo.load.sqlite.node import  NodeLoadItemSQL, NodeLoadGlobalSQL
#
# steelpy.f2uModel
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework
#
import numpy as np
#
#
class BasicLoadSQL(BasicLoadCase):
    __slots__ = ['db_file', '_component',
                 'gravity', '_nodes', '_beams']

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
        if isinstance(load_name, (int, np.integer)):
            load_name = int(load_name)
        
        try:
            index = self._labels.index(load_name)
            return LoadTypeSQL(load_name=load_name,
                               component=self._component, 
                               bd_file=self.db_file)
        except ValueError:
            raise IOError("load case not defined")

    #
    # -----------------------------------------------
    # SQL operations
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
    #
    # -----------------------------------------------
    # Loading Operations
    #
    def FER(self, elements):
        """
        Beams' Fixed End Reactions (FER)
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
    def FER_ENL(self):
        """
        Equivalent Beam Nodal Loads

        Return:
            Dataframe with global node loads
        """
        conn = create_connection(self.db_file)
        with conn:
            dfbeam = pull_ENL_FER(conn, self._component)
        return dfbeam
    #
    def FER_END(self):
        """
        Equivalent Beam Nodal Displacements

        Return:
            Dataframe with global node displacements
        """
        conn = create_connection(self.db_file)
        with conn:
            dfbeam = pull_END_FER(conn, self._component)
        return dfbeam
    #
    # -----------------------------------------------
    #
    def NF_global(self, plane):
        """
        Nodal Force (total)

        Return:
            Dataframe consisting of summation of
            node and elements (FER) nodal force
        """
        # FIXME: step needs to be sorted
        columns = [*plane.hforce]
        head = ['load_name', 'mesh_name',
                'load_id', 'load_level',
                'load_title','load_system',
                'node_name', 'node_index']
        # beam FER
        conn = create_connection(self.db_file)
        with conn:
            beam_fer = pull_ENL_FER(conn, self._component)
        #
        if not beam_fer.empty:
            # 3D sign correction hack
            #beam_fer['Mz'] *= -1
            Fn_df = beam_fer.groupby(head)[columns].sum()
            #beam_fer.reset_index(inplace=True)
        #
        # Node load
        node_force =self._nodes.force
        if node_force.empty:
            if beam_fer.empty:
                db = DBframework()
                return db.DataFrame()                
            #Fn_df = beam_grp
        else:
            node_grp = node_force.groupby(head)[columns].sum()
            if beam_fer.empty:
                Fn_df = node_grp
            else:
                Fn_df = node_grp.add(Fn_df, fill_value=0, axis='columns')
        #node_df = self._nodes.df
        #node_df = node_df.groupby(head)[columns].sum()
        #node_grp = node_df.groupby(head)
        #node_grp = node_grp.getgroup('node_name')
        #node_df.reset_index(inplace=True)
        #
        #dfnodal[columns] = dfnodal[columns].add(dfbeam[columns], fill_value=0)
        #head = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        #Fn_df = node_grp.add(beam_grp, fill_value=0, axis='columns')
        Fn_df.reset_index(inplace=True)
        #
        return Fn_df.reindex(columns=['load_name', 'mesh_name', 
                                      'load_id', 'load_level',
                                      'load_title','load_system',
                                      'node_name', 'node_index', 
                                      'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'])
    #
    def ND_global(self, plane):
        """
        Nodal Displacement (total)

        Return:
            Dataframe consisting of summation of
            node and elements (FER) nodal displacement
        """
        db = DBframework()
        # FIXME: step needs to be sorted
        #
        columns = [*plane.hdisp]
        head = ['load_name', 'mesh_name',
                'load_id', 'load_level',
                'load_title','load_system',
                'node_name', 'node_index']
        # beam FER
        conn = create_connection(self.db_file)
        with conn:
            beam_fer = pull_END_FER(conn, self._component)
        #
        #Fn_df = None
        if not beam_fer.empty:
            # 3D sign correction hack
            #beam_fer['rz'] *= -1            
            Fn_df = beam_fer.groupby(head)[columns].sum()
            Fn_df = Fn_df[(Fn_df != 0).any(axis=1)]
            if Fn_df.empty:
                beam_fer = db.DataFrame()
        #
        # Node displacement
        node_disp = self._nodes.displacement
        if node_disp.empty:
            if beam_fer.empty:
                return db.DataFrame()
            pass
        else:
            node_grp = node_disp.groupby(head)[columns].sum()
            #dfnodal = self._nodes.df
            #dfnodal = dfnodal.groupby(head)[columns].sum()
            #
            if beam_fer.empty:
                Fn_df = node_grp
            else:
                #dfnodal[columns] = dfnodal[columns].add(dfbeam[columns], fill_value=0)
                Fn_df = node_grp.add(Fn_df, fill_value=0, axis='columns')
        #
        Fn_df.reset_index(inplace=True)
        #
        return Fn_df.reindex(columns=['load_name', 'mesh_name', 
                                      'load_id', 'load_level',
                                      'load_title','load_system',
                                      'node_name', 'node_index', 
                                      'x', 'y', 'z', 'rx', 'ry', 'rz'])    
#
#
def pull_ENL_FER(conn, component: int):
    """
    Return:
        Equivalent Nodal Loads in dataframe form"""
    df = pull_FER_data(conn, component)
    df = df[['load_name', 'mesh_name', 
             'load_title', 'load_level',
             'load_id', 'load_system', 'load_comment',
             'element_name',
             'node_name', 'node_index', 
             'load_type',
             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz', 
             'Psi', 'B', 'Tw', 'step']]
    # 3D sign correction
    #df['Mz'] *= -1
    return df
#
def pull_END_FER(conn, component: int):
    """
    Return :
        Equivalent Nodal Displacement in dataframe form"""
    df = pull_FER_data(conn, component)
    df = df[['load_name', 'mesh_name',
               'load_title', 'load_level',
               'load_id', 'load_system', 'load_comment',
               'element_name',
               'node_name', 'node_index',
               'load_type',
               'x', 'y', 'z', 'rx', 'ry', 'rz', 'step']]
    # 3D sign correction
    #df['rz'] *= -1
    return df
#
def pull_FER_df(conn, component: int):
    """ Return Fix End Forces in dataframe form"""
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
    """ Return Fix End Forces from database in dataframe form"""
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
class LoadTypeSQL(BasicLoadType):
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
