#
# Copyright (c) 2009 steelpy
# 

# Python stdlib imports
from __future__ import annotations
#from dataclasses import dataclass
import time

from numpy.ma.extras import column_stack

#from typing import NamedTuple
#import re

# package imports
# steelpy.f2uModel
from steelpy.ufo.load.process.actions import SelfWeight
from steelpy.ufo.load.process.load_case import BasicLoadCase, BasicLoadRoot
from steelpy.ufo.load.process.node.utils import find_NodeLoad_item
from steelpy.ufo.load.process.beam.utils import find_BeamLoad_item
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
    __slots__ = ['db_file', '_mesh_id', '_name',
                 'gravity', '_nodes', '_beams']

    #
    def __init__(self, db_file:str, #plane: NamedTuple,
                 mesh_id: int, name:int|str) -> None:
        """
        """
        super().__init__(name)
        self._mesh_id = mesh_id
        self.db_file = db_file
        #
        self._nodes = NodeLoadGlobalSQL(mesh_id=self._mesh_id,
                                        db_file=self.db_file)
        #
        self._beams = BeamLoadGloabalSQL(mesh_id=self._mesh_id,
                                         db_file=self.db_file)
        #self._beams = BeamLoadItemSQL(load_name='*',
        #                              mesh_id=self._mesh_id,
        #                              db_file=self.db_file)        
        #
        conn = create_connection(self.db_file)
        with conn: 
            self._new_table(conn)        
    #
    @property
    def _labels(self):
        """ """
        query = ('basic', self._mesh_id, )
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
            return BasicLoadTypeSQL(load_name=load_name,
                                    mesh_id=self._mesh_id,
                                    #name=self._name,
                                    bd_file=self.db_file)
        except ValueError:
            raise IOError("load case not defined")

    #
    # -----------------------------------------------
    # SQL operations
    #
    def _push_load(self, conn, load_name:int|str, load_title:str):
        """ """
        query = (load_name,  self._mesh_id, "basic", load_title, 'ufo', None)
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
            dfbeam = pull_ENL_FER(conn, self._mesh_id)
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
            dfbeam = pull_END_FER(conn, self._mesh_id)
        return dfbeam
    #
    # -----------------------------------------------
    #
    def _Fnt(self):
        """
        Total Nodal Force (user input + FER)

        Return:
            Dataframe consisting of summation of
            node and elements (FER) nodal force
        """
        # FIXME: step needs to be sorted
        #columns = [*plane.hforce]
        columns = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        head = ['load_name', 'mesh_name',
                'load_id', 'load_level',
                'load_title','system',
                'node_name', 'node_index']
        # beam FER
        conn = create_connection(self.db_file)
        with conn:
            beam_fer = pull_ENL_FER(conn, self._mesh_id)
        #
        if not beam_fer.empty:
            # Select forces
            Fn_df = beam_fer.groupby(head)[columns].sum()
            #beam_fer.reset_index(inplace=True)
        #
        # Node load
        node_force =self._nodes.force
        if node_force.empty:
            if beam_fer.empty:
                #db = DBframework()
                #return db.DataFrame()
                return beam_fer
            #Fn_df = beam_grp
        else:
            node_grp = node_force.groupby(head)[columns].sum()
            if beam_fer.empty:
                Fn_df = node_grp
            else:
                Fn_df = node_grp.add(Fn_df, fill_value=0, axis='columns')
        #
        head += columns # ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        Fn_df.reset_index(inplace=True)
        return Fn_df.reindex(columns=head)
    #
    def _Dnt(self):
        """
        Total Nodal Displacement (user input + FER)

        Return:
            Dataframe consisting of summation of
            node and elements (FER) nodal displacement
        """
        #db = DBframework()
        # FIXME: step needs to be sorted
        #
        #columns = [*plane.hdisp]
        columns = ['x', 'y', 'z', 'rx', 'ry', 'rz']
        head = ['load_name', 'mesh_name',
                'load_id', 'load_level',
                'load_title','system',
                'node_name', 'node_index']
        # TODO: determine if FER-displacements are required
        # beam FER
        conn = create_connection(self.db_file)
        with conn:
            beam_fer = pull_END_FER(conn, self._mesh_id)
        #
        #Fn_df = None
        if not beam_fer.empty:
            # Select displacement
            Dn_df = beam_fer.groupby(head)[columns].sum()
            Dn_df = Dn_df[(Dn_df != 0).any(axis=1)]
            #if Dn_df.empty:
            #    beam_fer = db.DataFrame()
        #
        # Node displacement
        node_disp = self._nodes.displacement
        if node_disp.empty:
            if beam_fer.empty:
                return beam_fer
            pass
        else:
            node_grp = node_disp.groupby(head)[columns].sum()
            #dfnodal = self._nodes.df
            #dfnodal = dfnodal.groupby(head)[columns].sum()
            #
            if beam_fer.empty:
                Dn_df = node_grp
            else:
                #dfnodal[columns] = dfnodal[columns].add(dfbeam[columns], fill_value=0)
                Dn_df = node_grp.add(Dn_df, fill_value=0, axis='columns')
        #
        head += columns
        Dn_df.reset_index(inplace=True)
        return Dn_df.reindex(columns=head)
#
#
def pull_ENL_FER(conn, mesh_id: int):
    """
    Return:
        Equivalent Nodal Loads in dataframe form"""
    df = pull_FER_data(conn, mesh_id)
    df = df[['load_name', 'mesh_name', 
             'load_title', 'load_level',
             'load_id', 'system', 'comment',
             'element_name',
             'node_name', 'node_index', 
             'load_type',
             'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz', 
             'Psi', 'B', 'Tw', 'step']]
    # 3D sign correction
    #df['Mz'] *= -1
    return df
#
def pull_END_FER(conn, mesh_id: int):
    """
    Return :
        Equivalent Nodal Displacement in dataframe form"""
    df = pull_FER_data(conn, mesh_id)
    df = df[['load_name', 'mesh_name',
               'load_title', 'load_level',
               'load_id', 'system', 'comment',
               'element_name',
               'node_name', 'node_index',
               'load_type',
               'x', 'y', 'z', 'rx', 'ry', 'rz', 'step']]
    # 3D sign correction
    #df['rz'] *= -1
    return df
#
def pull_FER_df(conn, mesh_id: int):
    """ Return Fix End Forces in dataframe form"""
    df = pull_FER_data(conn, mesh_id)
    #
    return df[['load_name', 'mesh_name', 
               'load_title', 'load_level',
               'load_id', 'system', 'comment',
               'element_name',
               'node_name', 'node_index',
               'load_type',
               'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
               'x', 'y', 'z', 'rx', 'ry', 'rz']]
#
#
def pull_FER_data(conn, mesh_id: int):
    """ Return Fix End Forces from database in dataframe form"""
    query = (mesh_id, )
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
            'comment', 'system','load_type',
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
class BasicLoadTypeSQL(BasicLoadRoot):
    """
    """
    __slots__ = ['name', 'number', 'title', '_db_file',
                 '_node', '_beam', '_selfweight', '_mesh_id']

    def __init__(self, load_name: str|int,
                 mesh_id: int,
                 #name:int|str,
                 bd_file:str):
        """
        """
        super().__init__(load_name)
        self._db_file = bd_file
        self._mesh_id = mesh_id
        #
        self.number, self.title = self._load_spec(self.name)
        #
        self._node = NodeLoadItemSQL(load_name=load_name,
                                     mesh_id=mesh_id, 
                                     db_file=self._db_file)
        #
        self._beam = BeamLoadItemSQL(load_name=load_name,
                                     mesh_id=mesh_id, 
                                     db_file=self._db_file)
        #
        self._selfweight = SelfWeight()
    #
    @property
    def _labels(self):
        """ """
        query = (self._mesh_id, )
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
    def __setitem__(self, load_name: str | int,
                    properties: list[float]) -> None:
        """
        """
        try:
            self._labels.index(load_name)
            self.name = load_name
            self.number, self.title = self._load_spec(self.name)
            #self.title = properties[0]
            #self.number = properties[1]
        except ValueError:
            raise Exception(f'Load {load_name} already exist')
    
    
    def __getitem__(self, load_name: int|str):
        """
        node_name : node number
        """
        try:
            index = self._labels.index(load_name)
            return self
        except ValueError:
            raise KeyError(f'   *** Load {load_name} does not exist')
    #
    # -----------------------------------------------
    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += f"Load Name : {str(self.name):12s}  Number : {self.number:8.0f}  Title : {self.title}\n"
        # node load
        if self._node:
            output += f"--- Node\n"
            output += self._node.__str__()
        # beam line
        if self._beam:
            output += f"--- Beam \n"
            output += self._beam.__str__()
        #
        if self._selfweight:
            output += f"--- Gravity/Selfweight\n"
            output += self._selfweight.__str__()
        #
        #output += "\n"
        # print('---')
        return output
    # 
    # -----------------------------------------------    
    #
    def _load_spec(self, load_name: str|int):
        """ """
        query = (load_name, self._mesh_id)
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
        1 / 0
        query = (self._mesh_id, )
        table = 'SELECT Component.name \
                 FROM Component, Mesh \
                 WHERE Component.number = Mesh.component_id \
                 AND Mesh.number = ?;'
        conn = create_connection(self._db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            items = cur.fetchone()
        return items[0]
    #
    # -----------------------------------------------
    #
    def function(self, beam_name: str | int,
                 load_name: str | int, 
                 steps: int, Pa: float,
                 factor: float = 1.0)->DBframework.DataFrame:
        """
        beam_name:
        load_name:
        steps: beam steps where load function will be calculated
        factor: load factor

        Return:
            Beam's load functions
        """
        beam = self._beam[beam_name]
        bfunction = beam.function(steps=steps,
                                  Pa=Pa, factor=factor)
        #
        # Axial   [FP, blank, blank, Fu]
        # torsion [T, B, Psi, Phi, Tw]
        # Bending [V, M, theta, w]
        #
        # [Fx, Fy, Fz, Mx, My, Mz]
        # [V, M, w, theta]
        header = ['load_name', 'mesh_name',
                  'comment', 'load_type',
                  'load_level', 'system',
                  'element_name', 'length',
                  'axial', 'torsion', 'VM_inplane', 'VM_outplane']
        #
        #          'FP', 'blank1', 'blank2', 'Fu',
        #          'T', 'B', 'Psi', 'Phi', 'Tw',
        #          'Vy', 'Mz', 'theta_y', 'w_y',
        #          'Vz', 'My', 'theta_z', 'w_z']
        df = DBframework()
        bfunction = df.DataFrame(data=bfunction, columns=header, index=None)
        grp_bfunc = bfunction.groupby(['load_name', 'element_name', 'length'])
        #grp_bfunc.get_group(())
        bfunction = grp_bfunc[['axial', 'torsion', 'VM_inplane', 'VM_outplane']].sum()
        bfunction.reset_index(inplace=True)
        grp_bfunc = bfunction.groupby(['load_name', 'element_name'])
        grp_bfunc = grp_bfunc.get_group((load_name, beam_name, )).reset_index()
        return grp_bfunc
    #
    # -----------------------------------------------
    #
    #
    def node(self, values:tuple|list|None=None,
             df=None):
        """ Nodal load"""
        if values:
            # Input data for specific basic node load
            if isinstance(values, dict):
                columns = list(values.keys())
                header = {item: find_NodeLoad_item(item)
                          for item in columns}
                values ={header[key]: item
                         for key, item in values.items()}
                nodeid = values['node']
                if isinstance(nodeid, (list, tuple)):
                    db = DBframework()
                    dfnew = db.DataFrame(data=values)
                    dfnew['load'] = self.name
                    self._node.df = dfnew
                else:
                    self._node[nodeid] = values
            elif isinstance(values, (list, tuple, dict)):
                if isinstance(values[0], (list, tuple, dict)):
                    for item in values:
                        if isinstance(item, dict):
                            header = {item: find_NodeLoad_item(item)
                                      for item in item}
                            update = {header[key]: item
                                      for key, item in item.items()}                           
                            nodeid = update['node']
                            load = update
                        elif isinstance(item, (list, tuple)):
                            nodeid = item[0]
                            load = item[1:]
                        #
                        self._node[nodeid] = load
                else:
                    self._node[values[0]] = values[1:]
        #
        # dataframe input
        try:
            columns = list(df.columns)
            header = {key: find_NodeLoad_item(key)
                      for key in columns}
            df.rename(columns=header, inplace=True)
            #columns = list(df.columns)
            df['load'] = self.name
            self._node.df = df            
        except AttributeError:
            pass
        #
        return self._node
    #
    #
    def beam(self, values:tuple|list|dict|None=None,
             df=None):
        """ beam loading """
        if values:
            if isinstance(values, dict):
                columns = list(values.keys())
                header = {item: find_BeamLoad_item(item)
                          for item in columns}
                values ={header[key]: item
                         for key, item in values.items()}                
                beamid = values['beam']
                if isinstance(beamid, (list, tuple)):
                    db = DBframework()
                    dfnew = db.DataFrame(data=values)
                    dfnew['load'] = self.name
                    self._beam.df = dfnew
                else:
                    self._beam[beamid] = values
            elif isinstance(values, (list, tuple)):
                if isinstance(values[0], (list, tuple, dict)):
                    for item in values:
                        if isinstance(item, dict):
                            header = {item: find_BeamLoad_item(item)
                                      for item in item}
                            update = {header[key]: item
                                      for key, item in item.items()}
                            #load_name = update['load']                            
                            beamid = update['beam']
                            load =  update
                        elif isinstance(item, (list, tuple)):
                            beamid = item[0]
                            load =  item[1:]
                        #
                        self._beam[beamid] = load
                else:
                    self._beam[values[0]] = values[1:]
        # dataframe input
        try:
            columns = list(df.columns)
            header = {key: find_bload_item(key)
                      for key in columns}
            df.rename(columns=header, inplace=True)              
            #columns = list(df.columns)
            df['load'] = self.name
            self._beam.df = df
            return 
        except AttributeError:
            pass
        #
        return self._beam
    #    
    #
    # -----------------------------------------------
    #
    #
#
#
# ---------------------------------
#
