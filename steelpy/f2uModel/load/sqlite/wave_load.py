#
# Copyright (c) 2009-2023 fem2ufo
# 

# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
from dataclasses import dataclass
from typing import NamedTuple
#import re
from operator import itemgetter
from itertools import groupby

# package imports
# steelpy.f2uModel
#from steelpy.f2uModel.load.process.actions import SelfWeight
#from steelpy.f2uModel.load.process.basic_load import BasicLoadBasic, LoadTypeBasic
from steelpy.f2uModel.load.sqlite.beam import BeamLoadItemSQL,  BeamToNodeSQL
#from steelpy.f2uModel.load.sqlite.node import  NodeLoadItemSQL
#
# steelpy.f2uModel.load
from ..process.beam import BeamLoadItem #, BeamLoad
#
from steelpy.f2uModel.mesh.sqlite.beam import BeamItemSQL
# steelpy.f2uModel
from steelpy.f2uModel.mesh.process.process_sql import (create_connection, create_table,
                                                       get_load_data, check_element)
#
from steelpy.f2uModel.load.process.beam import LineBeam, BeamLoad
#
import pandas as pd
from steelpy.process.dataframe.main import DBframework
#
#

#
#@dataclass
class WaveLoadItemSQL(BeamLoadItem):
    """ """
    __slots__ = ['_wave','_bd_file', '_name', '_load']
    
    def __init__(self, load_name: int|str, bd_file:str):
        """ """
        super().__init__()
        #
        self._wave = []
        self._design_load = []
        self._bd_file = bd_file
        self._name = load_name
        #
        self._node_eq = BeamToNodeSQL(load_name=load_name, 
                                        bd_file=bd_file)
        #
        #self._load = BeamLoad()
        #
        conn = create_connection(self._bd_file)
        with conn:        
            self._create_table(conn)
            self._create_node_table(conn)
    #
    #
    def __setitem__(self, load_name: int|str,
                    wave_load: list) -> None:
        """
        """
        try:
            self._labels.index(load_name)
            raise Exception('wave load case {:} already exist'.format(load_name))
        except ValueError:
            self._labels.append(load_name)
            self._wave.append(wave_load[0])
            self._design_load.append(wave_load[1])
            #conn = create_connection(self._bd_file)
            #1 / 0
    #
    def __getitem__(self, load_name: int | str):
        """
        """
        try:
            index = self._labels.index(load_name)
        except ValueError:
            raise KeyError('Invalid load name : {:}'.format(load_name))
        #
        #conn = create_connection(self._bd_file)
        #1 / 0        
        return self._wave[index]

    #
    #
    def _create_table(self, conn) -> None:
        """ """
        _table_element_line = "CREATE TABLE IF NOT EXISTS tb_WaveLoad(\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                load_number INTEGER NOT NULL REFERENCES tb_Load(number),\
                                element_number INTEGER NOT NULL REFERENCES tb_Elements(number),\
                                title TEXT,\
                                system INTEGER NOT NULL,\
                                L_end1 DECIMAL,\
                                qx1 DECIMAL,\
                                qy1 DECIMAL,\
                                qz1 DECIMAL,\
                                qx1i DECIMAL,\
                                qy1i DECIMAL,\
                                qz1i DECIMAL,\
                                L_end2 DECIMAL,\
                                qx2 DECIMAL,\
                                qy2 DECIMAL,\
                                qz2 DECIMAL,\
                                qx2i DECIMAL,\
                                qy2i DECIMAL,\
                                qz2i DECIMAL,\
                                x DECIMAL,\
                                y DECIMAL,\
                                z DECIMAL);"
        #
        create_table(conn, _table_element_line)
    #
    def _push_load(self, conn, load: list):
        """ """
        #
        sql = 'INSERT INTO tb_WaveLoad(load_number, element_number,\
                                        title, system,\
                                        L_end1, qx1, qy1, qz1, qx1i, qy1i, qz1i,\
                                        L_end2, qx2, qy2, qz2, qx2i, qy2i, qz2i,\
                                        x, y, z)\
                                        VALUES(?,?,?,?,\
                                               ?,?,?,?,?,?,?,\
                                               ?,?,?,?,?,?,?,\
                                               ?,?,?)'
        cur = conn.cursor()
        cur.executemany(sql, load)
    #
    def get_wave_load(self, conn, load_name:int):
        """ """
        db = DBframework()      
        # get beam data
        #beam = check_element(conn, beam_name)
        #beam_number = beam[0]
        # beam line load
        load_data = get_load_data(conn, self._name, load_type='basic')
        load_number = load_data[0]
        #
        cur = conn.cursor()
        cur.execute("SELECT tb_Load.name, tb_Elements.name, tb_WaveLoad.*\
                    FROM tb_WaveLoad, tb_Elements, tb_Load\
                    WHERE tb_WaveLoad.load_number = {:}\
                    AND tb_WaveLoad.element_number = tb_Elements.number\
                    AND tb_WaveLoad.load_number = tb_Load.number;"
                    .format(load_number))
        rows = cur.fetchall()
        #
        cols = ['load_name', 'element_name', 'number', 'load_number', 'element_number', 
                'load_comment', 'load_system',
                'L_end1', 'qx1', 'qy1', 'qz1', 'qx1i', 'qy1i', 'qz1i',
                'L_end2', 'qx2', 'qy2', 'qz2', 'qx2i', 'qy2i', 'qz2i',
                'x', 'y', 'z']        
        df = db.DataFrame(data=rows, columns=cols)
        df = df[['load_name', 'load_number', 
                 'load_comment', 'load_system',
                 'element_name', 'element_number', 
                 'L_end1', 'qx1', 'qy1', 'qz1',
                 'L_end2', 'qx2', 'qy2', 'qz2',
                 'x', 'y', 'z']]
        return df
    #
    def get_elements(self,conn):
        """ """
        cur = conn.cursor()
        cur.execute ("SELECT * FROM tb_Elements;")
        row = cur.fetchall()
        return row
    #
    #
    #
    # -----------------------------------------------
    #
    #
    def _create_node_table(self, conn) -> None:
        """ """
        _table_element_point = "CREATE TABLE IF NOT EXISTS tb_LoadBeamToNode(\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                load_number INTEGER NOT NULL REFERENCES tb_Load(number),\
                                title TEXT,\
                                system TEXT,\
                                element_number INTEGER NOT NULL REFERENCES tb_Elements(number),\
                                node_number INTEGER NOT NULL REFERENCES tb_Nodes(number),\
                                fx DECIMAL,\
                                fy DECIMAL,\
                                fz DECIMAL,\
                                mx DECIMAL,\
                                my DECIMAL,\
                                mz DECIMAL,\
                                fxi DECIMAL,\
                                fyi DECIMAL,\
                                fzi DECIMAL,\
                                mxi DECIMAL,\
                                myi DECIMAL,\
                                mzi DECIMAL);"
        #
        create_table(conn, _table_element_point)    
    #
    #
    def _push_node_load(self, conn,
                        node_load:list[float]):
        """ """
        #print("-->")
        #beam = check_element(conn, beam_name)
        #beam_number = beam[0]
        #
        #node_name = node_load.pop(0)
        #node = check_nodes(conn, node_name)
        #node_number = node[0]
        #
        #load_data = get_load_data(conn, self._name, load_type='basic')
        #load_number = load_data[0]
        #
        #try:
        #    1 / load_system
        #    load_system = 'local'
        #except ZeroDivisionError:
        #    raise RuntimeError('node load in global system')
        #
        #project = (load_number, load_title, load_system,
        #           beam_number, node_number,
        #           *node_load,
        #           'NULL', 'NULL', 'NULL',
        #           'NULL', 'NULL', 'NULL')
        #
        sql = 'INSERT INTO tb_LoadBeamToNode(load_number, title, system, \
                                             element_number,\
                                             node_number, fx, fy, fz, mx, my, mz,\
                                             fxi, fyi, fzi, mxi, myi, mzi)\
                                            VALUES(?,?,?,?,?,\
                                                   ?,?,?,?,?,?,?,\
                                                   ?,?,?,?,?)'
        cur = conn.cursor()
        cur.executemany(sql, node_load)    
    #
    #
    # -----------------------------------------------
    #
    def process(self, load_name: int | str):
        """ """
        print('# Calculating wave beam forces')
        wave = self.__getitem__(load_name)
        wname = f'{wave.name}_{wave.title}'
        wload = wave.load()
        #
        df_bload = self.get_beam_load(wname, wload)
        #
        self.df(data=df_bload)
        #
        print('# End process')
    #
    #
    def get_beam_load(self, wname, wload):
        """ """
        #
        conn = create_connection(self._bd_file)
        with conn:
            labels = self.get_elements(conn)
        #
        df_bload = None
        for idx, key in enumerate(labels):
            beam =  BeamItemSQL(key[1], self._bd_file)
            df_load = wload.wave_force(beam=beam)
            df_load['load_name'] = self._name
            df_load['load_title'] = wname
            #df_load['load_title'] = df_load.apply(lambda row: f"{wname}_{round(row.x, 2)}_{round(row.y, 2)}_{round(row.z, 2)}", axis=1)
            df_load['load_system'] = 1
            try:
                1/idx
                df_bload = pd.concat([df_bload, df_load], ignore_index=True)
            except ZeroDivisionError:
                df_bload = df_load
        #
        return df_bload        
    #
    # -----------------------------------------------
    #
    def df(self, data):
        """ """
        #
        conn = create_connection(self._bd_file)
        #
        beams = data.groupby(['element_type']).get_group('beam')
        grpbeam = beams.groupby(['element_name', 'load_name'])
        #
        for key, item in grpbeam:
            subgrp = item.groupby('load_type')
            line = subgrp.get_group('line')
            #
            beam_name = key[0]
            with conn:  
                beam_data = check_element(conn, beam_name)
                load_data = get_load_data(conn, key[1], load_type='basic')
            #          
            line['load_number'] = load_data[0]
            line['beam_number'] = beam_data[0]
            line['qx0i'] = 'NULL'
            line['qy0i'] = 'NULL'
            line['qz0i'] = 'NULL'
            line['qx1i'] = 'NULL'
            line['qy1i'] = 'NULL'
            line['qz1i'] = 'NULL'
            #
            line = line[['load_number', 'beam_number',
                         'load_title', 'load_system', 
                         'L0', 'qx0', 'qy0', 'qz0', 'qx0i', 'qy0i', 'qz0i',
                         'L1', 'qx1', 'qy1', 'qz1', 'qx1i', 'qy1i', 'qz1i',
                         'x', 'y', 'z']].values
            #line
            with conn:
                self._push_load(conn, load=line)
        #print('===')
        #1/0
    #
    #
    def fer(self):
        """ Return Fix End Reactions (FER) global system"""
        #beams = self._f2u_beams
        conn = create_connection(self._bd_file)       
        for load_name in set(self._labels):
            with conn: 
                bldf = self.get_wave_load(conn, load_name=load_name)
            #
            # TODO : select coordinate for design load
            #
            grpm = bldf.groupby(['element_name'])
            #
            ipart = ['NULL', 'NULL', 'NULL','NULL', 'NULL', 'NULL']
            res = []
            for key, items in grpm:
                beam =  BeamItemSQL(key[0], self._bd_file)
                end_nodes = beam.connectivity
                #
                for item in items.itertuples():
                    data = [item.qx1, item.qy1, item.qz1,
                            item.qx2, item.qy2, item.qz2,
                            item.L_end1, item.L_end2,
                            item.element_name, item.load_comment,
                            item.load_name, item.load_system,
                            0, 'Line Load']
                    load = LineBeam._make(data)
                    gnload = load.fer(L=beam.L)
                    res.extend([[gnload[1], gnload[2], 1, beam.number,
                                 end_nodes[0], *gnload[4], *ipart],
                                [gnload[1], gnload[2], 1, beam.number,
                                 end_nodes[1], *gnload[5], *ipart]])
                #
                #
                #self._load.df(data=items)
                #res = self._load(beam=beam).fer()
                #
                #1 / 0
                #for gnload in res:
                #    self._node_eq[key[0]] = [[end_nodes[0], *gnload[4], gnload[1], gnload[2]],
                #                             [end_nodes[1], *gnload[5], gnload[1], gnload[2]]]
                #
            #
            with conn: 
                self._push_node_load(conn, res)
            #
            #print('---')
        print('--> get_end_forces')
        #1 / 0    
    #
    #
#