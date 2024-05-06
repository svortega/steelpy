#
# Copyright (c) 2009 steelpy
# 
#
# Python stdlib imports
from __future__ import annotations
#from array import array
from collections.abc import Mapping
from concurrent.futures import ProcessPoolExecutor
#from multiprocessing import Pool
#from collections import defaultdict
from dataclasses import dataclass
from typing import NamedTuple
#import re
import os
#from operator import itemgetter
#from itertools import groupby
from math import dist
from operator import sub, add

# package imports
from steelpy.ufo.load.process.wave_load import WaveData, MetoceanLoad
from steelpy.ufo.load.sqlite.utils import get_load_data
from steelpy.ufo.load.process.beam.beam import LineBeam
from steelpy.ufo.mesh.sqlite.beam import BeamItemSQL, BeamSQL
from steelpy.ufo.mesh.sqlite.utils import (get_connectivity,
                                           get_element_data,
                                           check_element)
from steelpy.sections.sqlite.utils import get_section
from steelpy.ufo.mesh.sqlite.nodes import pull_node
#
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.math.operations import linspace, trnsload
#
#
import pandas as pd
from steelpy.utils.dataframe.main import DBframework
#
#
#
#
#
#
#@dataclass
class WaveLoadItemSQLXX(Mapping): #WaveLoadItem
    """ """
    __slots__ = ['_db_file', '_name',
                 '_load', '_plane', '_mg', 
                 '_seastate', '_design_load', '_criterion']

    def __init__(self, load_name: str|int,
                 plane: NamedTuple, db_file:str):
        """ """
        #super().__init__()
        #
        self._name = load_name
        self._plane = plane
        self._db_file = db_file
        #
        #self._node_eq = BeamToNodeSQL(load_name=load_name,
        #                              db_file=self._db_file)
        self._seastate:list = []
        self._design_load:list = []
        self._criterion:list = []
        self._mg:list = []
        #
        conn = create_connection(self._db_file)
        with conn:        
            self._create_table(conn)
            #self._create_node_table(conn)
    #
    @property
    def _labels(self):
        """ """
        table = 'SELECT Load.name FROM Load'
        conn = create_connection(self._db_file)
        with conn:        
            cur = conn.cursor()
            cur.execute(table)
            items = cur.fetchall()
        return [item[0] for item in items]
    #
    #
    # -----------------------------------------------
    #
    def __setitem__(self, load_name: int|str,
                    wave_load: list) -> None:
        """
        criterion : wave_comb, mg, design_load, criterion
        """
        #1 / 0
        try:
            self._labels.index(load_name)
            self._seastate.append(wave_load[0])
            self._mg.append(wave_load[1])
            self._design_load.append(wave_load[2])
            self._criterion.append(wave_load[3])            
        except ValueError:
            raise Exception(f'wave load case {load_name} missing')
    
    def __getitem__(self, load_name: int | str):
        """
        """
        try:
            index = self._labels.index(load_name)
        except ValueError:
            raise KeyError('Invalid load name : {:}'.format(load_name))
        #
        1 / 0
        return WaveData(name=self._labels[index], 
                        seastate=self._seastate[index],
                        design_load=self._design_load[index],
                        criterion=self._criterion[index])     
    #    
    #
    #
    # -----------------------------------------------
    #
    def __contains__(self, value) -> bool:
        return value in self._labels

    def __len__(self) -> int:
        return len(self._labels)

    def __iter__(self):
        """
        """
        items = set(self._labels)
        return iter(items)
    #
    # -----------------------------------------------   
    #    
    #
    # -----------------------------------------------
    #    
    #
    def _create_table(self, conn) -> None:
        """ """
        _table_element_line = "CREATE TABLE IF NOT EXISTS WaveLoad(\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                load_id INTEGER NOT NULL REFERENCES Load(number),\
                                element_id INTEGER NOT NULL REFERENCES Element(number),\
                                title TEXT,\
                                system INTEGER NOT NULL,\
                                L0 DECIMAL,\
                                qx0 DECIMAL,\
                                qy0 DECIMAL,\
                                qz0 DECIMAL,\
                                qx0i DECIMAL,\
                                qy0i DECIMAL,\
                                qz0i DECIMAL,\
                                L1 DECIMAL,\
                                qx1 DECIMAL,\
                                qy1 DECIMAL,\
                                qz1 DECIMAL,\
                                qx1i DECIMAL,\
                                qy1i DECIMAL,\
                                qz1i DECIMAL,\
                                BS DECIMAL,\
                                OTM DECIMAL,\
                                x DECIMAL,\
                                y DECIMAL,\
                                z DECIMAL);"
        #
        create_table(conn, _table_element_line)
    #
    def _push_load(self, conn, load: list):
        """ """
        #
        sql = 'INSERT INTO WaveLoad(load_id, element_id,\
                                        title, system,\
                                        L0, qx0, qy0, qz0, qx0i, qy0i, qz0i,\
                                        L1, qx1, qy1, qz1, qx1i, qy1i, qz1i,\
                                        BS, OTM,\
                                        x, y, z)\
                                        VALUES(?,?,?,?,\
                                               ?,?,?,?,?,?,?,\
                                               ?,?,?,?,?,?,?,\
                                               ?,?,?,?,?)'
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
        load_data = get_load_data(conn, self._name, load_level='basic')
        load_id = load_data[0]
        #
        cur = conn.cursor()
        cur.execute("SELECT Load.name, Element.name, WaveLoad.*\
                    FROM WaveLoad, Element, Load\
                    WHERE WaveLoad.load_id = {:}\
                    AND WaveLoad.element_id = Element.number\
                    AND WaveLoad.load_id = Load.number;"
                    .format(load_id))
        rows = cur.fetchall()
        #
        cols = ['load_name', 'element_name', 'number', 'load_id', 'element_id', 
                'load_comment', 'load_system',
                'L0', 'qx0', 'qy0', 'qz0', 'qx0i', 'qy0i', 'qz0i',
                'L1', 'qx1', 'qy1', 'qz1', 'qx1i', 'qy1i', 'qz1i',
                'BS', 'OTM',
                'x', 'y', 'z']        
        df = db.DataFrame(data=rows, columns=cols)
        df = df[['load_name', 'load_id', 
                 'load_comment', 'load_system',
                 'element_name', 'element_id', 
                 'L0', 'qx0', 'qy0', 'qz0',
                 'L1', 'qx1', 'qy1', 'qz1',
                 'BS', 'OTM',
                 'x', 'y', 'z']]
        return df
    #
    def get_elements(self,conn):
        """ """
        cur = conn.cursor()
        cur.execute ("SELECT * FROM Element;")
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
        1 / 0
        _table_element_point = "CREATE TABLE IF NOT EXISTS LoadBeamToNode(\
                                number INTEGER PRIMARY KEY NOT NULL,\
                                load_id INTEGER NOT NULL REFERENCES Load(number),\
                                title TEXT,\
                                system TEXT,\
                                element_id INTEGER NOT NULL REFERENCES Element(number),\
                                node_id INTEGER NOT NULL REFERENCES Node(number),\
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
        1 / 0
        sql = 'INSERT INTO LoadBeamToNode(load_id, title, system, \
                                             element_id,\
                                             node_id, fx, fy, fz, mx, my, mz,\
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
        try:
            index = self._labels.index(load_name)
        except ValueError:
            raise KeyError('Invalid load name : {:}'.format(load_name))        
        print('# Calculating wave beam forces')
        #index = self._labels.index(load_name)
        design_load = self._design_load[index]
        seastate = self._seastate[index]
        #wave = self.__getitem__(load_name)
        sname = f'{seastate.name}_{seastate.title}'
        sload = seastate.load()
        #
        df_bload = self.get_beam_load(sname, sload)
        #
        #grpbeam = df_bload.groupby(['element_type', 'element_name'])
        #
        self.df(data=df_bload)
        #
        print('# End process')
    #
    #
    def get_beam_load(self, wname, wload):
        """ """
        #
        conn = create_connection(self._db_file)
        with conn:
            labels = self.get_elements(conn)
        #
        df_bload = None
        for idx, key in enumerate(labels):
            beam = BeamItemSQL(key[1],
                               plane=self._plane,
                               db_file=self._db_file)
            Fwave = wload.Fwave(beam=beam)
            df_load = Fwave.df
            df_load['load_name'] = self._name
            df_load['load_title'] = wname
            #df_load['load_title'] = df_load.apply(lambda row: f"{wname}_{round(row.x, 2)}_{round(row.y, 2)}_{round(row.z, 2)}", axis=1)
            df_load['load_system'] = 'local'
            try:
                1/idx
                df_bload = pd.concat([df_bload, df_load], ignore_index=True)
            except ZeroDivisionError:
                df_bload = df_load
            #
            # process to select wave point based on user request
            #
            #Fx, Fy, OTM = udl.span_loading()
            #indmax = Fx.argmax(dim='length').values
            #vmax = Fx.idxmax(dim='length').values
            #print('')
            #print('Total combined force [kN-m]')
            #print(f'Fx ={np.max(Fx) / 1000: 1.3e}, Fy={np.max(Fy) / 1000: 1.3e}, OTM={np.max(OTM)/1000: 1.3e}')
            #print('---')
            #
            #Fx.sel(x=0).plot.line(hue='z')
            #plt.show()            
            #
        #1 / 0
        return df_bload
    #
    # -----------------------------------------------
    #
    def df(self, data):
        """ """
        #
        conn = create_connection(self._db_file)
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
            line['load_id'] = load_data[0]
            line['beam_number'] = beam_data[0]
            line['qx0i'] = 'NULL'
            line['qy0i'] = 'NULL'
            line['qz0i'] = 'NULL'
            line['qx1i'] = 'NULL'
            line['qy1i'] = 'NULL'
            line['qz1i'] = 'NULL'
            #
            line = line[['load_id', 'beam_number',
                         'load_title', 'load_system', 
                         'L0', 'qx0', 'qy0', 'qz0', 'qx0i', 'qy0i', 'qz0i',
                         'L1', 'qx1', 'qy1', 'qz1', 'qx1i', 'qy1i', 'qz1i',
                         'BS', 'OTM', 
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
        conn = create_connection(self._db_file)       
        for load_name in set(self._labels):
            with conn: 
                bldf = self.get_wave_load(conn, load_name=load_name)
            #
            # ------------------------------------------
            # TODO : select coordinate for design load
            #
            waveinput = self.__getitem__(load_name=load_name)
            design_load = waveinput.design_load.split("_")
            criteria = waveinput.criterion
            #
            value_type = design_load[0]
            load_type = design_load[1]
            #
            #
            #grpm = bldf.groupby(['element_name','y'])
            grpwave = bldf.groupby(['element_name','y'])[['BS', 'OTM']].sum()
            #
            if criteria == 'global':
                # FIXME
                1 / 0
            else: # default 
                grpm = grpwave[load_type].abs().groupby('element_name')
            #
            if value_type.lower() == 'max':
                idx = grpm.idxmax()
            else:
                1 / 0
            #
            bldf.set_index(['element_name','y'], inplace=True)
            grpm = bldf.loc[idx]
            #
            #
            # ----------------------------------------
            #
            1 / 0
            ipart = ['NULL', 'NULL', 'NULL','NULL', 'NULL', 'NULL']
            res = []
            for item in grpm.itertuples():
                element_name = item.Index[0]
                beam =  BeamItemSQL(element_name,
                                    plane=self._plane, 
                                    db_file=self._db_file)
                node1, node2 = beam.nodes
                #
                #for item in items.itertuples():
                #
                data = [item.qx0, item.qy0, item.qz0,    # q_inplane [qx, qy, qz, qt]
                        item.qx1, item.qy1, item.qz1,    # q_outplane [qx, qy, qz, qt]
                        item.L0, item.L1,                # L0, L1,
                        element_name, item.load_comment, # name, title, load_name, component_name
                        item.load_name, item.load_system, # system, load_complex,
                        1, 'Line Load']
                load = LineBeam._make(data)
                gnload = load.fer_beam(L=beam.L)
                # load local system to global 
                gnload = [*gnload[4], *gnload[5]]
                lnload = trnsload(gnload, beam.T3D())
                #
                res.extend([[item.load_id, item.load_comment, 
                             'global', beam.number,
                             node1.number, *lnload[:6], *ipart],
                            [item.load_id, item.load_comment, 
                             'global', beam.number,
                             node2.number, *lnload[6:], *ipart]])
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
        #print('--> get_end_forces')
        #1 / 0    
    #
    #
    # -----------------------------------------------
    #
    def beam_load(self, steps:int = 10):
        """ """
        1 / 0
        conn = create_connection(self._db_file)
        #beamfun = defaultdict(list)
        beamfun = []
        for load_name in set(self._labels):
            with conn: 
                bldf = self.get_wave_load(conn, load_name=load_name)
            #
            # TODO : select coordinate for design load
            #
            grpm = bldf.groupby(['element_name'])
            #
            for key, items in grpm:
                beam = BeamItemSQL(key[0],
                                   plane=self._plane, 
                                   db_file=self._db_file)
                mat = beam.material
                sec = beam.section.properties()
                Lsteps = linspace(start=0, stop=beam.L, num=steps+1, endpoint=True)
                #
                for item in items.itertuples():
                    data = [item.qx0, item.qy0, item.qz0,  # q_inplane [qx, qy, qz, qt]
                            item.qx1, item.qy1, item.qz1,  # q_outplane [qx, qy, qz, qt]
                            item.L0, item.L1,              # L0, L1,
                            item.element_name, item.load_comment, # name, title, load_name, component_name
                            item.load_name, item.load_system,     # system, load_complex,
                            1, 'Line Load']
                    bitem = LineBeam._make(data)
                    lout = bitem.Fx(x=Lsteps, L=beam.L,
                                    E=mat.E, G=mat.G, 
                                    Iy=sec.Iy, Iz=sec.Iy,
                                    J=sec.J, Cw=sec.Cw, Area=sec.area)
                    #beamfun[key[0]].extend(lout)
                    beamfun.extend(lout)
                #
        #print('---')
        return beamfun
#
#
class MetoceanLoadSQL(MetoceanLoad):
    __slots__ = ['db_file', '_plane', '_component']
    #
    def __init__(self, db_file:str,
                 plane: NamedTuple, component: int) -> None:
        """
        """
        super().__init__()
        #
        self.db_file = db_file
        self._plane = plane
        self._component = component
    #
    #
    #
    # -----------------------------------------------
    #
    def __setitem__(self, name: int|str, condition) -> None:
        """
        """
        try:
            index = self._labels.index(name)
            raise IOError(f"Load name {name} exist")
        except ValueError:
            self._labels.append(name)
            self._condition.append(condition)
        #
        #cases= ((name, condition.name, "basic", 'metocean'), )
        cases= ((name, self._component, "basic", condition.name), )
        conn = create_connection(self.db_file)
        with conn:
            self._push_basicload(conn, cases)
        #
    #
    def __getitem__(self, name:int|str) :
        """
        """
        try:
            index = self._labels.index(name)
            condition = self._condition[index]
        except ValueError:
            raise IOError(f"Load name {name} not found")
        #
        # get load data
        #conn = create_connection(self._db_file)
        #with conn:  
        #    node = get_node(conn, node_name)
        1 / 0
    #
    #
    #@property
    #def condition(self):
    #    """ """
    #    return self._condition
    #
    
    #@condition.setter
    #def condition(self, value):
    #    """ """
    #    cases = []
    #    for key, item in value.items():
    #        self._condition[key] = item
    #        cases.append((key, f'MET_{key}', "basic", 'metocean'))
    #    #
    #    conn = create_connection(self.db_file)
    #    with conn:        
    #        self._push_basicload(conn, tuple(cases))
    #    #print('-->')
    #
    # -----------------------------------------------
    #
    def _push_basicload(self, conn, cases: tuple):
        """ """
        #project = (load_name, load_title, "basic", 'metocean')
        table = 'INSERT INTO Load(name, component_id, level, title) \
                 VALUES(?,?,?,?)'
        cur = conn.cursor()
        cur.executemany(table, cases)
    #
    #
    #
    # -----------------------------------------------
    # Process
    #
    def process(self):
        """ """
        if self._condition:
            print('--- Metocean Process')
            #
            # ------------------------------------------
            #
            beams = BeamSQL(db_file=self.db_file,
                            component=self._component,
                            plane=self._plane)
            #
            # ------------------------------------------
            conn = create_connection(self.db_file)
            for item in self._condition:
                #sname =  f'MET_{key}'
                #sload = item.load()
                #self.get_beam_load(wnumber=item.number,
                #                   title=sname,
                #                   wload=sload)
                #
                df_bload = item.get_beam_load(beams)
                #
                # ------------------------------------------
                # push data in database
                #
                with conn:
                    df_bload.to_sql('LoadBeamLine', conn,
                                    #index_label=header,
                                    if_exists='append', index=False)
                #res.Fwave()
        #print('-->')
        #1 / 0
    #
    #
    def get_beam_loadXX(self, wnumber, title, wload):
        """ """
        #
        conn = create_connection(self.db_file)
        with conn:
            labels = self.get_elements(conn)
            labels = [item[1] for item in labels]
        #
        #
        dftemp = []
        #
        # ------------------------------------------
        # TODO: Multiprocess --> sqlite issue
        #
        #cpuno = os.cpu_count() - 1
        #
        #def myfunc(beam):
        #    beam = BeamItemSQL(name,
        #                       plane=self._plane,
        #                       db_file=self.db_file)
        #    Fwave = wload.Fwave(beam=beam)
        #    return Fwave.solve()
        #
        #beams = [BeamItemWave(name,
        #                      db_file=self.db_file)
        #         for name in labels]
        #
        #with ProcessPoolExecutor(max_workers=cpuno) as executor:
        #    for r in executor.map(myfunc, beams):
        #        dftemp.append(r)
        #
        #engine = Engine(wload=wload)
        #                plane=self._plane,
        #                wload=wload)
        #
        #dftemp = engine.run(beams)
        #
        #with ProcessPoolExecutor(max_workers=cpuno) as executor:
        #    for r in executor.map(engine, beams):
        #        dftemp.append(r)
        #
        #try:
        #    pool = Pool(cpuno)
        #    engine = Engine(db_file=self.db_file,
        #                    plane=self._plane,
        #                    wload=wload)
        #    dftemp = pool.map(engine, beams)
        #finally:
        #    pool.close()
        #    pool.join()
        #
        # ------------------------------------------
        #
        #for beam in beams:
        for name in labels:
            # get beam 
            beam = BeamItemSQL(name,
                               plane=self._plane,
                               db_file=self.db_file)
            # solve beam forces
            Fwave = wload.Fwave(beam=beam)
            dftemp.extend(Fwave.solve())
            #
            #print('-->')
        #
        # ------------------------------------------
        # create database
        #
        header = ['element_type', 'element_name', 'element_id', 
                  'type',
                  'qx0', 'qy0', 'qz0', 'qx1', 'qy1', 'qz1',
                  'L0', 'L1', 'BS', 'OTM', 
                  'x', 'y', 'z']        
        #
        df = DBframework()
        df_bload = df.DataFrame(data=dftemp, columns=header, index=None)        
        #
        header = ['load_id', 'element_id',
                  'title', 'system', 'type',
                  'L0', 'qx0', 'qy0', 'qz0',
                  'L1', 'qx1', 'qy1', 'qz1',
                  'BS', 'OTM', 'x', 'y', 'z']        
        #
        df_bload['load_id'] = wnumber
        df_bload['title'] = title        
        df_bload['system'] = 'local'
        #df_bload.rename(columns={'load_type': 'type'}, inplace=True)
        df_bload =  df_bload[header]
        #
        #
        # ------------------------------------------
        # select data in database        
        #
        #
        #
        # ------------------------------------------
        # push data in database        
        #
        with conn:
            df_bload.to_sql('LoadBeamLine', conn,
                            index_label=header, 
                            if_exists='append', index=False)
        1 / 0
        return df_bload
    #
    # -----------------------------------------------
    #
    #def get_elements(self,conn):
    #    """ """
    #    cur = conn.cursor()
    #    cur.execute ("SELECT * FROM Element;")
    #    row = cur.fetchall()
    #    return row
    #
    #    
#
#
#
class HydroLoadSQLXX(Mapping):
    
    __slots__ = ['db_file', '_plane', '_condition', '_labels']
    
    def __init__(self, plane: NamedTuple, db_file:str) -> None:
        """
        """
        self._plane = plane
        self.db_file = db_file
        self._condition:dict = {}
        self._labels:list = []
    #
    # -----------------------------------------------
    #
    def __setitem__(self, load_name: int|str,
                    condition) -> None:
        """
        criterion : wave_comb, mg, design_load, criterion
        """
        try:
            self._condition[load_name]
            #index = self._labels.index(load_name)
            raise IOError('Invalid load name : {:}'.format(load_name))
        except KeyError:
            self._labels.append(load_name)
            self._condition[load_name] = condition
    
    def __getitem__(self, load_name: int | str):
        """
        """
        #try:
        #    index = self._labels.index(load_name)
        #except ValueError:
        #    raise KeyError('Invalid load name : {:}'.format(load_name))
        #
        #1 / 0
        #print('-->')
        return self._condition[load_name]
    #
    # -----------------------------------------------
    #
    def __contains__(self, value) -> bool:
        return value in self._labels

    def __len__(self) -> int:
        return len(self._labels)

    def __iter__(self):
        """
        """
        return iter(self._labels)
    #
    # -----------------------------------------------
    # Process
    #
    def process(self):
        """ """
        print('-->')
        for key, item in self.items():
            sname =  f'{key}_{item.name}'
            sload = item.load()
            self.get_beam_load(wnumber=item.number, title=sname,
                               wload=sload)
            #res.Fwave()
        print('-->')
        #1 / 0
    #
    #
    def get_beam_load(self, wnumber, title, wload):
        """ """
        #
        conn = create_connection(self.db_file)
        with conn:
            labels = self.get_elements(conn)
            labels = [item[1] for item in labels]
        #
        #
        dftemp = []
        #
        # ------------------------------------------
        # TODO: Multiprocess --> sqlite issue
        #
        #cpuno = os.cpu_count() - 1
        #
        #def myfunc(beam):
        #    beam = BeamItemSQL(name,
        #                       plane=self._plane,
        #                       db_file=self.db_file)
        #    Fwave = wload.Fwave(beam=beam)
        #    return Fwave.solve()
        #
        #beams = [BeamItemWave(name,
        #                      db_file=self.db_file)
        #         for name in labels]
        #
        #with ProcessPoolExecutor(max_workers=cpuno) as executor:
        #    for r in executor.map(myfunc, beams):
        #        dftemp.append(r)
        #
        #engine = Engine(wload=wload)
        #                plane=self._plane,
        #                wload=wload)
        #
        #dftemp = engine.run(beams)
        #
        #with ProcessPoolExecutor(max_workers=cpuno) as executor:
        #    for r in executor.map(engine, beams):
        #        dftemp.append(r)
        #
        #try:
        #    pool = Pool(cpuno)
        #    engine = Engine(db_file=self.db_file,
        #                    plane=self._plane,
        #                    wload=wload)
        #    dftemp = pool.map(engine, beams)
        #finally:
        #    pool.close()
        #    pool.join()
        #
        # ------------------------------------------
        #
        #for beam in beams:
        for name in labels:
            # get beam 
            beam = BeamItemSQL(name,
                               plane=self._plane,
                               db_file=self.db_file)
            # solve beam forces
            Fwave = wload.Fwave(beam=beam)
            dftemp.extend(Fwave.solve())
            #
            #df_load = Fwave.df
            #df_load['load_id'] = wnumber
            #df_load['title'] = title
            #df_load['load_title'] = df_load.apply(lambda row: f"{wname}_{round(row.x, 2)}_{round(row.y, 2)}_{round(row.z, 2)}", axis=1)
            #df_load['system'] = 'local'
            #df_load.rename(columns={'load_type': 'type'}, inplace=True)
            #df_load.drop(columns=['element_name'], inplace=True, axis=1)
            #df_load =  df_load[header]
            #try:
            #    1/idx
            #    df_bload = pd.concat([df_bload, df_load], ignore_index=True)
            #except ZeroDivisionError:
            #    df_bload = df_load
            #
            # process to select wave point based on user request
            #
            #Fx, Fy, OTM = udl.span_loading()
            #indmax = Fx.argmax(dim='length').values
            #vmax = Fx.idxmax(dim='length').values
            #print('')
            #print('Total combined force [kN-m]')
            #print(f'Fx ={np.max(Fx) / 1000: 1.3e}, Fy={np.max(Fy) / 1000: 1.3e}, OTM={np.max(OTM)/1000: 1.3e}')
            #print('---')
            #
            #Fx.sel(x=0).plot.line(hue='z')
            #plt.show()
            #
            #print('-->')
        #
        # ------------------------------------------
        # create database
        #
        header = ['element_type', 'element_name', 'element_id', 
                  'type',
                  'qx0', 'qy0', 'qz0', 'qx1', 'qy1', 'qz1',
                  'L0', 'L1', 'BS', 'OTM', 
                  'x', 'y', 'z']        
        #
        df = DBframework()
        df_bload = df.DataFrame(data=dftemp, columns=header, index=None)        
        #
        header = ['load_id', 'element_id',
                  'title', 'system', 'type',
                  'L0', 'qx0', 'qy0', 'qz0',
                  'L1', 'qx1', 'qy1', 'qz1',
                  'BS', 'OTM', 'x', 'y', 'z']        
        #
        df_bload['load_id'] = wnumber
        df_bload['title'] = title        
        df_bload['system'] = 'local'
        #df_bload.rename(columns={'load_type': 'type'}, inplace=True)
        df_bload =  df_bload[header]
        #
        #
        # ------------------------------------------
        # select data in database        
        #
        #
        #
        # ------------------------------------------
        # push data in database        
        #
        with conn:
            df_bload.to_sql('LoadBeamLine', conn,
                            index_label=header, 
                            if_exists='append', index=False)
        1 / 0
        return df_bload
    #
    # -----------------------------------------------
    #
    def get_elements(self,conn):
        """ """
        cur = conn.cursor()
        cur.execute ("SELECT * FROM Element;")
        row = cur.fetchall()
        return row
    #
    #
#
#
@dataclass
class Engine:
    def __init__(self, wload):
        #self.db_file = db_file
        self.wload = wload
        #self.plane = plane
    
    #def __call__(self, beam):
    #    #beam = BeamItemSQL(name,
    #    #                   plane=self.plane,
    #    #                   db_file=self.db_file)
    #    Fwave = self.wload.Fwave(beam=beam)
    #    return Fwave.solve()
    #
    def myfuct(self, beam):
        """ """
        Fwave = self.wload.Fwave(beam=beam)
        return Fwave.solve()        
    #
    def run(self, beams):
        """ """
        cpuno = os.cpu_count() - 1
        #
        dftemp = []
        with ProcessPoolExecutor(max_workers=cpuno) as executor:
            for r in executor.map(self.myfuct, beams):
                dftemp.append(r)
        #
        return dftemp
#
#
#
@dataclass
class BeamItemWave:
    """ """
    __slots__ = ['name', 'number', 'type', 'db_file', '_component', 
                 'nodes', 'unit_vector', 'section', 'connectivity']
    
    def __init__(self, name:int|str,
                 component: int, 
                 db_file:str) -> None:
        """
        """
        self.name = name
        self.type: str = "beam"
        self._component = component
        self.db_file = db_file
        #
        self.connectivity = self._connectivity()
        self.number = self._number()
        self.nodes = self._nodes()
        self.unit_vector =  self._unit_vector()
        self.section = self._section()
    #
    #
    def _connectivity(self) -> list:
        """
        """
        conn = create_connection(self.db_file)
        with conn:
            connodes = get_connectivity(conn, self.name)
        return connodes
    #
    def _number(self) -> int:
        """ """
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn,
                                    element_name=self.name,
                                    element_type=self.type,
                                    component=self._component)
        return data[1]    
    #
    def _nodes(self) -> list:
        """
        """
        _nodes = []
        conn = create_connection(self.db_file)
        for _node in self.connectivity:
            with conn:
                _nodes.append(pull_node(conn, node_name=_node))
        return _nodes
    #
    def _unit_vector(self) -> list[ float ]:
        """
        """
        nodes = self.connectivity
        conn = create_connection(self.db_file)
        with conn:
            node1 = pull_node(conn, node_name=nodes[0])
            node2 = pull_node(conn, node_name=nodes[1])
        # direction cosines
        L =  dist(node1[:3], node2[:3])
        uv = list(map(sub, node2[:3], node1[:3]))
        return [item / L for item in uv]       
    #
    def _section(self) -> list:
        """
        """
        #BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        conn = create_connection(self.db_file)
        with conn:
            data = get_element_data(conn,
                                    element_name=self.name,
                                    element_type=self.type,
                                    component=self._component)
            
            sect =  get_section(conn, data[5])
        #
        return sect
    