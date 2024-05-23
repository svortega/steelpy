#
# Copyright (c) 2009 steelpy
# 
#
# Python stdlib imports
from __future__ import annotations
#from array import array
#from collections.abc import Mapping
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
from steelpy.ufo.load.process.wave_load import MetoceanLoad # WaveData, 
#from steelpy.ufo.load.sqlite.utils import get_load_data
#from steelpy.ufo.load.process.beam.beam import LineBeam
from steelpy.ufo.mesh.sqlite.beam import BeamSQL # BeamItemSQL, 
from steelpy.ufo.mesh.sqlite.utils import (get_connectivity,
                                           get_element_data) # check_element)
from steelpy.sections.sqlite.utils import get_section
from steelpy.ufo.mesh.sqlite.nodes import pull_node
#
from steelpy.utils.sqlite.utils import create_connection #, create_table
#from steelpy.utils.math.operations import linspace, trnsload
#
#
#import pandas as pd
#from steelpy.utils.dataframe.main import DBframework
#
#
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
                #
                df_bload = item.beam_load(beams)
                #
                # ------------------------------------------
                # push data in database
                #
                with conn:
                    df_bload.to_sql('LoadBeamLine', conn,
                                    if_exists='append', index=False)
                #res.Fwave()
        #print('-->')
        #1 / 0
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
    