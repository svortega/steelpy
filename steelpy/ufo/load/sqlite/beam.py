#
# Copyright (c) 2009 steelpy
# 
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Mapping
from typing import NamedTuple
import re
from math import isclose
#
# package imports
#
from steelpy.ufo.mesh.sqlite.beam import BeamItemSQL, BeamSQL
from steelpy.ufo.mesh.sqlite.utils import check_element #, check_nodes

from steelpy.ufo.load.process.beam.beam import LineBeam, PointBeam
from steelpy.ufo.load.process.beam.main import (BeamTypeBasic,
                                                BeamLineBasic,
                                                BeamPointBasic)
from steelpy.ufo.load.sqlite.utils import get_load_data

from steelpy.utils.math.operations import linstep
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework
#
#
#
# ---------------------------------
#
class BeamSQLMaster(Mapping):
    
    def __init__(self) -> None:
        """
        """
        self._system_flag: int = 0   
    #
    # -----------------------------------------------
    #
    def __len__(self) -> int:
        return len(self._labels)
    #
    def __contains__(self, value) -> bool:
        return value in self._labels
    #
    def __iter__(self):
        """
        """
        return iter(self._labels) 
    #  
#
#
def get_load_basics(conn, component: int):
    """ """
    query = (component, )
    #with conn:
    cur = conn.cursor()
    #
    table = "SELECT Node.name, Node.number FROM Node \
            WHERE Node.component_id = ? ;"
    cur.execute(table, query)
    nodes = cur.fetchall()
    nodes = {item[0]:item[1] for item in nodes}
    #
    table = "SELECT Element.name, Element.number FROM Element \
             WHERE Element.component_id = ? ;"
    cur.execute(table, query)
    elements = cur.fetchall()
    elements = {item[0]:item[1] for item in elements}            
    #
    query = ('basic', component, )
    table = "SELECT Load.name, Load.number FROM Load \
             WHERE Load.level = ? \
             AND Load.component_id = ? ;"
    cur.execute(table, query)
    basic = cur.fetchall()
    basic = {item[0]:item[1] for item in basic}
    #
    return nodes, elements, basic
#
#
# ---------------------------------
#
class BeamLoadItemSQL(BeamSQLMaster):
    __slots__ = ['_name',  '_load', '_plane', '_component', 
                 '_plane', '_node_eq', 'db_file', '_system_flag', 
                 '_beam']

    def __init__(self, load_name: str|int, plane: NamedTuple,
                 component: int, 
                 db_file: str) -> None:
        """
        """
        #
        self._name = load_name
        self._plane = plane
        self._component = component
        #self._system_flag: int = 0
        super().__init__()
        #
        self.db_file = db_file
        # create table
        conn = create_connection(self.db_file)
        with conn:
            self._create_table(conn)
        #
        #
        self._load = BeamLoadTypeSQL(load_name=self._name,
                                     component=self._component,
                                     plane=self._plane, 
                                     db_file=self.db_file)
    #
    @property
    def _labels(self):
        """ """
        query = (self._name, self._component)
        #
        #if isinstance(self._name, str):
        #    load_name = f"AND Load.name = '{self._name}' "
        #else:
        #    load_name = f"AND Load.name = {self._name}"        
        #
        # point
        table = "SELECT Element.name \
                FROM Element, LoadBeamPoint, Load \
                WHERE LoadBeamPoint.element_id = Element.number \
                AND Load.number = LoadBeamPoint.load_id \
                AND Load.name = ? AND Load.component_id = ? ;"
        #table += load_name
        #
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            rows = cur.fetchall()
        labels = [item[0] for item in rows]
        #
        # line
        table = "SELECT Element.name \
                FROM Element, LoadBeamLine, Load \
                WHERE LoadBeamLine.element_id = Element.number \
                AND Load.number = LoadBeamLine.load_id \
                AND Load.name = ? AND Load.component_id = ? ;"
        #table += load_name
        #
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            rows = cur.fetchall()
        labels.extend([item[0] for item in rows])
        #
        return labels
    #    
    # -----------------------------------------------
    #
    #
    def __setitem__(self, beam_name: int|str,
                    beam_load: list) -> None:
        """
        """
        #
        if isinstance(beam_load, dict):
            load_type = beam_load['type']
        else:
            load_type = beam_load[0]
        #
        if re.match(r"\b(point|mass)\b", load_type, re.IGNORECASE):
            self._load._point[beam_name] = beam_load
            
        elif re.match(r"\b(line|u(niform(ly)?)?d(istributed)?l(oad)?)\b",
                      load_type, re.IGNORECASE):
            self._load._line[beam_name] = beam_load
        
        else:
            raise IOError(f'Beam lod type {beam_load[0]} not implemented')    
    #
    def __getitem__(self, beam_name: int | str):
        """
        """
        conn = create_connection(self.db_file)
        with conn:  
            #beam = check_element(conn, beam_name)
            beam =  BeamItemSQL(beam_name=beam_name,
                                plane=self._plane,
                                component=self._component, 
                                db_file=self.db_file)
        #try:
        #memb_type = beam.type # beam[3]
        if beam.type != 'beam':
            raise ValueError(f"element {beam_name} type {beam.type} not valid")
        #except TypeError:
        #    raise IOError(f"beam {beam_name} not found")
        #
        #if not beam_name in self._labels:
        #if not beam_name in  self._labels:
        #self._labels.append(beam_name)
        #return self._load(beam=beam)
        #return BeamLoadTypeSQL(load_name=self._name,
        #                       beam=beam,
        #                       component=self._component, 
        #                       db_file=self.db_file)
        #
        return self._load(beam=beam)
    #
    # -----------------------------------------------
    #
    #
    def __str__(self, units: str = "si") -> str:
        """ """
        conn = create_connection(self.db_file)
        with conn:         
            line, point = get_beam_load(conn,
                                        load_name=self._name,
                                        beam_name="*",
                                        component=self._component)
        #
        output = ""
        for item in line:
            output += item.__str__()
        #
        for item in point:
            output += item.__str__()
        #
        return output
    #
    #
    # -----------------------------------------------
    #
    def fer2(self, beams, load_name: str|int):
        """
        Beam reacition global system according to boundaries
        """
        conn = create_connection(self.db_file)
        with conn:         
            line, point = get_beam_load(conn,
                                        load_name=load_name,
                                        beam_name='*',
                                        component=self._component)
        #
        b2n = []
        global_system = 0
        # line loadreactions
        for bload in line:
            beam = beams[bload.name]
            node1, node2 = beam.nodes
            #print(f'line load {key}')
            lnload = bload.fer_beam(beam=beam,
                                    system=global_system)
                                    #Pdelta=Pdelta)
            # local to global system
            #gnload = [*res[7], *res[8]]
            #lnload = trnsload(gnload, beam.T3D())
            #
            load7 = zerofilter(lnload[7])
            load8 = zerofilter(lnload[8])
            #
            b2n.append([bload.load_name,
                        bload.title,
                        global_system, 
                        beam.number,
                        node1.number,
                        load7,
                        node2.number,
                        load8])
        # point load
        for bload in point:
            beam = beams[bload.name]
            node1, node2 = beam.nodes
            #print(f'point load {key}')
            lnload = bload.fer_beam(beam=beam,
                                    system=global_system)
                                    #Pdelta=Pdelta)
            #gnload = [*res[7], *res[8]]
            #lnload = trnsload(gnload, beam.T3D())
            #
            load7 = zerofilter(lnload[7])
            load8 = zerofilter(lnload[8])
            #
            b2n.append([bload.load_name,
                        bload.title,
                        global_system, 
                        beam.number,
                        node1.number,
                        load7,
                        node2.number,
                        load8])
        #
        return b2n    
    #
    def fer(self, beams) -> None:
        """ Push Fix End Reactions (FER) global system in sqlite """
        conn = create_connection(self.db_file)
        with conn:
            load_data = get_load_data(conn, self._name,
                                      load_level='basic',
                                      component=self._component)
            load_id = load_data[0]
        #
        #ipart = [None, None, None, None, None, None]
        res =[]
        #res1 = self._load.fer(beams)
        res1 = self.fer2(beams, self._name)
        #
        for gnload in res1:
            try:
                1 / gnload[2]
                raise RuntimeError('node load in local system')
            except ZeroDivisionError:
                load_system = 'global'
                #
            #res.extend([[load_id, gnload[1], load_system, gnload[3], gnload[4], *gnload[5], *ipart],
            #            [load_id, gnload[1], load_system, gnload[3], gnload[6], *gnload[7], *ipart]])
            #
            # [fx, fy, fz, mx, my, mz, rx, ry,rz, theta', theta'',theta''']
            res.extend([[load_id, gnload[3], gnload[4],
                         gnload[1], load_system, 'load', *gnload[5]],
                        [load_id, gnload[3], gnload[6],
                         gnload[1], load_system, 'load', *gnload[7]]])            
        #
        #print('--> get_end_forces')
        if res:
            with conn:  
                push_FER(conn, node_load=res)
    #
    # -----------------------------------------------
    #
    def _create_table(self, conn) -> None:
        """ """
        # -------------------------------------
        # line 
        table = "CREATE TABLE IF NOT EXISTS LoadBeamLine(\
                number INTEGER PRIMARY KEY NOT NULL,\
                load_id INTEGER NOT NULL REFERENCES Load(number),\
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
                BS DECIMAL,\
                OTM DECIMAL,\
                x DECIMAL,\
                y DECIMAL,\
                z DECIMAL);"
        create_table(conn, table)
        # -------------------------------------
        # point
        table = "CREATE TABLE IF NOT EXISTS LoadBeamPoint(\
                number INTEGER PRIMARY KEY NOT NULL,\
                load_id INTEGER NOT NULL REFERENCES Load(number),\
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
                rz DECIMAL);"
        create_table(conn, table)
        # -------------------------------------
        # FER
        table = "CREATE TABLE IF NOT EXISTS LoadBeamFER(\
                number INTEGER PRIMARY KEY NOT NULL,\
                load_id INTEGER NOT NULL REFERENCES Load(number),\
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
                Tw DECIMAL);"
        create_table(conn, table)
    #
    #
    # -----------------------------------------------
    #
    @property
    def df(self):
        """ beam load df"""
        print(' beam load')
        conn = create_connection(self.db_file)
        with conn:
            bload = get_beam_load(conn,
                                  beam_name='*',
                                  load_name=self._name,
                                  component=self._component)
        bload
        1 / 0
    
    @df.setter
    def df(self, df):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            nodes, elements, basic = get_load_basics(conn, self._component)
            #print('--')
        #
        df['load_id'] = df['name'].apply(lambda x: basic[x])
        df['element_id'] =  df['beam'].apply(lambda x: elements[x])
        #df['node_id'] = df['node'].apply(lambda x: nodes[x])
        df['system'] = 'local'
        df['type'] = df['type'].apply(lambda x: x.lower())
        #
        bgroup = df.groupby('type')
        try:
            bpoint = bgroup.get_group('point')
            #pheader = ['name', 'beam', 'title', 'type', 
            #           'L0', 'fx', 'fy', 'fz', 'mx', 'my', 'mz']
            #bpoint = bpoint[pheader]
            #bpoint['type'] = 'force'
            bpoint[['x', 'y', 'z', 'rx', 'ry', 'rz']] = None
            #1 / 0
            #bload._beam._point.df = bpoint
            #
            header = ['load_id', 'element_id',
                      'title', 'system', 'type', 'L0', 
                      'fx', 'fy', 'fz', 'mx', 'my', 'mz',
                      'x', 'y', 'z', 'rx', 'ry', 'rz']        
            #
            bconn = bpoint[header].copy()
            bconn.replace(to_replace=[''], value=[float(0)], inplace=True)
            bconn['title'] = bpoint['title']
            #
            with conn:
                bconn.to_sql('LoadBeamPoint', conn,
                             index_label=header, 
                             if_exists='append', index=False)            
        except KeyError:
            pass
        #
        try:
            bline = bgroup.get_group('line')
            #header = ['name', 'beam', 'title', 'type', 
            #          'L0', 'qx0', 'qy0', 'qz0', 'L1','qx1', 'qy1', 'qz1']
            #bline = bline[lheader]
            #bline['type'] = 'load'
            bline[['BS', 'OTM', 'x', 'y', 'z']] = None
            #1 / 0
            #bload._beam._line.df = bline
            #        
            #
            header = ['load_id', 'element_id',
                      'title', 'system', 'type',
                      'L0', 'qx0', 'qy0', 'qz0', 'qt0',
                      'L1', 'qx1', 'qy1', 'qz1', 'qt1',
                      'BS', 'OTM', 'x', 'y', 'z']
            #
            bconn = bline[header].copy()
            bconn[['BS', 'OTM', 'x', 'y', 'z']] = None
            bconn.replace(to_replace=[''], value=[float(0)], inplace=True)
            bconn['title'] = bline['title']
            #
            with conn:
                bconn.to_sql('LoadBeamLine', conn,
                             index_label=header,
                             if_exists='append', index=False)
        except KeyError:
            pass
        #
        try:
            bmass = bgroup.get_group('mass')
            1 / 0
        except KeyError:
            pass        
    #
    #
    def dfX(self):
        """ """
        #
        conn = create_connection(self.db_file)
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
            #
            #if not beam_name in self._labels:
            self._labels.append(beam_name)           
            #
            line['load_id'] = load_data[0]
            line['beam_number'] = beam_data[0]
            line[['BS', 'OTM', 'x', 'y', 'z']] = None
            #line['qz1i'] = 'NULL'
            #
            line = line[['load_id', 'beam_number',
                         'load_title', 'load_system', 
                         'L0', 'qx0', 'qy0', 'qz0',
                         'L1', 'qx1', 'qy1', 'qz1',
                         'BS', 'OTM', 'x', 'y', 'z']].values
            #line
            with conn:
                self._load._line._push_load(conn, beam_name=beam_name, load=line)
        #print('===')
        #1/0
    #
    # -----------------------------------------------
    #
    def load_function(self, Fb, steps:int,
                      factor:float):
        """ """
        #print('-->')
        beams = BeamSQL(db_file=self.db_file,
                        component=self._component,
                        plane=self._plane)
        #
        #print(' beam load')
        conn = create_connection(self.db_file)
        with conn:
            line, point = get_beam_load(conn,
                                        beam_name='*',
                                        load_name=self._name,
                                        component=self._component)
        #
        #Fbgrp = Fb.groupby('element_name')
        #for bname, bitem in Fbgrp:
        #    bload = self.__getitem__(beam_name=bname)
        #    for key, items in bload._line.items():
        #        print(key)
        #
        loadfun = []
        # line load
        for item in line:
            beam = beams[item.name]
            nodes = beam.connectivity
            mat = beam.material
            sec = beam.section.properties(poisson=mat.poisson)
            #
            Lsteps = linstep(d=beam.section.geometry.d,
                             L=beam.L, steps=steps)
            #
            Fbeam = Fb.loc[item.name]
            Fbeam.set_index('node_name', inplace=True)
            Pdelta = Fbeam.loc[nodes[1]]
            #
            #for bitem in items:
            lout = item.Fx(x=Lsteps, L=beam.L,
                            E=mat.E, G=mat.G,
                            Iy=sec.Iy, Iz=sec.Iz,
                            J=sec.J, Cw=sec.Cw,
                            Area=sec.area,
                            Asy=sec.Asy, Asz=sec.Asz,
                            P=Pdelta.Fx, factor=factor)
            #
            #lout2 = [[Pdelta.load_name, *item]
            #         for item in lout]
            #
            loadfun.extend([[Pdelta.load_name, *item]
                            for item in lout])
        #
        # point load
        for item in point:
            beam = beams[item.name]
            nodes = beam.connectivity
            mat = beam.material
            sec = beam.section.properties(poisson=mat.poisson)
            #
            Lsteps = linstep(d=beam.section.geometry.d,
                             L=beam.L, steps=steps)
            #
            Fbeam = Fb.loc[item.name]
            Fbeam.set_index('node_name', inplace=True)
            Pdelta = Fbeam.loc[nodes[1]]
            #
            #for bitem in items:
            lout = item.Fx(x=Lsteps, L=beam.L,
                            E=mat.E, G=mat.G,
                            Iy=sec.Iy, Iz=sec.Iz,
                            J=sec.J, Cw=sec.Cw,
                            Area=sec.area,
                            Asy=sec.Asy, Asz=sec.Asz,
                            P=Pdelta.Fx, factor=factor)
            #
            #loadfun.extend(lout)
            loadfun.extend([[Pdelta.load_name, *item]
                            for item in lout])
        #
        #1 / 0
        return loadfun
#
#
def zerofilter(items: list):
    """remove near zero items"""
    return [0 if isclose(abs(item), 0.0, abs_tol=1e-09)
            else item
            for item in items]
#
#
def get_beam_load(conn, beam_name:int|str,
                  load_name:int|str, component: int):
    """ """
    line = get_line_load(conn, beam_name,
                         load_name, component)
    point = get_point_load(conn, beam_name,
                           load_name, component)
    return line, point
#
#
# ---------------------------------
#
#
@dataclass
class BeamLoadTypeSQL(BeamTypeBasic):
    __slots__ = ['_system_flag', '_db_file', 
                 '_line', '_point', '_name',
                 '_beam', '_component']

    def __init__(self, load_name: str|int,
                 component: int,
                 plane, 
                 db_file: str): #, beams
        """
        """
        super().__init__()
        self._db_file = db_file
        self._name = load_name
        self._component = component
        #
        self._line = BeamDistributedSQL(load_name=self._name,
                                        component=component,
                                        plane=plane, 
                                        db_file=self._db_file)
        
        self._point = BeamPointSQL(load_name=self._name,
                                   component=component,
                                   plane=plane, 
                                   db_file=self._db_file)
        #
        #self._beam =  beam
        #self._beam_name = beam.name
    #
    #
    #def __call__(self, beam):
    #    """ """
    #    self._beam = beam
    #    return self    
    #
    # ------------------
#
#
#
# ---------------------------------
#
# TODO: Why this global?
class BeamLoadGloabalSQL(BeamSQLMaster):
    """ """
    __slots__ = ['_db_file', '_component', '_plane']
    
    def __init__(self,  component: int, plane, db_file: str) -> None:
        """
        """
        #super().__init__()
        self._db_file =  db_file
        self._component = component
        self._plane = plane
        
        #
        self._line = BeamDistributedSQL(load_name='*',
                                        component=component,
                                        plane=plane, 
                                        db_file=self._db_file)
        self._point = BeamPointSQL(load_name='*',
                                   component=component,
                                   plane=plane, 
                                   db_file=self._db_file)        
    #
    # ------------------
    #
    def __setitem__(self, beam_name: int|str,
                    beam_load: list) -> None:
        """
        """
        1 / 0
        conn = create_connection(self._db_file)
        with conn:
            self._beam =  BeamItemSQL(beam_name,
                                      plane=self._plane,
                                      component=self._component,
                                      db_file=self._db_file)
            #beam = check_element(conn, beam_name)        
        #try:
        #    self._beam.name
        #except (TypeError, IndexError):
        #    raise IOError(f"beam {beam_name} not found")
        #
        # TODO: check if _beam_id affects something else
        #self._beam_id = beam_name
        #self._labels.append(beam_name)
        #
        load_type = beam_load.pop(0)
        #
        if re.match(r"\b(point|node|mass)\b", str(load_type), re.IGNORECASE):
            #self._load(beam).point =  beam_load[1:]
            #self._point[beam_name] = beam_load[1:]
            if isinstance(beam_load, dict):
                load = self._get_point(beam_load)
                load.insert(0, 'force')
                #self._point[beam_name] = load
                1 / 0
            else:
                if isinstance(beam_load[0], (list, tuple)):
                    for item in beam_load:
                        load = self._get_point(item)
                        #load.insert(0, 'force')
                        #self._point[beam_name] = load
                        with conn:
                            push_point_load(conn,
                                            load_name=self._name,
                                            beam_name=beam_name,
                                            point_type='force', 
                                            point_load=load)                        
                else:
                    load =  self._get_point(beam_load)
                    #load.insert(0, 'force')
                    #self._point[beam_name] = load
                    with conn:
                        push_point_load(conn,
                                        load_name=self._name,
                                        beam_name=beam_name,
                                        point_type='force', 
                                        point_load=load)
            
        elif re.match(r"\b(line|udl)\b", str(load_type), re.IGNORECASE):
            #self._load(beam).line = beam_load[1:]
            if isinstance(beam_load, dict):
                1 / 0
                load = self._get_line(beam_load)
                load.insert(0, 'load')
                self._line[beam_name] = load
            else:
                if isinstance(beam_load[0], (list, tuple)):
                    #load = []
                    for item in beam_load:
                        #load.append(self._get_line(item))
                        #load[-1].insert(0, 'load')
                        #self._line[beam_name] = load
                        load = self._get_line(item)
                        with conn:
                            push_line_load(conn,
                                           load_name=self._name,
                                           beam_name=beam_name,
                                           load_type='load', 
                                           udl=load)
                else:
                    load =  self._get_line(beam_load)
                    #load.insert(0, 'load')
                    #self._line[beam_name] = load
                    with conn:
                        push_line_load(conn,
                                       load_name=self._name,
                                       beam_name=beam_name,
                                       load_type='load', 
                                       udl=load)                    
            #
            #self._line[beam_name] = beam_load[1:]
            #with conn:
            #    push_line_load(conn,
            #                   load_name=self._name,
            #                   beam_name=beam_name,
            #                   load_title=None,
            #                   load_system=None,
            #                   load_type='load', 
            #                   udl=load)
            
        else:
            raise IOError(f'Beam lod type {beam_load[0]} not implemented')

    #
    def __getitem__(self, beam_name: int | str):
        """
        """
        conn = create_connection(self._db_file)
        with conn:  
            #beam = check_element(conn, beam_name)
            beam =  BeamItemSQL(beam_name=beam_name,
                                plane=self._plane,
                                component=self._component,
                                db_file=self._db_file)
        #try:
        #memb_type = beam.type # beam[3]
        if beam.type != 'beam':
            raise ValueError(f"element {beam_name} type {beam.type} not valid")
        #except TypeError:
        #    raise IOError(f"beam {beam_name} not found")
        #
        #return BeamLoadTypeSQL(load_name=self._name,
        #                       beam=beam, 
        #                       db_file=self._db_file)
        1 / 0
    #
    # ------------------
    #
    @property
    def point(self):
        """ Concentrated force """
        return self._point
    #
    # ------------------
    #
    @property
    def line(self):
        """
        Linear Varying Load (lvl) - Non Uniformly Distributed Load
        
        value : [qx1, qy1, qz1, qt1,
                 qx2, qy2, qz2, qt2,
                 L1, L2]
    
                        |
             q0         | q1
        o------|        |----------o
        |                          |
        +  L0  +        +    L1    +

        """
        return self._line
    #
    #
    # -----------------------------------------------
    #
    def load_function(self, steps:int,
                      Pa:float, factor:float):
        """ """
        #print('-->')
        beams = BeamSQL(db_file=self._db_file,
                        component=self._component,
                        plane=self._plane)        
        #1 / 0
        #beamfun = defaultdict(list)
        loadfun = []
        # line load
        for key, items in self._line.items():
            beam = beams[key]
            nodes = beam.connectivity
            mat = beam.material
            sec = beam.section.properties(poisson=mat.poisson)
            #
            Lsteps = linstep(d=beam.section.geometry.d,
                             L=beam.L, steps=steps)
            #Lsteps = linspace(start=0, stop=beam.L, num=steps+1, endpoint=True)
            #
            for bitem in items:
                #header = [bitem.load_name, bitem.component_name, *nodes]
                lout = bitem.Fx(x=Lsteps, L=beam.L,
                                E=mat.E, G=mat.G, 
                                Iy=sec.Iy, Iz=sec.Iz,
                                J=sec.J, Cw=sec.Cw,
                                Area=sec.area,
                                Asy=sec.Asy, Asz=sec.Asz,
                                P=Pa, factor=factor)
                #beamfun[key].extend(lout)
                #beamfun[key].append(lout)
                #
                loadfun.extend(lout)
                #loadfun.append([key, bitem.name, lout])
                #loadfun.extend([[bitem.name, 'local', key, *step]
                #                for step in lout])
                #for step in lout:
                #    # load_name, load_title, [Fx, Fy, Fz, Mx, My, Mz]
                #    loadfun.append([bitem.name, 'local', key, *step])
        # point load
        for key, items in self._point.items():
            beam = beams[key]
            mat = beam.material
            sec = beam.section.properties(poisson=mat.poisson)
            #
            Lsteps = linstep(d=beam.section.geometry.d,
                             L=beam.L, steps=steps)
            #Lsteps = linspace(start=0, stop=beam.L, num=steps+1, endpoint=True)
            #
            for bitem in items:
                lout = bitem.Fx(x=Lsteps, L=beam.L,
                                E=mat.E, G=mat.G, 
                                Iy=sec.Iy, Iz=sec.Iz,
                                J=sec.J, Cw=sec.Cw,
                                Area=sec.area,
                                Asy=sec.Asy, Asz=sec.Asz,
                                P=Pa, factor=factor)
                #beamfun[key].append([bitem.name, lout])
                #beamfun[key].extend(lout)
                #beamfun[key].append(lout)
                #
                loadfun.extend(lout)
                #loadfun.append([key, bitem.name, lout])
                #loadfun.extend([[bitem.name, 'local', key, *step]
                #                for step in lout])                
                #for step in lout:
                #    loadfun.append([bitem.name, 'local', key, *step])
        #print('---')
        return loadfun # beamfun #         
    #    
    #
    # -----------------------------------------------
    #
    @property
    def df(self):
        """ beam load df"""
        #print(' beam load')
        #
        line = self._line.df
        line[['L0', 'qx0', 'qy0', 'qz0', 'qt0',
              'L1', 'qx1', 'qy1', 'qz1', 'qt1']] = None        
        line['load_type'] = 'line'
        #
        point = self._point.df
        point[['L0',
               'Fx', 'Fy', 'Fz',
               'Mx', 'My', 'Mz']] = None
        point['load_type'] = 'point'
        #
        db = DBframework()
        df_beam = db.concat([line, point], ignore_index=True, sort=False)
        #df_beam.rename(columns={'L0': 'a', 'L1': 'b'},
        #               inplace=True)
        #
        return df_beam
#
#
# ---------------------------------
#
#
class BeamDistributedSQL(BeamLineBasic): 
    __slots__ = ['_db_file', '_name', '_beam', 
                 '_plane', '_component']
    
    def __init__(self, load_name: str|int,
                 component: int, plane, 
                 db_file: str) -> None:
        """
        """
        super().__init__()
        self._db_file = db_file
        self._plane = plane
        self._name = load_name
        self._component = component
        #
        # create table
        #conn = create_connection(self._db_file)
        #with conn:
        #    self._create_table(conn)         
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._component, )
        table = "SELECT Element.name \
                 FROM Element, LoadBeamLine \
                 WHERE LoadBeamLine.element_id = Element.number \
                 AND Element.component_id = ? ;"
        #
        # Extract data from sqlite
        conn = create_connection(self._db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            rows = cur.fetchall()
        #
        labels = set([item[0] for item in rows])
        return list(labels)    
    #    
    # -----------------------------------------------
    #
    def __setitem__(self, beam_name: int|str, line_load: list) -> None:
        """
        Linear Varying Load (lvl) - Non Uniformly Distributed Load
        value : [qx1, qy1, qz1, qt1,
                 qx2, qy2, qz2, qt2,
                 L1, L2]
    
                        |
             q0         | q1
        o------|        |----------o
        |                          |
        +  L0  +        +    L1    +
        """ 
        #
        #
        conn = create_connection(self._db_file)
        with conn:
            self._beam =  BeamItemSQL(beam_name,
                                      plane=self._plane,
                                      component=self._component, 
                                      db_file=self._db_file)
               
        #try:
        #    self._beam.name
        #except (TypeError, IndexError):
        #    raise IOError(f"beam {beam_name} not found")        
        #
        load_source = 'force' #line_load.pop(0)
        # clean load input
        line_load =  self._get_line(line_load)        
        #
        # get load data
        # set element load
        #self._labels.append(beam_name)
        #title = line_load.pop()
        #self._title.append(title)
        #system = line_load.pop() #line_load[8]
        #
        # push to SQL
        with conn:
            push_line_load(conn,
                           load_name=self._name,
                           beam_name=beam_name,
                           component=self._component,
                           load_type=load_source,
                           udl=line_load)
        #print("-->")
    #
    def __getitem__(self, beam_name:int|str) -> list:
        """
        """
        conn = create_connection(self._db_file)      
        with conn:
            udl = get_line_load(conn,
                                beam_name=beam_name,
                                load_name='*',
                                component=self._component)
        return udl
    #
    # -----------------------------------------------
    #
    def __str__(self, load_name: int|str) -> str:
        """ """
        output = ""
        conn = create_connection(self._db_file) 
        bload = get_line_load(conn, beam_name="*", load_name=load_name)
        #
        for item in bload:
            output += item.__str__()
        #print('---')
        return output
    #
    #
    # -----------------------------------------------
    #
    #def _create_table(self, conn) -> None:
    #    """ """
    #    table = "CREATE TABLE IF NOT EXISTS LoadBeamLine(\
    #                number INTEGER PRIMARY KEY NOT NULL,\
    #                load_id INTEGER NOT NULL REFERENCES Load(number),\
    #                element_id INTEGER NOT NULL REFERENCES Element(number),\
    #                title TEXT,\
    #                system INTEGER NOT NULL,\
    #                type TEXT NOT NULL,\
    #                L0 DECIMAL,\
    #                qx0 DECIMAL,\
    #                qy0 DECIMAL,\
    #                qz0 DECIMAL,\
    #                qt0 DECIMAL,\
    #                L1 DECIMAL,\
    #                qx1 DECIMAL,\
    #                qy1 DECIMAL,\
    #                qz1 DECIMAL,\
    #                qt1 DECIMAL,\
    #                BS DECIMAL,\
    #                OTM DECIMAL,\
    #                x DECIMAL,\
    #                y DECIMAL,\
    #                z DECIMAL);"
    #    #
    #    create_table(conn, table)    
    #
    def _push_load(self, conn, beam_name: str, load: list):
        """ """
        1 / 0
        self._labels.append(beam_name)
        #
        sql = 'INSERT INTO LoadBeamLine(load_id, element_id,\
                                        title, system,\
                                        L0, qx0, qy0, qz0, qt0,\
                                        L1, qx1, qy1, qz1, qt1,\
                                        BS, OTM, x, y, z)\
                                        VALUES(?,?,?,?,\
                                                ?,?,?,?,?,?,?,\
                                                ?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.executemany(sql, load)     
    #
    # -----------------------------------------------
    #
    @property
    def df(self):
        """ """
        db = DBframework()
        conn = create_connection(self._db_file)
        # 
        with conn:
            query = (self._component, )
            table = "SELECT Load.*, \
                    Element.name, \
                    LoadBeamLine.title, LoadBeamLine.system,\
                    LoadBeamLine.L0, LoadBeamLine.qx0, LoadBeamLine.qy0, \
                    LoadBeamLine.qz0, LoadBeamLine.qt0,\
                    LoadBeamLine.L1, LoadBeamLine.qx1, LoadBeamLine.qy1, \
                    LoadBeamLine.qz1, LoadBeamLine.qt1 \
                    FROM Load, Element, LoadBeamLine \
                    WHERE LoadBeamLine.load_id = Load.number\
                    AND LoadBeamLine.element_id = Element.number \
                    AND Load.component_id = ? ;"
            #
            cur = conn.cursor()
            cur.execute(table, query)            
            rows = cur.fetchall()
        #
        #beam_load = []
        #for row in rows:
        #    data = [*row[:2]]
        #
        cols = ['load_id','load_name', 'load_title', 'load_level', 'load_type',
                'element_name',
                'load_comment', 'load_system',
                'L0', 'qx0', 'qy0', 'qz0', 'qt0',
                'L1', 'qx1', 'qy1', 'qz1', 'qt1']
        df = db.DataFrame(data=rows, columns=cols)
        #
        df = df[['load_name', 'load_level', 'load_id', 'load_system', 'load_comment',
                 'element_name',
                'L0', 'qx0', 'qy0', 'qz0', 'qt0',
                'L1', 'qx1', 'qy1', 'qz1', 'qt1']]
        #       
        #print('--->')
        return df
    
    @df.setter
    def df(self, df):
        """ """
        conn = create_connection(self._db_file)
        with conn:
            nodes, elements, basic = get_load_basics(conn, self._component)
        #
        df['load_id'] = df['name'].apply(lambda x: basic[x])
        df['element_id'] =  df['beam'].apply(lambda x: elements[x])
        #df['node_id'] = df['node'].apply(lambda x: nodes[x])
        df['system'] = 'local'
        df['type'] = df['type'].apply(lambda x: x.lower())         
        #
        header = ['load_id', 'element_id',
                  'title', 'system', 'type',
                  'L0', 'qx0', 'qy0', 'qz0', 'qt0',
                  'L1', 'qx1', 'qy1', 'qz1', 'qt1',
                  'BS', 'OTM', 'x', 'y', 'z']
        #
        bconn = df[header].copy()
        bconn.replace(to_replace=[''], value=[float(0)], inplace=True)
        bconn['title'] = df['title']
        #
        with conn:
            bconn.to_sql('LoadBeamLine', conn,
                         index_label=header,
                         if_exists='append', index=False)
        #print('--->')
#
#
#
def push_line_load(conn, load_name: str|int,
                   beam_name:int|str,
                   component: int,
                   load_type: str, 
                   udl:list[float],
                   load_system: str = "local"):
    """ """  
    #
    # Beam check
    beam = check_element(conn, beam_name,
                         component=component)
    beam_number = beam[0]
    #
    # Load check
    load_data = get_load_data(conn, load_name,
                              load_level='basic',
                              component=component)
    try:
        load_id = load_data[0]
    except TypeError:
        raise IOError(f"Load {load_name} not found")    
    #
    load_title = udl.pop()
    system = udl.pop()
    try:
        1 / system
    except ZeroDivisionError:
        raise RuntimeWarning('load in global system')
    #
    if load_type in ['wave']:
        raise NotImplemented
    
    else:
        project = (load_id, beam_number,
                   load_title, load_system,
                   load_type, 
                   udl[8], *udl[:4],   # L0, qx0, qy0, qz0, qt0
                   udl[9], *udl[4:8],  # L1, qx1, qy1, qz1, qt1
                   None, None,         # BS, OTM,
                   None, None, None,)  # x, y, z
    #
    sql = 'INSERT INTO LoadBeamLine(load_id, element_id,\
                                        title, system, type,\
                                        L0, qx0, qy0, qz0, qt0,\
                                        L1, qx1, qy1, qz1, qt1,\
                                        BS, OTM, x, y, z)\
                                        VALUES(?,?,?,?,?,\
                                               ?,?,?,?,?,?,?,\
                                               ?,?,?,?,?,?,?,?)'
    #
    with conn:
        cur = conn.cursor()
        cur.execute(sql, project)
        #cur.executemany(sql, udl)
#
#
def get_line_load(conn, beam_name:int|str,
                  load_name:int|str, component: int):
    """ """
    query = [component]
    table = "SELECT Load.name, Element.name, \
             LoadBeamLine.*, Component.name \
             FROM Load, Element, LoadBeamLine, Component \
             WHERE LoadBeamLine.load_id = Load.number \
             AND LoadBeamLine.element_id = Element.number \
             AND Component.number = ? "
    #
    # get beam data
    if beam_name in ['*', '']:
        pass
    else:
        query.extend([beam_name])
        table += f"AND  Element.name = ? "
    #
    # beam line load
    if load_name in ['*', '']:
        pass
    else:
        query.extend([load_name])
        table += f"AND Load.name = ? "
    #
    table += " ;"
    #
    # Extract data from sqlite
    with conn:
        cur = conn.cursor()
        cur.execute(table, tuple(query))
        rows = cur.fetchall()
    #
    #
    beam_line = [] # defaultdict(list)
    for row in rows:
        #if row[7] in ['wave']:
        #    raise NotImplemented
        #else:
        beam_line.append(LineBeam(row[1], row[5],    # name, title,
                                  row[0], row[-1],   # load_name, component_name
                                  row[6],            # system,
                                  #
                                  *row[9:13],        # q_inplane [qx, qy, qz, qt]
                                  *row[14:18],       # q_outplane [qx, qy, qz, qt]
                                  row[8], row[12]))  # L0, L1
    return beam_line
#
#
#
class BeamPointSQL(BeamPointBasic):
    __slots__ = ['_db_file', '_name', '_beam', 
                 '_component', '_plane']

    def __init__(self, load_name: str|int,
                 component: int, plane,
                 db_file: str) -> None:
        """
        """
        super().__init__()
        self._db_file = db_file
        self._plane = plane
        self._name = load_name
        self._component = component
        #
        # create table
        #conn = create_connection(self._db_file)
        #with conn:
        #    self._create_table(conn)
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._component, )
        table = "SELECT Element.name \
                 FROM Element, LoadBeamPoint \
                 WHERE LoadBeamPoint.element_id = Element.number \
                 AND Element.component_id = ? ;"
        #
        # Extract data from sqlite
        conn = create_connection(self._db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            rows = cur.fetchall()
        #
        labels = set([item[0] for item in rows])
        return list(labels)
    #    
    # -----------------------------------------------
    #
    def __setitem__(self, beam_name: int|str, point_load: list) -> None:
        """
        """
        #
        conn = create_connection(self._db_file)
        with conn:
            self._beam =  BeamItemSQL(beam_name,
                                      plane=self._plane,
                                      component=self._component, 
                                      db_file=self._db_file)
               
        #try:
        #    self._beam.name
        #except (TypeError, IndexError):
        #    raise IOError(f"beam {beam_name} not found")        
        #
        point_type = 'force' #point_load.pop(0)
        point_load =  self._get_point(point_load)
        # get load data
        # set element load
        #self._labels.append(beam_name)
        #title = point_load.pop()
        #self._title.append(title)
        #system = point_load.pop() #line_load[8]        
        # 
        # push to SQL
        conn = create_connection(self._db_file)
        with conn:
            push_point_load(conn,
                            load_name=self._name,
                            beam_name=beam_name,
                            component=self._component,
                            point_type=point_type,
                            point_load=point_load)
        # print("-->")
    #
    def __getitem__(self, beam_name:int|str)-> list:
        """
        """
        #bd_file = self._db_file
        # get beam load
        conn = create_connection(self._db_file)
        with conn:
            pl = get_point_load(conn, beam_name=beam_name,
                                load_name='*',
                                component=self._component)
        return pl
    #
    # -----------------------------------------------
    #
    def __str__(self, load_name: int|str) -> str:
        """ """
        output = ""
        conn = create_connection(self.db_file)
        bload = get_point_load(conn, beam_name="*",
                               component=self._component, 
                               load_name=load_name)
        #
        for item in bload:
            output += item.__str__()
        #print('---')
        return output
    #
    # -----------------------------------------------
    #
    #def _create_table(self, conn) -> None:
    #    """ """
    #    table = "CREATE TABLE IF NOT EXISTS LoadBeamPoint(\
    #                number INTEGER PRIMARY KEY NOT NULL,\
    #                load_id INTEGER NOT NULL REFERENCES Load(number),\
    #                element_id INTEGER NOT NULL REFERENCES Element(number),\
    #                title TEXT,\
    #                system INTEGER NOT NULL,\
    #                type TEXT NOT NULL,\
    #                L0 DECIMAL,\
    #                fx DECIMAL,\
    #                fy DECIMAL,\
    #                fz DECIMAL,\
    #                mx DECIMAL,\
    #                my DECIMAL,\
    #                mz DECIMAL,\
    #                x DECIMAL,\
    #                y DECIMAL,\
    #                z DECIMAL,\
    #                rx DECIMAL,\
    #                ry DECIMAL,\
    #                rz DECIMAL);"
    #    #
    #    create_table(conn, table)
    #
    # -----------------------------------------------
    #
    @property
    def df(self):
        """ """
        db = DBframework()
        conn = create_connection(self._db_file)
        #
        with conn:
            query = (self._component, )
            table = "SELECT Load.*, \
                        Element.name, \
                        LoadBeamPoint.title, LoadBeamPoint.system,\
                        LoadBeamPoint.L0, LoadBeamPoint.fx, LoadBeamPoint.fy, LoadBeamPoint.fz, \
                        LoadBeamPoint.mx, LoadBeamPoint.my, LoadBeamPoint.mz \
                        FROM Load, Element, LoadBeamPoint \
                        WHERE LoadBeamPoint.load_id = Load.number\
                        AND LoadBeamPoint.element_id = Element.number \
                        AND Load.component_id = ? ;"
            #
            cur = conn.cursor()
            cur.execute(table, query)            
            rows = cur.fetchall()
        #
        cols = ['load_id','load_name', 'load_title', 'load_level', 'load_type', 
                'element_name',
                'load_comment', 'load_system',
                'L0', 'Fx', 'Fy', 'Fz',
                'Mx', 'My', 'Mz']
        df = db.DataFrame(data=rows, columns=cols)
        #
        df = df[['load_name', 'load_level', 'load_id', 'load_system', 'load_comment',
                 'element_name',
                 'L0', 'Fx', 'Fy', 'Fz',
                 'Mx', 'My', 'Mz']]
        #       
        return df
    
    @df.setter
    def df(self, df):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            nodes, elements, basic = get_load_basics(conn, self._component)
        #
        df['load_id'] = df['name'].apply(lambda x: basic[x])
        df['element_id'] =  df['beam'].apply(lambda x: elements[x])
        #df['node_id'] = df['node'].apply(lambda x: nodes[x])
        df['system'] = 'local'
        df['type'] = df['type'].apply(lambda x: x.lower())        
        #
        header = ['load_id', 'element_id',
                  'title', 'system', 'type', 'L0', 
                  'fx', 'fy', 'fz', 'mx', 'my', 'mz',
                  'x', 'y', 'z', 'rx', 'ry', 'rz']        
        #
        bconn = df[header].copy()
        bconn.replace(to_replace=[''], value=[float(0)], inplace=True)
        bconn['title'] = df['title']
        #
        with conn:
            bconn.to_sql('LoadBeamPoint', conn,
                         index_label=header, 
                         if_exists='append', index=False)
        #print('--->')
#
#
def push_point_load(conn, load_name: str|int,
                    beam_name:int|str,
                    component: int, 
                    point_type: str, 
                    point_load:list[float],
                    load_system: str = "local"):
    """ """
    # Beam check 
    beam = check_element(conn, beam_name,
                         component=component)
    beam_number = beam[0]
    #
    # Load check
    load_data = get_load_data(conn, load_name,
                              load_level='basic',
                              component=component)
    try:
        load_id = load_data[0]
    except TypeError:
        raise IOError(f"Load {load_name} not found")  
    #
    #
    load_title = point_load.pop()
    system = point_load.pop()
    try:
        1 / system
    except ZeroDivisionError:
        raise RuntimeWarning('load in global system')    
    #
    #
    if re.match(r"\b(force|point)\b", point_type, re.IGNORECASE):
        project = (load_id, beam_number,
                   load_title, load_system, 'force', 
                   point_load[6],   # L0
                   *point_load[:6], # fx, fy, fz, mx, my, mz
                   None, None, None,
                   None, None, None,)
    else:
        raise NotImplemented
    #
    sql = 'INSERT INTO LoadBeamPoint(load_id, element_id,\
                                    title, system, type, \
                                    L0, fx, fy, fz, mx, my, mz,\
                                    x, y, z, rx, ry, rz)\
                                    VALUES(?,?,?,?,?,\
                                           ?,?,?,?,?,?,?,\
                                           ?,?,?,?,?,?)'
    with conn:
        cur = conn.cursor()
        cur.execute(sql, project)
#
#
def get_point_load(conn, beam_name:int|str,
                   load_name:int|str,
                   component: int):
    """ """
    query = [component]
    #
    table = "SELECT Load.name, Element.name, \
            LoadBeamPoint.*, Component.name \
            FROM Load, Element, LoadBeamPoint, Component \
            WHERE LoadBeamPoint.load_id = Load.number \
            AND LoadBeamPoint.element_id = Element.number \
            AND Component.number = ? "
    #
    #
    # get beam data
    if beam_name in ['*', '']:
        pass
    else:
        query.extend([beam_name])
        table += f"AND  Element.name = ? "        
    #
    #
    # beam line load
    if load_name in ['*', '']:
        pass
    else:
        query.extend([load_name])
        table += f"AND Load.name = ? "
    #
    table += " ;"
    #
    # Extract data from sqlite
    with conn:
        cur = conn.cursor()
        cur.execute(table, tuple(query))
        rows = cur.fetchall()
    #
    #
    beam_point = []
    for row in rows:
        if row[7] in ['mass']:
            raise NotImplemented
        else:
            beam_point.append(PointBeam(row[1], row[5],  # name, title, 
                                        row[0], row[-1], # load_name, component_name
                                        row[6],          # system,
                                        #
                                        *row[9:15],      # fx, fy, fz, mx, my, mz
                                        row[8]))         # L0
    return beam_point
#
#
def push_FER(conn, node_load:list):
    """
    update sql with beam's Fixed End Reactions
    """
    #
    sql = 'INSERT INTO LoadBeamFER(load_id, element_id, node_id, \
                                   title, system, type, \
                                   fx, fy, fz, mx, my, mz,\
                                   x, y, z, rx, ry, rz,\
                                   Psi, B, Tw) \
                        VALUES(?,?,?,?,?,\
                               ?,?,?,?,?,?,?,\
                               ?,?,?,?,?,?,\
                               ?,?,?)'
    cur = conn.cursor()
    cur.executemany(sql, node_load)
#
#

