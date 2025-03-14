#
# Copyright (c) 2009 steelpy
# 
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
#from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
#from typing import NamedTuple
import re
from math import isclose
import os
#
# package imports
#
from steelpy.ufo.mesh.sqlite.beam import BeamItemSQL, BeamSQL
from steelpy.ufo.mesh.sqlite.utils import check_element #, check_nodes

from steelpy.ufo.load.process.beam.operations import LineBeam, PointBeam

from steelpy.ufo.load.process.beam.main import (BeamTypeBasic,
                                                BeamLineBasic,
                                                BeamPointBasic,
                                                BeamLoadBasic)

from steelpy.ufo.load.process.beam.utils import (get_BeamLoad_df,
                                                 find_BeamLoad_item,
                                                 find_load_type)

from steelpy.ufo.load.sqlite.utils import pull_basic, get_load_basics

from steelpy.utils.math.operations import linstep
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework
#
#
#
#
# ---------------------------------
#
def fer_b2n(beam, bload, global_system:int):
    """"""
    node1, node2 = beam.nodes
    # print(f'line load {key}')
    lnload = bload.fer_beam(beam=beam,
                            system=global_system)
    #
    load7 = zerofilter(lnload[7])
    load8 = zerofilter(lnload[8])
    #
    return [bload.load_name,
            bload.comment,
            global_system,
            beam.number,
            node1.number,
            load7,
            node2.number,
            load8,
            bload.load_step]
#
#
# ---------------------------------
#
#
class BeamLoadItemSQL(BeamLoadBasic):
    __slots__ = ['_name',  '_load', '_mesh_id', 
                 '_node_eq', 'db_file', '_system_flag', 
                 '_beam']

    def __init__(self, load_name: str|int, #plane: NamedTuple,
                 mesh_id: int, 
                 db_file: str) -> None:
        """
        """
        super().__init__()
        #
        self._name = load_name
        self._mesh_id = mesh_id
        #
        self.db_file = db_file
        # create table
        conn = create_connection(self.db_file)
        with conn:
            self._new_table(conn)
        #
        #
        self._load = BeamLoadTypeSQL(load_name=self._name,
                                     mesh_id=self._mesh_id,
                                     db_file=self.db_file)
    #
    @property
    def _labels(self):
        """ """
        query = (self._name, self._mesh_id)
        #
        # point
        table = "SELECT Element.name \
                FROM Element, LoadBeamPoint, LoadBasic, Load \
                WHERE LoadBeamPoint.element_id = Element.number \
                AND LoadBasic.number = LoadBeamPoint.basic_id \
                AND LoadBasic.load_id = Load.number \
                AND Load.name = ? \
                AND Load.mesh_id = ? ;"
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
                FROM Element, LoadBeamLine, LoadBasic, Load \
                WHERE LoadBeamLine.element_id = Element.number \
                AND LoadBasic.number = LoadBeamLine.basic_id \
                AND LoadBasic.load_id = Load.number \
                AND Load.name = ? AND Load.mesh_id = ? ;"
        #table += load_name
        #
        conn = create_connection(self.db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(table, query)
            rows = cur.fetchall()
        labels.extend([item[0] for item in rows])
        #
        return list(set(labels))
    #    
    # -----------------------------------------------
    #
    #
    def __setitem__(self, name: int|str,
                    load: list|tuple|dict) -> None:
        """
        """
        if isinstance(load, dict):
            columns = list(load.keys())
            header = {item:find_BeamLoad_item(item) for item in columns}
            for key, item in header.items():
                if item in ['type']:
                    load_type = find_load_type(load[key])
                    break
        else:
            load_type = find_load_type(load[0])
        #
        if re.match(r"\b(point|mass)\b", load_type, re.IGNORECASE):
            self._load._point[name] = load
        elif re.match(r"\b(line|u(niform(ly)?)?d(istributed)?l(oad)?)\b",
                      load_type, re.IGNORECASE):
            self._load._line[name] = load
        else:
            raise NotImplementedError(f'Beam lod {load_type} not implemented')    
    #
    def __getitem__(self, beam_name: int | str):
        """
        """
        conn = create_connection(self.db_file)
        with conn:  
            beam =  BeamItemSQL(beam_name=beam_name,
                                mesh_id=self._mesh_id, 
                                db_file=self.db_file)
        #
        if beam.type != 'beam':
            raise ValueError(f"element {beam_name} type {beam.type} not valid")
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
                                        mesh_id=self._mesh_id)
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
    #
    def fer(self, beams) -> None:
        """ Push Fix End Reactions (FER) global system in sqlite """
        conn = create_connection(self.db_file)
        with conn:
            basic = pull_basic(conn,
                               load_name=self._name)
            line, point = get_beam_load(conn,
                                        load_name=self._name,
                                        beam_name='*',
                                        mesh_id=self._mesh_id)
        try:
            basic_id =  basic[0]
        except IndexError:
            raise IOError(f"Load {self._name} not found")
        #
        # ------------------------------------------
        #
        psystem = os.name
        if psystem == 'posix':
            cpuno = max(1, os.cpu_count() - 1)
            executor = ProcessPoolExecutor(max_workers=cpuno)
        else:
            executor = ThreadPoolExecutor()
        #
        # ------------------------------------------
        #
        global_system = 0
        b2n = []
        with executor:
            for item in line:
                beam = beams[item.name]
                b2n.append(executor.submit(fer_b2n, beam, item, global_system).result())
            #
            for item in point:
                beam = beams[item.name]
                b2n.append(executor.submit(fer_b2n, beam, item, global_system).result())
        #
        fer_res =[]
        for gnload in b2n:
            try:
                1 / gnload[2]
                raise RuntimeError('node load in local system')
            except ZeroDivisionError:
                load_system = 'global'
            #
            # [fx, fy, fz, mx, my, mz, rx, ry,rz, Psi, B, Tw, step]
            fer_res.extend([[basic_id, gnload[3], gnload[4],
                             gnload[1], load_system, 'reaction', *gnload[5], gnload[8]],
                           [basic_id, gnload[3], gnload[6],
                            gnload[1], load_system, 'reaction', *gnload[7], gnload[8]]])
        #
        # FER results to database
        if fer_res:
            with conn:
                push_FER(conn, node_load=fer_res)
    #
    # -----------------------------------------------
    #
    def _new_table(self, conn) -> None:
        """ """
        # -------------------------------------
        # line 
        table = "CREATE TABLE IF NOT EXISTS LoadBeamLine(\
                number INTEGER PRIMARY KEY NOT NULL,\
                basic_id INTEGER NOT NULL REFERENCES LoadBasic(number),\
                element_id INTEGER NOT NULL REFERENCES Element(number),\
                comment TEXT,\
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
                comment TEXT,\
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
                comment TEXT,\
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
                                  mesh_id=self._mesh_id)
        bload
        1 / 0
    
    @df.setter
    def df(self, df):
        """ """
        dfdict = get_BeamLoad_df(df, system=self._system_flag)
        #
        conn = create_connection(self.db_file)
        with conn:
            nodes, elements, basic = get_load_basics(conn, self._mesh_id)
        #
        for key, item in  dfdict.items():
            bconn = item.copy()
            bconn['basic_id'] = bconn['name'].apply(lambda x: basic[x])
            bconn['element_id'] =  bconn['beam'].apply(lambda x: elements[x])
            bconn.replace(to_replace=[''], value=[float(0)], inplace=True)
            #
            header = ['basic_id', 'element_id',
                      'comment', 'system', 'type', 'L0',
                      'fx', 'fy', 'fz', 'mx', 'my', 'mz',
                      'x', 'y', 'z', 'rx', 'ry', 'rz']
            #
            if key in ['point']:
                bconn[['x', 'y', 'z', 'rx', 'ry', 'rz']] = None
                with conn:
                    bconn[header].to_sql('LoadBeamPoint', conn,
                                         index_label=header, 
                                         if_exists='append', index=False)
            elif key in ['line']:
                header = ['basic_id', 'element_id',
                          'comment', 'system', 'type',
                          'L0', 'qx0', 'qy0', 'qz0', 'qt0',
                          'L1', 'qx1', 'qy1', 'qz1', 'qt1',
                          'step']
                with conn:
                    bconn[header].to_sql('LoadBeamLine', conn,
                                         index_label=header,
                                         if_exists='append', index=False)                
            
            elif key in ['mass']:
                bconn[['fx', 'fy', 'fz', 'mx', 'my', 'mz']] = None
                with conn:
                    bconn[header].to_sql('LoadBeamPoint', conn,
                                         index_label=header, 
                                         if_exists='append', index=False)
            
            else:
                raise IOError(f'beam load {key} not valid')
    #
    #   
    #
    # -----------------------------------------------
    #
    def load_functionXX(self, Fb, steps:int,
                      factor:float):
        """ """
        #print('-->')
        beams = BeamSQL(db_file=self.db_file,
                        mesh_id=self._mesh_id)
        #
        conn = create_connection(self.db_file)
        with conn:
            line, point = get_beam_load(conn,
                                        beam_name='*',
                                        load_name=self._name,
                                        mesh_id=self._mesh_id)
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
            Fx = round(Pdelta.Fx, 6)
            if isclose(a=abs(Fx), b=0.0):
                Fx = 0.0
            #for bitem in items:
            lout = item.Fx(x=Lsteps, L=beam.L,
                            E=mat.E, G=mat.G,
                            Iy=sec.Iy, Iz=sec.Iz,
                            J=sec.J, Cw=sec.Cw,
                            Area=sec.area,
                            Asy=sec.Asy, Asz=sec.Asz,
                            P=Fx, factor=factor)
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
            lout = item.Fx(x=Lsteps, L=beam.L,
                            E=mat.E, G=mat.G,
                            Iy=sec.Iy, Iz=sec.Iz,
                            J=sec.J, Cw=sec.Cw,
                            Area=sec.area,
                            Asy=sec.Asy, Asz=sec.Asz,
                            P=Pdelta.Fx, factor=factor)
            #
            loadfun.extend([[Pdelta.load_name, *item]
                            for item in lout])
        #
        #1 / 0
        return loadfun
#
# ---------------------------------
#
def zerofilter(items: list):
    """remove near zero items"""
    return [0 if isclose(abs(item), 0.0, abs_tol=1e-09)
            else item
            for item in items]
#
#
def get_beam_load(conn, beam_name:int|str,
                  load_name:int|str, mesh_id: int):
    """ """
    line = get_line_load(conn, beam_name,
                         load_name, mesh_id)
    point = get_point_load(conn, beam_name,
                           load_name, mesh_id)
    return line, point
#
# ---------------------------------
#
#
@dataclass
class BeamLoadTypeSQL(BeamTypeBasic):
    __slots__ = ['_system_flag', '_db_file', 
                 '_line', '_point', '_name',
                 '_beam', '_mesh_id']

    def __init__(self, load_name: str|int,
                 mesh_id: int,
                 db_file: str):
        """
        """
        super().__init__()
        self._db_file = db_file
        self._name = load_name
        self._mesh_id = mesh_id
        #
        self._line = BeamDistributedSQL(load_name=self._name,
                                        mesh_id=mesh_id,
                                        db_file=self._db_file)
        
        self._point = BeamPointSQL(load_name=self._name,
                                   mesh_id=mesh_id,
                                   db_file=self._db_file)
    #
    #
    #
    # -----------------------------------------------
    #
    def function(self, steps:int,
                 Pa:float, factor:float)->list[list]:
        """
        Return:
               load_function = [load_name, mesh_name,
                                load_comment, load_type,
                                load_level, load_system,
                                element_name, length,
                                axial, torsion, VM_inplane, VM_outplane]
        """
        #print('-->')
        beam = self._beam
        Lb = beam.L
        mat = beam.material
        sec = beam.section.properties
        Asy, Asz = beam.section.As(poisson=mat.poisson)
        geometry = beam.section.geometry
        #
        Lsteps = linstep(d=geometry.d,
                         L=Lb, steps=steps)
        #
        load_function = []
        # line load
        line_load = self._line[self._beam_id]
        for bitem in line_load:
            lout = bitem.Fx(x=Lsteps, L=Lb,
                            E=mat.E, G=mat.G, 
                            Iy=sec.Iy, Iz=sec.Iz,
                            J=sec.J, Cw=sec.Cw,
                            Area=sec.area,
                            Asy=Asy, Asz=Asz,
                            P=Pa, factor=factor)
            load_function.extend(lout)
        # point load
        point_load = self._point[self._beam_id]
        for bitem in point_load:
            lout = bitem.Fx(x=Lsteps, L=Lb,
                            E=mat.E, G=mat.G, 
                            Iy=sec.Iy, Iz=sec.Iz,
                            J=sec.J, Cw=sec.Cw,
                            Area=sec.area,
                            Asy=Asy, Asz=Asz,
                            P=Pa, factor=factor)
            load_function.extend(lout)
        #print('---')
        return load_function
    #
    # ------------------
#
#
#
# ---------------------------------
#
# TODO: Why this global?
class BeamLoadGloabalSQL(BeamLoadBasic):
    """ """
    __slots__ = ['_db_file', '_mesh_id',
                 '_line', '_point', '_load']
    
    def __init__(self,  mesh_id: int, db_file: str) -> None:
        """
        """
        super().__init__()
        #
        self._db_file =  db_file
        self._mesh_id = mesh_id
        #
        #self._line = BeamDistributedSQL(load_name='*',
        #                                mesh_id=mesh_id,
        #                                db_file=self._db_file)
        #
        #self._point = BeamPointSQL(load_name='*',
        #                           mesh_id=mesh_id,
        #                           db_file=self._db_file)
        #
        self._load = BeamLoadTypeSQL(load_name='*',
                                     mesh_id=self._mesh_id,
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
                                      mesh_id=self._mesh_id,
                                      db_file=self._db_file)
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
        #1 / 0
        conn = create_connection(self._db_file)
        with conn:  
            beam =  BeamItemSQL(beam_name=beam_name,
                                mesh_id=self._mesh_id,
                                db_file=self._db_file)
        #
        if beam.type != 'beam':
            raise ValueError(f"element {beam_name} type {beam.type} not valid")
        
        return self._load(beam=beam)
    #
    # ------------------
    #
    @property
    def point(self):
        """ Concentrated force """
        return self._load._point
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
        return self._load._line
    #
    #
    # -----------------------------------------------
    #
    def load_function(self, steps:int,
                      Pa:float, factor:float):
        """ """
        #print('-->')
        beams = BeamSQL(db_file=self._db_file,
                        mesh_id=self._mesh_id)
                        #plane=self._plane)        
        1 / 0
        loadfun = []
        # line load
        for key, items in self._load._line.items():
            beam = beams[key]
            mat = beam.material
            sec = beam.section.properties(poisson=mat.poisson)
            #
            Lsteps = linstep(d=beam.section.geometry.d,
                             L=beam.L, steps=steps)
            #
            for bitem in items:
                lout = bitem.Fx(x=Lsteps, L=beam.L,
                                E=mat.E, G=mat.G, 
                                Iy=sec.Iy, Iz=sec.Iz,
                                J=sec.J, Cw=sec.Cw,
                                Area=sec.area,
                                Asy=sec.Asy, Asz=sec.Asz,
                                P=Pa, factor=factor)
                loadfun.extend(lout)
        # point load
        for key, items in self._load._point.items():
            beam = beams[key]
            mat = beam.material
            sec = beam.section.properties(poisson=mat.poisson)
            #
            Lsteps = linstep(d=beam.section.geometry.d,
                             L=beam.L, steps=steps)
            #
            for bitem in items:
                lout = bitem.Fx(x=Lsteps, L=beam.L,
                                E=mat.E, G=mat.G, 
                                Iy=sec.Iy, Iz=sec.Iz,
                                J=sec.J, Cw=sec.Cw,
                                Area=sec.area,
                                Asy=sec.Asy, Asz=sec.Asz,
                                P=Pa, factor=factor)
                loadfun.extend(lout)
        #print('---')
        return loadfun
    #    
    #
    # -----------------------------------------------
    #
    @property
    def df(self):
        """ beam load df"""
        #print(' beam load')
        #
        line = self._load._line.df
        line[['L0', 'qx0', 'qy0', 'qz0', 'qt0',
              'L1', 'qx1', 'qy1', 'qz1', 'qt1']] = None        
        line['load_type'] = 'line'
        #
        point = self._load._point.df
        point[['L0',
               'Fx', 'Fy', 'Fz',
               'Mx', 'My', 'Mz']] = None
        point['load_type'] = 'point'
        #
        db = DBframework()
        df_beam = db.concat([line, point], ignore_index=True, sort=False)
        #
        return df_beam
#
#
# ---------------------------------
#
class BeamDistributedSQL(BeamLineBasic): 
    __slots__ = ['_db_file', '_name', '_beam', 
                 '_mesh_id']
    
    def __init__(self, load_name: str|int,
                 mesh_id: int, 
                 db_file: str) -> None:
        """
        """
        super().__init__()
        self._db_file = db_file
        self._name = load_name
        self._mesh_id = mesh_id        
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._mesh_id, )
        table = "SELECT Element.name \
                 FROM Element, LoadBeamLine \
                 WHERE LoadBeamLine.element_id = Element.number \
                 AND Element.mesh_id = ? ;"
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
    def __setitem__(self, beam_name: int|str,
                    line_load: list|dict) -> None:
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
                                      mesh_id=self._mesh_id, 
                                      db_file=self._db_file)
        #
        load_source = 'force' #line_load.pop(0)
        # clean load input
        line_load =  self._get_line(line_load)
        #
        # push to SQL
        with conn:
            push_line_load(conn,
                           load_name=self._name,
                           beam_name=beam_name,
                           mesh_id=self._mesh_id,
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
                                mesh_id=self._mesh_id)
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
    def _push_load(self, conn, beam_name: str, load: list):
        """ """
        1 / 0
        self._labels.append(beam_name)
        #
        sql = 'INSERT INTO LoadBeamLine(basic_id, element_id,\
                                        title, system,\
                                        L0, qx0, qy0, qz0, qt0,\
                                        L1, qx1, qy1, qz1, qt1,\
                                        BS, OTM)\
                                        VALUES(?,?,?,?,\
                                               ?,?,?,?,?,?,?,\
                                               ?,?,?,?,?)'
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
            query = (self._mesh_id, )
            table = "SELECT Load.number, Load.name, Load.level,\
                    Load.title, LoadBeamLine.step,\
                    Element.name, \
                    LoadBeamLine.comment, LoadBeamLine.system,\
                    LoadBeamLine.L0, LoadBeamLine.qx0, LoadBeamLine.qy0, \
                    LoadBeamLine.qz0, LoadBeamLine.qt0,\
                    LoadBeamLine.L1, LoadBeamLine.qx1, LoadBeamLine.qy1, \
                    LoadBeamLine.qz1, LoadBeamLine.qt1 \
                    FROM Load, Element, LoadBeamLine, LoadBasic \
                    WHERE LoadBeamLine.basic_id = LoadBasic.number\
                    AND LoadBeamLine.element_id = Element.number \
                    AND LoadBasic.load_id = Load.number \
                    AND Load.mesh_id = ? ;"
            #
            cur = conn.cursor()
            cur.execute(table, query)            
            rows = cur.fetchall()
        #
        # FIXME: step tb sorted
        cols = ['load_id','load_name', 'load_level', 
                'load_title', 'step',
                'element_name',
                'load_comment', 'load_system',
                'L0', 'qx0', 'qy0', 'qz0', 'qt0',
                'L1', 'qx1', 'qy1', 'qz1', 'qt1']
        df = db.DataFrame(data=rows, columns=cols)
        #
        df = df[['load_name', 'load_level', 'load_id',
                 'load_system', 'load_comment',
                 'element_name',
                'L0', 'qx0', 'qy0', 'qz0', 'qt0',
                'L1', 'qx1', 'qy1', 'qz1', 'qt1']]
        #       
        return df
    
    @df.setter
    def df(self, df):
        """ """
        1 / 0 # FIXME : basic load id
        conn = create_connection(self._db_file)
        with conn:
            nodes, elements, basic = get_load_basics(conn, self._mesh_id)
        #
        df['basic_id'] = df['name'].apply(lambda x: basic[x])
        df['element_id'] =  df['beam'].apply(lambda x: elements[x])
        #df['node_id'] = df['node'].apply(lambda x: nodes[x])
        df['system'] = 'local'
        df['type'] = df['type'].apply(lambda x: x.lower())         
        #
        header = ['basic_id', 'element_id',
                  'comment', 'system', 'type',
                  'L0', 'qx0', 'qy0', 'qz0', 'qt0',
                  'L1', 'qx1', 'qy1', 'qz1', 'qt1',
                  'BS', 'OTM']
        #
        bconn = df[header].copy()
        bconn.replace(to_replace=[''], value=[float(0)], inplace=True)
        #bconn['title'] = df['title']
        #
        with conn:
            bconn.to_sql('LoadBeamLine', conn,
                         index_label=header,
                         if_exists='append', index=False)
        #print('--->')
#
#
def push_line_load(conn, load_name: str|int,
                   beam_name:int|str,
                   mesh_id: int,
                   load_type: str, 
                   udl:list[float],
                   load_system: str = "local"):
    """ """  
    #
    # Beam check
    beam = check_element(conn, beam_name,
                         mesh_id=mesh_id)
    beam_id = beam[0]
    #
    # Load check
    basic = pull_basic(conn, load_name=load_name)
    try:
        basic_id =  basic[0]
    except IndexError:
        raise IOError(f"Load {load_name} not found")
    #
    #udl2 = list(zip(*udl))
    #
    load_comment= udl.pop()
    system = udl.pop()
    try:
        1 / system
    except ZeroDivisionError:
        raise RuntimeWarning('load in global system')
    #
    #
    project = (basic_id, beam_id,
               load_comment, load_system,
               load_type, 
               udl[8], *udl[:4],   # L0, qx0, qy0, qz0, qt0
               udl[9], *udl[4:8],  # L1, qx1, qy1, qz1, qt1
               None)               # step
               #None, None, None,)  # x, y, z
    #
    sql = 'INSERT INTO LoadBeamLine(basic_id, element_id,\
                                    comment, system, type,\
                                    L0, qx0, qy0, qz0, qt0,\
                                    L1, qx1, qy1, qz1, qt1,\
                                    step)\
                                    VALUES(?,?,?,?,?,\
                                           ?,?,?,?,?,?,?,\
                                           ?,?,?,?)'
    #
    with conn:
        cur = conn.cursor()
        cur.execute(sql, project)
        #cur.executemany(sql, project)
#
#
def get_line_load(conn, beam_name:int|str,
                  load_name:int|str, mesh_id: int):
    """ """
    query = [mesh_id]
    table = "SELECT Element.name, LoadBeamLine.comment,\
             Load.name, Mesh.name, LoadBeamLine.system,\
             LoadBeamLine.qx0, LoadBeamLine.qy0, LoadBeamLine.qz0, LoadBeamLine.qt0,\
             LoadBeamLine.qx1, LoadBeamLine.qy1, LoadBeamLine.qz1, LoadBeamLine.qt1,\
             LoadBeamLine.L0, LoadBeamLine.L1,\
             LoadBeamLine.step, LoadBeamLine.type \
             FROM Load, Element, LoadBasic, LoadBeamLine, Mesh \
             WHERE LoadBeamLine.basic_id = LoadBasic.number \
             AND LoadBeamLine.element_id = Element.number \
             AND LoadBasic.load_id = Load.number \
             AND Mesh.number = ? "
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
    beam_line = []
    for row in rows:
        beam_line.append(LineBeam(row[0], row[1],    # name, comment,
                                  row[2], row[3],   # load_name, mesh_name
                                  row[4],            # system, 
                                  #
                                  *row[5:9],        # q_inplane [qx, qy, qz, qt]
                                  *row[9:13],       # q_outplane [qx, qy, qz, qt]
                                  row[13], row[14],   # L0, L1
                                  row[15]))          # load_step
    return beam_line
#
# ---------------------------------
#
class BeamPointSQL(BeamPointBasic):
    __slots__ = ['_db_file', '_name', '_beam', 
                 '_mesh_id']

    def __init__(self, load_name: str|int,
                 mesh_id: int,
                 db_file: str) -> None:
        """
        """
        super().__init__()
        self._db_file = db_file
        self._name = load_name
        self._mesh_id = mesh_id
    #
    #
    @property
    def _labels(self):
        """ """
        query = (self._mesh_id, )
        table = "SELECT Element.name \
                 FROM Element, LoadBeamPoint \
                 WHERE LoadBeamPoint.element_id = Element.number \
                 AND Element.mesh_id = ? ;"
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
        conn = create_connection(self._db_file)
        with conn:
            self._beam =  BeamItemSQL(beam_name,
                                      mesh_id=self._mesh_id, 
                                      db_file=self._db_file)
                      
        #
        point_type = 'force'
        point_load =  self._get_point(point_load)       
        # 
        # push to SQL
        conn = create_connection(self._db_file)
        with conn:
            push_point_load(conn,
                            load_name=self._name,
                            beam_name=beam_name,
                            mesh_id=self._mesh_id,
                            load_type=point_type,
                            point_load=point_load)
        # print("-->")
    #
    def __getitem__(self, beam_name:int|str)-> list:
        """
        """
        # get beam load
        conn = create_connection(self._db_file)
        with conn:
            pl = get_point_load(conn, beam_name=beam_name,
                                load_name='*',
                                mesh_id=self._mesh_id)
        return pl
    #
    # -----------------------------------------------
    #
    def __str__(self, load_name: int|str) -> str:
        """ """
        output = ""
        conn = create_connection(self.db_file)
        bload = get_point_load(conn, beam_name="*",
                               mesh_id=self._mesh_id, 
                               load_name=load_name)
        #
        for item in bload:
            output += item.__str__()
        #print('---')
        return output
    #
    # -----------------------------------------------
    #
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
            query = (self._mesh_id, )
            table = "SELECT Load.number, Load.name, Load.level, Load.title, \
                        Element.name, \
                        LoadBeamPoint.comment, LoadBeamPoint.system,\
                        LoadBeamPoint.L0, LoadBeamPoint.fx, LoadBeamPoint.fy, LoadBeamPoint.fz, \
                        LoadBeamPoint.mx, LoadBeamPoint.my, LoadBeamPoint.mz \
                    FROM Load, Element, LoadBasic, LoadBeamPoint \
                    WHERE LoadBeamPoint.basic_id = LoadBasic.number\
                        AND LoadBeamPoint.element_id = Element.number \
                        AND LoadBasic.load_id = Load.number \
                        AND Load.mesh_id = ? ;"
            #
            cur = conn.cursor()
            cur.execute(table, query)            
            rows = cur.fetchall()
        #
        cols = ['load_id','load_name', 'load_level', 'load_title', 
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
        1 / 0 # FIXME : basic load id
        conn = create_connection(self.db_file)
        with conn:
            nodes, elements, basic = get_load_basics(conn, self._mesh_id)
        #
        df['basic_id'] = df['name'].apply(lambda x: basic[x])
        df['element_id'] =  df['beam'].apply(lambda x: elements[x])
        #df['node_id'] = df['node'].apply(lambda x: nodes[x])
        df['system'] = 'local'
        df['type'] = df['type'].apply(lambda x: x.lower())        
        #
        header = ['basic_id', 'element_id',
                  'comment', 'system', 'type', 'L0',
                  'fx', 'fy', 'fz', 'mx', 'my', 'mz',
                  'x', 'y', 'z', 'rx', 'ry', 'rz']        
        #
        bconn = df[header].copy()
        bconn.replace(to_replace=[''], value=[float(0)], inplace=True)
        #bconn['title'] = df['title']
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
                    mesh_id: int, 
                    load_type: str, 
                    point_load:list[float],
                    load_system: str = "local"):
    """ """
    # Beam check 
    beam = check_element(conn, beam_name,
                         mesh_id=mesh_id)
    beam_id = beam[0]
    #
    # Load check 
    basic = pull_basic(conn, load_name=load_name)
    try:
        basic_id =  basic[0]
    except IndexError:
        raise IOError(f"Load {load_name} not found")
    #
    load_comment = point_load.pop()
    system = point_load.pop()
    try:
        1 / system
    except ZeroDivisionError:
        raise RuntimeWarning('load in global system')    
    #  
    #
    if re.match(r"\b(force|point)\b", load_type, re.IGNORECASE):
        project = (basic_id, beam_id,
                   load_comment, load_system, 'force',
                   point_load[6],    # L0
                   *point_load[:6],  # fx, fy, fz, mx, my, mz
                   None, None, None, # x, y, z, 
                   None, None, None, # rx, ry, rz,
                   None)             # step
    else:
        raise NotImplemented
    #
    sql = 'INSERT INTO LoadBeamPoint(basic_id, element_id,\
                                     comment, system, type, \
                                     L0, fx, fy, fz, mx, my, mz,\
                                     x, y, z, rx, ry, rz,\
                                     step)\
                                    VALUES(?,?,?,?,?,\
                                           ?,?,?,?,?,?,?,\
                                           ?,?,?,?,?,?,?)'
    with conn:
        cur = conn.cursor()
        cur.execute(sql, project)
#
#
def get_point_load(conn, beam_name:int|str,
                   load_name:int|str,
                   mesh_id: int):
    """ """
    query = [mesh_id]
    table = "SELECT Element.name, LoadBeamPoint.comment,\
            Load.name, Mesh.name, LoadBeamPoint.system,\
            LoadBeamPoint.fx, LoadBeamPoint.fy, LoadBeamPoint.fz, \
            LoadBeamPoint.mx, LoadBeamPoint.my, LoadBeamPoint.mz, \
            LoadBeamPoint.L0, LoadBeamPoint.step, LoadBeamPoint.type \
            FROM Load, Element, LoadBeamPoint, LoadBasic, Mesh \
            WHERE LoadBeamPoint.basic_id = LoadBasic.number \
            AND LoadBeamPoint.element_id = Element.number \
            AND LoadBasic.load_id = Load.number \
            AND Mesh.number = ? "
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
        if row[13] in ['mass']:
            raise NotImplemented
        else:
            beam_point.append(PointBeam(row[0], row[1],   # name, comment,
                                        row[2], row[3],   # load_name, mesh_id_name
                                        row[4],           # system,
                                        #
                                        *row[5:11],        # fx, fy, fz, mx, my, mz
                                        row[11], row[12])) # L0, load_step
    return beam_point
#
# ---------------------------------
#
def push_FER(conn, node_load:list):
    """
    update sql with beam's Fixed End Reactions
    """
    #
    sql = 'INSERT INTO LoadBeamFER(basic_id, element_id, node_id, \
                                   comment, system, type, \
                                   fx, fy, fz, mx, my, mz,\
                                   x, y, z, rx, ry, rz,\
                                   Psi, B, Tw, step) \
                        VALUES(?,?,?,?,?,\
                               ?,?,?,?,?,?,?,\
                               ?,?,?,?,?,?,\
                               ?,?,?,?)'
    cur = conn.cursor()
    cur.executemany(sql, node_load)
#
#

