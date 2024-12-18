#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime as dt
import time

#
# package imports
from steelpy.trave.postprocess.sql.Un import UnSQL
from steelpy.trave.postprocess.sql.beam import BeamResSQL
from steelpy.trave.postprocess.sql.node import NodeResSQL
from steelpy.trave.postprocess.utils.operations import MainPostProcess
from steelpy.utils.sqlite.utils import create_table, create_connection

from steelpy.ufo.mesh.sqlite.utils import (pull_node_mesh,
                                           pull_element_mesh,
                                           pull_load_mesh,
                                           pull_results_mesh)
#
from steelpy.utils.dataframe.main import DBframework


#
class PostProcessSQL:
    __slots__ = ['_mesh', '_process', '_Un', '_Pdelta',
                 'db_file', '_result_name', '_results']
    
    def __init__(self, mesh, result_name:int|str,
                 db_file: str) -> None:
        """
        """
        self._mesh = mesh
        self._result_name = result_name
        self.db_file = db_file
        #
        self._Un = UnSQL(mesh=self._mesh,
                         result_name=result_name,
                         db_file=self.db_file)
        #
        self._process = MainProcessSQL(mesh=self._mesh,
                                       result_name=result_name,
                                       db_file=self.db_file)
        #
        self._results = ResultSQL(mesh=self._mesh,
                                  result_name=result_name,
                                  db_file=self.db_file)
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
    def run(self, beam_steps: int= 10):
        """ """
        print("** Postprocessing")
        #
        Pdelta: bool = self._Pdelta
        beam_force, Qn = self._process.solve(self._Un,
                                             beam_steps,
                                             Pdelta=Pdelta)
        # load comb update
        if not Pdelta:
            # combination
            load_comb = self._mesh._load.combination()
            lcomb = load_comb.to_basic()
            beam_force, Qn = self._process._add_comb(beam_force, Qn, lcomb)
        # Node run
        self._process.node_reactions(Qn)
        # element run
        self._process.element_run(beam_force)
        # element stress run
        self._process.solve_stress(beam_force=beam_force)
        #
        #print('-->')
#
#
@dataclass
class MainProcessSQL(MainPostProcess):
    __slots__ = ['_plane', '_mesh', '_result_name', '_db_file']

    def __init__(self, mesh, result_name:int|str,
                 db_file: str) -> None: # name: str,
        """
        """
        super().__init__(mesh, result_name)
        self._db_file = db_file
        #self._result_name = result_name
        # TODO : link results
        #conn = create_connection(self._db_file)
        #with conn:
        #    self._new_table(conn)
        #    self._push_result(conn)
    #
    #
    # -----------------------------------------------------------
    # Operations
    # update load
    def _add_comb(self, beam_force, Qn, lcomb):
        """
        Update load with combinations

        Qn:
        beam_force:
        lcomb:
        """
        beam_force = update_member_df(member=beam_force,
                                      lcomb=lcomb,
                                      item='length',
                                      values=['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                                              'x', 'y', 'z', 'rx', 'ry', 'rz',
                                              'Psi', 'B', 'Tw'])
        #
        # End node force combination
        Qn = update_member_df(member=Qn,
                              lcomb=lcomb,
                              item='node_name',
                              values=self._plane.hforce)
        #
        return beam_force, Qn
    #
    # Node
    def node_reactions(self, Qn:DBframework.DataFrame)-> None:
        """ Push node results to database"""
        Qndf = Qn.copy()
        mesh_name=self._mesh._name
        conn = create_connection(self._db_file)
        with conn:
            ## node reactions
            nforce, nreac = self.get_reactions(Qndf)
            #
            item = 'element_name'
            element_id =  pull_element_mesh(conn, mesh_name)
            update_name2id_mesh(items=element_id,
                                Qn=Qndf,
                                component=item,
                                new_name='element_id')
            #
            load_id = pull_load_mesh(conn, mesh_name)
            item = 'load_name'
            update_name2id_mesh(items=load_id,
                                Qn=Qndf,
                                component=item,
                                new_name='load_id')

            update_name2id_mesh(items=load_id,
                                Qn=nforce,
                                component=item,
                                new_name='load_id')

            update_name2id_mesh(items=load_id,
                                Qn=nreac,
                                component=item,
                                new_name='load_id')
            #
            result_id = pull_results_mesh(conn, mesh_name)
            item = 'result_name'
            update_name2id_mesh(items=result_id,
                                Qn=Qndf,
                                component=item,
                                new_name='result_id')

            update_name2id_mesh(items=result_id,
                                Qn=nforce,
                                component=item,
                                new_name='result_id')

            update_name2id_mesh(items=result_id,
                                Qn=nreac,
                                component=item,
                                new_name='result_id')
            #
            # node reactions
            #
            node_id = pull_node_mesh(conn, mesh_name)
            item = 'node_name'
            update_name2id_mesh(items=node_id,
                                Qn=Qndf,
                                component=item,
                                new_name='node_id')
            #
            update_name2id_mesh(items=node_id,
                                Qn=nforce,
                                component=item,
                                new_name='node_id')
            #
            update_name2id_mesh(items=node_id,
                                Qn=nreac,
                                component=item,
                                new_name='node_id')
        #
        Qndf.rename(columns={'load_system': 'system'}, inplace=True)
        nforce.rename(columns={'load_system': 'system'}, inplace=True)
        nreac.rename(columns={'load_system': 'system'}, inplace=True)
        #
        # Push to sql
        header = ['result_id', 'load_id',
                  'element_id', 'node_id',
                  'system', *self._plane.hforce]
        #
        conn = create_connection(self._db_file)
        with conn:
            Qndf[header].to_sql('ResultMemberEndForce', conn,
                                index_label=header,
                                if_exists='append',
                                index=False)
            #
            header.remove('element_id')
            nforce[header].to_sql('ResultNodeForce', conn,
                                  index_label=header,
                                  if_exists='append',
                                  index=False)
            #
            nreac[header].to_sql('ResultNodeReaction', conn,
                                 index_label=header,
                                 if_exists='append',
                                 index=False)
        #
    #
    # Elements
    def element_run(self, beam_force:DBframework.DataFrame)-> DBframework.DataFrame:
        """ """
        beamf = beam_force.copy()
        mesh_name=self._mesh._name
        conn = create_connection(self._db_file)
        with conn:
            item = 'element_name'
            element_id =  pull_element_mesh(conn, mesh_name)
            update_name2id_mesh(items=element_id,
                                Qn=beamf,
                                component=item,
                                new_name='element_id')
            #
            load_id = pull_load_mesh(conn, mesh_name)
            item = 'load_name'
            update_name2id_mesh(items=load_id,
                                Qn=beamf,
                                component=item,
                                new_name='load_id')
            #
            result_id = pull_results_mesh(conn, mesh_name)
            item = 'result_name'
            update_name2id_mesh(items=result_id,
                                Qn=beamf,
                                component=item,
                                new_name='result_id')
        #
        beamf.rename(columns={'load_system': 'system'}, inplace=True)
        #
        df_Qbeam = beamf[['result_id', 'load_id', #'mesh_name',
                          'element_id', 'length', 'system',
                          'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz',
                          'Psi', 'B', 'Tw']]
        #
        df_Dbeam = beamf[['result_id', 'load_id',
                          'element_id', 'length', 'system',
                          'x', 'y', 'z', 'rx', 'ry', 'rz']]
        # Push to sql
        conn = create_connection(self._db_file)
        with conn:
            df_Qbeam.to_sql('ResultBeamForce', conn,
                            if_exists='append', index=False)
            #
            df_Dbeam.to_sql('ResultBeamDeflection', conn,
                            if_exists='append', index=False)
    #
    def get_reactions(self, node_force:DBframework.DataFrame)-> (DBframework.DataFrame, DBframework.DataFrame):
        """ get nodal reactions """
        #
        header = ['mesh_name', 'result_name',
                  'load_name', 'load_level',
                  'node_name', 'system']
        #
        nforce = (node_force.groupby(header,
                                     as_index=False)
                   [self._plane.hforce].sum())
        #
        # TODO : maybe there is a better pandas way that the code below..
        supports = self._mesh._boundaries.support()
        nname = list(supports)
        nrsupp = nforce.loc[nforce['node_name'].isin(nname)]
        db = DBframework()
        nrsupp = db.DataFrame(data=nrsupp)
        return nforce, nrsupp
    #
    def solve_stress(self, beam_force:DBframework.DataFrame):
        """get elements stress"""
        #print("** Beam stress calculation")
        start_time = time.time()
        #
        beamfdf = beam_force.copy()
        mesh_name=self._mesh._name
        elements = self._mesh._elements
        memb_grp = beamfdf.groupby(['mesh_name', 'result_name',
                                   'load_name',  'load_level',
                                   'system', 'element_name'])
        #
        for key, memb_item in memb_grp:
            member = elements[key[-1].tolist()]
            section = member.section
            material = member.material
            beam_stress = section.stress(beam_result=memb_item,
                                         E=material.E,
                                         G=material.G,
                                         poisson=material.poisson)
            conn = create_connection(self._db_file)
            with conn:
                element_id =  pull_element_mesh(conn, mesh_name)
                item = 'element_name'
                update_name2id_mesh(items=element_id,
                                    Qn=beam_stress,
                                    component=item,
                                    new_name='element_id')
                #
                load_id = pull_load_mesh(conn, mesh_name)
                item = 'load_name'
                update_name2id_mesh(items=load_id,
                                    Qn=beam_stress,
                                    component=item, 
                                    new_name='load_id')
                #
                result_id = pull_results_mesh(conn, mesh_name)
                item = 'result_name'
                update_name2id_mesh(items=result_id,
                                    Qn=beam_stress,
                                    component=item,
                                    new_name='result_id')
                #
                beam_stress.rename(columns={'load_system': 'system'}, inplace=True)
                #
                # Push stress data
                header = ['result_id', 'load_id',
                          'element_id', 'length', 'system',
                          'stress_point', 'y', 'z',
                          'tau_x', 'tau_y', 'tau_z',
                          'sigma_x', 'sigma_y', 'sigma_z']

                beam_stress[header].to_sql('ResultBeamStress', conn,
                                           if_exists='append', index=False)
        #
        uptime = time.time() - start_time
        print(f"** Beam Stress Calculation: {uptime:1.4e} sec")
    #
    #
    # -----------------------------------------------------------
    # SQL ops
    #
    def _new_table(self, conn) -> None:
        """ """
        table = "CREATE TABLE IF NOT EXISTS Result (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    type TEXT,\
                    mesh_id INTEGER REFERENCES Mesh(number),\
                    units TEXT NOT NULL,\
                    plane TEXT NOT NULL,\
                    date TEXT);"
        #
        create_table(conn, table)
    #
    #
    def _push_result(self, conn) -> None:
        """ """
        table = 'INSERT INTO Result(name, type, mesh_id, \
                                    units, plane, date)\
                                    VALUES(?,?,?,?,?,?)'
        #
        time = dt.now().strftime('%Y-%m-%d')
        plane = '3D'
        if self._plane.plane2D:
            plane = '2D'
        data = (self._name, None, self._mesh._name,  # name, type, mesh_name
                'si', plane, time)                   # units, plane, date
        # push
        cur = conn.cursor()
        cur.execute(table, data)
#
#
# -----------------------------------------------------------
#
#
@dataclass
class ResultSQL:
    __slots__ = ['_mesh', '_node', '_beam', 'shell',
                 '_db_file','_plane']
    
    def __init__(self, mesh, result_name:int|str,
                 db_file: str) -> None:
        """
        """
        self._mesh = mesh
        self._db_file = db_file
        self._node = NodeResSQL(mesh=self._mesh,
                                result_name=result_name,
                                db_file=self._db_file)
        self._beam = BeamResSQL(mesh=self._mesh,
                                result_name=result_name,
                                db_file=self._db_file)
        self._plane = mesh._plane
        #
        conn = create_connection(self._db_file)
        with conn:
            self._new_table(conn)
    #
    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += self._node.__str__()
        output += self._beam.__str__()
        return output
    #
    # -----------------------------------------------------------
    # SQL ops
    #
    def _new_table(self, conn) -> None:
        """ """
        #
        # ---------------------------------------------------
        # Node
        # ---------------------------------------------------
        #
        # Node end forces
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultNodeForce (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        result_id INTEGER NOT NULL REFERENCES Result(number),\
                        load_id INTEGER NOT NULL REFERENCES Load(number),\
                        node_id INTEGER NOT NULL REFERENCES Node(number),\
                        system TEXT NOT NULL,\
                        Fx DECIMAL,\
                        Fy DECIMAL,\
                        Fz DECIMAL,\
                        Mx DECIMAL,\
                        My DECIMAL,\
                        Mz DECIMAL);"
        create_table(conn, table_nodes)
        #
        # Node Reactions
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultNodeReaction (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        result_id INTEGER NOT NULL REFERENCES Result(number),\
                        load_id INTEGER NOT NULL REFERENCES Load(number),\
                        node_id INTEGER NOT NULL REFERENCES Node(number),\
                        system TEXT NOT NULL,\
                        Fx DECIMAL,\
                        Fy DECIMAL,\
                        Fz DECIMAL,\
                        Mx DECIMAL,\
                        My DECIMAL,\
                        Mz DECIMAL);"
        create_table(conn, table_nodes)
        #
        # ---------------------------------------------------
        # Member global end forces
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultMemberEndForce (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        result_id INTEGER NOT NULL REFERENCES Result(number),\
                        load_id INTEGER NOT NULL REFERENCES Load(number),\
                        element_id INTEGER NOT NULL REFERENCES Element(number),\
                        node_id INTEGER NOT NULL REFERENCES Node(number),\
                        system TEXT NOT NULL,\
                        Fx DECIMAL,\
                        Fy DECIMAL,\
                        Fz DECIMAL,\
                        Mx DECIMAL,\
                        My DECIMAL,\
                        Mz DECIMAL);"
        create_table(conn, table_nodes)
        #
        # ---------------------------------------------------
        # Beam
        # ---------------------------------------------------
        #
        # Beam local deflection
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultBeamDeflection (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        result_id INTEGER NOT NULL REFERENCES Result(number),\
                        load_id INTEGER NOT NULL REFERENCES Load(number),\
                        element_id INTEGER NOT NULL REFERENCES Element(number),\
                        length DECIMAL NOT NULL,\
                        system TEXT NOT NULL,\
                        x DECIMAL,\
                        y DECIMAL,\
                        z DECIMAL,\
                        rx DECIMAL,\
                        ry DECIMAL,\
                        rz DECIMAL);"
        create_table(conn, table_nodes)
        #
        # Beam Local internal forces
        table_nodes = ("CREATE TABLE IF NOT EXISTS ResultBeamForce (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        result_id INTEGER NOT NULL REFERENCES Result(number),\
                        load_id INTEGER NOT NULL REFERENCES Load(number),\
                        element_id INTEGER NOT NULL REFERENCES Element(number),\
                        length DECIMAL NOT NULL,\
                        system TEXT NOT NULL,\
                        Fx DECIMAL,\
                        Fy DECIMAL,\
                        Fz DECIMAL,\
                        Mx DECIMAL,\
                        My DECIMAL,\
                        Mz DECIMAL,\
                        psi DECIMAL,\
                        B DECIMAL,\
                        Tw DECIMAL);")
        create_table(conn, table_nodes)
        #
        # Beam local stress
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultBeamStress (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        result_id INTEGER NOT NULL REFERENCES Result(number),\
                        load_id INTEGER NOT NULL REFERENCES Load(number),\
                        element_id INTEGER NOT NULL REFERENCES Element(number),\
                        length DECIMAL NOT NULL,\
                        system TEXT NOT NULL,\
                        stress_point INTEGER NOT NULL,\
                        y DECIMAL,\
                        z DECIMAL,\
                        tau_x DECIMAL,\
                        tau_y DECIMAL,\
                        tau_z DECIMAL,\
                        sigma_x DECIMAL,\
                        sigma_y DECIMAL,\
                        sigma_z DECIMAL);"
        create_table(conn, table_nodes)
        #
        # ---------------------------------------------------
        # Shell
        # ---------------------------------------------------
        #
    #
    #
    # -----------------------------------------------------------
    #
    #def reactions(self):
    #    """ """
    #    pass
    #
    def nodes(self):
        """node results"""
        return self._node
    #
    def beam(self):
        """node results"""
        return self._beam
    #
    # -----------------------------------------------------------
    #
    #
    #def results(self, df_beamf, df_Qn,
    #            Pdelta: bool):
    #    """
    #    beam_steps : Integration points beam element (10 default)
    #    """
    #    #print("** Postprocessing")
    #    #
    #    #
    #    #1 / 0
    #    #return self
    #
    # -----------------------------------------------------------
    #
    # -----------------------------------------------------------
    #
#
# -----------------------------------------------------------
# Combination utils
# -----------------------------------------------------------
#
#
def update_member_df(member, lcomb,
                     item: str,
                     values:list[str]):
    """
    Update node displacements to include lcomb
    """
    db = DBframework()
    # group basic load by name
    grp = member.groupby('load_name')
    combgrp = lcomb.groupby('load_name')
    for key, combfactors in combgrp:
        for row in combfactors.itertuples():
            comb = grp.get_group(row.basic_load).copy()
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
        dftemp = dftemp.groupby(['mesh_name', 'result_name', 
                                 'load_name', 'load_level', 
                                 'element_name', item, 'system'],
                                  as_index=False)[values].sum()
        #test
        member = db.concat([member, dftemp], ignore_index=True)
    except UnboundLocalError:
        pass
    #
    return member
#
#
def update_name2id(conn, function, Qn,
                   item:int|str, item2:str,
                   new_name:str):
    """ """
    Qn[item] = [function(conn,
                         name=name,
                         mesh_name=mesh_name)
                for name, mesh_name in zip(Qn[item], Qn[item2])]

    Qn.rename(columns={item: new_name}, inplace=True)
    #return Qn
#
#
def update_name2id_mesh(items, Qn,
                        component:str|int,
                        new_name:str):
    """ """
    items = {item[1]: item[0] for item in items}
    Qn[component] = [items[name] for name in Qn[component]]
    Qn.rename(columns={component: new_name}, inplace=True)
    #return Qn
#