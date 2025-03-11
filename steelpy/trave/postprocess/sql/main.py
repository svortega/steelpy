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
from steelpy.trave.postprocess.sql.beam import BeamResultSQL
from steelpy.trave.postprocess.sql.node import (NodeResultSQL,
                                                get_displacement,
                                                get_reactions)
from steelpy.trave.postprocess.utils.main import MainPostProcess
#from steelpy.trave.postprocess.utils.node import NodeDeflection
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
        self._results = ResultSQL(mesh=self._mesh,
                                  result_name=result_name,
                                  db_file=self.db_file)
        #
        self._process = MainPostProcess(mesh=self._mesh,
                                        result_name=result_name,
                                           #db_file=self.db_file,
                                        results=self._results)
    #
    # ---------------------------------
    # Results
    # ---------------------------------
    #
    @property
    def Un(self)-> DBframework:
        """ """
        #units:str='si'
        conn = create_connection(self.db_file)
        with conn:         
            df = get_displacement(conn,
                                  result_name=self._result_name,
                                  mesh_name=self._mesh._name)
            
        #return self._results._node.displacement()
        #return NodeDeflection(df, units=units)
        return df

    @Un.setter
    def Un(self, Un:DBframework)-> None:
        """ """
        header = ['result_id', 'node_id',
                  'load_id', 'system',
                  *self._mesh._plane.hdisp]
        #          #'x', 'y', 'z', 'rx', 'ry', 'rz']
        conn = create_connection(self.db_file)
        with conn:
            results_id = pull_results_mesh(conn, mesh_name=self._mesh._name)
            results_id = {item[1]:item[0] for item in results_id}
            Un['result_id'] = [results_id[item] for item in Un['result_name']]
            #
            node_id = pull_node_mesh(conn, mesh_name=self._mesh._name)
            node_id ={item[1]: item[0] for item in node_id}
            Un['node_id'] = [node_id[item] for item in Un['node_name']]
            #
            load_id = pull_load_mesh(conn, mesh_name=self._mesh._name)
            load_id = {item[1]: item[0] for item in load_id}
            Un['load_id'] = [load_id[item] for item in Un['load_name']]
            #
            Un.rename(columns={'load_system':'system'}, inplace=True)
            Un[header].to_sql('ResultNodeDisplacement', conn,
                               index_label=header,
                               if_exists='append', index=False)
        #
        #print('nodal disp results saved')
    #
    #
    @property
    def Rn(self)-> DBframework:
        """ """
        conn = create_connection(self.db_file)
        with conn:
            df = get_reactions(conn,
                               result_name=self._result_name,
                               mesh_name=self._mesh._name)
        return df
    
    @Rn.setter
    def Rn(self, Rn:DBframework)-> None:
        """ """
        self._results._push_node_reaction(Rn)
    #
    #def _pull_displacement(self, conn):
    #    """ """
    #    table = ("SELECT Mesh.name, Result.name, Load.name, Load.Level,\
    #                     Node.name, ResultNodeDisplacement.system,\
    #                     ResultNodeDisplacement.x, ResultNodeDisplacement.y, ResultNodeDisplacement.z,\
    #                     ResultNodeDisplacement.rx, ResultNodeDisplacement.ry, ResultNodeDisplacement.rz\
    #             FROM ResultNodeDisplacement, Load, Mesh, Node, Result \
    #             WHERE Load.number = ResultNodeDisplacement.load_id \
    #             AND Node.number = ResultNodeDisplacement.node_id \
    #             AND Result.number = ResultNodeDisplacement.result_id \
    #             AND Result.mesh_id = Mesh.number ;")
    #    #
    #    cur = conn.cursor()
    #    cur.execute(table)
    #    rows = cur.fetchall()
    #    return rows
    #
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
        self._node = NodeResultSQL(mesh=self._mesh,
                                   result_name=result_name,
                                   db_file=self._db_file)
        self._beam = BeamResultSQL(mesh=self._mesh,
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
        #table_nodes = "CREATE TABLE IF NOT EXISTS ResultNodeForce (\
        #                number INTEGER PRIMARY KEY NOT NULL,\
        #                result_id INTEGER NOT NULL REFERENCES Result(number),\
        #                load_id INTEGER NOT NULL REFERENCES Load(number),\
        #                node_id INTEGER NOT NULL REFERENCES Node(number),\
        #                system TEXT NOT NULL,\
        #                Fx DECIMAL,\
        #                Fy DECIMAL,\
        #                Fz DECIMAL,\
        #                Mx DECIMAL,\
        #                My DECIMAL,\
        #                Mz DECIMAL);"
        #create_table(conn, table_nodes)
        #
        # Node displacement
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultNodeDisplacement (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        result_id INTEGER NOT NULL REFERENCES Result(number),\
                        node_id NOT NULL REFERENCES Node(number),\
                        load_id INTEGER NOT NULL REFERENCES Load(number),\
                        system TEXT NOT NULL,\
                        x DECIMAL NOT NULL,\
                        y DECIMAL NOT NULL,\
                        z DECIMAL,\
                        rx DECIMAL,\
                        ry DECIMAL,\
                        rz DECIMAL NOT NULL);"
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
        #table_nodes = "CREATE TABLE IF NOT EXISTS ResultMemberEndForce (\
        #                number INTEGER PRIMARY KEY NOT NULL,\
        #                result_id INTEGER NOT NULL REFERENCES Result(number),\
        #                load_id INTEGER NOT NULL REFERENCES Load(number),\
        #                element_id INTEGER NOT NULL REFERENCES Element(number),\
        #                node_id INTEGER NOT NULL REFERENCES Node(number),\
        #                system TEXT NOT NULL,\
        #                Fx DECIMAL,\
        #                Fy DECIMAL,\
        #                Fz DECIMAL,\
        #                Mx DECIMAL,\
        #                My DECIMAL,\
        #                Mz DECIMAL);"
        #create_table(conn, table_nodes)
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
    def _new_tableXX(self, conn) -> None:
        """ """
        1 / 0
        # Result table
        table = "CREATE TABLE IF NOT EXISTS Result (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    type TEXT,\
                    mesh_id INTEGER REFERENCES Mesh(number),\
                    units TEXT NOT NULL,\
                    date TEXT);"
        create_table(conn, table)
        # Steps table
        table = "CREATE TABLE IF NOT EXISTS ResultStep (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    result_id INTEGER NOT NULL REFERENCES Result(number),\
                    load_id INTEGER NOT NULL REFERENCES Load(number),\
                    type TEXT NOT NULL,\
                    step NOT NULL);"
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
        data = (self._name, None, self._mesh._name, 'si', time)
        #
        cur = conn.cursor()
        cur.execute(table, data)
    #
    # -----------------------------------------------------------
    #
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
    #
    # -----------------------------------------------------------
    #
    # -----------------------------------------------------------
    #
    # Node
    def _push_node_reaction(self, Qn:DBframework.DataFrame)-> None:
        """ Push node results to database"""
        #Qndf = Qn.copy()
        values = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        #header = ['mesh_name', 'result_name',
        #          'load_name', 'load_level',
        #          'node_name', 'system',
        #          *values]
        #1/0
        nreac = Qn.copy()  #[header]
        mesh_name=self._mesh._name
        conn = create_connection(self._db_file)
        with conn:
            ## node reactions
            #nforce, nreac = self.get_reactions(Qndf)
            #
            #item = 'element_name'
            #element_id =  pull_element_mesh(conn, mesh_name)
            #update_name2id_mesh(items=element_id,
            #                    Qn=Qndf,
            #                    component=item,
            #                    new_name='element_id')
            #
            load_id = pull_load_mesh(conn, mesh_name)
            item = 'load_name'
            #update_name2id_mesh(items=load_id,
            #                    Qn=Qndf,
            #                    component=item,
            #                    new_name='load_id')

            #update_name2id_mesh(items=load_id,
            #                    Qn=nforce,
            #                    component=item,
            #                    new_name='load_id')

            update_name2id_mesh(items=load_id,
                                Qn=nreac,
                                component=item,
                                new_name='load_id')
            #
            result_id = pull_results_mesh(conn, mesh_name)
            item = 'result_name'
            #update_name2id_mesh(items=result_id,
            #                    Qn=Qndf,
            #                    component=item,
            #                    new_name='result_id')

            #update_name2id_mesh(items=result_id,
            #                    Qn=nforce,
            #                    component=item,
            #                    new_name='result_id')

            update_name2id_mesh(items=result_id,
                                Qn=nreac,
                                component=item,
                                new_name='result_id')
            #
            # node reactions
            #
            node_id = pull_node_mesh(conn, mesh_name)
            item = 'node_name'
            #update_name2id_mesh(items=node_id,
            #                    Qn=Qndf,
            #                    component=item,
            #                    new_name='node_id')
            #
            #update_name2id_mesh(items=node_id,
            #                    Qn=nforce,
            #                    component=item,
            #                    new_name='node_id')
            #
            update_name2id_mesh(items=node_id,
                                Qn=nreac,
                                component=item,
                                new_name='node_id')
        #
        #Qndf.rename(columns={'load_system': 'system'}, inplace=True)
        #nforce.rename(columns={'load_system': 'system'}, inplace=True)
        #nreac.rename(columns={'load_system': 'system'}, inplace=True)
        #
        # Push to sql
        header = ['result_id', 'load_id',
                  'element_id', 'node_id',
                  'system', #*values]
                  *self._plane.hforce]
        #
        conn = create_connection(self._db_file)
        with conn:
            #Qndf[header].to_sql('ResultMemberEndForce', conn,
            #                    index_label=header,
            #                    if_exists='append',
            #                    index=False)
            #
            header.remove('element_id')
            #nforce[header].to_sql('ResultNodeForce', conn,
            #                      index_label=header,
            #                      if_exists='append',
            #                      index=False)
            #
            nreac[header].to_sql('ResultNodeReaction', conn,
                                 index_label=header,
                                 if_exists='append',
                                 index=False)
        #
    #
    # Elements
    def _push_beam_result(self, beam_force:DBframework.DataFrame)-> None:
        """ """
        #1 / 0
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
    #
    def _push_beam_stress(self, beam_stress:DBframework.DataFrame)-> None:
        """
        """
        mesh_name = self._mesh._name
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
# -----------------------------------------------------------
# Combination utils
# -----------------------------------------------------------
#
#
def update_ndf(self, dfnode, dfcomb,
                values :list[str]|None = None)-> DBframework.DataFrame:
    """
    Update node displacements to include lcomb
    """
    # set default 3D model
    if not values: 
        values = ['x', 'y', 'z', 'rx', 'ry', 'rz']
    #
    db = DBframework()
    # group basic load by name
    ndgrp = dfnode.groupby('load_name')
    # get combinations with basic loads
    combgrp = dfcomb.groupby('load_name')
    dftemp = []
    for key, combfactors in combgrp:
        for row in combfactors.itertuples():
            comb = ndgrp.get_group(row.basic_load).copy()
            comb.loc[:, values] *= row.factor
            comb['load_level'] = 'combination'
            comb['load_name'] = row.load_name
            comb['load_id'] = row.load_id
            comb['load_title'] = row.load_title
            #
            dftemp.append(comb)
    #
    dftemp = db.concat(dftemp, ignore_index=True)
    #
    try:
        dftemp = dftemp.groupby(['load_name', 'load_id' ,'load_level',
                                 'load_title', 'load_system' ,'node_name'],
                                as_index=False)[values].sum()
        dfnode = db.concat([dfnode, dftemp], ignore_index=True)
    except UnboundLocalError:
        pass
    #
    return dfnode
#

def update_member_dfXX(member, combination,
                     item: list[str],
                     values:list[str]):
    """
    Update node displacements to include combination
    """
    db = DBframework()
    # group basic load by name
    grp = member.groupby('load_name')
    comb_grp = combination.groupby('load_name')
    temp = []
    for key, comb_factors in comb_grp:
        for row in comb_factors.itertuples():
            comb = grp.get_group(row.basic_load).copy()
            comb.loc[:, values] *= row.factor
            comb['load_level'] = 'combination'
            comb['load_name'] = row.load_name
            comb['load_id'] = row.load_id
            comb['load_title'] = row.load_title
            temp.append(comb)
    #
    try:
        temp = db.concat(temp, ignore_index=True)
        temp = temp.groupby(item, as_index=False)[values].sum()
        member = db.concat([member, temp], ignore_index=True)
    except (UnboundLocalError, ValueError):
        pass
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