#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from dataclasses import dataclass

#
# package imports
from steelpy.utils.sqlite.utils import create_table, create_connection
from steelpy.ufo.mesh.sqlite.utils import (pull_node_mesh,
                                           pull_results_mesh,
                                           pull_load_mesh)
from steelpy.utils.dataframe.main import DBframework
#
#
#
#
@dataclass
class UnSQL:
    """ """
    __slots__ = ['_mesh', '_db_file', '_result_name']
    def __init__(self, mesh, result_name:int|str,
                 db_file:str) -> None:
        self._mesh = mesh
        self._db_file = db_file
        self._result_name = result_name
        #
        conn = create_connection(self._db_file)
        with conn:
            self._new_table(conn)
    #
    # ---------------------------------
    #
    @property
    def df(self):
        """ """
        conn = create_connection(self._db_file)
        with conn:
            ndata = self._pull_node_displacement(conn)
        # dataframe setup
        cols = ['mesh_name', 'result_name',
                'load_name',  'load_level',
                'node_name', 'system',
                'x', 'y', 'z', 'rx', 'ry', 'rz']
        db = DBframework()
        Un = db.DataFrame(data=ndata, columns=cols)
        Un = Un.astype({'x': 'float64', 'y': 'float64', 'z': 'float64',
                        'rx': 'float64', 'ry': 'float64', 'rz': 'float64'}).fillna(value=0.0)
        return Un

    @df.setter
    def df(self, Un:DBframework):
        """ """
        header = ['result_id', 'node_id',
                  'load_id', 'system',
                  *self._mesh._plane.hdisp]
        #          #'x', 'y', 'z', 'rx', 'ry', 'rz']
        conn = create_connection(self._db_file)
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
    def _update_ndf(self, dfnode, dfcomb,
                    values :list[str]|None = None):
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
        #
        combgrp = dfcomb.groupby('load_name')
        for key, combfactors in combgrp:
            for row in combfactors.itertuples():
                comb = ndgrp.get_group(row.basic_load).copy()
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
            # check = dftemp.groupby(['node_name', 'c']).sum().reset_index()
            dftemp = dftemp.groupby(['load_name', 'load_id' ,'load_level',
                                     'load_title', 'load_system' ,'node_name'],
                                    as_index=False)[values].sum()
            # test
            dfnode = db.concat([dfnode, dftemp], ignore_index=True)
        except UnboundLocalError:
            pass
        #
        return dfnode # , memb_comb
    #
    # ---------------------------------
    # SQL ops
    #
    def _new_table(self, conn) -> None:
        """ """
        # Node displacement solution
        table_nodes = "CREATE TABLE IF NOT EXISTS ResultNodeDisplacement (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        result_id INTEGER NOT NULL REFERENCES Result(number),\
                        node_id NOT NULL REFERENCES Node(number),\
                        load_id INTEGER NOT NULL REFERENCES Load(number),\
                        system TEXT NOT NULL,\
                        x DECIMAL,\
                        y DECIMAL,\
                        z DECIMAL,\
                        rx DECIMAL,\
                        ry DECIMAL,\
                        rz DECIMAL);"
        #
        create_table(conn, table_nodes)
    #
    #
    def _pull_node_displacement(self, conn):
        """ """
        table = ("SELECT Mesh.name, Result.name, Load.name, Load.Level, \
                         Node.name, ResultNodeDisplacement.system, \
                         ResultNodeDisplacement.x, ResultNodeDisplacement.y, ResultNodeDisplacement.z,\
                         ResultNodeDisplacement.rx, ResultNodeDisplacement.ry, ResultNodeDisplacement.rz \
                 FROM ResultNodeDisplacement, Load, Mesh, Node, Result \
                 WHERE Load.number = ResultNodeDisplacement.load_id \
                 AND Node.number = ResultNodeDisplacement.node_id \
                 AND Result.number = ResultNodeDisplacement.result_id \
                 AND Result.mesh_id = Mesh.number ;")
        #
        cur = conn.cursor()
        cur.execute(table)
        rows = cur.fetchall()
        return rows
#
#