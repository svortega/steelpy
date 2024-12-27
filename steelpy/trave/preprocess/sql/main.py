#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from datetime import datetime as dt

# package imports
#
#from steelpy.trave.beam.main import Beam
#from steelpy.utils.dataframe.main import DBframework
from steelpy.ufo.mesh.sqlite.utils import pull_mesh_id
from steelpy.utils.sqlite.utils import create_connection, create_table
#
class TraveItemSQL:
    """
    A program for static & dynamic analysis
    of 3-d framed structures
    """
    __slots__ = ['_mesh', 'db_file',
                 '_log'] # '_postprocess', '_name',

    def __init__(self, mesh, #name: str,
                 sql_file: str,
                 log: bool = False) -> None:
        """
        """
        #self._name = name
        self._mesh = mesh
        self.db_file = sql_file
        self._log = log
        #
        # create results table
        conn = create_connection(self.db_file)
        with conn:
            self._new_table(conn)
    #
    # --------------------------------------------
    # SQL ops
    #
    def _new_table(self, conn) -> None:
        """ """
        table = "CREATE TABLE IF NOT EXISTS Result (\
                    number INTEGER PRIMARY KEY NOT NULL,\
                    name NOT NULL,\
                    mesh_id INTEGER NOT NULL REFERENCES Mesh(number), \
                    type TEXT NOT NULL, \
                    p_delta BOOL NOT NULL, \
                    ineleastic BOOL NOT NULL, \
                    plane TEXT NOT NULL, \
                    date TEXT NOT NULL,\
                    units TEXT NOT NULL,\
                    title TEXT);"
        create_table(conn, table)
    #
    #
    def _push_data(self, conn,
                   name: int |str,
                   mesh_id: int,
                   analysis_type: str,
                   Pdelta: bool,
                   ineleastic:bool,
                   plane: str,
                   title: str |None = None):
        """ """
        table = 'INSERT INTO Result(name, mesh_id, \
                                    type, p_delta, ineleastic,\
                                    plane, date, units, title)\
                            VALUES(?,?,?,?,?,?,?,?,?)'
        #
        date = dt.now().strftime('%Y-%m-%d')
        data = (name, mesh_id,
                analysis_type, Pdelta, ineleastic,
                plane, date, 'si', title)
        # push
        cur = conn.cursor()
        out = cur.execute(table, data)
        return out.lastrowid
    #
    # --------------------------------------------
    #
    def _push_analysis(self, name:str,
                       analysis_type:str,
                       plane:str,
                       Pdelta: bool,
                       ineleastic: bool):
        """
        push static data
        Pdelta : Second order (True/False)
        ineleastic : non-linear (True/False)
        """
        # update table
        conn = create_connection(self.db_file)
        with conn:
            #component = self._mesh._component
            #mesh_name = self._mesh._name
            #mesh_id = pull_mesh_id(conn,
            #                       name=mesh_name,
            #                       component=component)
            #
            mesh_id = self._mesh._id
            result_id = self._push_data(conn, name, mesh_id,
                                        analysis_type, Pdelta,
                                        ineleastic, plane)
        return result_id
#
#