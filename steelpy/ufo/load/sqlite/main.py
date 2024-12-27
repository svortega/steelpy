#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from typing import NamedTuple

# package imports
from steelpy.ufo.load.process.main import MasterLoad
from steelpy.ufo.mesh.sqlite.node import NodeSQL
from steelpy.ufo.mesh.sqlite.element import ElementsSQL
from steelpy.ufo.load.sqlite.wave_load import MetoceanLoadSQL
#
from steelpy.ufo.load.sqlite.load_case import BasicLoadSQL
from steelpy.ufo.load.sqlite.combination import LoadCombSQL
from steelpy.utils.sqlite.utils import create_connection, create_table
#
#
#
#
class MeshLoad(MasterLoad):
    """
    """
    __slots__ = ['_df_nodal', '_mesh_id',
                 '_hydro', '_basic', '_combination',
                 '_elements', '_nodes', '_boundaries',
                 '_db_file']

    def __init__(self,
                 mesh_type: str,
                 mesh_id: str|int,
                 name: str|str,
                 db_file: str | None = None):
        """
        """
        super().__init__(name=name)
        self._mesh_id = mesh_id
        self._db_file = db_file
        #
        # mesh
        self._nodes = NodeSQL(db_system=mesh_type,
                              mesh_id=mesh_id,
                              db_file=db_file,
                              name=name)
        #
        self._elements = ElementsSQL(db_system=mesh_type,
                                     mesh_id=mesh_id,
                                     db_file=db_file,
                                     name=name)
        #
        #
        self._basic = BasicLoadSQL(db_file,
                                   mesh_id=mesh_id,
                                   name=name)
        #
        self._combination = LoadCombSQL(db_file,
                                        mesh_id=mesh_id,
                                        name=name)
        #
        self._hydro = MetoceanLoadSQL(db_file=db_file,
                                      mesh_id=mesh_id,
                                      name=name)
        #
        conn = create_connection(self._db_file)
        with conn:
            self._new_table(conn)
    #
    #
    #def __str__(self) -> str:
    #    """ """
    #    output = "\n"
    #    output += self._basic.__str__()
    #    output += self._combination.__str__()
    #    return output
    #
    # ----------------------------
    # sql operations
    def _new_table(self, conn):
        """ """
        # Basic load
        table = "CREATE TABLE IF NOT EXISTS Load(\
                number INTEGER PRIMARY KEY NOT NULL,\
                name NOT NULL,\
                mesh_id INTEGER NOT NULL REFERENCES Mesh(number), \
                level TEXT NOT NULL,\
                title TEXT,\
                input_type TEXT, \
                input_file TEXT);"
        create_table(conn, table)
        # Step Load
        table = "CREATE TABLE IF NOT EXISTS LoadStep(\
                number INTEGER PRIMARY KEY NOT NULL,\
                load_id INTEGER NOT NULL REFERENCES Load(number), \
                step_type TEXT NOT NULL,\
                step NOT NULL);"
        create_table(conn, table)
    #
    #
    #
    #
#